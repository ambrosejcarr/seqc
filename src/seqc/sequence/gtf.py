import re
import fileinput
import string
from collections import defaultdict
from intervaltree import IntervalTree
from copy import deepcopy
import pandas as pd
from seqc import reader
from seqc.sparse_frame import SparseFrame


def first(iterable):
    return next(iter(iterable))


class Record:
    """
    Simple namespace object that makes the fields of a GTF record available. Subclassed
    to create records specific to exons, transcripts, and genes
    """

    __slots__ = ['_fields', '_attribute']

    _del_letters = string.ascii_letters.encode()
    _del_non_letters = ''.join(set(string.printable).difference(string.ascii_letters))\
        .encode()

    def __init__(self, fields: list):

        self._fields = fields
        self._attribute = {}

    def __repr__(self) -> str:
        return '<Record: %s>' % bytes(self).decode()

    def __bytes__(self) -> bytes:
        return b'\t'.join(self._fields)

    def _parse_attribute(self) -> None:
        for field in self._fields[8].rstrip(b';\n').split(b';'):
            key, *value = field.strip().split()
            self._attribute[key] = b' '.join(value).strip(b'"')

    def __hash__(self) -> int:
        """concatenate strand, start, end, and chromosome and hash the resulting bytes"""
        return hash(self._fields[6] + self._fields[3] + self._fields[4] + self._fields[0])

    @property
    def seqname(self) -> bytes:
        return self._fields[0]

    @property
    def chromosome(self) -> bytes:
        return self._fields[0]  # synonym for seqname

    @property
    def source(self) -> bytes:
        return self._fields[1]

    @property
    def feature(self) -> bytes:
        return self._fields[2]

    @property
    def start(self) -> int:
        return int(self._fields[3])

    @property
    def end(self) -> int:
        return int(self._fields[4])

    @property
    def score(self) -> bytes:
        return self._fields[5]

    @property
    def strand(self) -> bytes:
        return self._fields[6]

    @property
    def frame(self) -> bytes:
        return self._fields[7]

    @property
    def size(self) -> int:
        return self.end - self.start

    @property
    def fields(self) -> list:
        return self._fields

    def attribute(self, item):
        """
        access an item from the attribute field of a GTF file.

        :param item: item to access
        :return: value of item
        """
        try:
            return self._attribute[item]
        except KeyError:
            if not self._attribute:
                self._parse_attribute()
                return self._attribute[item]
            else:
                raise KeyError('%s is not a stored attribute of this gtf record' %
                               repr(item))

    @property
    def integer_gene_id(self) -> int:
        """ENSEMBL gene id without the organism specific prefix, encoded as an integer"""
        return int(self.attribute(b'gene_id').split(b'.')[0]
                   .translate(None, self._del_letters))

    @property
    def organism_prefix(self) -> bytes:
        """Organism prefix of ENSEMBL gene id (e.g. ENSG for human, ENSMUSG)"""
        return self.attribute(b'gene_id').translate(None, self._del_non_letters)

    @property
    def string_gene_id(self) -> bytes:
        """ENSEMBL gene id, including organism prefix."""
        return self.attribute(b'gene_id')

    @staticmethod
    def int2str_gene_id(integer_id: int, organism_prefix: bytes) -> bytes:
        """
        converts an integer gene id (suffix) to a string gene id (including organism-
        specific suffix)
        :param organism_prefix: bytes
        :param integer_id: int
        """
        bytestring = str(integer_id).encode()
        diff = 11 - len(bytestring)
        return organism_prefix + (b'0' * diff) + bytestring

    def __eq__(self, other):
        """equivalent to testing if start, end, chrom and strand are the same."""
        return hash(self) == hash(other)

    def __ne__(self, other):
        return not self.__eq__(other)


class Exon(Record):

    def __repr__(self) -> str:
        return '<Exon: %s>' % bytes(self).decode()


class Gene(Record):
    """
    Gene record object. In addition to base properties of the Record object, provides:

    :property exons: set, exon objects associated with this gene
    :method genomic_intervals: iterator, yields tuples of sequence that are covered by
      this gene.

    """

    __slots__ = ['_exons']

    def __init__(self, fields: list):
        super().__init__(fields)
        self._exons = set()

    @property
    def exons(self) -> set:
        return self._exons

    def genomic_intervals(self):
        """
        Iterator, yields tuples of genomic sequence that codes for this gene. In cases
        where multiple exons overlap, the union of the intervals is returned.
        """
        assert self.exons
        ivs = sorted(((e.start, e.end) for e in self.exons))
        ivs = sorted(ivs)
        saved = list(ivs[0])
        for st, en in ivs:
            if st <= saved[1]:
                saved[1] = max(saved[1], en)
            else:
                yield tuple(saved)
                saved[0] = st
                saved[1] = en
        yield tuple(saved)


class GeneIntervals:
    """
    Encodes genomic ranges in an Intervaltree

    :method translate: translates a genomic coordinate on a stranded chromosome into the
      gene identifier that occupies that location (if any exists)

    """

    def __init__(self, gtf: str):
        interval_tree = defaultdict(IntervalTree)
        reader_ = Reader(gtf)
        for gene in reader_.iter_genes():
            if gene.exons:
                for start, end in gene.genomic_intervals():
                    try:
                        interval_tree[(gene.chromosome, gene.strand)].addi(
                                start, end, gene.integer_gene_id)
                    except ValueError:
                        if start == end:  # ensure this is the reason the error was raised
                            continue
                        else:
                            raise

        self._interval_tree = interval_tree

    def translate(self, chromosome: bytes, strand: bytes, position: int):
        """
        translates a genomic coordinate on a stranded chromosome into the gene identifier
        that occupies that location. If no gene ids or an ambiguous location with more
        than one gene id is identified, returns None.

        :param chromosome: bytes, chromosome id
        :param strand: bytes, [b'+', b'-'], the strand
        :param position: int, position on the chromosome
        :return: int, gene id.
        """
        ivs = self._interval_tree[(chromosome, strand)].search(position)
        if len(ivs) == 1:
            return first(ivs).data
        else:
            return None


class Reader(reader.Reader):
    """
    SubClass of reader.Reader, returns an Reader with several specialized iterator
    methods.

    :method __iter__: Iterator over all non-header records in gtf; yields Record objects.
    :method iter_genes: Iterator over all genes in gtf; yields Gene objects.
    """

    def __iter__(self):
        """return an iterator over all non-header records in gtf"""
        hook = fileinput.hook_compressed
        with fileinput.input(self._files, openhook=hook, mode='rb') as f:

            # get rid of header lines
            file_iterator = iter(f)
            first_record = next(file_iterator)
            while first_record.startswith(b'#'):
                first_record = next(file_iterator)
                continue
            yield first_record.split(b'\t')

            for record in file_iterator:  # now, run to exhaustion
                yield record.split(b'\t')

    def iter_genes(self):
        """iterate over all the records for each gene in passed gtf"""

        records = iter(self)

        # get the first gene record
        record = next(records)
        while record[2] != b'gene':
            record = next(records)

        # aggregate exons for each gene
        gene = Gene(record)
        record = next(records)
        while record:
            if record[2] == b'exon':
                gene.exons.add(Exon(record))
            elif record[2] == b'gene':
                yield gene
                gene = Gene(record)
            record = next(records)

        # yield the final gene record
        yield gene


def create_phix_annotation(phix_fasta):
    """
    Several tools in this package require ENSEMBL formatting for .gtf files. However,
    the PhiX genome provided by NCBI does not come with a .gtf file. This tool creates
    a companion gtf file for the phiX fasta file

    :param phix_fasta: str, name of the phix fasta file.
    """
    import numpy as np

    with open(phix_fasta, 'r') as f:
        header = f.readline()  # phiX has only one chromosome
        data = f.readlines()

    # concatenate data
    contig = ''
    for line in data:
        contig += line.strip()

    # get chromosome
    chromosome = header.split()[0].strip('>')
    source = 'seqc'
    score = '.'
    frame = '.'
    gene_meta = 'gene_id "PHIXG00{NUM}"; gene_name "PHIX{NAME!s}";'
    exon_meta = ('gene_id "PHIXG00{NUM}"; gene_name "PHIX{NAME!s}"; '
                 'exon_id "PHIX{NAME!s}";')

    # SEQC truncates genes at 1000b from the end of each transcript. However, phiX DNA
    # that is spiked into an experiment is not subject to library construction. Thus,
    # we will create artificial transcripts for phiX that ensure that all of the DNA is
    # correctly identified.

    length = len(contig)
    transcript_starts = np.arange(length // 1000 + 1) * 1000
    transcript_ends = np.array([min(s + 1000, length) for s in transcript_starts])

    phix_gtf = phix_fasta.replace('.fa', '.gtf')

    with open(phix_gtf, 'w') as f:
        for i, (s, e) in enumerate(zip(transcript_starts, transcript_ends)):
            # add forward strand gene
            gene = [chromosome, source, 'gene', str(s), str(e), score, '+', frame,
                    gene_meta.format(NUM=str(i + 1) * 9, NAME=i + 1)]
            f.write('\t'.join(gene) + '\n')
            exon = [chromosome, source, 'exon', str(s), str(e), score, '+', frame,
                    exon_meta.format(NUM=str(i + 1) * 9, NAME=i + 1)]
            f.write('\t'.join(exon) + '\n')
            # add reverse strand gene
            gene = [chromosome, source, 'gene', str(s), str(e), score, '-', frame,
                    gene_meta.format(NUM=str(i + 1) * 9, NAME=i + 1)]
            f.write('\t'.join(gene) + '\n')
            exon = [chromosome, source, 'exon', str(s), str(e), score, '-', frame,
                    exon_meta.format(NUM=str(i + 1) * 9, NAME=i + 1)]
            f.write('\t'.join(exon) + '\n')


def create_gene_id_to_official_gene_symbol_map(gtf: str):
    """
    create_gene_id_to_official_gene_symbol_map: map integer ENSEMBL ids to
    official gene symbols.

    :param gtf: str, filename of gtf file from which to create the map.
    """
    pattern = re.compile(
        r'(^.*?gene_id "[^0-9]*)([0-9]*)(\.?.*?gene_name ")(.*?)(".*?$)')
    gene_id_map = defaultdict(set)
    with open(gtf, 'r') as f:
        for line in f:
            match = re.match(pattern, line)
            if match:
                gene_id_map[int(match.group(2))].add(match.group(4).upper())
    return gene_id_map


def ensembl_gene_id_to_official_gene_symbol(
        data: pd.DataFrame or SparseFrame, gtf=None, gene_id_map=None):
    """convert data containing ensembl gene ids into an index of gene symbols

    :param data: cells x genes DataFrame or SparseFrame
    :param gtf: str, filename of a gtf file
    :param gene_id_map: gene_id_map constructed from
      Experiment.create_gene_id_to_official_gene_symbol_map. If converting multiple
      objects, it is much faster to only construct the map a single time.
    :return: data: pd.DataFrame or SparseFrame object
    """
    if gene_id_map is None:
        if gtf is None:
            raise ValueError('User must pass either GTF or a gene_id_map object')
        gene_id_map = create_gene_id_to_official_gene_symbol_map(gtf)
    data = deepcopy(data)
    data.columns = ['-'.join(gene_id_map[i]) for i in data.columns]

    return data
