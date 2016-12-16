import re
import fileinput
import string
from collections import defaultdict
from seqc import reader
from intervaltree import IntervalTree


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


def first(iterable):
    """simple function to return first object in an iterable"""
    return next(iter(iterable))


class GeneIntervals:
    """
    Encodes genomic ranges in an Intervaltree

    :method translate: translates a genomic coordinate on a stranded chromosome into the
      gene identifier that occupies that location (if any exists)

    """

    def __init__(self, gtf: str, max_transcript_length=1000):
        """Construct a dictionary containing genomic intervals that map to genes. Allows
        the translation of alignment coordinates (chromosome, strand, position) to
        gene identifiers.

        The returned object is accessible by the translate() method, also defined for this
        object.

        :param gtf: annotation file in GTF format. Can be gz or bz2 compressed
        :param max_transcript_length: Indicates that for an alignment to be valid, it
          must align within this (spliced) distances from a transcript termination site.
          This logic stems from 3' sequencing methods, wherein the distribution of
          transcript sizes indicates that the majority of non-erroneous fragments of
          mRNA molecules should align within this region.
        """
        self._chromosomes_to_genes = self.construct_translator(gtf, max_transcript_length)

    @staticmethod
    def iterate_adjusted_exons(exons, strand, max_transcript_length):
        """
        :param list exons: a list of exon records from a .gtf file
        :param strand: the strand of the exons, options = ['+', '-']
        :param int max_transcript_length: maximum allowable distance of an alignment from
          a transcription termination site
        :yield (int, int): tuple of (exon start, exon end). If the transcript exceeds
          maximum transcript length, then this list will be truncated, and the final
          exon may be shortened.
        """
        for exon in exons:  # closest exon to TTS is first (see Reader.iter_transcripts)
            start, end = int(exon[3]), int(exon[4])
            size = end - start
            if size >= max_transcript_length:
                if strand == '+':
                    yield end - max_transcript_length, end
                else:
                    yield start, start + max_transcript_length
                break  # we've exhausted the allowable transcript length
            else:
                yield start, end
                max_transcript_length -= size

    # todo implement me.
    @staticmethod
    def _remove_overlapping_intervals(dictionary):
        """
        The intervaltree is a dictionary-based structure, so automatically stores a
        unique set of intervals.

        we could parse the tree interval by interval by calling overlap queries using
        each interval. Then, we could update based on the result.

        removing duplicates within the same gene is very easy; we simply replace each
        duplicate with their union, as all alignments in these intervals correspond to a
        unique assignment.

        removing overlaps between exons of different genes can be done using symmetric
        difference, because we want to eliminate overlaps that are ambiguous between more
        than one gene

        What follows is a skeleton that could be used to parse intervals. I suspect using
        bisection may be the quickest way to find out which category (above) an interval
        is in. However, optimization of this function is not critical, as it need only
        be called once per run.

        :param dict dictionary: result of self.construct_translator()
        :return dict: dictionary whose intervals correspond to unique gene assignments
        """

        def parse_intervals(query, other, tree):
            """private function to parse intervals, should side-effect to modify
            interval_tree

            :param query: the interval tested for overlaps
            :param set other: set of Interval objects which overlapped the query
              interval
            :param tree: the tree which should be modified
            :return None: side-effecting function
            """
            raise NotImplementedError

        # get each intervaltree # todo can be parallelized
        for chromosome in dictionary:
            for strand, interval_tree in dictionary[chromosome].items():
                for iv in interval_tree:
                    intervals = interval_tree.search(iv.begin, iv.end)
                    if intervals:
                        parse_intervals(iv, intervals, interval_tree)

    def construct_translator(self, gtf, max_transcript_length):
        """Construct a dictionary containing genomic intervals that map to genes. Allows
        the translation of alignment coordinates (chromosome, strand, position) to
        gene identifiers.

        The returned object is accessible by the translate() method, also defined for this
        object.

        :param gtf: annotation file in GTF format. Can be gz or bz2 compressed
        :param max_transcript_length: Indicates that for an alignment to be valid, it
          must align within this (spliced) distances from a transcript termination site.
          This logic stems from 3' sequencing methods, wherein the distribution of
          transcript sizes indicates that the majority of non-erroneous fragments of
          mRNA molecules should align within this region.

        :return dict: {chromosome: {strand: {position: gene}}}, nested dictionary of
          chromosome -> strand -> position which returns a gene.
        """
        results_dictionary = defaultdict(dict)
        for (tx_chromosome, tx_strand, gene_id), exons in Reader(gtf).iter_transcripts():
            for start, end in self.iterate_adjusted_exons(
                    exons, tx_strand, max_transcript_length):
                if start == end:
                    continue  # zero-length exons apparently occur in the gtf
                try:
                    results_dictionary[tx_chromosome][tx_strand].addi(
                        start, end, gene_id)
                except KeyError:
                    results_dictionary[tx_chromosome][tx_strand] = IntervalTree()
                    results_dictionary[tx_chromosome][tx_strand].addi(
                        start, end, gene_id)
        return dict(results_dictionary)

    def translate(self, chromosome, strand, pos):
        """translates a chromosome, position, and strand into a gene identifier

        Uses the IntervalTree data structure to rapidly search for the corresponding
        identifier.

        :param bytes chromosome: chromosome for this alignment
        :param bytes strand: strand for this alignment (one of ['+', '-'])
        :param int pos: position of the alignment within the chromosome
        :return int|None: Returns either an integer gene_id if a unique gene was found
          at the specified position, or None otherwise
        """
        # todo remove duplicate exons during construction to save time
        try:
            result = set(x.data for x in
                         self._chromosomes_to_genes[chromosome][strand][pos])
            if len(result) == 1:
                return first(result)  # just right
            else:
                return None  # too many genes
        except KeyError:
            return None  # no gene


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
        with fileinput.input(self._files, openhook=hook, mode='r') as f:

            # get rid of header lines
            file_iterator = iter(f)
            first_record = next(file_iterator)
            while first_record.startswith('#'):
                first_record = next(file_iterator)

            yield first_record.split('\t')  # avoid loss of first non-comment line

            for record in file_iterator:  # now, run to exhaustion
                yield record.split('\t')

    @staticmethod
    def strip_gene_num(attribute_str):
        try:
            gene_start = attribute_str.index('gene_id')
        except ValueError:
            raise ValueError(
                'Gene_id field is missing in annotations file: {}'.format(attribute_str))

        try:
            gene_end = attribute_str.index('";', gene_start)
        except ValueError:
            raise ValueError(
                'no "; in gene_id attribute, gtf file might be corrupted: {}'.format(
                    attribute_str))

        try:
            id_start = attribute_str.index('0', gene_start)
        except ValueError:
            raise ValueError(
                'Corrupt gene_id field in annotations file - {}'.format(attribute_str))

        # ignore the gene version, which is located after a decimal in some gtf files
        try:
            gene_end = attribute_str.index('.', id_start, gene_end)
        except ValueError:
            pass

        return int(attribute_str[id_start:gene_end])

    def iter_transcripts(self):
        """Iterate over transcripts in a gtf file, returning a transcripts's chromosome
        strand, gene_id, and a list of tab-split exon records

        :yield (str, str, int), [[str], ...]: (chromosome, strand, id), [exon_records]
        """
        iterator = iter(self)
        record = next(iterator)

        # skip to first transcript record and store chromosome and strand
        while record[2] != 'transcript':
            record = next(iterator)
        transcript_chromosome = record[0]
        transcript_strand = record[6]
        transcript_gene_id = self.strip_gene_num(record[8])

        # Notes on the following loop:
        # ----------------------------
        # - gtf files are ordered gene -> transcript1 -> exons1 -> transcript2 -> exons2
        #
        # - therefore, the below loop will begin by parsing the exons for the first
        #   transcript identified above, since those are next in sequence.
        #
        # - Then, if the gene had only one transcript, the else clause will skip the next
        #   gene record, and will arrive at the following transcript record. Otherwise,
        #   a second transcript of the same gene will directly follow.
        #
        # - The loop yields the chromosome, strand and gene_id from the first transcript,
        #   stores the new transcript values, and zeros out the exon list. Then the loop
        #   continues.
        #
        # - Downstream, we want to process the exons from the closest exon to the TTS to
        #   the most distant. The minus strand is already in this order, but we need to
        #   invert the plus strand exons to obtain this ordering.

        exons = []
        for record in iterator:
            if record[2] == 'exon':
                exons.append(record)
            elif record[2] == 'transcript':
                # we want exons in inverse order, - is already in inverse order.
                if transcript_strand == '+':
                    exons = exons[::-1]
                yield (
                    (transcript_chromosome, transcript_strand, transcript_gene_id), exons)
                exons = []
                transcript_chromosome = record[0]
                transcript_strand = record[6]
                transcript_gene_id = self.strip_gene_num(record[8])

        # yield the final transcript
        if transcript_strand == '+':
            exons = exons[::-1]
        yield (transcript_chromosome, transcript_strand, transcript_gene_id), exons


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
            # Skip comment lines
            if line.startswith('#'):
                continue

            fields = line.split('\t')  # speed-up, only run regex on gene lines
            if fields[2] != 'gene':
                continue

            match = re.match(pattern, line)  # run regex
            if match:
                gene_id_map[int(match.group(2))].add(match.group(4).upper())
    return gene_id_map


def ensembl_gene_id_to_official_gene_symbol(ids, gene_id_map):
    """convert data containing ensembl gene ids into an index of gene symbols

    :param Iterable ids: an iterable containing integer ids to be converted
    :param gene_id_map: gene_id_map constructed from
      Experiment.create_gene_id_to_official_gene_symbol_map. If converting multiple
      objects, it is much faster to only construct the map a single time.
    :return list: converted ids
    """
    return ['-'.join(gene_id_map[i]) for i in ids]
