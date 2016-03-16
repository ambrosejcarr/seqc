import fileinput
import string
from collections import defaultdict

from intervaltree import IntervalTree

import seqc


def first(iterable):
    return next(iter(iterable))


class Record:

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
        return int(self.attribute(b'gene_id').split(b'.')[0]
                   .translate(None, self._del_letters))

    @property
    def organism_prefix(self) -> bytes:
        return self.attribute(b'gene_id').translate(None, self._del_non_letters)

    @property
    def string_gene_id(self) -> bytes:
        return self.attribute(b'gene_id')

    @staticmethod
    def int2str_gene_id(integer_id: int, organism_prefix: bytes) -> bytes:
        """converts an integer gene id to a string gene id
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

    __slots__ = ['_exons']

    def __init__(self, fields: list):
        super().__init__(fields)
        self._exons = set()

    @property
    def exons(self) -> set:
        return self._exons

    def genomic_intervals(self):
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

    def __init__(self, gtf: str):
        interval_tree = defaultdict(IntervalTree)
        reader = Reader(gtf)
        for gene in reader.iter_genes():
            if gene.exons:
                for start, end in gene.genomic_intervals():
                    try:
                        interval_tree[(gene.chromosome, gene.strand)].addi(
                                start, end, gene.integer_gene_id)
                    except ValueError:
                        if start == end: # ensure this is the reason the error was raised
                            continue
                        else:
                            raise

        self._interval_tree = interval_tree

    def translate(self, chromosome: bytes, strand: bytes, position: int):
            ivs = self._interval_tree[(chromosome, strand)].search(position)
            if len(ivs) == 1:
                return first(ivs).data
            else:
                return None


class Reader(seqc.reader.Reader):

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
        """iterate over all the records for each gene in passed gtf

        :yields: Gene
        """

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
