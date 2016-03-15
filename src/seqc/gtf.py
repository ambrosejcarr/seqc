import fileinput
import string
from collections.abc import Iterator
from copy import deepcopy
from functools import lru_cache
from operator import attrgetter
from random import randint
from sys import maxsize
import numpy as np
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
        return '<Record: %s>' % str(self)

    def __str__(self) -> str:
        return b'\t'.join(self._fields).decode()

    def __bytes__(self) -> bytes:
        return b'\t'.join(self._fields)

    def __getitem__(self, slice_object):
        """
        Implements slicing of Record objects. In instances where slices go over the
        bounds of the original object, the item will return only up the the end/beginning
        of the original object.

        usage:
        ------
        if record.start = 10, record.end = 2010:
        record[-1000:] -> returns a new record object with only the last 1000 bases of
         the original record: record.start = 1010, record.end = 2010
        record[500:1000] -> returns a new record object with only the 500 to 1000th bases:
         record.start = 500, record.end = 1000

        if record.start = 50, record.end = 500:
        record[-1000:] -> returns a new record object with only the last 1000 bases of
         the original record: record.start = 50, record.end = 500
        record[500:1000] -> raises IndexError, no item falls within the object.

        returns:
        --------
        Record object with new coordinates
        """
        fields = deepcopy(self._fields)

        # if exon is minus stranded, reverse slicing direction

        # -1000: -> :1000
        # 10:50 -> -10:-50

        if self.strand == b'-':
            try:
                start = -slice_object.stop
            except TypeError:
                start = None
            try:
                stop = -slice_object.start
            except TypeError:
                stop = None
        else:
            start = slice_object.start
            stop = slice_object.stop

        try:
            if start:
                if start > 0:
                    fields[3] = str(self.start + start).encode()
                else:
                    fields[3] = str(self.end + start).encode()

            if stop:
                if stop > 0:
                    fields[4] = str(self.start + stop).encode()
                else:
                    fields[4] = str(self.end + stop).encode()
        except AttributeError:
            raise TypeError('slice_object must be %s' % repr(slice))

        assert(fields[3] < fields[4])

        # todo
        # this might break; trying to find a way to make getitem inheritable to subclasses
        return Record(fields)

    def _parse_attribute(self) -> None:
        for field in self._fields[8].rstrip(b';\n').split(b';'):
            key, *value = field.strip().split()
            self._attribute[key] = b' '.join(value).strip(b'"')

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
    def gene_name(self) -> bytes:
        return self.attribute(b'gene_name')

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
    def int2str_gene_id(integer_id, organism_prefix) -> bytes:
        """converts an integer gene id to a string gene id
        :param integer_id:  convert integer geneid to string
        """
        bytestring = str(integer_id).encode()
        diff = 11 - len(bytestring)
        return organism_prefix + (b'0' * diff) + bytestring

    def __eq__(self, other):
        return (self.start == other.start and self.end == other.end and
                self.strand == other.strand and self.chromosome == other.chromosome)

    def __ne__(self, other):
        return not self.__eq__(other)

    def __hash__(self):
        return id(self)


class Exon(Record):

    def __repr__(self) -> str:
        return '<Exon: %s>' % str(self)


class Transcript(Record):

    __slots__ = ['_exons']

    def __init__(self, transcript: list, exons: list):

        # define transcript
        super().__init__(transcript)

        # sorting necessary to get - strand transcripts' UTRs in correct order
        if transcript[6] == b'+':
            self._exons = exons
        else:
            self._exons = sorted(exons, key=attrgetter('start'), reverse=True)

    def __iter__(self) -> Iterator:
        return iter(self.exons)

    def __len__(self) -> int:
        return len(self.exons)

    def __repr__(self) -> str:
        return '<Transcript (composed of %d exon(s) and UTR(s)): %s>' % \
               (len(self), str(self)[:-1])  # remove newline from str(self) for repr

    def __getitem__(self, slice_object: slice):
        """
        slice the exons stored in self, utilizing the exon slicing routine. If no exons
        result from the passed slice, returns None."""

        # get absolute size of exons
        exons = deepcopy(self.exons)
        diffs = list((e.end - e.start) for e in exons)
        total_length = sum(diffs)

        # convert any relative (e.g. x[-1000:]) to absolute positions
        if slice_object.start:
            if slice_object.start > 0:
                start = slice_object.start
            else:
                start = total_length + slice_object.start
        else:
            start = 0

        if slice_object.stop:
            if slice_object.stop > 0:
                stop = slice_object.stop
            else:
                stop = total_length + slice_object.stop
        else:
            stop = maxsize

        # move through exons, trimming according to the slicer
        delete = []
        position = 0
        i = -1
        exon_iterator = iter(exons)
        if start:
            while position <= start:
                try:
                    exon = next(exon_iterator)
                    i += 1
                except StopIteration:
                    break
                if position + diffs[i] <= start:
                    delete.append(i)
                elif start <= position + diffs[i] < stop:
                    exons[i] = exon[start - position:]
                elif start < stop < position + diffs[i]:
                    exons[i] = exon[start - position: stop - position]
                else:
                    raise ValueError('Invalid absolute slicer! %d:%d, %s, %d' %
                                     (start, stop, repr(slice_object), self.size))
                position += diffs[i]
        if stop:
            while position < stop:
                try:
                    exon = next(exon_iterator)
                    i += 1
                except StopIteration:
                    break

                if position + diffs[i] < stop:
                    pass
                else:
                    exons[i] = exon[:stop - position]
                position += diffs[i]

        while True:  # add all remaining exons to the delete queue
            try:
                _ = next(exon_iterator)
                i += 1
                delete.append(i)
            except StopIteration:
                break

        # delete exons
        for i in delete[::-1]:
            del exons[i]

        if not exons:
            return None

        # get new start/stop for transcript
        tx_fields = deepcopy(self._fields)
        tx_fields[3] = str(min(e.start for e in exons)).encode()
        tx_fields[4] = str(max(e.end for e in exons)).encode()

        return Transcript(tx_fields, exons)

    @classmethod
    def from_records(cls, transcript, *exons):
        exons = [Exon(record) for record in exons]
        return cls(transcript, exons)

    @property
    def exons(self) -> list:
        return self._exons

    @property
    def size(self) -> int:
        return sum(exon.end - exon.start for exon in self)

    @property
    def TTS(self) -> int:
        return self.end


class Gene(Record):

    __slots__ = ['_transcripts']

    def __init__(self, gene: list, transcripts: list):

        # define gene
        super().__init__(gene)
        self._transcripts = transcripts

    def __iter__(self) -> Iterator:
        return iter(self.transcripts)

    def __len__(self) -> int:
        return len(self.transcripts)

    def __repr__(self) -> str:
        return '<Gene (composed of %d transcript(s)): %s>' % (len(self), str(self))

    def __getitem__(self, slice_):
        """
        slice each transcript and merge them together to get the valid sections of the
        gene record. If a slice would return no data (e.g. gene[:-10000], but gene is
        100 bases in length), returns None.
        """
        tx = []
        for t in self.transcripts:
            sliced = t[slice_.start:slice_.stop]
            if sliced:
                tx.append(sliced)

        # some slices return no data; in this case, return None
        if not tx:
            return None

        gene_fields = deepcopy(self.fields)
        gene_fields[3] = str(min(t.start for t in tx)).encode()
        gene_fields[4] = str(max(t.end for t in tx)).encode()
        return Gene(gene_fields, tx)

    @classmethod
    def from_records(cls, gene, *transcript_records):
        transcripts = []
        try:
            current = [transcript_records[0]]
        except IndexError:
            print('No transcript records detected: %s; Invalid gene.' %
                  repr(transcript_records))
            raise
        for record in transcript_records[1:]:
            if record[2] in [b'exon', b'UTR']:
                current.append(record)
            elif record[2] == b'transcript':
                transcripts.append(Transcript.from_records(*current))
                current = [record]
            else:
                pass  # CDS, stop_codon are also in GTF, but we don't use these records

        # append the last transcript
        if len(current) > 1:
            transcripts.append(Transcript.from_records(*current))

        return cls(gene, transcripts)

    @property
    def transcripts(self) -> list:
        return self._transcripts

    @property
    def transcript_sizes(self) -> list:
        return [tx.size for tx in self]

    @property
    def size(self):
        raise NotImplementedError  # gene size is the superset of all merged transcripts

    @staticmethod
    def _merge_intervals(ivs: (int, int)) -> Iterator:
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

    @lru_cache(4)  # don't recalculate same slice twice
    def intervals(self, start=None, stop=None) -> tuple:
        """
        return the set of genomic intervals defined by gene, given constraints posed by
        slice. if the slice would return no data, returns None.

        e.g. if slice= gene[-1000:], return the set of intervals formed by the union of
        the final 1000 bases of all transcripts of gene.

        args:
        -----
        slice_: the positions of each transcript to consider. Must be a contiguous section
         of the transcript

        returns:
        --------
        intervals: iterator of merged (start, end) tuples in chromosome coordinates
        :param stop:
        :param start:

        """
        if start or stop:
            gene = self[start:stop]
        else:
            gene = self

        if not gene:
            return tuple()

        ivs = []
        for tx in gene:
            ivs.extend((e.start, e.end) for e in tx if e.end - e.start > 0)

        return tuple(self._merge_intervals(ivs))


class Annotation:

    def __init__(self, gtf: str, fasta=None, max_insert_size=1000):
        """
        returns a complete genomic annotation from gtf and fasta files
        """

        # create gtf objects
        # self._len = 0  # tracker for length
        gtf_rd = Reader(gtf)

        # organize by gene
        self._genes = {}
        for gene_records in gtf_rd.iter_gene_sets():
            if len(gene_records) < 3:
                continue  # not a valid transcript; 3 records needed: gene, tx, exon.
            g = Gene.from_records(*gene_records)
            self._genes[g.integer_gene_id] = g

        # create fasta object
        if fasta:
            self._fasta = seqc.sequence.fasta.Genome.from_file(fasta)

        # create empty interval tree
        self.slice_genes(start=-max_insert_size, stop=None)
        self._interval_tree = None
        self.create_interval_tree()

    @property
    def genes(self) -> dict:
        return self._genes

    def __len__(self) -> int:
        """return the number of genes in the annotation"""
        return len(self.genes)

    def __repr__(self) -> str:
        return "<Annotation composed of %d genes>" % len(self)

    def __getitem__(self, id_):
        """
        both referencing with self[chromosome] and self[chromosome][gene] will return
        their respective objects due to heirarchical organization of dictionary
        """
        return self.genes[id_]

    def keys(self) -> Iterator:
        return self.genes.keys()

    def items(self) -> Iterator:
        return self.genes.items()

    def values(self) -> Iterator:
        return self.genes.values()

    def __iter__(self) -> Iterator:
        """iterate over genes"""
        for gene in self.keys():
            yield gene

    def random_genes(self, n):
        """generate n random genes
        :param n: number of genes to generate
        """
        indices = np.random.randint(0, len(self.genes), (n,))
        return map(lambda i: self.genes[i].integer_gene_id, indices)

    def random_sequences(self, n, sequence_length, start=None,
                         stop=None, filters=None) -> (list, list):
        """
        generate n random sequences; note that < n sequences may be generated. Check the
        length of the output sequence

        args:
        -----
        n: number of sequences to generate
        return_genes (True): return a second array with the gene ids
        start (None): if not None, restrict generation of sequences to begin at this
         position in the gene
        stop (None): if not None, restrict generation of sequences to end at this position
         in the gene

        returns:
        --------
        seqs: list of n nucleotide sequences of sequence_length
        genes (if return_genes is True): list of genes that seqs were drawn from

        """

        if not self._fasta:
            raise RuntimeError('rebuild GTF object with a passed fasta file to create '
                               'random sequences. self.fasta == None')

        if filters is None:
            all_genes = list(self.values())
        else:
            all_genes = []
            for gene in self.values():
                if not any(gene.attribute(k) == v for k, v in filters.items()):
                    all_genes.append(gene)
        n_genes = len(all_genes) - 1
        s = 0
        genes = []
        seqs = []
        while s < n:
            gene = all_genes[randint(0, n_genes)]
            ivs = gene.intervals(start, stop)
            if not ivs:
                continue
            i = randint(0, len(ivs) - 1)
            istart, istop = ivs[i]
            iend = istop - sequence_length
            if istart < iend:
                pos = randint(istart, iend - 1)
                sequence = self._fasta[gene.chromosome][pos:pos + sequence_length]
                if gene.strand == b'-':
                    sequence = seqc.sequence.revcomp(sequence)
                elif gene.strand == b'+':
                    pass
                else:
                    raise ValueError('strand must be + or -.')
                seqs.append(sequence)
                genes.append(gene)
                s += 1

        return seqs, genes

    def translate(self, strand: bytes, chromosome: bytes, position: int,
                  filters=None) -> list:
        # todo
        # some genes multimap for reasons noted in worklog on jan 13; eliminate for now
        try:
            ivs = self._interval_tree[(chromosome, strand)].search(position)
            if len(ivs) == 1:
                return first(ivs).data
            else:
                return None
        except KeyError:
            return None
        except TypeError:
            self.create_interval_tree(filters=filters)
            try:
                ivs = self._interval_tree[(chromosome, strand)].search(position)
                if len(ivs) == 1:
                    return first(ivs).data
                else:
                    return None
                # old version
                # return [iv.data for iv in
                #         self._interval_tree[(chromosome, strand)].search(position)]
            except KeyError:
                return None

    def create_interval_tree(self, filters: dict=None) -> None:
        """create a tree that maps genomic intervals to their corresponding gene ids

        args:
        -----
        filters: remove these types of genes (e.g. gene_type: miRNA)
        """

        # map integer ENSEMBL ids to intervals
        data = {}

        if filters is None:
            candidate_genes = self.values()
        else:
            candidate_genes = []
            for gene in self.values():
                if not any(gene.attribute(k) == v for k, v in filters.items()):
                    candidate_genes.append(gene)

        for gene in candidate_genes:
            for iv in gene.intervals():
                try:
                    id_ = gene.integer_gene_id
                    data[(gene.seqname, gene.strand)].addi(iv[0], iv[1], id_)
                except KeyError:
                    data[(gene.seqname, gene.strand)] = IntervalTree()
                    id_ = gene.integer_gene_id
                    data[(gene.seqname, gene.strand)].addi(iv[0], iv[1], id_)

        self._interval_tree = data

    def slice_genes(self, start=None, stop=None):
        """reduce the size of the annotation, slicing the genes; makes changes in-place"""
        for id_, gene in self.items():
            self._genes[id_] = gene[start:stop]


class Reader(seqc.reader.Reader):

    def __iter__(self) -> Iterator:
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

    def iter_exons(self) -> Iterator:
        """iterate over all exons (and UTRs) in gtf

        yields:
        -------
        iterator of Record objects

        """
        for record in self:
            if record[2] in [b'exon', b'UTR']:
                yield Record(record)

    def iter_genes(self) -> Iterator:
        """iterate over all genes in gtf

        yields:
        -------
        iterator of Record objects

        """
        for record in self:
            if record[2] == b'gene':
                yield Record(record)

    def iter_transcripts(self) -> Iterator:
        """iterate over all transcripts in gtf

        yields:
        -------
        iterator of Record objects

        """
        for record in self:
            if record[2] == b'transcript':
                yield Record(record)

    def iter_gene_sets(self) -> Iterator:
        """iterate over all the records for each gene in passed gtf

        yields:
        -------
        iterator of lists of Records objects
        """

        records = iter(self)

        # get the first gene record
        record = next(records)
        while record[2] != b'gene':
            record = next(records)

        # aggregate all records until the next gene record pops up
        gene = [record]
        record = next(records)
        while record:
            if record[2] != b'gene':
                gene.append(record)
            else:
                yield gene
                gene = [record]
            record = next(records)

        # yield the final gene record
        yield gene
