import seqc
import gzip
import bz2
from copy import deepcopy
from collections import Counter, defaultdict
from sys import maxsize
from operator import attrgetter, itemgetter
from functools import lru_cache, partial
from collections.abc import Iterator
from random import randint
import numpy as np
import string

_del_letters = string.ascii_letters.encode()
_del_non_letters = ''.join(set(string.printable).difference(string.ascii_letters)).encode()


class Record:

    __slots__ = ['_fields', '_attribute']

    def __init__(self, fields):
        """

        caution: these objects must be created quickly in large numbers, so type checking
         has been omitted.

        args:
        -----
        fields: list of str gtf record fields
        """

        self._fields = fields
        self._attribute = {}

    def __repr__(self):
        return '<Record: %s>' % str(self)

    def __str__(self):
        return b'\t'.join(self._fields).decode()

    def __bytes__(self):
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

        return Exon(fields)

    def _parse_attribute(self):
        for field in self._fields[8].rstrip(b';\n').split(b';'):
            key, value = field.strip().split()
            self._attribute[key] = value.strip(b'"')

    @property
    def seqname(self):
        return self._fields[0]

    @property
    def chromosome(self):
        return self._fields[0]  # synonym for seqname

    @property
    def source(self):
        return self._fields[1]

    @property
    def feature(self):
        return self._fields[2]

    @property
    def start(self):
        return int(self._fields[3])

    @property
    def end(self):
        return int(self._fields[4])

    @property
    def score(self):
        return self._fields[5]

    @property
    def strand(self):
        return self._fields[6]

    @property
    def frame(self):
        return self._fields[7]

    @property
    def size(self):
        return self.end - self.start

    @property
    def fields(self):
        return self._fields

    def attribute(self, item):
        try:
            return self._attribute[item]
        except KeyError:
            if not self._attribute:
                self._parse_attribute()  # todo implement
                return self._attribute[item]
            else:
                raise KeyError('%s is not a stored attribute of this gtf record' %
                               repr(item))

    @property
    def gene_name(self):
        return self.attribute(b'gene_name')

    @property
    def integer_gene_id(self):
        return int(self.attribute(b'gene_id').split(b'.')[0]
                   .translate(None, _del_letters))

    @property
    def organism_prefix(self):
        return self.attribute(b'gene_id').translate(None, _del_non_letters)

    @property
    def string_gene_id(self):
        return self.attribute(b'gene_id')

    @staticmethod
    def int2str_gene_id(integer_id, organism_prefix):
        """converts an integer gene id to a string gene id"""
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

    def __repr__(self):
        return '<Exon: %s>' % str(self)


class UTR(Record):

    def __repr__(self):
        return '<UTR: %s>' % str(self)


class Transcript(Record):

    __slots__ = ['_exons']

    def __init__(self, transcript, exons):

        # define transcript
        super().__init__(transcript)

        # sorting necessary to get - strand transcripts' UTRs in correct order
        if transcript[6] == b'+':
            self._exons = exons
        else:
            self._exons = sorted(exons, key=attrgetter('start'), reverse=True)

    def __iter__(self):
        return iter(self.exons)

    def __len__(self):
        return len(self.exons)

    def __repr__(self):
        return '<Transcript (composed of %d exon(s) and UTR(s)): %s>' % \
               (len(self), str(self))

    def __getitem__(self, slice_object):
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
                    exon = next(exon_iterator)
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
    def exons(self):
        return self._exons

    @property
    def size(self):
        return sum(exon.end - exon.start for exon in self)


class Gene(Record):

    __slots__ = ['_transcripts']

    def __init__(self, gene, transcripts):

        # define gene
        super().__init__(gene)
        self._transcripts = transcripts

    def __iter__(self):
        return iter(self.transcripts)

    def __len__(self):
        return len(self.transcripts)

    def __repr__(self):
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
        return Gene(self.fields, tx)

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
    def transcripts(self):
        return self._transcripts

    @property
    def transcript_sizes(self):
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
    def intervals(self, start=None, stop=None):
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

        """
        if start or stop:
            gene = self[start:stop]
        else:
            gene = self

        if not gene:
            return tuple()

        ivs = []
        for tx in gene:
            ivs.extend((e.start, e.end) for e in tx)

        return tuple(self._merge_intervals(ivs))


class Annotation:

    def __init__(self, gtf, fasta):
        """
        returns a complete genomic annotation from gtf and fasta files
        """

        # create gtf objects
        self._len = 0  # tracker for length
        gtf_rd = Reader(gtf)

        # organize by chromosome and gene
        self._chromosomes = defaultdict(dict)
        self._genes = []
        for gene_records in gtf_rd.iter_gene_sets():
            if len(gene_records) < 3:
                continue  # not a valid transcript; 3 records needed: gene, tx, exon.
            g = Gene.from_records(*gene_records)
            self._genes.append(g)
            self._chromosomes[g.chromosome][g.attribute(b'gene_id')] = g
            self._len += 1

        # create fasta object
        self._fasta = seqc.fasta.Fasta.from_file(fasta)

    def __len__(self):
        """return the number of genes in the annotation"""
        return self._len

    def __repr__(self):
        return "<Annotation: %d genes over %d chromosomes: %s>" % (
            len(self), len(self.chromosomes), repr(self.chromosomes.keys()))

    def __getitem__(self, item):
        """
        both referencing with self[chromosome] and self[chromosome][gene] will return
        their respective objects due to heirarchical organization of dictionary
        """
        return self.chromosomes[item]

    def __iter__(self):
        """iterate over genes"""
        for chromosome in self.chromosomes:
            for gene in self[chromosome]:
                yield self[chromosome][gene]

    @property
    def chromosomes(self):
        return self._chromosomes

    @property
    def genes(self):
        return self._genes

    @property
    def fasta(self):
        return self._fasta

    def random_genes(self, n):
        """generate n random genes"""
        indices = np.random.randint(0, len(self.genes), (n,))
        return map(lambda i: self.genes[i].string_gene_id, indices)

    def random_sequences(self, n, sequence_length, return_genes=True, start=None,
                         stop=None):
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

        n_genes = len(self.genes) - 1
        s = 0
        genes = []
        seqs = []
        while s < n:
            gene = self.genes[randint(0, n_genes)]
            ivs = gene.intervals(start, stop)
            if not ivs:
                continue
            i = randint(0, len(ivs) - 1)
            istart, istop = ivs[i]
            iend = istop - sequence_length
            if istart < iend:
                pos = randint(istart, iend - 1)
                seqs.append(self.fasta[gene.chromosome][pos:pos + sequence_length])
                genes.append(gene)
                s += 1

        if return_genes:
            return seqs, genes
        else:
            return seqs


class Reader:

    def __init__(self, gtf):

        seqc.util.check_type(gtf, str, 'gtf')
        seqc.util.check_file(gtf, 'gtf')

        self.gtf = gtf
        try:
            gtf_iterator = iter(self)
            next(gtf_iterator)
        except:
            raise ValueError('gtf is an invalid GTF file. Please check file formatting.')

    def _open(self):
        """
        seamlessly open self._gtf, whether gzipped or uncompressed

        returns:
        --------
        fobj: open file object
        """
        if self.gtf.endswith('.gz'):
            fobj = gzip.open(self.gtf, 'rb')
        elif self.gtf.endswith('.bz2'):
            fobj = bz2.open(self.gtf, 'rb')
        else:
            fobj = open(self.gtf, 'rb')
        return fobj

    def __iter__(self):
        """return an iterator over all non-header records in gtf"""
        fobj = self._open()
        records = iter(fobj)

        try:
            # get rid of headers
            record = next(records)
            while record.startswith(b'#'):
                record = next(records)

            # yield all records
            while record:
                yield record.split(b'\t')
                record = next(records)
        finally:
            fobj.close()

    def iter_exons(self):
        """iterate over all exons (and UTRs) in gtf

        yields:
        -------
        iterator of Record objects

        """
        for record in self:
            if record.feature == 'exon' or record.feature == 'UTR':
                yield record

    def iter_genes(self):
        """iterate over all genes in gtf

        yields:
        -------
        iterator of Record objects

        """
        for record in self:
            if record.feature == 'gene':
                yield record

    def iter_transcripts(self):
        """iterate over all transcripts in gtf

        yields:
        -------
        iterator of Record objects

        """
        for record in self:
            if record.feature == 'transcript':
                yield record

    def iter_gene_sets(self):
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

