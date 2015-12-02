__author__ = 'ambrose'

import gzip
from collections import namedtuple, defaultdict
from collections.abc import Iterator, Iterable
from more_itertools import peekable
from itertools import chain
from operator import itemgetter, attrgetter
import numpy as np
import io
import seqc
import pickle

# create a simple mock-class for returned gtf records
Record = namedtuple('Record', ['seqname', 'source', 'feature', 'start', 'end',
                               'score', 'strand', 'frame', 'attribute'])


# same as Record except 'start' and 'end' are replaced by 'intervals', a list of
#   [(start, end)] coordinates
MultiRecord = namedtuple('MultiRecord', ['seqname', 'source', 'feature', 'intervals',
                                         'score', 'strand', 'frame', 'attribute'])


class Sample:
    """Sample from intervals specified by a gtf file"""

    def __init__(self, gtf):
        seqc.util.check_type(gtf, str, 'gtf')
        self._gtf = gtf
        self._gene_data = None
        self._tx_data = None
        self._exon_data = None

    def sample_genes(self, fragment_length: int=1000) -> np.array:

class Reader:
    """GTF reader optimized for utility and simplicity"""

    def __init__(self, gtf):

        seqc.util.check_type(gtf, str, 'gtf')
        seqc.util.check_file(gtf, 'gtf')

        self._gtf = gtf
        try:
            gtf_iterator = iter(self)
            next(gtf_iterator)
        except:
            raise ValueError('gtf is an invalid GTF file. Please check file formatting.')

    @property
    def gtf(self):
        return self._gtf

    def _open(self) -> io.TextIOBase:
        """
        seamlessly open self._gtf, whether gzipped or uncompressed

        returns:
        --------
        fobj: open file object
        """
        if self._gtf.endswith('.gz'):
            fobj = gzip.open(self._gtf, 'rt')
        else:
            fobj = open(self._gtf)
        return fobj

    def __iter__(self):
        """return an iterator over all non-header records in gtf"""
        fobj = self._open()
        try:
            for line in fobj:
                if line.startswith('#'):
                    continue
                line = line.strip().split('\t')
                meta = line[-1]
                dmeta = {}
                for field in meta.rstrip(';').split(';'):
                    field = field.strip().split()
                    dmeta[field[0]] = field[1].strip('"')
                line[-1] = dmeta
                yield Record(*line)
        finally:
            fobj.close()

    def iter_exons(self) -> Iterator:
        """iterate over all exons (and UTRs) in gtf

        yields:
        -------
        iterator of Record objects

        """
        for record in self:
            if record.feature == 'exon' or record.feature == 'UTR':
                yield record

    def iter_genes(self) -> Iterator:
        """iterate over all genes in gtf

        yields:
        -------
        iterator of Record objects

        """
        for record in self:
            if record.feature == 'gene':
                yield record

    def iter_transcripts(self) -> Iterator:
        """iterate over all transcripts in gtf

        yields:
        -------
        iterator of Record objects

        """
        for record in self:
            if record.feature == 'transcript':
                yield record

    def iter_exon_sets(self) -> Iterator:
        """iterate over the sets of exons for each transcript in self.gtf

        yields:
        -------
        iterator of lists of Records objects
        """
        tx = 'transcript_id'
        exon_iterator = peekable(self.iter_exons())

        # get the first exon
        records = [next(exon_iterator)]

        for exon in exon_iterator:
            if exon.attribute[tx] == records[0].attribute[tx]:
                records.append(exon)
            else:
                yield records
                records = [exon]
        yield records

    @staticmethod
    def _take_final_n(exons: list, n: int, strand: str) -> MultiRecord:
        """
        return a MultiRecord containing all intervals within n bases of the TTS for the
        exons of a given transcript

        note: if the input are the exons of a single transcript from a properly formatted
        gtf file, then the sorting is superfluous, and the set operation should do
        nothing, however both seem necessary for ENSEMBL gtf files.

        args:
        -----
        n: the number of bases from the TTS to consider when returning intervals
        strand: whether the intervals are found on the + or - strand. Options: '+', '-'

        returns:
        MultiRecord object whose start and end fields are lists of intervals that
          fall within n bases of the TTS.

        """

        if n <= 0:
            raise ValueError('n must be a positive integer')

        # convert exons into intervals and eliminate duplicates
        intervals = set((int(r.start), int(r.end)) for r in exons)

        output_intervals = []
        template = exons[0]

        # process exons
        intervals = sorted(intervals, key=itemgetter(1, 0), reverse=True)
        for start, end in intervals:

            # add complete exon
            if end - start < n:
                output_intervals.append((start, end))
                n -= end - start  # decrement n

            # exon is larger than remaining bases; add part of exon
            else:
                if strand is '+':
                    output_intervals.append((end - n, end))
                elif strand is '-':
                    output_intervals.append((start, start + n))
                else:
                    raise ValueError('strand must be "+" or "-", not %s' % repr(strand))

                return MultiRecord(
                    template.seqname, template.source, template.feature, output_intervals,
                    template.score, template.strand, template.frame, template.attribute)

        # exons exhausted, transcript was smaller than n; return entire list
        return MultiRecord(
            template.seqname, template.source, template.feature, output_intervals,
            template.score, template.strand, template.frame, template.attribute)

    @staticmethod
    def _merge_transcripts(transcripts: list) -> MultiRecord:
        """
        returns a MultiRecord whose start and stop fields are tuples containing all
        alignable intervals in a gene. If any intervals overlap, they are merged. If run
        on the results of self.take_final_n(), then the result is all alignable sequence
        for a given gene within n bases of all TTSs for that gene.

        args:
        -----
        transcripts: iterable of MultiRecord objects. this function is normally
          called on the output of self.take_final_n() or self.merge_exons().


        returns:
        --------
        record: merged Record object

        """

        def merge_intervals(ivs: list) -> Iterator:
            saved = list(ivs[0])
            for st, en in ivs:
                if st <= saved[1]:
                    saved[1] = max(saved[1], en)
                else:
                    yield tuple(saved)
                    saved[0] = st
                    saved[1] = en
            yield tuple(saved)

        # get a template for the output MultiRecord; we're going to change the intervals
        template = vars(transcripts[0])

        # sort the complete list of intervals from all transcripts
        sorted_intervals = sorted(chain(*(tx.intervals for tx in transcripts)))

        # merge any overlapping intervals
        template['intervals'] = list(merge_intervals(sorted_intervals))

        return MultiRecord(**template)

    def iter_genes_final_nbases(self, n: int=1000) -> Iterator:
        """
        Iterator of MultiRecords containing all intervals within n bases of a gene's TTSs.

        args:
        -----
        n: the number of bases from the TTS to consider when returning intervals

        returns:
        --------
        iterator: iterator of records containing the final n bases of each gene

        """

        seqc.util.check_type(n, int, 'n')

        exon_set_iterator = self.iter_exon_sets()

        # create the first transcript
        first_exon_set = next(exon_set_iterator)
        transcripts = [self._take_final_n(first_exon_set, n,
                                          strand=first_exon_set[0].strand)]

        for exons in self.iter_exon_sets():
            tx = self._take_final_n(exons, n, exons[0].strand)
            if tx.attribute['gene_id'] == transcripts[0].attribute['gene_id']:
                transcripts.append(tx)
            else:  # new gene; merge and yield the last one
                gene = self._merge_transcripts(transcripts)
                yield gene
                transcripts = [tx]

        # yield the final gene
        gene = self._merge_transcripts(transcripts)
        yield gene

    def scid_to_gene(self, save=''):
        # seems like there could be a bug with scid=0: it has a huge number of genes.
        gmap = defaultdict(set)
        for record in self.iter_transcripts():
            gmap[int(record.attribute['scseq_id'].strip('SC'))].add(
                record.attribute['gene_name'])

        # merge features that are not unique into single strings
        for k, v in gmap.items():
            gmap[k] = '-'.join(v)

        if save:
            with open(save, 'wb') as f:
                pickle.dump(gmap, f)

        return gmap

    def get_phix_id(self):
        for record in self:
            if record.seqname is 'NC_001422.1':
                if record.feature == 'transcript':
                    return int(record.attribute['scseq_id'].strip('SC'))

    def get_mitochondrial_ids(self):
        mt_ids = set()
        for record in self:
            if record.feature == 'transcript':
                if record.attribute['transcript_id'].startswith('MT-'):
                    mt_ids.add(int(record.attribute['scseq_id'].strip('SC')))
        return mt_ids
