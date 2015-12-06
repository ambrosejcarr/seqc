__author__ = 'ambrose'

import os
import seqc
import unittest
import nose2
from nose2.tools import params
from operator import itemgetter
from seqc.test.test_unit import config


class GTFReaderTakeFinalNTest(unittest.TestCase):

    test_dir = 'test_seqc/'
    gtf = 'test_seqc/test_data.gtf'

    @classmethod
    def setUpClass(cls):
        """mock-up a test_data file from the first 100 lines of chr19_mm38.gtf

        list genes: awk '{if($3=="gene"){print $0}}' test_seqc/test_data.gtf

        result: ENSMUSG00000100969.1, ENSMUSG00000093983.1, and ENSMUSG00000024831.8
        """

        if not os.path.isdir(cls.test_dir):
            os.mkdir(cls.test_dir)

        i = 0
        with open(config.gtf, 'rb') as fin:
            with open(cls.gtf, 'wb') as fout:
                while i < 100:
                    fout.write(fin.readline())
                    i += 1

    @params(250, 500, 1000)
    def test_take_final_n_and_iter_genes_final_nbases(self, n):

        # native code calls take_final_n given a set of exons, and a number of bases.
        # the exons are generated as follows, with self replaced by a .gtf reader
        # object. The following code is taken from gtf.Reader.iter_genes_final_nbases():
        # ------------------------------------------------------------------------------

        rd = seqc.gtf.Reader(self.gtf)
        seqc.util.check_type(n, int, 'n')

        exon_set_iterator = rd.iter_exon_sets()

        # create the first transcript
        first_exon_set = next(exon_set_iterator)

        exon_sets = list(rd.iter_exon_sets())

        # run some checks on the exon inputs
        # ------------------------------------------------------------------------------

        # make sure correct number of exon sets are present
        self.assertEqual(9, len(exon_sets))

        for exons in exon_sets:
            # check that all exons have the same strand
            strands = set(e.strand for e in exons)
            self.assertEqual(len(strands), 1)

            # check that all exons have the same transcript id
            tx_ids = set(e.attribute['transcript_id'] for e in exons)
            self.assertEqual(len(tx_ids), 1)

            # check that all exons have a larger end than start
            intervals = sum(1 if int(e.start) < int(e.end) else 0 for e in exons)
            self.assertEqual(intervals, len(exons))

            # strand is derived from exons[0].strand

        # set inputs
        exons = exon_sets[0]
        strand = exons[0].strand  # all strands are identical within set, checked above.

        # test_data _take_final_nbases() given this input
        # ------------------------------------------------------------------------------

        if n <= 0:
            raise ValueError('n must be a positive integer')

        # convert exons into intervals and eliminate duplicates
        intervals = set((int(r.start), int(r.end)) for r in exons)
        self.assertEqual(len(intervals), len(exons), 'Transcript input is sorted and '
                                                     'should contain no duplicate exons')

        output_intervals = []
        template = exons[0]

        # process exons
        intervals = sorted(intervals, key=itemgetter(1, 0), reverse=True)

        for start, end in intervals:
            self.assertTrue(end > start)

            # add complete exon
            if end - start < n:
                output_intervals.append((start, end))
                n -= end - start  # decrement n

            # exon is larger than remaining bases; add part of exon
            else:
                if strand is '+':
                    output_intervals.append((end - n, end))
                    self.assertTrue(end > end - n)
                elif strand is '-':
                    self.assertTrue(start + n > start)
                    output_intervals.append((start, start + n))
                else:
                    raise ValueError('strand must be "+" or "-", not %s' % repr(strand))

                res = seqc.gtf.MultiRecord(
                    template.seqname, template.source, template.feature, output_intervals,
                    template.score, template.strand, template.frame, template.attribute)

                self.assertTrue(all(iv[1] > iv[0] for iv in output_intervals))

        # exons exhausted, transcript was smaller than n; return entire list
        res = seqc.gtf.MultiRecord(
            template.seqname, template.source, template.feature, output_intervals,
            template.score, template.strand, template.frame, template.attribute)

        self.assertTrue(all(iv[1] > iv[0] for iv in output_intervals))


if __name__ == "__main__":
    nose2.main()