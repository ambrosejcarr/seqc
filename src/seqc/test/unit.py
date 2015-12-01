__author__ = 'ambrose'


from numpy.lib.recfunctions import append_fields
import nose2
from nose2.tools import params
from more_itertools import first
import unittest
import os
import pickle
import tables as tb
import seqc
from io import StringIO
import numpy as np
import xml.dom.minidom
import random
import shutil
# from memory_profiler import profile, memory_usage
# from itertools import product
# from xml.etree import ElementTree as ET
# from subprocess import Popen, PIPE
# import cProfile
# from pstats import Stats
# import hashlib
# import logging
# import socket
# import re
# import pandas as pd
# import gzip
# import ftplib


################################ STATE CONFIGURATION ####################################

# noinspection PyPep8Naming
class config:
    """
    Dummy class to hold state configuration variables that is importable into other
    scopes for testing
    """

    # this is the universal data dir for these tests
    seqc_dir = '/'.join(seqc.__file__.split('/')[:-3]) + '/'

    # set of dtypes
    data_types = ['in_drop', 'drop_seq']

    # universal file patterns
    samfile_pattern = seqc_dir + 'test_data/%s/.seqc_test/Aligned.out.sam'
    forward_pattern = seqc_dir + 'test_data/%s/fastq/test_seqc_r1.fastq'
    reverse_pattern = seqc_dir + 'test_data/%s/fastq/test_seqc_r2.fastq'
    barcode_pattern = seqc_dir + 'test_data/%s/barcodes/barcodes.p'
    barcode_link_pattern = 's3://dplab-data/barcodes/%s/barcodes.p'
    h5_name_pattern = seqc_dir + 'test_data/%s/test_seqc.h5'

    # universal index files
    gtf = seqc_dir + 'test_data/genome/annotations.gtf'
    fasta = seqc_dir + 'test_data/genome/mm38_chr19.fa'
    index = seqc_dir + 'test_data/genome/'
    index_link = 's3://dplab-data/genome/mm38_chr19/'

    # config parameters
    n_threads = 7


############################### SETUP FILE GENERATION ###################################

def check_fastq(data_type: str):
    """check if the required fastq files are present, else generate them"""

    # replace any dashes with underscores
    data_type = data_type.replace('-', '_')

    # get forward and reverse file ids
    forward = config.forward_pattern % data_type
    reverse = config.reverse_pattern % data_type

    if not all(os.path.isfile(f) for f in [forward, reverse]):
        gen_func = getattr(seqc.fastq.GenerateFastq, data_type)
        gen_func(10000, config.seqc_dir + 'test_data/%s/fastq/test_seqc', config.fasta,
                 config.gtf)

    return forward, reverse


def check_sam(data_type: str):
    """check if the required sam files are present, else generate them"""

    # replace any dashes with underscores
    data_type = data_type.replace('-', '_')

    # get sam pattern
    samfile = config.samfile_pattern % data_type

    # generation params
    n = 10000
    prefix = config.seqc_dir + 'test_data/%s/fastq/test_seqc'

    if not os.path.isfile(samfile):
        gen_func = getattr(seqc.sam.GenerateSam, data_type)
        gen_func(n=n, prefix=prefix, fasta=config.fasta, gtf=config.gtf, index=config.index)

    return samfile


def check_index():
    """ensure that there is an index present. If not, download it."""
    gfiles = ['Genome', 'SA', 'SAindex', 'annotations.gtf']
    if not os.path.isdir(config.index) or not all(os.path.isfile(config.index + f) for f
                                                  in gfiles):
        index_bucket = 'dplab-data'
        index_prefix = 'genomes/mm38_chr19/'
        seqc.io.S3.download_files(
            bucket=index_bucket, key_prefix=index_prefix, output_prefix=config.index,
            cut_dirs=2)


##################################### UNIT TESTS ########################################


### FASTQ TESTS ###

class FastqRevcompTest(unittest.TestCase):

    def test_revcomp_palindrome(self):
        s = 'ACGT'
        rc = seqc.fastq.revcomp(s)
        self.assertEqual(s, rc)

    def test_revcomp_asymetric(self):
        s = 'ACCGGTT'
        rc = seqc.fastq.revcomp(s)
        self.assertEqual(rc, 'AACCGGT')

    def test_revcomp_raises_value_error(self):
        s = 'AAGAGAV'  # V should raise
        self.assertRaises(ValueError, seqc.fastq.revcomp, s)

    def test_revcomp_wrong_input_type(self):
        s = list('AAGAGA')
        self.assertRaises(TypeError, seqc.fastq.revcomp, s)


class FastqTruncateSequenceLengthTest(unittest.TestCase):

    def setUp(self):
        self.fastq_file = seqc.fastq.GenerateFastq.simple_fastq(10, 100)
        self.n = 50
        self.fname = 'fastq_truncateSequenceLengthTest.fastq'

    def test_truncate_sequence_length_wrong_input_raises_type_error(self):
        self.assertRaises(TypeError, seqc.fastq.truncate_sequence_length, self.fastq_file,
                          '10', self.fname)
        self.assertRaises(TypeError, seqc.fastq.truncate_sequence_length, self.fastq_file,
                          self.n, (self.fname,))
        self.assertRaises(TypeError, seqc.fastq.truncate_sequence_length,
                          [self.fastq_file], self.n, self.fname)

    def test_truncate_sequence_length_produces_correct_length(self):
        seqc.fastq.truncate_sequence_length(self.fastq_file, self.n, self.fname)
        with open(self.fname, 'r') as f:
            for r in seqc.fastq.iter_records(f):
                self.assertEqual(len(r[1]), self.n + 1)  # + 1 for '\n'

    def test_truncate_sequence_length_produces_well_formed_fastq(self):
        seqc.fastq.truncate_sequence_length(self.fastq_file, self.n, self.fname)
        with open(self.fname, 'r') as f:
            for r in seqc.fastq.iter_records(f):
                self.assertEqual(r[0][0], '@')
                self.assertEqual(r[2][0], '+')

        # make sure final quality string ends with '\n'
        self.assertEqual(r[3][-1], '\n')

    def tearDown(self):
        if os.path.isfile(self.fname):
            os.remove(self.fname)


class FastqEstimateSequenceLengthTest(unittest.TestCase):

    def setUp(self):
        self.fname = 'fastq_estimateSequenceLengthTest.fastq'

    def write_simple_fastq(self, n, length):
        fastq_data = seqc.fastq.GenerateFastq.simple_fastq(n, length).read()
        with open(self.fname, 'w') as f:
            f.write(fastq_data)

    def test_estimate_sequence_length_wrong_input_raises_type_error(self):
        self.write_simple_fastq(10, 100)
        self.assertRaises(TypeError, seqc.fastq.estimate_sequence_length,
                          open(self.fname, 'rb'))

    def test_estimate_sequence_length_small_input_file(self):
        """should produce the exact correct result"""
        self.write_simple_fastq(10, 95)
        mean, std, (counts, freqs) = seqc.fastq.estimate_sequence_length(self.fname)
        self.assertEqual(mean, 95)
        self.assertEqual(std, 0)
        self.assertEqual(np.array([95]), counts)
        self.assertEqual(np.array([10]), freqs)

    def tearDown(self):
        if os.path.isfile(self.fname):
            os.remove(self.fname)


### IO TESTS ###

class IoS3Test(unittest.TestCase):

    def setUp(self):
        self.bucket = 'dplab-data'
        self.testfile_key = 'genomes/mm38/chrStart.txt'
        self.testfile_download_name = 'test_s3_download.txt'

    def test_wrong_input_type_raises_type_error(self):
        self.assertRaises(TypeError, seqc.io.S3.download_file, self.bucket, 10,
                          self.testfile_download_name)
        self.assertRaises(TypeError, seqc.io.S3.download_file, self.bucket,
                          self.testfile_key, StringIO('1010'))

    def test_download_incorrect_filepath_raises_file_not_found_error(self):
        self.assertRaises(FileNotFoundError, seqc.io.S3.download_file, self.bucket,
                          'foobar', self.testfile_download_name, overwrite=False)

    def test_download_files_gets_expected_data(self):
        pass  # todo implement

    def test_download_files_does_not_overwrite(self):
        pass  # todo implement; can use os.stat to check modification time.


### GTF TESTS ###

class GTFReaderTest(unittest.TestCase):

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
                while i < 203:
                    fout.write(fin.readline())
                    i += 1

    def test_reader_input_is_not_a_file_raises_file_not_found_error(self):
        self.assertRaises(FileNotFoundError, seqc.gtf.Reader, 'not_a_real_file.gtf')

    def test_reader_input_is_not_a_string_raises_type_error(self):
        self.assertRaises(TypeError, seqc.gtf.Reader, ('not_a_real_file.gtf',))

    def test_reader_input_is_not_a_gtf_file_raises_value_error(self):
        self.assertRaises(ValueError, seqc.gtf.Reader, config.fasta)

    def test_reader_input_is_valid_gtf_file_does_not_raise(self):
        seqc.gtf.Reader(self.gtf)

    def test_reader_iterate_does_not_return_headers(self):

        # get original gtf data
        with open(self.gtf, 'rb') as f:
            data = f.read()

        # prepend some headers to it.
        with open(self.gtf, 'wb') as f:
            f.write(b'# this is a header line\n')
            f.write(b'# and so is this\n')
            f.write(data)

        for record in seqc.gtf.Reader(self.gtf):
            pass  # simply iterating is test_data enough; any header would cause iter to fail

    def test_reader_iter_exons_returns_all_exons(self):
        # 120 lines are exons or UTRs;
        # awk '{if($3=="exon"||$3=="UTR"){print $0}}' test_seqc/test_data.gtf | wc -l
        rd = seqc.gtf.Reader(self.gtf)
        n_exons = sum(1 for exon in rd.iter_exons())
        self.assertEqual(n_exons, 120)

    def test_reader_iter_transcripts_returns_all_transcripts(self):
        # 17 lines are transcripts;
        # awk '{if($3=="transcript"){print $0}}' test_seqc/test_data.gtf | wc -l
        rd = seqc.gtf.Reader(self.gtf)
        n_transcripts = sum(1 for tx in rd.iter_transcripts())
        self.assertEqual(n_transcripts, 17)

    def test_reader_iter_transcripts_returns_all_genes(self):
        # 4 lines are genes;
        # awk '{if($3=="gene"){print $0}}' test_seqc/test_data.gtf | wc -l
        rd = seqc.gtf.Reader(self.gtf)
        n_genes = sum(1 for gene in rd.iter_genes())
        self.assertEqual(n_genes, 4)

    def test_reader_iter_exon_sets_returns_correct_number_sets(self):
        # if there are 17 transcripts, there should be 17 exon sets
        rd = seqc.gtf.Reader(self.gtf)
        n_sets = sum(1 for exon_set in rd.iter_exon_sets())
        self.assertEqual(n_sets, 17)

        for exon_set in rd.iter_exon_sets():
            # test_data that the resulting iterable is a list of Record objects
            self.assertIsInstance(exon_set, list)

            # test_data that all list items are Record objects
            self.assertTrue(all(isinstance(exon, seqc.gtf.Record) for exon in exon_set),
                            'Exon set should return lists of Record objects')

            # test_data that all record objects have the same transcript id
            unique_transcripts = set(e.attribute['transcript_id'] for e in exon_set)
            self.assertEqual(1, len(unique_transcripts))

    def test_reader__take_final_n_negative_or_zero_n_raises(self):
        rd = seqc.gtf.Reader(self.gtf)
        exons = first(rd.iter_exon_sets())
        self.assertRaises(ValueError, rd._take_final_n, exons, -1, exons[0].strand)
        self.assertRaises(ValueError, rd._take_final_n, exons, 0, exons[0].strand)

    def test_gtf_file_all_interval_ends_larger_than_starts(self):
        rd = seqc.gtf.Reader(self.gtf)
        for i, record in enumerate(rd):
            self.assertTrue(int(record.end) > int(record.start),
                            'potentially malformed record #%d: %s' % (i, repr(record)))

    @unittest.skip('This does not appear to be true for ENSEMBL GTF files.')
    def test_gtf_file_no_duplicate_exons_within_transcripts(self):
        rd = seqc.gtf.Reader(self.gtf)
        for exons in rd.iter_exon_sets():
            all_intervals = [(int(r.start), int(r.end)) for r in exons]
            unique_intervals = set(all_intervals)
            self.assertEqual(len(all_intervals), len(unique_intervals))

    def test_reader__take_final_n_returns_correct_number_bases(self):
        rd = seqc.gtf.Reader(self.gtf)

        def sum_unique_intervals(multi_record: seqc.gtf.MultiRecord):
            n = 0
            for iv in multi_record.intervals:
                # check that the interval is positive
                assert iv[1] > iv[0], '%d <= %d' % (iv[1], iv[0])
                n += iv[1] - iv[0]
            return n

        # go through the exon sets for each gene, calculate the total sequence size for
        # each transcript
        tx_sizes = []
        for exons in rd.iter_exon_sets():
            seen = set()
            size = 0
            for exon in exons:
                iv = (int(exon.start), int(exon.end))
                if not iv in seen:
                    seen.add(iv)
                    size += iv[1] - iv[0]
            tx_sizes.append(size)

        # check that intervals return the correct interval sizes
        for size, exons in zip(tx_sizes, rd.iter_exon_sets()):
            # very small
            quarter = rd._take_final_n(exons, n=5, strand=exons[0].strand)
            self.assertEqual(sum_unique_intervals(quarter), min(size, 5))

            # normal
            half = rd._take_final_n(exons, n=500, strand=exons[0].strand)
            self.assertEqual(sum_unique_intervals(half), min(size, 500))

            # very large
            full = rd._take_final_n(exons, n=100000, strand=exons[0].strand)
            self.assertEqual(sum_unique_intervals(full), min(size, 100000))

        # check that the final positions are the ones being returned, not the start
        # for both small and large interval allowances
        for exons in rd.iter_exon_sets():
            last_interval = (int(exons[-1].start), int(exons[-1].end))

            # 1e10 > whole transcriptome (large)
            tx = rd._take_final_n(exons, int(1e10), exons[0].strand)
            self.assertIn(last_interval, tx.intervals)

            # 500 will be smaller than some exons (small)
            tx = rd._take_final_n(exons, 500, exons[0].strand)
            start, end = int(exons[-1].start), int(exons[-1].end)
            if exons[0].strand == '+':
                last_interval = (start if end - start < 500 else end - 500, end)
                self.assertIn(last_interval, tx.intervals)
            elif exons[0].strand == '-':
                last_interval = (start, end if end - start < 500 else start + 500, repr(exons))
            else:
                raise ValueError

    @params(10, 1000, int(1e9))
    def test_reader_merge_transcripts_final_n(self, n):

        # basic run boilerplate
        rd = seqc.gtf.Reader(self.gtf)
        txs = [rd._take_final_n(exons, n, exons[0].strand)
               for exons in rd.iter_exon_sets()]
        gene = rd._merge_transcripts(txs)

        # test_data that output is sorted
        self.assertEqual(gene.intervals, sorted(gene.intervals))

        # test_data that merged size is greater than or equal to the largest individual
        # transcript size.
        interval_size = lambda t: sum(e - s for s, e in t.intervals)
        max_tx_size = max(interval_size(tx) for tx in txs)
        self.assertLessEqual(max_tx_size, interval_size(gene))

    def test_reader_convert_scid_to_gene(self):
        pass  # todo implement

    def test_reader_get_phix_id(self):
        pass  # todo implement

    def test_reader_get_mitochondrial_ids(self):
        pass  # todo implement

    # @classmethod
    # def tearDownClass(cls):
    #     if os.path.isdir(cls.test_dir):
    #         shutil.rmtree(cls.test_dir)
    #     pass


### CORE TESTS ###

class CoreFixOutputPathsTest(unittest.TestCase):

    def setUp(self):
        self.target_prefix = os.path.expanduser('~') + '/testing_seqc/test_seqc'
        self.target_dir = os.path.expanduser('~') + '/testing_seqc/test_seqc/'


    def test_fix_output_paths_wrong_input_raises(self):
        self.assertRaises(TypeError, seqc.core.fix_output_paths, ('~/any_directory/',))

    def test_fix_output_paths_correctly_returns_abspaths(self):

        self.assertFalse(os.path.isdir(self.target_dir))

        # change to user directory
        cwd = os.getcwd()
        os.chdir(os.path.expanduser('~'))
        test_paths = ('testing_seqc/test_seqc', 'testing_seqc/../testing_seqc/test_seqc',
                      './testing_seqc/test_seqc', '~/testing_seqc/test_seqc')
        for p in test_paths:
            prefix, dir_ = seqc.core.fix_output_paths(p)
            self.assertEqual(dir_, self.target_dir)
            self.assertEqual(prefix, self.target_prefix)

        # get back to current working directory
        os.chdir(cwd)

    def tearDown(self):
        if os.path.isdir(self.target_dir):
            shutil.rmtree(self.target_dir)


class CoreCheckIndex(unittest.TestCase):

    def setUp(self):
        self.test_dir = 'test_seqc/'

    def test_check_index_wrong_input_raises(self):
        self.assertRaises(TypeError, seqc.core.check_index, ('index',), self.test_dir)

    @unittest.skip('Test downloads a large file, takes minutes')
    def test_check_index_downloads_when_not_present_and_creates_directory(self):

        # make sure these files are present in the downloaded directory
        critical_index_files = ['SA', 'SAindex', 'Genome', 'annotations.gtf',
                                'p_coalignment.pckl']

        # make sure directory deoesn't already exist
        self.assertFalse(os.path.isdir('test_seqc/'))

        # download index and check integrity
        index = seqc.core.check_index(config.index_link, 'test_seqc/')
        self.assertTrue(os.path.isdir('test_seqc'))
        for f in critical_index_files:
            self.assertTrue(os.path.isfile(f))

        # make sure correct directory was passed as index
        self.assertEqual(os.path.abspath('test_seqc/index/'), index)

    @unittest.skip('hangs; need to debug this')
    def test_download_non_index_file_raises_file_not_found_error(self):
        self.assertFalse(os.path.isdir('test_seqc/'))
        self.assertRaises(FileNotFoundError, seqc.core.check_index,
                          's3://dplab-data/seqc/test_seqc/', self.test_dir)

    def tearDown(self):
        if os.path.isdir(self.test_dir):
            shutil.rmtree(self.test_dir)


class CoreCheckLoadBarcodesTest(unittest.TestCase):
    """nose2 test_data generator"""

    def setUp(self):
        self.output_dir = 'test_seqc/'

    def test_load_barcodes_wrong_input_type_raises(self):
        self.assertRaises(TypeError, seqc.core.check_and_load_barcodes,
                          10, self.output_dir)

    @params(*config.data_types)
    @unittest.skip('Downloads from s3, takes time, bandwidth')
    def test_load_barcodes_from_aws_link(self, dtype):
        barcode_link = config.barcode_link_pattern % dtype
        self.assertFalse(os.path.isdir('test_seqc/'))
        barcodes = seqc.core.check_and_load_barcodes(barcode_link, self.output_dir)
        self.assertIsInstance(barcodes, seqc.barcodes.CellBarcodes)

    @params(*config.data_types)
    @unittest.skip('This function takes a long time to load the serialized object. Think '
                   'about improving its run time.')
    def test_load_barcodes_local(self, dtype):
        if dtype == 'drop_seq':
            return  # None should be passed, not a link that does not target a file obj.

        barcode_file = config.barcode_pattern % dtype
        barcodes = seqc.core.check_and_load_barcodes(barcode_file)
        barcodes = self.assertIsInstance(barcodes, seqc.barcodes.CellBarcodes)

    def test_load_pickle_containing_non_barcode_data(self):
        pass # todo implement


### CONVERT FEATURES TESTS ###

class ConvertFeaturesConvertGeneCoordinatesTest(unittest.TestCase):

    test_dir = 'test_seqc/'
    gtf = 'test_seqc/test_data.gtf'

    @classmethod
    def setUpClass(cls):

        # create the test_dir if it doesn't exist
        if not os.path.isdir(cls.test_dir):
            os.mkdir(cls.test_dir)

        # create a very small gtf file for testing
        i = 0
        with open(config.gtf, 'rb') as fin:
            with open(cls.gtf, 'wb') as fout:
                while i < 100:
                    fout.write(fin.readline())
                    i += 1

    @classmethod
    def tearDownClass(cls):
        if os.path.isdir(cls.test_dir):
            shutil.rmtree(cls.test_dir)

    def test_convert_features_empty_input_raises(self):
        self.assertRaises(ValueError, seqc.convert_features.ConvertGeneCoordinates,
                          {})

    def test_convert_features_wrong_input_type_raises(self):
        self.assertRaises(TypeError, seqc.convert_features.ConvertGeneCoordinates,
                          {'chr1': 'this_should_be_intervaltree_not_string'})

    def test_convert_features_from_gtf(self):
        cgc = seqc.convert_features.ConvertGeneCoordinates.from_gtf(self.gtf)

        # test_data that the method returns the correct object type
        print(type(cgc))
        self.assertIsInstance(cgc, seqc.convert_features.ConvertGeneCoordinates)

        # make sure that the object gets all of the genes in the test_data file
        pass  # todo implement


    def test_convert_features_reads_all_gene_records(self):
        pass  # todo implement



# todo all fastq functions below estimate_sequence_length are missing tests()

# todo keep
@unittest.skip('')
class TestJaggedArray(unittest.TestCase):

    def generate_input_iterable(self, n):
        for i in range(n):
            yield [random.randint(0, 5) for _ in range(random.randint(0, 5))]

    def test_jagged_array(self):
        n = int(1e3)
        data = list(self.generate_input_iterable(n))
        data_size = sum(len(i) for i in data)
        jarr = seqc.arrays.JaggedArray.from_iterable(data_size, data)
        self.assertTrue(jarr._data.dtype == np.uint32)
        self.assertTrue(jarr._data.shape == (data_size,))
        print(jarr[10])



# todo keep
@unittest.skip('')
class TestThreeBitInDrop(unittest.TestCase):

    def setUp(self):
        self.tb = seqc.three_bit.ThreeBit.default_processors('in-drop')
        # note: don't test_data with palindromic sequences, these can produce unexpected
        # errors, given the rotation occuring in the method.
        c11 = 'TAAAAAAA'
        c12 = 'TAAAAAAAA'
        c13 = 'TAAAAAAAAA'
        c14 = 'TAAAAAAAAAA'
        self.c1s = [c11, c12, c13, c14]
        self.c2 = 'GCCCCCCC'
        self.rmt = 'ACGGCT'
        self.example_seqs = [(c1 + 'GAGTGATTGCTTGTGACGCCTT' + self.c2 +
                             self.rmt + 'TTTAT') for c1 in self.c1s]

    def test_bin2str_inverts_str2bin(self):
        s = seqc.three_bit.ThreeBit.str2bin(self.example_seqs[0])
        seq = seqc.three_bit.ThreeBit.bin2str(s)
        self.assertTrue(seq == self.example_seqs[0])

    def test_bin2str_inverts_str2bin_with_ints2int(self):
        b1 = seqc.three_bit.ThreeBit.str2bin(self.c1s[0])
        b2 = seqc.three_bit.ThreeBit.str2bin(self.c1s[1])
        b12 = seqc.three_bit.ThreeBit.ints2int([b1, b2])
        self.assertTrue(b12 == seqc.three_bit.ThreeBit.str2bin(self.c1s[0] + self.c1s[1]))

    def test_str2bin_converts_forwards_not_backwards(self):
        s = seqc.three_bit.ThreeBit.str2bin('ACGT')
        self.assertEqual(s, 0b100110101011)
        self.assertNotEqual(s, 0b011101110100)

    def test_ints2int_and_cells(self):
        s1 = 'AAAAAAAA'
        s2 = 'TTTTTTTT'
        seqs = [seqc.three_bit.ThreeBit.str2bin(s) for s in [s1, s2]]
        c = seqc.three_bit.ThreeBit.ints2int(seqs)
        cells = self.tb.split_cell(c)
        c1, c2 = [seqc.three_bit.ThreeBit.bin2str(c) for c in cells]
        self.assertTrue(s1 == c1)
        self.assertTrue(s2 == c2)

    def test_extract_cell_spacer_extraction(self):

        results = [3446, 2478, 2869, 1894]
        for i, strseq in enumerate(self.example_seqs):
            # note this tracks code in extract cell; test_data will not reflect function
            # if extract cell is changed (not an ideal test_data)
            s = seqc.three_bit.ThreeBit.str2bin(strseq)
            bitlen = s.bit_length()

            # correct for leading T-nucleotide (011) which gets trimmed
            if bitlen % 3:
                bitlen += 1

            # extract spacer localizer
            tetramer = s >> (bitlen - 3 * 28) & 0o7777
            self.assertEqual(tetramer, results[i])

    def test_process_forward_read(self):
        for i, strseq in enumerate(self.example_seqs):
            cell, rmt, n_poly_t = self.tb.process_forward_sequence(strseq)

            # check cell
            str_cell = seqc.three_bit.ThreeBit.bin2str(cell)
            self.assertEqual(self.c1s[i] + self.c2, str_cell)

            # check rmt
            str_rmt = self.tb.bin2str(rmt)
            self.assertEqual(self.rmt, str_rmt)

            # check num_poly_t
            self.assertEqual(n_poly_t, 4)  # poly-t has a single 'A' contaminant

    def test_process_forward_read_truncated_input(self):
        full_len = 'TAAAAAAA' + 'GAGTGATTGCTTGTGACGCCTT' + 'GCCCCCCC' + 'ACGTAC' + 'TTTTT'
        no_poly_t = 'TAAAAAAA' + 'GAGTGATTGCTTGTGACGCCTT' + 'GCCCCCCC' + 'ACGTAC'
        short_rmt = 'TAAAAAAA' + 'GAGTGATTGCTTGTGACGCCTT' + 'GCCCCCCC' + 'ACG'
        cell_only = 'TAAAAAAA' + 'GAGTGATTGCTTGTGACGCCTT' + 'GCCCCCCC'
        no_data = 'TAAAAAAA' + 'GAGTGATTGCTTGTGACGCCTT' + 'GCCCCC'

        cell, rmt, n_poly_t = self.tb.process_forward_sequence(full_len)
        self.assertEqual(n_poly_t, 5)
        self.assertEqual(rmt, self.tb.str2bin('ACGTAC'))
        self.assertEqual(cell, self.tb.str2bin('TAAAAAAAGCCCCCCC'))

        cell, rmt, n_poly_t = self.tb.process_forward_sequence(no_poly_t)
        self.assertEqual(n_poly_t, 0)
        self.assertEqual(rmt, self.tb.str2bin('ACGTAC'))
        self.assertEqual(cell, self.tb.str2bin('TAAAAAAAGCCCCCCC'))

        cell, rmt, n_poly_t = self.tb.process_forward_sequence(short_rmt)
        self.assertEqual(n_poly_t, 0)
        self.assertEqual(rmt, 0)
        self.assertEqual(cell, self.tb.str2bin('TAAAAAAAGCCCCCCC'))

        cell, rmt, n_poly_t = self.tb.process_forward_sequence(cell_only)
        self.assertEqual(n_poly_t, 0)
        self.assertEqual(rmt, 0)
        self.assertEqual(cell, self.tb.str2bin('TAAAAAAAGCCCCCCC'))

        cell, rmt, n_poly_t = self.tb.process_forward_sequence(no_data)
        self.assertEqual(n_poly_t, 0)
        self.assertEqual(rmt, 0)
        self.assertEqual(cell, 0)


# todo keep
@unittest.skip('')
class TestThreeBitGeneral(unittest.TestCase):

    def test_3bit_mars_seq(self):
        self.assertRaises(NotImplementedError, seqc.three_bit.ThreeBit.default_processors,
                          'mars-seq')

    def test_3bit_drop_seq(self):
        cell = 'ACGTACGTACGT'
        rmt = 'AACCGGTT'
        poly_t = ''
        sample_seq = cell + rmt
        tb = seqc.three_bit.ThreeBit.default_processors('drop-seq')
        icell, irmt, num_poly_t = tb.process_forward_sequence(sample_seq)
        str_cell = tb.bin2str(icell)
        str_rmt = tb.bin2str(irmt)
        self.assertEqual(cell, str_cell)
        self.assertEqual(rmt, str_rmt)
        self.assertEqual(poly_t.count('T'), num_poly_t)

    def test_3bit_cel_seq(self):
        cell = 'ACGTAC'
        rmt = 'AACCGG'
        poly_t = 'TTTATTTAT'
        sample_seq = cell + rmt + poly_t
        tb = seqc.three_bit.ThreeBit.default_processors('cel-seq')
        icell, irmt, num_poly_t = tb.process_forward_sequence(sample_seq)
        str_cell = tb.bin2str(icell)
        str_rmt = tb.bin2str(irmt)
        self.assertEqual(cell, str_cell)
        self.assertEqual(rmt, str_rmt)
        self.assertEqual(poly_t.count('T'), num_poly_t)

    def test_3bit_avo_seq(self):
        cell = 'ACGTACGT'
        rmt = 'AACC'
        # adding N's breaks the poly-t somehow?
        poly_t = 'TAATTTNTTTT'
        sample_seq = cell + rmt + poly_t
        tb = seqc.three_bit.ThreeBit.default_processors('avo-seq')
        icell, irmt, num_poly_t = tb.process_forward_sequence(sample_seq)
        str_cell = tb.bin2str(icell)
        str_rmt = tb.bin2str(irmt)
        self.assertEqual(cell, str_cell)
        self.assertEqual(rmt, str_rmt)
        self.assertEqual(poly_t.count('T'), num_poly_t)

    def test_3bit_avo_seq_truncated_input(self):
        full_length = 'ACGTACGT' + 'AACC' 'TAATTTNTTTT'
        no_poly_t = 'ACGTACGT' + 'AACC'
        truncated_rmt = 'ACGTACGT' + 'AAC'
        no_data = 'ACGTAC'

        tb = seqc.three_bit.ThreeBit.default_processors('avo-seq')

        # results
        cell = tb.str2bin('ACGTACGT')
        rmt = tb.str2bin('AACC')
        n_poly_t = 8

        self.assertEqual(tb.process_forward_sequence(full_length), (cell, rmt, n_poly_t))
        self.assertEqual(tb.process_forward_sequence(no_poly_t), (cell, rmt, 0))
        self.assertEqual(tb.process_forward_sequence(truncated_rmt), (cell, 0, 0))
        self.assertEqual(tb.process_forward_sequence(no_data), (0, 0, 0))

    def test_gc_content(self):
        test_string = 'TGCGCAAAAG'
        expected_result = 0.5
        bin_string = seqc.three_bit.ThreeBit.str2bin(test_string)
        result = seqc.three_bit.ThreeBit.gc_content(bin_string)
        self.assertEqual(result, expected_result)


# todo keep; speed up
@unittest.skip('')
class TestThreeBitCellBarcodes(unittest.TestCase):

    def setUp(self):
        self.barcode_dir = ('/'.join(seqc.__file__.split('/')[:-2]) +
                            '/data/in_drop/barcodes/')

    @unittest.skip('creation of the cb_3bit.p object is slow')
    def test_create_in_drop_cb(self):

        # create barcode files
        barcode_files = [self.barcode_dir + 'cb1.txt', self.barcode_dir + 'cb2.txt']
        cb = seqc.barcodes.CellBarcodes(*barcode_files)
        with open(self.barcode_dir + 'cb_3bit.p', 'wb') as f:
            pickle.dump(cb, f)

    def test_created_in_drop_cb(self):
        with open(self.barcode_dir + 'cb_3bit.p', 'rb') as f:
            cb = pickle.load(f)

        # now check and make sure the object works
        with open(self.barcode_dir + 'cb1.txt') as f:
            cb1s = [l.strip() for l in f.readlines()]
        with open(self.barcode_dir + 'cb2.txt') as f:
            cb2s = [l.strip() for l in f.readlines()]

        # get a random cb1 and a random cb2
        one = random.choice(cb1s)
        two = random.choice(cb2s)
        umi = 'AACCGG'

        # convert to 3bit
        tb = seqc.three_bit.ThreeBit.default_processors('in-drop')
        one3bit = tb.str2bin(one)
        two3bit = tb.str2bin(two)
        cell = tb.ints2int([one3bit, two3bit])

        # check that this does the same thing as getting the cell fom _extract_cell
        seq2 = one + tb._w1 + two + umi + 'T' * 5
        cell2, rmt, n_poly_t = tb.process_forward_sequence(seq2)
        self.assertEqual(cell, cell2)

        # check that this sequence is in the cell barcodes
        self.assertTrue(cb.perfect_match(cell))
        self.assertTrue(cb.close_match(cell))

    def test_in_drop_error_code_construction_and_recognition(self):
        cb1 = StringIO('AAA\nTTT\n')
        cb2 = StringIO('CCC\nGGG\n')
        cb = seqc.barcodes.CellBarcodes(cb1, cb2)
        perfect_code = seqc.three_bit.ThreeBit.str2bin('AAACCC')
        error_code = seqc.three_bit.ThreeBit.str2bin('NTTCCC')
        self.assertTrue(cb.perfect_match(perfect_code))
        self.assertFalse(cb.perfect_match(error_code))
        self.assertTrue(cb.close_match(perfect_code))
        self.assertTrue(cb.close_match(error_code))

    def test_in_drop_error_code_error_type_identification(self):
        cb1 = StringIO('AAA\nTTT\n')
        cb2 = StringIO('CCC\nGGG\n')
        cb = seqc.barcodes.CellBarcodes(cb1, cb2)
        error_code = seqc.three_bit.ThreeBit.str2bin('NTTCCC')

        errors = cb.map_errors(error_code)
        self.assertEqual({'TN'}, set(errors))


# todo import these tests from scseq/seqdb
@unittest.skip('')
class TestIndexGeneration(unittest.TestCase):
    pass


# todo
# contains memory and time profiling examples for other unit tests; final is a
# relevant test_data
@unittest.skip('')
class TestSamToReadArray(unittest.TestCase):

    def setUp(self):
        self.data_dir = '/'.join(seqc.__file__.split('/')[:-2]) + '/data/'
        self.in_drop_samfile = (
            self.data_dir + 'test_seqc_merge_in_drop_fastq/Aligned.out.sam')
        self.drop_seq_samfile = (
            self.data_dir + 'test_seqc_merge_drop_seq_fastq/Aligned.out.sam')
        self.gtf = self.data_dir + 'genome/mm38_chr19/annotations.gtf'
        self.ftable, self.fpositions = seqc.sam.construct_feature_table(
            self.gtf, fragment_len=1000)
        self.h5name = self.data_dir + 'test_seqc_merge_in_drop_fastq/ra.h5'

    @unittest.skip('')
    def test_construct_read_array(self):
        a = seqc.arrays.ReadArray.from_samfile(
            self.in_drop_samfile, self.fpositions, self.ftable)
        array_lengths = (len(a.data), len(a.features), len(a.positions))
        self.assertTrue(len(set(array_lengths)) == 1)  # all arrays should have same #rows

    @unittest.skip('')
    def test_generate_indices_for_error_correction(self):
        a = seqc.arrays.ReadArray.from_samfile(
            self.in_drop_samfile, self.fpositions, self.ftable)
        _ = (len(a.data), len(a.features), len(a.positions))
        _ = a.group_for_error_correction(required_poly_t=4)

    @unittest.skip('')
    def test_save_and_load_readarray(self):
        ra = seqc.arrays.ReadArray.from_samfile(
            self.in_drop_samfile, self.fpositions, self.ftable)
        ra.save_h5(self.h5name)
        ra2 = seqc.arrays.ReadArray.from_h5(self.h5name)

        self.assertTrue(np.array_equiv(ra.data, ra2.data))
        self.assertTrue(np.array_equiv(ra._features._data, ra2._features._data))
        self.assertTrue(np.array_equiv(ra._features._index, ra2._features._index))
        self.assertTrue(np.array_equiv(ra._positions._data, ra2._positions._data))
        self.assertTrue(np.array_equiv(ra._positions._index, ra2._positions._index))

    # def memory_usage(self, sizes, directory, index):
    #
    #     # generate fastq files
    #     for s in sizes:
    #         fastq_filename = directory + str(s)
    #         seqc.fastq.GenerateFastq.in_drop(s, fastq_filename)
    #     forward = [directory + str(s) + '_r1.fastq' for s in sizes]
    #     reverse = [directory + str(s) + '_r2.fastq' for s in sizes]
    #
    #     # get cell barcodes
    #     with open(self.data_dir + 'in_drop/barcodes/cb_3bit.p', 'rb') as f:
    #         cb = pickle.load(f)
    #
    #     # merge the data
    #     merged_files = []
    #     for f, r in zip(forward, reverse):
    #         fq, _ = seqc.fastq.merge_fastq([f], [r], 'in-drop', directory, cb)
    #         merged_files.append(fq)
    #
    #     # align the data
    #     samfiles = seqc.align.STAR.align_multiple_files(
    #         merged_files, index, 7, directory
    #     )
    #
    #     ft, fp = seqc.convert_features.construct_feature_table(self.gtf, 1000)
    #
    #     logfile = open(directory + 'report.txt', 'w+')
    #
    #     h5files = [directory + '%dm.h5' % s for s in sizes]
    #
    #     # define and run functions
    #     def test_memory_usage(idx):
    #         seqc.memory_usage((seqc.arrays.ReadArray.from_samfile, (samfiles[idx], ft, fp)))
    #         # arr.save_h5(h5files[index])
    #
    #     for i in range(len(samfiles)):
    #         test_memory_usage(i)
    #
    #     logfile.close()
    #
    # @unittest.skip('')
    # def test_ra_memory_usage_small(self):
    #     index = self.data_dir + 'genome/mm38_chr19/'
    #     working_directory = self.data_dir + 'test_ra_memory_usage/'
    #     if not os.path.isdir(working_directory):
    #         os.makedirs(working_directory)
    #     self.memory_usage([int(1e6)], working_directory, index)
    #
    # @unittest.skip('')
    # def test_profile_mem_usage(self):
    #     samfile = ('/Users/ambrose/PycharmProjects/SEQC/src/data/test_ra_memory_usage/'
    #                'merged_temp/Aligned.out.sam')
    #     ft, fp = seqc.convert_features.construct_feature_table(self.gtf, 1000)
    #     usage = seqc.memory_usage((seqc.arrays.ReadArray.from_samfile, (samfile, ft, fp)))
    #     print(np.array(usage))
    #
    # @unittest.skip('')
    # def test_ra_memory_usage(self):
    #     data_dir = '/'.join(seqc.__file__.split('/')[:-2]) + '/data/'
    #
    #     # generate fastq data for testing
    #     files = [data_dir + 'test_memory_usage/' + prefix for prefix in
    #              ['1m', '2m', '4m', '8m', '16m']]
    #     index = data_dir + 'genome/mm38_chr19/'
    #     seqc.generate_in_drop_fastq_data(int(1e6), files[0])
    #     seqc.generate_in_drop_fastq_data(int(2e6), files[1])
    #     seqc.generate_in_drop_fastq_data(int(4e6), files[2])
    #     seqc.generate_in_drop_fastq_data(int(8e6), files[3])
    #     seqc.generate_in_drop_fastq_data(int(16e6), files[4])
    #
    #     # align the data, generating sam files
    #     samfiles = seqc.align.STAR.align_multiple_files(
    #         files, index, 7, data_dir + 'test_memory_usage/')
    #
    #     ft, fp = seqc.convert_features.construct_feature_table(self.gtf, 1000)
    #
    #     logfile = open(data_dir + 'test_memory_usage/report.txt', 'w+')
    #
    #     h5files = [f + '.h5' for f in files]
    #
    #     @profile(stream=logfile)
    #     def test_1m():
    #         arr = seqc.arrays.ReadArray.from_samfile(samfiles[0], ft, fp)
    #         arr.save_h5(h5files[0])
    #
    #     @profile(stream=logfile)
    #     def test_2m():
    #         arr = seqc.arrays.ReadArray.from_samfile(samfiles[1], ft, fp)
    #         arr.save_h5(h5files[1])
    #
    #     @profile(stream=logfile)
    #     def test_4m():
    #         arr = seqc.arrays.ReadArray.from_samfile(samfiles[2], ft, fp)
    #         arr.save_h5(h5files[2])
    #
    #     @profile(stream=logfile)
    #     def test_8m():
    #         arr = seqc.arrays.ReadArray.from_samfile(samfiles[3], ft, fp)
    #         arr.save_h5(h5files[3])
    #
    #     @profile(stream=logfile)
    #     def test_16m():
    #         arr = seqc.arrays.ReadArray.from_samfile(samfiles[4], ft, fp)
    #         arr.save_h5(h5files[4])
    #
    #     for f in [test_1m, test_2m, test_4m, test_8m, test_16m]:
    #         f()
    #
    #     logfile.close()

    @unittest.skip('')
    def test_create_ra_object(self):
        samfile = ('/Users/ambrose/PycharmProjects/SEQC/src/data/test_ra_memory_usage/'
                   'merged_temp/Aligned.out.sam')
        fc = seqc.convert_features.ConvertFeatureCoordinates.fromCfg.gtf(self.gtf, 1000)
        res = seqc.arrays.ReadArray.from_samfile(samfile, fc)

    # def test_profile_counts_creation(self):
    #     ra = arrays.ReadArray.from_h5(self.h5name)
    #     pr = cProfile.Profile()
    #     pr.enable()
    #     ra.to_sparse_counts(collapse_molecules=True, n_poly_t_required=0)
    #     pr.disable()
    #     p = Stats(pr)
    #     p.strip_dirs()
    #     p.sort_stats('cumtime')
    #     p.print_stats()


# todo rewrite for standard locations
@unittest.skip('')
class TestResolveAlignments(unittest.TestCase):
    """this may need more tests for more sophisticated input data."""

    def setUp(self):
        self.data_dir = '/'.join(seqc.__file__.split('/')[:-2]) + '/data/'
        self.expectations = self.data_dir + 'genome/mm38_chr19/p_coalignment.pckl'

    @unittest.skip('only run this if new ra disambiguation input is needed')
    def test_generate_disambiguation_read_array(self):
        cell_barcodes = self.data_dir + 'in_drop/barcodes/cb_3bit.p'
        save = self.data_dir + 'in_drop/disambiguation_ra_input.p'
        _ = seqc.generate_in_drop_read_array(self.expectations, cell_barcodes, 1000, 10,
                                        save=save)

    # @unittest.skip('')
    def test_ra_disambiguate(self):
        expectations = self.data_dir + 'genome/mm38_chr19/p_coalignment_array.p'
        arr_pickle = self.data_dir + 'in_drop/disambiguation_ra_input.p'
        with open(arr_pickle, 'rb') as f:
            ra = pickle.load(f)
        ra.resolve_alignments(expectations)
        # 4 == complete disambiguation.
        vals, counts = np.unique(ra.data['disambiguation_results'], return_counts=True)
        self.assertTrue(np.array_equal(vals, np.array([2, 4])))


@unittest.skip('')
class TestGeneTable(unittest.TestCase):

    def setUp(self):
        self.genome_dir = '/'.join(seqc.__file__.split('/')[:-2]) + '/data/genome/'

    def test_gene_table(self):
        gt = seqc.convert_features.GeneTable(config.gtf)
        chromosome = 'chr19'
        start = 4007800
        end = 4007900
        strand = '-'
        genes = gt.coordinates_to_gene_ids(chromosome, start, end, strand)


# todo rewrite; keep
@unittest.skip('')
class TestSamToCount(unittest.TestCase):

    def setUp(self):
        self.data_dir = '/'.join(seqc.__file__.split('/')[:-2]) + '/data/'

    def test_sam_to_count_single_file(self):
        samfile = self.data_dir + 'test_seqc_merge_drop_seq_fastq/Aligned.out.sam'
        gtf = self.data_dir + '/genome/mm38_chr19/annotations.gtf'
        coo, gene_index, ci = seqc.sam.to_count_single_file(samfile, gtf, 100)
        print(repr(coo))
        print(len(gene_index))

        samfile = self.data_dir + 'test_seqc_merge_in_drop_fastq/Aligned.out.sam'
        gtf = self.data_dir + '/genome/mm38_chr19/annotations.gtf'
        coo, gene_index, ci = seqc.sam.to_count_single_file(samfile, gtf, 100)
        print(repr(coo))
        print(len(gene_index))

    def test_sam_to_count_multiple_files(self):
        samfile = self.data_dir + 'test_seqc_merge_drop_seq_fastq/Aligned.out.sam'
        gtf = self.data_dir + '/genome/mm38_chr19/annotations.gtf'
        coo, gene_index, ci = seqc.sam.to_count_multiple_files([samfile, samfile], gtf, 100)
        print(repr(coo))
        print(len(gene_index))

        samfile = self.data_dir + 'test_seqc_merge_in_drop_fastq/Aligned.out.sam'
        gtf = self.data_dir + '/genome/mm38_chr19/annotations.gtf'
        coo, gene_index, ci = seqc.sam.to_count_multiple_files([samfile, samfile], gtf, 100)
        print(repr(coo))
        print(len(gene_index))


@unittest.skip('')
class TestS3Client(unittest.TestCase):

    def setUp(self):
        # create a file
        self.txt1 = 'testing aws s3\n'
        self.txt2 = 'testing aws s3 a second time\n'
        self.txt3 = 'testing a nested file\n'
        os.mkdir('.test_aws_s3')
        self.file1 = '.test_aws_s3/test_aws_s3_file1.txt'
        self.file2 = '.test_aws_s3/test_aws_s3_file2.txt'
        os.mkdir('.test_aws_s3/test2')
        self.file3 = '.test_aws_s3/test2/subfile.txt'
        with open(self.file1, 'w') as f:
            f.write(self.txt1)
        with open(self.file2, 'w') as f:
            f.write(self.txt2)
        with open(self.file3, 'w') as f:
            f.write(self.txt3)

    @unittest.skip('')
    def test_aws_s3_single_upload(self):
        bucket = 'dplab-home'
        key1 = 'ajc2205/testing_aws/'
        seqc.io.S3.upload_file(self.file1, bucket=bucket, key=key1)

        download_key = key1 + self.file1.split('/')[-1]
        fout = '.test_aws_s3/recovered_file1.txt'

        seqc.io.S3.download_file(bucket=bucket, key=download_key, fout=fout)
        with open(fout, 'r') as f:
            self.assertEqual(f.read(), self.txt1)

        seqc.io.S3.remove_file(bucket, download_key)

        os.remove('.test_aws_s3/recovered_file1.txt')

    def test_aws_s3_multiple_upload(self):
        # ResourceWarnings are generated by an interaction of io.S3.list() and
        # the unittesting suite. These should not occur in normal usage.
        bucket = 'dplab-home'
        key_prefix = 'ajc2205/testing_aws/'
        file_prefix = '.test_aws_s3/*'

        # upload file1, file2, and file3
        seqc.io.S3.upload_files(file_prefix, bucket, key_prefix)

        aws_dir = set(seqc.io.S3.list(bucket, key_prefix))
        aws_files = {key_prefix + f for f in
                     ['test_aws_s3_file1.txt', 'test_aws_s3_file2.txt',
                      'test2/subfile.txt']}
        self.assertEqual(aws_files, aws_dir)

        # download the files again
        output_prefix = '.test_aws_s3/download_test/'
        seqc.io.S3.download_files(bucket, key_prefix, output_prefix)
        all_downloaded_files = []
        for path, subdirs, files in os.walk(output_prefix):
            for name in files:
                all_downloaded_files.append(os.path.join(path, name))

        # remove the files that we uploaded
        seqc.io.S3.remove_files(bucket, key_prefix)

    def tearDown(self):
        shutil.rmtree('.test_aws_s3/')


@unittest.skip('')
class TestGenerateSRA(unittest.TestCase):

    def setUp(self):
        self.data_dir = '/'.join(seqc.__file__.split('/')[:-2]) + '/data/'
        self.sra_dir = self.data_dir + '.test_dump_sra/'
        if not os.path.isdir(self.sra_dir):
            os.makedirs(self.sra_dir)

    @unittest.skip("doesn't work")
    def test_md5sum(self):
        res = seqc.sra.SRAGenerator.md5sum(self.data_dir + 'in_drop/sample_data_r1.fastq')
        print(res)

    @unittest.skip('')
    def test_generate_experiment_xml_file(self):

        seqc.sra.SRAGenerator.make_experiment_xml('mm38', 'DUMMY_SAMPLE',
                                         fout_stem=self.sra_dir + 'EXP_XML')
        xml_data = xml.dom.minidom.parse(self.sra_dir + 'EXP_XML_experiment.xml')
        pretty = xml_data.toprettyxml()
        print('\n', pretty)

    @unittest.skip('')
    def test_generate_run_xml_file(self):
        forward_fastq = ['testfile_r1.fastq']
        reverse_fastq = ['testfile_r2.fastq']
        forward_checksum = [str(912301292)]
        reverse_checksum = [str(102931209)]
        data_block = str(1)
        fout_stem = self.sra_dir + 'RUN_XML'
        experiment_alias = 'EXP_XML_experiment.xml'
        seqc.sra.SRAGenerator.make_run_xml(
            forward_fastq, experiment_alias, forward_checksum, fout_stem, data_block,
            reverse_fastq, reverse_checksum)
        xml_data = xml.dom.minidom.parse(self.sra_dir + 'RUN_XML_run.xml')
        pretty = xml_data.toprettyxml()
        print('\n', pretty)

    # @unittest.skip('')
    def test_generate_sra_paired_end(self):
        data_dir = '/'.join(seqc.__file__.split('/')[:-2]) + '/data/'
        forward_fastq = [data_dir + 'in_drop/sample_data_r1.fastq']
        forward_md5 = ['DUMMY' for _ in forward_fastq]
        reverse_fastq = [data_dir + 'in_drop/sample_data_r2.fastq']
        reverse_md5 = ['DUMMY' for _ in reverse_fastq]
        file_directory = data_dir + 'in_drop/'
        experiment_alias = 'in_drop_test_experiment_2014'
        reference = 'mm38'
        platform_name = 'ILLUMINA'
        instrument = 'Illumina HiSeq 2500'
        single_or_paired_end = 'PAIRED'
        fout_stem = self.sra_dir + experiment_alias
        exp_xml = seqc.sra.SRAGenerator.make_experiment_xml(
            reference, experiment_alias, platform_name, single_or_paired_end, instrument,
            fout_stem)
        run_xml = seqc.sra.SRAGenerator.make_run_xml(
            forward_fastq, experiment_alias, forward_md5, fout_stem, data_block=None,
            reverse_fastq_files=reverse_fastq, reverse_checksum_results=reverse_md5)
        output_path = data_dir + 'test_generate_sra'
        seqc.sra.SRAGenerator.fastq_load(file_directory, run_xml, exp_xml, output_path)


# @unittest.skip('')
# class TestProcessSingleFileSCSEQExperiment(unittest.TestCase):
#
#     def setUp(self):
#         self.forward, self.reverse = check_fastq('in_drop')
#         self.s3_bucket = 'dplab-home'
#         self.s3_key = 'ajc2205/test_in_drop.npz'
#
#     @unittest.skip('')
#     def test_process_single_file_no_sra_download(self):
#
#         # set some variables
#         index_bucket = None
#         index_key = None
#
#         experiment_name = 'test_in_drop'
#         s3_bucket = self.s3_bucket
#         s3_key = self.s3_key
#         cell_barcodes = ('/Users/ambrose/PycharmProjects/SEQC/src/data/in_drop/barcodes/'
#                          'in_drop_barcodes.p')
#
#         # set the index
#         if not config.index:  # download the index
#             index_dir = working_directory + 'index/'
#             S3.download_files(bucket=index_bucket, key_prefix=index_key,
#                               output_prefix=index_dir, no_cut_dirs=True)
#             index = index_dir + index_key.lstrip('/')
#         if not os.path.isdir(index):
#             raise FileNotFoundError('Index does not lead to a directory')
#
#         # merge fastq files
#         merged_fastq, _ = fastq.merge_fastq(
#             self.forward, self.reverse, 'in-drop', self.working_directory, cell_barcodes)
#
#         # align the data
#         sam_file = STAR.align(
#             merged_fastq, index, n_threads, working_directory, reverse_fastq_file=None)
#
#         # create the matrix
#         gtf_file = index + 'annotations.gtf'
#         coo, rowind, colind = sam_to_count_single_file(sam_file, gtf_file)
#
#         numpy_archive = experiment_name + '.npz'
#         with open(numpy_archive, 'wb') as f:
#             np.savez(f, mat=coo, row=rowind, col=colind)
#
#         # upload the matrix to amazon s3
#         S3.upload_file(numpy_archive, s3_bucket, s3_key)
#
#     def test_process_multiple_file_no_sra_download(self):
#         # set some variables
#         index = self.index
#         working_directory = self.working_directory
#         index_bucket = None
#         index_key = None
#         S3 = io.S3
#         STAR = align.STAR
#         n_threads = 7
#         sam_to_count_multiple_files = qc.sam_to_count_multiple_files
#         experiment_name = 'test_in_drop'
#         s3_bucket = self.s3_bucket
#         s3_key = self.s3_key
#         cell_barcodes = config.barcode_pattern % dtype
#
#         # potential issue: reverse should never map..
#         forward = [self.forward[0]] * 3
#         reverse = [self.reverse[0]] * 3
#
#         # set the index
#         if not index:  # download the index
#             index_dir = working_directory + 'index/'
#             S3.download_files(bucket=index_bucket, key_prefix=index_key,
#                               output_prefix=index_dir, no_cut_dirs=True)
#             index = index_dir + index_key.lstrip('/')
#         if not os.path.isdir(index):
#             raise FileNotFoundError('Index does not lead to a directory')
#
#         # align the data
#         sam_files = STAR.align_multiple_files(
#             forward, index, n_threads, working_directory, reverse_fastq_files=reverse)
#
#         # create the matrix
#         gtf_file = index + 'annotations.gtf'
#         coo, rowind, colind = sam_to_count_multiple_files(sam_files, gtf_file)
#
#         numpy_archive = experiment_name + '.npz'
#         with open(numpy_archive, 'wb') as f:
#             np.savez(f, mat=coo, row=rowind, col=colind)
#
#         # upload the matrix to amazon s3
#         S3.upload_file(numpy_archive, s3_bucket, s3_key)


@unittest.skip('')
class TestGroupForErrorCorrection(unittest.TestCase):

    def setUp(self):
        # create a dummy array with 10 records
        dtype = seqc.arrays.ReadArray._dtype
        dummy_data = np.zeros((20,), dtype=dtype)

        # we'll use these cell barcodes; no palindromes
        s2b = seqc.three_bit.ThreeBit.str2bin
        cb11 = s2b('CACGGACA')
        cb21 = s2b('GTGTGGTT')
        cb12 = s2b('TTTCCTTT')
        cb22 = s2b('AGAAAGGAAA')

        # we'll use these rmts, again no palindromes:
        rmt1 = s2b('AGGTTC')
        rmt2 = s2b('GGATAC')

        # convert everything to ints the same way as in a fastq
        # in fastq, we first convert each to an int, and concatenate
        # the integers using ints2int; we'll create two cells this way
        i2i = seqc.three_bit.ThreeBit.ints2int
        c1 = i2i([cb11, cb21])
        c2 = i2i([cb12, cb22])

        # set defaults for the other fields
        def create_record(cell, rmt):
            n_poly_t = 4
            valid_cell = True
            dust_score = 0
            rev_quality = 40
            fwd_quality = 40
            is_aligned = True
            alignment_score = 100
            return (cell, rmt, n_poly_t, valid_cell, dust_score, rev_quality,
                    fwd_quality, is_aligned, alignment_score)

        # position isn't used, so we'll ignore it, but we need features
        f1 = [10, 11]
        f2 = [5]

        # create an array with 20 records: 10 each for f1 and f2,
        # for each, 4 are c2, 6 are c2, each are evenly split among rmts
        features = [f1] * 10 + [f2] * 10
        positions = [[0]] * 20
        data = [create_record(c1, rmt1)] * 2
        data += [create_record(c1, rmt2)] * 2
        data += [create_record(c2, rmt1)] * 3
        data += [create_record(c2, rmt2)] * 3
        data += [create_record(c1, rmt1)] * 2
        data += [create_record(c1, rmt2)] * 2
        data += [create_record(c2, rmt1)] * 3
        data += [create_record(c2, rmt2)] * 3

        for i, r in enumerate(data):
            dummy_data[i] = r

        positions = seqc.arrays.JaggedArray.from_iterable(positions)
        features = seqc.arrays.JaggedArray.from_iterable(features)

        self.ra = seqc.arrays.ReadArray(dummy_data, features, positions)

    @unittest.skip('')
    def test_ints2int(self):
        i2i = seqc.three_bit.ThreeBit.ints2int
        for record in self.ra.data:
            print(type(record['cell'].astype(int)))
            print(i2i([record['cell'].astype(int), record['rmt'].astype(int)]))

    # @unittest.skip('')
    def test_group_for_error_correction(self):
        grp = self.ra.group_for_error_correction()
        print(self.ra.data)
        print(grp)
        print(grp.keys())
        b2s = seqc.three_bit.ThreeBit.bin2str
        for feature in grp.keys():
            for seq in grp[feature].keys():
                if seq:
                    print(seq)
                    print(b2s(seq))
                else:
                    print('Returned None')


# todo keep; rewrite to use check_sam()
@unittest.skip('')
class TestParallelConstructSam(unittest.TestCase):

    @unittest.skip('')
    def test_parallel_construct_sam(self):
        h5_name = 'test_data.h5'
        n_processes = 7
        chunk_size = 10000
        seqc.sam.to_h5(self.samfile, h5_name, n_processes, chunk_size,
                       self.gtf, fragment_length=1000)

    @unittest.skip('')
    def test_writing_non_parallel(self):
        h5_name = 'test_data.h5'
        n_processes = 7
        chunk_size = 100000
        # seqc.sam.to_h5(self.samfile, h5_name, n_processes, chunk_size,
        #                self.gtf, fragment_length=1000)

        filesize = os.stat(self.samfile).st_size

        # get a bunch of records to check average size of a record
        with open(self.samfile) as f:
            records = []
            for line in f:
                if line.startswith('@'):
                    continue
                records.append(line)
                if len(records) > 1000:
                    break
        average_record_size = np.mean([len(r) for r in records])

        # estimate expected rows for the h5 database
        expectedrows = filesize / average_record_size

        # create a feature_converter object to convert genomic -> transcriptomic features
        fc = seqc.convert_features.ConvertFeatureCoordinates.fromCfg.gtf(self.gtf, 1000)

        # open all our tables
        blosc5 = tb.Filters(complevel=5, complib='blosc')
        h5f = tb.open_file(h5_name, mode='w', filters=blosc5,
                           title='Data for ReadArray constructed from %s' % self.samfile)
        try:
            # description = tb.descr_from_dtype(_dtype)
            a = tb.UInt32Atom()
            # create the tables and arrays needed to store data
            h5f.create_table(h5f.root, 'data', seqc.sam._DataTable, 'ReadArray.data',
                             filters=blosc5, expectedrows=expectedrows)
            h5f.create_earray(h5f.root, 'index', a, (0,), filters=blosc5,
                              expectedrows=expectedrows)
            h5f.create_earray(h5f.root, 'features', a, (0,), filters=blosc5,
                              expectedrows=expectedrows)
            h5f.create_earray(h5f.root, 'positions', a, (0,), filters=blosc5,
                              expectedrows=expectedrows)

            itersam = seqc.sam._iterate_chunk(self.samfile, n=chunk_size)
            for chunk in itersam:
                processed = seqc.sam._process_chunk(chunk, fc)
                seqc.sam._write_chunk(processed, h5f)
            # read_kwargs = dict(samfile=self.samfile, n=chunk_size)
            # process_kwargs = dict(feature_converter=fc)
            # write_kwargs = dict(h5_file=h5f)
            # seqc.parallel.process_parallel(n_processes, sam._iterate_chunk,
            #                                sam._process_chunk,
            #                                sam._write_chunk, read_kwargs=read_kwargs,
            #                                process_kwargs=process_kwargs,
            #                                write_kwargs=write_kwargs)
        finally:
            h5f.close()

    @unittest.skip('')
    def test_writing_non_parallel_writeobj(self):
        h5_name = 'test_data.h5'
        chunk_size = 1000000

        filesize = os.stat(self.samfile).st_size

        # get a bunch of records to check average size of a record
        with open(self.samfile) as f:
            records = []
            for line in f:
                if line.startswith('@'):
                    continue
                records.append(line)
                if len(records) > 1000:
                    break
        average_record_size = np.mean([len(r) for r in records])

        # estimate expected rows for the h5 database
        expectedrows = filesize / average_record_size

        # create a feature_converter object to convert genomic -> transcriptomic features
        fc = seqc.convert_features.ConvertFeatureCoordinates.fromCfg.gtf(self.gtf, 1000)

        # open all our tables
        h5f = seqc.sam.ReadArrayH5Writer(h5_name)
        h5f.create(expectedrows)
        itersam = seqc.sam._iterate_chunk(self.samfile, n=chunk_size)
        for chunk in itersam:
            processed = seqc.sam._process_chunk(chunk, fc)
            h5f.write(processed)
        h5f.close()

    def test_writing_parallel_writeobj(self):
        h5_name = 'test_data.h5'
        n_processes = 7
        chunk_size = 100000
        seqc.sam.to_h5(self.samfile, h5_name, n_processes, chunk_size,
                       self.gtf, fragment_length=1000)
        nlines = self.samfile


@unittest.skip('')
class TestCounting(unittest.TestCase):

    def setUp(self):

        # design a UniqueReadArray to test_data validity of method
        dtype = seqc.arrays.ReadArray._dtype

        # several fields don't matter; only using cell, rmt, features
        def create_records(rmt, cell, feature, n):
            n_poly_t = 5
            valid_cell = True
            dust_score = 0
            rev_quality = 40
            fwd_quality = 40
            is_aligned = True
            alignment_score = 100
            records = [(cell, rmt, n_poly_t, valid_cell, dust_score, rev_quality,
                        fwd_quality, is_aligned, alignment_score) for _ in range(n)]
            positions = np.zeros(n)
            features_ = np.array([feature] * n)
            return records, features_, positions

        def create_unique_read_array(rmts, cells, features):
            assert len(rmts) == len(cells) == len(features)
            # start with an easy case where the method should get things right
            # use these ids

            n = len(rmts)
            i = 0
            j = 0
            udata = np.zeros((n,), dtype=dtype)
            ufeatures = np.zeros(n, dtype=np.int32)
            upositions = np.zeros(n, dtype=np.int32)
            for (r, c, f) in zip(rmts, cells, features):
                recs, feats, posn = create_records(r, c, f, 1)
                for v in recs:
                    udata[i] = v
                    i += 1
                ufeatures[j: j + len(feats)] = feats
                upositions[j: j + len(posn)] = posn
                j += len(feats)

            # these are all unique; should fail when required support > 0 and return empty
            # molecule and read arrays
            return seqc.arrays.UniqueReadArray(udata, ufeatures, upositions)

        rmts = [10] * 9 + [6] * 9 + [111] * 9
        cells = [1, 2, 3] * 9
        features = [7, 7, 7, 8, 8, 8, 9, 9, 9] * 3
        self.simple_unique = create_unique_read_array(rmts, cells, features)
        self.simple_duplicate = create_unique_read_array(rmts * 2, cells * 2, features * 2)

    def test_print_data(self):
        self.simple_unique._sort()
        data = append_fields(self.simple_unique.data, 'features', self.simple_unique.features)
        # print(data[self.simple_unique._sorted][['cell', 'features', 'rmt']])

    def test_counting_filter_too_high(self):

        # test_data that sorting is working
        # seems to be sorting cell, feature, rmt
        self.assertRaises(ValueError, self.simple_unique.to_experiment, 1)

    def test_counting_simple(self):
        exp = self.simple_unique.to_experiment(0)
        rdata = np.array(exp.reads.counts.todense())
        mdata = np.array(exp.molecules.counts.todense())
        # print(pd.DataFrame(rdata, exp.reads.index, exp.reads.columns))
        # print(pd.DataFrame(mdata, exp.molecules.index, exp.molecules.columns))

    def test_counting_doublets(self):
        exp = self.simple_duplicate.to_experiment(0)
        rdata = np.array(exp.reads.counts.todense())
        mdata = np.array(exp.molecules.counts.todense())
        # print(pd.DataFrame(rdata, exp.reads.index, exp.reads.columns))
        # print(pd.DataFrame(mdata, exp.molecules.index, exp.molecules.columns))


@unittest.skip('')
class TestUniqueArrayCreation(unittest.TestCase):
    """
    suspicion -> unique array creation is breaking something in the pipeline;
    could be an earlier method but lets check upstream so we know how each part is
    working.

    method: find the smallest real dataset that breaks it and use that to test_data.
    """

    def test_num_unique_samfile(self):
        fc = seqc.convert_features.ConvertFeatureCoordinates.fromCfg.gtf(config.gtf, 1000)
        # todo fix this to use check_sam()
        ra = seqc.arrays.ReadArray.from_samfile(config.samfile, fc)

        # verify that the reads aren't all unique
        # n_unique = 21984; n = 29999
        n_unique = sum([1 for f in ra.features if f.shape == (1,)])
        self.assertTrue(n_unique > 0, 'Zero unique reads.')
        self.assertTrue(n_unique < ra.shape[0], 'All reads were unique.')

        # determine if JaggedArray to_unique() is generating the right result
        unique = ra.features.to_unique()
        self.assertTrue(unique.shape[0] == n_unique, '%d != %d' %
                        (unique.shape[0], n_unique))

        # determine if ReadArray to_unique() is generating the right result
        fbool = ra.features.is_unique()  # not yet tested explicitly, but implicitly working based on jagged.to_unique()
        data = ra._data[fbool]
        features = ra.features.to_unique(fbool)
        positions = ra.positions.to_unique(fbool)
        ua = seqc.arrays.UniqueReadArray(data, features, positions)
        self.assertTrue(ua.shape[0] == n_unique, 'Incorrect number of unique reads: '
                                                 '%d != %d' % (ua.shape[0], n_unique))

        # check if other filters are causing it to fail
        ua2 = ra.to_unique(3)
        self.assertTrue(ua2.shape == ua.shape, 'Incongruent number of unique reads: '
                                               '%d != %d' % (ua2.shape[0], ua.shape[0]))

        self.assertTrue(ua2.shape[0] == 21984)  # double check that the number is right


if __name__ == '__main__':
    nose2.main()
