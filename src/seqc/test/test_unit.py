__author__ = "Ambrose J. Carr"

from copy import deepcopy
import nose2
import unittest
import os
import pickle
import seqc
import numpy as np
import xml.dom.minidom
import random
import shutil
import re
from nose2.tools import params
from more_itertools import first
from itertools import islice, permutations
from functools import partial
from operator import attrgetter
from numpy.lib.recfunctions import append_fields
import pandas as pd
from io import StringIO
from collections import defaultdict

# todo future test list:
# tests for downloading s3 inputs: barcodes, index, fastq
# tests for downloading from basespace
# test that error on remote instance should result in an email
# test sshutils
# test cluster_utils
# test sam reader (multialignments)


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
    samfile_pattern = seqc_dir + 'test_data/%s/sam/alignments.sam'
    forward_pattern = seqc_dir + 'test_data/%s/fastq/test_seqc_r1.fastq'
    reverse_pattern = seqc_dir + 'test_data/%s/fastq/test_seqc_r2.fastq'
    merged_pattern = seqc_dir + 'test_data/%s/fastq/merged.fastq'
    barcode_serial_pattern = seqc_dir + 'test_data/%s/barcodes/barcodes.p'
    barcode_partial_serial_pattern = seqc_dir + 'test_data/%s/barcodes/cb_partial.p'
    barcode_prefix_pattern = seqc_dir + 'test_data/%s/barcodes/'
    barcode_serialized_link_pattern = 's3://dplab-data/barcodes/%s/serial/barcodes.p'
    barcode_partial_serial_link_pattern = ('s3://dplab-data/barcodes/%s/serial/'
                                           'cb_partial.p')
    barcode_files_link_prefix_pattern = 's3://dplab-data/barcodes/%s/flat/'
    h5_name_pattern = seqc_dir + 'test_data/%s/h5/test_seqc.h5'
    gene_map_pattern = seqc_dir + 'test_data/%s/h5/test_seqc_gene_id_map.p'
    experiment_pattern = seqc_dir + 'test_data/%s/results/test_seqc.npz'

    # universal index files
    gtf = seqc_dir + 'test_data/genome/annotations.gtf'
    fasta = seqc_dir + 'test_data/genome/mm38_chr19.fa'
    index = seqc_dir + 'test_data/genome/'
    index_link = 's3://dplab-data/genome/mm38_chr19/'

    # config parameters
    n_threads = 7

    # s3 test file download links
    s3_forward_fastq_pattern = 's3://dplab-data/seqc/test_seqc/%s/fastq/forward/'
    s3_reverse_fastq_pattern = 's3://dplab-data/seqc/test_seqc/%s/fastq/reverse/'
    s3_merged_fastq_pattern = 's3://dplab-data/seqc/test_seqc/%s/fastq/merged/merged.fastq'
    s3_sam_pattern = 's3://dplab-data/seqc/test_seqc/%s/sam/alignments.sam'


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


def check_barcodes(data_type: str):
    """ensure that all necessary barcodes are present and downloaded"""

    # make directory if not present
    barcode_dir = config.barcode_prefix_pattern % data_type
    if not os.path.isdir(barcode_dir):
        os.makedirs(barcode_dir)

    if data_type == 'drop_seq':
        return  # drop-seq has no barcodes

    # download flat files if not present
    if not any(f.endswith('.txt') and 'partial' not in f
               for f in os.listdir(barcode_dir)):
        bucket, key_prefix = seqc.io.S3.split_link(
            config.barcode_files_link_prefix_pattern % data_type)
        seqc.io.S3.download_files(bucket, key_prefix, output_prefix=barcode_dir,
                                  cut_dirs=3)

    # check for serialized barcode file; download if missing.
    serial_barcodes = config.barcode_serial_pattern % data_type
    if not os.path.isfile(serial_barcodes):
        bucket, key = seqc.io.S3.split_link(
            config.barcode_serialized_link_pattern % data_type)
        seqc.io.S3.download_file(bucket, key, fout=serial_barcodes)

    # get all downloaded flat text files; create partial versions for rapid testing
    n_barcodes = 10
    flat_files = [barcode_dir + f for f in os.listdir(barcode_dir) if f.endswith('.txt')
                  and 'partial' not in f]
    for i, fin in enumerate(flat_files):
        with open(fin) as barcode_file:
            data = islice(barcode_file, n_barcodes)
            with open(barcode_dir + 'cb%d_partial.txt' % (i + 1), 'w') as fout:
                for line in data:
                    fout.write(line)

    # create a serialized partial file
    partial_files = [barcode_dir + f for f in os.listdir(barcode_dir) if
                     f.endswith('_partial.txt')]
    cb = seqc.barcodes.CellBarcodes.from_files(*partial_files)
    cb.pickle(barcode_dir + 'cb_partial.p')
    assert (config.barcode_partial_serial_pattern % data_type == barcode_dir +
            'cb_partial.p')


def check_fastq(data_type: str) -> (str, str):
    """check if the required fastq files are present, else generate them
    """

    # replace any dashes with underscores
    data_type = data_type.replace('-', '_')

    check_barcodes(data_type)

    # get forward and reverse file ids
    forward = config.forward_pattern % data_type
    reverse = config.reverse_pattern % data_type
    if data_type == 'drop_seq':
        barcodes = seqc.barcodes.DropSeqCellBarcodes()
    else:
        barcodes = config.barcode_partial_serial_pattern % data_type

    prefix = config.seqc_dir + 'test_data/%s/fastq/test_seqc' % data_type
    # if not os.path.isdir(prefix):
    #     os.makedirs(prefix)

    if not all(os.path.isfile(f) for f in [forward, reverse]):
        gen_func = getattr(seqc.fastq.GenerateFastq, data_type)
        loc, unloc, unalign = 3000, 3000, 3000
        seqlen, fraglen, rep = 100, 1000, 3
        gen_func(loc, unloc, unalign, prefix, config.fasta, config.gtf, config.index,
                 replicates=rep, n_threads=config.n_threads, fragment_length=fraglen,
                 sequence_length=seqlen, barcodes=barcodes)

    return forward, reverse


def check_merged_fastq(data_type: str) -> str:
    """check if the required merged fastq files are present, else generate them"""
    forward, reverse = check_fastq(data_type)

    if not os.path.isfile(config.merged_pattern % data_type):
        forward, reverse = [forward], [reverse]
        exp_type = data_type.replace('_', '-')
        n_processors = 4
        output_dir = config.seqc_dir + 'test_data/%s/fastq' % data_type
        if not os.path.isdir(output_dir):
            os.makedirs(output_dir)
        if data_type == 'drop_seq':
            cb = seqc.barcodes.DropSeqCellBarcodes()
        else:
            cb = seqc.barcodes.CellBarcodes.from_pickle(
                config.barcode_partial_serial_pattern % data_type)

        # create merged file
        merged = seqc.fastq.merge_fastq(
            forward, reverse, exp_type, output_dir, cb, n_processors)
        msg = ('output_file "%s" does not match the expected pattern: "%s"' %
               (merged, config.merged_pattern % data_type))
        assert merged == config.merged_pattern % data_type, msg

    return config.merged_pattern % data_type


def check_sam(data_type: str) -> str:
    """check if the required sam files are present, else generate them"""

    # replace any dashes with underscores
    data_type = data_type.replace('-', '_')

    merged = check_merged_fastq(data_type)

    # generation params
    samfile = config.samfile_pattern % data_type
    sam_dir = '/' + '/'.join(samfile.strip('/').split('/')[:-1])
    if not os.path.isdir(sam_dir):
        os.makedirs(sam_dir)

    if not os.path.isfile(samfile):
        gen_func = seqc.sam.GenerateSam.from_merged_fastq
        gen_func(samfile, merged, config.index, n_threads=config.n_threads)

    assert os.path.isfile(samfile)

    return samfile


def check_h5(data_type: str) -> str:

    samfile = check_sam(data_type)

    expected_h5 = config.h5_name_pattern % data_type
    h5_dir = '/'.join(expected_h5.split('/')[:-1])
    if not os.path.isdir(h5_dir):
        os.makedirs(h5_dir)

    if not os.path.isfile(config.h5_name_pattern % data_type):
        h5 = seqc.sam.to_h5(samfile, config.h5_name_pattern % data_type, 4, int(1e7),
                            config.gtf)
        assert os.path.isfile(config.h5_name_pattern % data_type)
        assert h5 == config.h5_name_pattern % data_type
        assert os.path.isfile(config.gene_map_pattern % data_type)
    else:
        h5 = config.h5_name_pattern % data_type

    return h5


def check_experiment(data_type: str) -> str:

    h5_file = check_h5(data_type)

    experiment_file = config.experiment_pattern % data_type
    exp_dir = '/'.join(experiment_file.split('/')[:-1])
    if not os.path.isdir(exp_dir):
        os.makedirs(exp_dir)

    if not os.path.isfile(experiment_file):
        ra = seqc.ReadArray.from_h5(h5_file)
        if data_type == 'drop_seq':
            exp = seqc.Experiment.from_read_array(ra, 0, 2)
        else:
            exp = seqc.Experiment.from_read_array(ra, 3, 2)
        exp.to_npz(experiment_file)

    return experiment_file


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


class FastqGenerateTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.localized = 1000
        cls.unlocalized = 1000
        cls.unalignable = 1000
        cls.sequence_length = 100
        cls.fragment_length = 1000
        cls.n_threads=7

    def test_reverse_three_prime(self):
        fasta = config.fasta
        gtf = config.gtf
        res = seqc.fastq.GenerateFastq._reverse_three_prime(
            5, read_length=30, fasta=fasta, gtf=gtf)
        self.assertTrue(True) # this test only tests that the data runs.

    def test_generate_reverse_complex(self):
        """generate a set of reads, only some of which align to expected locations"""

        res = seqc.fastq.GenerateFastq._reverse_split(
            self.localized, self.unlocalized, self.unalignable, config.index,
            config.fasta, config.gtf, self.fragment_length,
            self.sequence_length, self.n_threads)

    @params(*config.data_types)
    def test_generate_typed_fastq(self, data_type):
        """"""
        if data_type == 'drop_seq':
            barcodes = seqc.barcodes.DropSeqCellBarcodes()
        else:
            barcodes = config.barcode_partial_serial_pattern % data_type

        fgen = getattr(seqc.fastq.GenerateFastq, data_type)
        res = fgen(
            self.localized, self.unlocalized, self.unalignable, './', config.fasta,
            config.gtf, config.index, replicates=3, n_threads=self.n_threads,
            fragment_length=self.fragment_length, sequence_length=self.sequence_length,
            barcodes=barcodes)


class FastqMergeTest(unittest.TestCase):
    test_dir = 'test_seqc_fastq_merge/'

    @classmethod
    def setUpClass(cls):
        for t in config.data_types:
            check_barcodes(t)
            check_fastq(t)
        if not os.path.isdir(cls.test_dir):
            os.makedirs(cls.test_dir)

    # @unittest.skip('~6s; multiprocessing has some basic time-lag built in')
    @params(*config.data_types)
    def test_merge_fastq(self, data_type):
        # merge fastq should take lists of files and should work on both single and
        # multiple length lists
        forward = [config.forward_pattern % data_type] * 2
        reverse = [config.reverse_pattern % data_type]  # intentionally wrong input
        exp_type = data_type.replace('_', '-')
        output_dir = self.test_dir
        if data_type == 'drop_seq':
            cb = seqc.barcodes.DropSeqCellBarcodes()
        else:
            cb = seqc.barcodes.CellBarcodes.from_pickle(
                config.barcode_partial_serial_pattern % data_type)
        n_processes = 4

        # If input lists are different sizes, merge_fastq should raise a ValueError
        self.assertRaises(ValueError, seqc.fastq.merge_fastq, forward, reverse,
                          exp_type, output_dir, cb, n_processes)

        # fix input and create merged file
        reverse *= 2

        merged = seqc.fastq.merge_fastq(forward, reverse, exp_type, output_dir, cb,
                                        n_processes)

        # should return a string file id that maps to a file object
        self.assertIsInstance(merged, str)

        # read the merged file and run some additional tests:
        with open(merged) as f:
            data = f.readlines()
        merged_size = len(data)

        # get the size of the original file:
        with open(forward[0]) as f:
            original_size = 2 * sum(1 for _ in f)

        # no data is lost, original and merged files should have same amount of data
        self.assertEqual(original_size, merged_size)


        # merged file should be a valid fastq object (final character == '\n')
        self.assertEqual(data[-1][-1], '\n')

        # check that records have the right number of fields
        name = first(seqc.fastq.iter_records(data))[0]
        self.assertEqual(len(name.split(':')), 7)  # this is not a great test

        # check integrity of merged fastq
        n = int(len(data) / 4)
        cells = np.zeros(n, dtype=np.uint64)
        rmts = np.zeros(n, dtype=np.uint32)
        n_poly_t = np.zeros(n, dtype=np.uint16)
        valid_cells = np.zeros(n, dtype=np.bool)

        for i, record in enumerate(seqc.fastq.iter_records(data)):
            cell, rmt, npt, valid_cell, *unused = record[0][1:].split(':')
            cells[i] = int(cell)
            rmts[i] = int(rmt)
            n_poly_t[i] = int(npt)
            valid_cells[i] = int(valid_cell)

        # test merged fastq integrity
        self.assertTrue(np.all(cells))
        self.assertTrue(np.all(rmts))
        if data_type == 'in_drop':
            self.assertTrue(np.all(n_poly_t >= 2))
        elif data_type == 'drop_seq':
            self.assertTrue(~np.all(n_poly_t))
        self.assertTrue(np.all(valid_cells))

    @classmethod
    def tearDownClass(cls):
        if os.path.isdir(cls.test_dir):
            shutil.rmtree(cls.test_dir)


class IoS3ClientTest(unittest.TestCase):
    def setUp(self):
        # create a file
        self.bucket = 'dplab-data'
        self.testfile_key = 'genomes/mm38/chrStart.txt'
        self.testfile_download_name = 'test_s3_download.txt'
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

    def test_wrong_input_type_raises_type_error(self):
        self.assertRaises(TypeError, seqc.io.S3.download_file, self.bucket, 10,
                          self.testfile_download_name)
        self.assertRaises(TypeError, seqc.io.S3.download_file, self.bucket,
                          self.testfile_key, StringIO('1010'))

    def test_download_incorrect_filepath_raises_file_not_found_error(self):
        self.assertRaises(FileNotFoundError, seqc.io.S3.download_file, self.bucket,
                          'foobar', self.testfile_download_name, overwrite=False)

    @unittest.skip('slow')
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

    @unittest.skip('slow')
    def test_aws_s3_multiple_upload(self):
        # ResourceWarnings are generated by an interaction of io.S3.list() and
        # the unittesting suite. These should not occur in normal usage.
        bucket = 'dplab-home'
        key_prefix = 'ajc2205/testing_aws/'
        file_prefix = '.test_aws_s3/*'

        # upload file1, file2, and file3
        seqc.io.S3.upload_files(file_prefix, bucket, key_prefix)

        aws_dir = set(seqc.io.S3.listdir(bucket, key_prefix))
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

    def test_reader_iter_genes_returns_all_genes(self):
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
                last_interval = (
                    start, end if end - start < 500 else start + 500, repr(exons))
            else:
                raise ValueError

    def test_reader_iter_genes_final_n_returns_correct_number_genes(self):
        rd = seqc.gtf.Reader(self.gtf)
        self.assertEqual(len(list(rd.iter_genes_final_nbases(1000))), 4)

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

    @classmethod
    def tearDownClass(cls):
        if os.path.isdir(cls.test_dir):
            shutil.rmtree(cls.test_dir)


class GTFTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        check_index()
        cls.gtf = config.gtf
        cls._gene_id_pattern = b'(.*?gene_id ")(\D+)(\d+)(\.?)(\d*)(.*)'
        cls.fasta = config.fasta

    def test_reader_init(self):
        seqc.gtf_new.Reader(self.gtf)

    def test_reader_iter(self):
        rd = seqc.gtf_new.Reader(self.gtf)
        first(rd)  # check iterable

    def test_reader_iter_genes(self):
        rd = seqc.gtf_new.Reader(self.gtf)
        for i in rd.iter_gene_sets():
            break

        # make sure the gene_id is identical for each value
        gene_ids = set()
        for r in i:
            try:
                gene_ids.add(re.match(self._gene_id_pattern, r[-1]).groups()[2])
            except AttributeError:
                print(r)
                raise
        self.assertEqual(len(gene_ids), 1)

    def test_exon(self):
        rd = seqc.gtf_new.Reader(self.gtf)

        gene_records = first(rd.iter_gene_sets())

        # get an exon from this gene record
        for record in gene_records:
            if record[2] != b'exon':
                continue
            exon = seqc.gtf_new.Exon(record)
            break

        # test that all the properties are working
        self.assertIsInstance(exon.seqname, bytes)
        self.assertIsInstance(exon.chromosome, bytes)
        self.assertIsInstance(exon.source, bytes)
        self.assertIsInstance(exon.feature, bytes)
        self.assertIsInstance(exon.start, int)
        self.assertIsInstance(exon.end, int)
        self.assertIsInstance(exon.score, bytes)
        self.assertIn(exon.strand, [b'+', b'-'])
        self.assertIsInstance(exon.frame, bytes)

        # test that attribute parsing works
        try:
            self.assertIsInstance(exon.attribute(b'exon_id'), bytes)
        except:
            print(exon._fields[-1])
            raise

        # test repr & str
        self.assertIsInstance(repr(exon), str)
        self.assertIsInstance(str(exon), str)

        # test slicing
        plus_exon = deepcopy(exon)
        plus_exon._fields[6] = b'+'
        minus_exon = deepcopy(exon)
        minus_exon._fields[6] = b'-'
        size = exon.end - exon.start

        # one of these tests is yielding incorrect output!
        self.assertEqual(plus_exon[:round(size * .25)].end,
                         plus_exon.start + round(size * .25))
        self.assertEqual(plus_exon[-round(size * .25):].start,
                         plus_exon.end - round(size * .25))
        self.assertEqual(plus_exon[10:50].end, plus_exon.start + 50)
        self.assertEqual(plus_exon[10:50].start, plus_exon.start + 10)

        self.assertEqual(minus_exon[:round(size * .25)].end, minus_exon.end)
        self.assertEqual(minus_exon[:round(size * .25)].start,
                         minus_exon.end - round(size * .25))
        self.assertEqual(minus_exon[-round(size * .25):].start, minus_exon.start)
        self.assertEqual(minus_exon[-round(size * .25):].end, minus_exon.start +
                         round(size * .25))
        self.assertEqual(minus_exon[10:50].end, minus_exon.end - 10)
        self.assertEqual(minus_exon[10:50].start, minus_exon.end - 50)

    def test_transcript(self):
        rd = seqc.gtf_new.Reader(self.gtf)
        gene_records = first(rd.iter_gene_sets())

        iter_gene = iter(gene_records)

        # get the transcript record
        transcript_record = next(iter_gene)
        while transcript_record[2] != b'transcript':
            transcript_record = next(iter_gene)

        # now get all the UTRs
        exons = []
        record = next(iter_gene)  # get the first non-transcript record
        while record[2] not in [b'transcript', b'gene']:
            if record[2] in [b'exon', b'UTR']:
                exons.append(record)
            record = next(iter_gene)

        transcript = seqc.gtf_new.Transcript.from_records(transcript_record, *exons)

        # create a minus and plus transcript for testing
        if transcript._fields[6] == b'+':
            plus_transcript = transcript

            minus_transcript = deepcopy(transcript)
            minus_transcript._fields[6] = b'-'
            minus_transcript._exons = sorted(minus_transcript._exons,
                                             key=attrgetter('start'), reverse=True)
            for exon in minus_transcript:
                exon._fields[6] = b'-'
        else:
            minus_transcript = transcript
            plus_transcript = deepcopy(transcript)
            plus_transcript._fields[6] = b'+'
            plus_transcript._exons = sorted(plus_transcript._exons,
                                            key=attrgetter('start'))
            for exon in plus_transcript:
                exon._fields[6] = b'+'

        # test that all the properties are working
        self.assertIsInstance(transcript.seqname, bytes)
        self.assertIsInstance(transcript.chromosome, bytes)
        self.assertIsInstance(transcript.source, bytes)
        self.assertIsInstance(transcript.feature, bytes)
        self.assertIsInstance(transcript.start, int)
        self.assertIsInstance(transcript.end, int)
        self.assertIsInstance(transcript.score, bytes)
        self.assertIn(transcript.strand, [b'+', b'-'])
        self.assertIsInstance(transcript.frame, bytes)

        # test transcript slicing
        self.assertEqual(plus_transcript.size, minus_transcript.size)
        subset = round(plus_transcript.size * .25)

        self.assertEqual(plus_transcript[:subset].size, subset)
        self.assertEqual(plus_transcript[-subset:].size, subset)
        self.assertEqual(plus_transcript.start, plus_transcript[:subset].start)
        self.assertEqual(plus_transcript[-subset:].end, plus_transcript.end)

        self.assertEqual(minus_transcript[-subset:].size, subset)
        self.assertEqual(minus_transcript[:subset].size, subset)
        self.assertEqual(minus_transcript.end, minus_transcript[:subset].end)
        self.assertEqual(minus_transcript.start, minus_transcript[-subset:].start)

        self.assertEqual(plus_transcript[:subset].size, minus_transcript[:subset].size)
        self.assertEqual(plus_transcript[:-subset].size, minus_transcript[:-subset].size)
        self.assertEqual(plus_transcript[-subset:].size, minus_transcript[-subset:].size)
        self.assertEqual(plus_transcript[subset:].size, minus_transcript[subset:].size)
        self.assertEqual(plus_transcript[1:5].size, 4)
        self.assertEqual(minus_transcript[1:5].size, 4)
        self.assertEqual(plus_transcript[1:5].start, plus_transcript.start + 1)
        self.assertEqual(plus_transcript[1:5].end, plus_transcript.start + 5)
        self.assertEqual(minus_transcript[1:5].end, minus_transcript.end - 1)
        self.assertEqual(minus_transcript[1:5].start, minus_transcript.end - 5)

    def test_gene(self):
        rd = seqc.gtf_new.Reader(self.gtf)
        minus_gene = None
        plus_gene = None
        gene_iterator = rd.iter_gene_sets()
        while not all((minus_gene, plus_gene)):
            gene_records = list(next(gene_iterator))
            gene = seqc.gtf_new.Gene.from_records(*gene_records)
            if len(gene) > 1 and gene.strand == b'+':
                plus_gene = gene
            if len(gene) > 1 and gene.strand == b'-':
                minus_gene = gene

        plus_subset = int(np.mean([tx.size * .75 for tx in plus_gene]))
        minus_subset = int(np.mean([tx.size * .75 for tx in minus_gene]))

        for gene in (plus_gene, minus_gene):
            subset = int(np.mean([tx.size * .75 for tx in gene]))
            sliced_gene = gene[-subset:]
            self.assertTrue(all(tx.size == min(subset, original.size) for
                                tx, original in zip(sliced_gene, gene)))

        # test merge intervals
        res = list(plus_gene._merge_intervals([(1, 5), (4, 9)]))
        self.assertEqual(len(res), 1)
        self.assertEqual(res[0], (1, 9))

        res = list(plus_gene._merge_intervals([(1, 5), (4, 9), (-1, 3)]))
        self.assertEqual(len(res), 1)
        self.assertEqual(res[0], (-1, 9))

        res = list(plus_gene._merge_intervals([(4, 9), (-1, 3)]))
        self.assertEqual(len(res), 2)
        self.assertEqual(res, [(-1, 3), (4, 9)])

        # test intervals
        for gene in (plus_gene, minus_gene):
            subset = int(np.mean([tx.size * .75 for tx in gene]))
            res = gene.intervals(-1000)

    def test_gene_slicing(self):

        anno = seqc.gtf_new.Annotation(self.gtf, self.fasta)

        # test slicing
        sliced = [g[-1000:] for g in anno.genes]
        self.assertEqual(sum(1 for s in sliced if s.start and s.end), len(anno.genes))

        # test intervals
        ivs = [g.intervals(-1000, None) for g in anno.genes]
        self.assertTrue(all(len(iv) != 0 for iv in ivs))


    def test_annotation(self):

        anno = seqc.gtf_new.Annotation(self.gtf, self.fasta)

        # test properties
        self.assertIsInstance(len(anno), int)
        self.assertIsInstance(repr(anno), str)
        self.assertIsInstance(anno.chromosomes, dict)

        # test __getitem__
        sample_gene = anno.genes[0]
        chromosome = sample_gene.chromosome
        gene_id = sample_gene.string_gene_id
        self.assertIsInstance(anno[chromosome][gene_id], seqc.gtf_new.Gene)

        prefix = sample_gene.organism_prefix
        int_id = sample_gene.integer_gene_id
        string_id = sample_gene.string_gene_id.split(b'.')[0]
        self.assertEqual(string_id, seqc.gtf_new.Record.int2str_gene_id(int_id, prefix))

    def test_random_sequences(self):

        anno = seqc.gtf_new.Annotation(self.gtf, self.fasta)

        n_seqs = 1000
        seqlen = 98

        # check that we can generate sequences inside expected region
        seqs, genes = anno.random_sequences(n_seqs, seqlen, True, -1000, None)
        self.assertEqual(len(seqs), n_seqs)
        self.assertTrue(all(len(s) == seqlen for s in seqs))
        self.assertEqual(len(genes), n_seqs)

        # check that we can generate sequences outside of expected region
        seqs, genes = anno.random_sequences(n_seqs, seqlen, True, None, -1000)
        self.assertEqual(len(seqs), n_seqs)
        self.assertTrue(all(len(s) == seqlen for s in seqs))
        self.assertEqual(len(genes), n_seqs)

        # finally, want to be able to generate giberish sequences; this is in align


class TestFastaReader(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.fasta = config.fasta

    def test_fasta_reader(self):
        rd = seqc.fasta.Reader(self.fasta)
        for name, desc, seq in rd:
            pass
            # print(name)
            # print(desc)
            # print(seq[:20] + b'...')


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


class CoreCheckIndexTest(unittest.TestCase):
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

    output_dir = 'test_seqc/'

    def test_load_barcodes_wrong_input_type_raises(self):
        self.assertRaises(TypeError, seqc.core.check_barcodes,
                          10, self.output_dir)

    @params(*config.data_types)
    @unittest.skip('Downloads from s3, takes time, bandwidth')
    def test_load_barcodes_from_aws_link(self, dtype):
        barcode_link = config.barcode_serialized_link_pattern % dtype
        self.assertFalse(os.path.isdir('test_seqc/'))
        barcodes = seqc.core.check_barcodes(barcode_link, self.output_dir)
        self.assertIsInstance(barcodes, seqc.barcodes.CellBarcodes)

    @params(*config.data_types)
    @unittest.skip('This function takes a long time to load the serialized object. Think '
                   'about improving its run time.')
    def test_load_barcodes_local(self, dtype):
        if dtype == 'drop_seq':
            return  # None should be passed, not a link that does not target a file obj.

        barcode_file = config.barcode_serial_pattern % dtype
        barcodes = seqc.core.check_barcodes(barcode_file)
        barcodes = self.assertIsInstance(barcodes, seqc.barcodes.CellBarcodes)

    def test_load_pickle_containing_non_barcode_data(self):
        pass  # todo implement

    @classmethod
    def tearDownClass(cls):
        if os.path.isdir(cls.output_dir):
            shutil.rmtree(cls.output_dir)


class ConvertFeaturesConvertGeneCoordinatesTest(unittest.TestCase):
    test_dir = 'test_seqc_convert_features/'
    gtf = 'test_seqc_convert_features/test_data.gtf'

    @classmethod
    def setUpClass(cls):

        # create the test_dir if it doesn't exist
        if not os.path.isdir(cls.test_dir):
            os.mkdir(cls.test_dir)

        check_index()

        # create a very small gtf file for testing
        i = 0
        with open(config.gtf, 'rb') as fin:
            with open(cls.gtf, 'wb') as fout:
                while i < 203:
                    fout.write(fin.readline())
                    i += 1

    @classmethod
    def tearDownClass(cls):
        if os.path.isdir(cls.test_dir):
            shutil.rmtree(cls.test_dir)

    def test_convert_gene_features_empty_input_raises(self):
        self.assertRaises(ValueError, seqc.convert_features.ConvertGeneCoordinates,
                          {}, {}, {})

    def test_convert_gene_features_wrong_input_type_raises(self):
        self.assertRaises(TypeError, seqc.convert_features.ConvertGeneCoordinates,
                          {'chr1': 'this_should_be_intervaltree_not_string'}, {})

    def test_convert_gene_features_from_gtf_reads_all_genes(self):
        cgc = seqc.convert_features.ConvertGeneCoordinates.from_gtf(self.gtf)

        # test_data that the method returns the correct object type
        self.assertIsInstance(cgc, seqc.convert_features.ConvertGeneCoordinates)

        # make sure the interval tree is properly constructed; should have 4 genes
        genes = set()
        for key in cgc._data:
            for iv in cgc._data[key]:
                genes.add(iv.data)
        self.assertEqual(4, len(genes))

    def test_convert_gene_features_translate_returns_non_empty_list(self):
        sm = seqc.gtf.Sample(self.gtf)
        vals = sm.sample_final_n(10)
        self.assertTrue(all(isinstance(v, tuple) for v in vals))
        cgc = seqc.convert_features.ConvertGeneCoordinates.from_gtf(self.gtf)

        # test that method always returns list
        translated = [cgc.translate(*v) for v in vals]
        self.assertTrue(all(isinstance(t, list) and t for t in translated))

        # test that returned id can be reversed by id_map conversion to a string id
        for t in translated:
            str_id = cgc.int2str_id(t)
            self.assertIsInstance(str_id, str)
            int_id = cgc.str2int_id(str_id)
            self.assertIsInstance(int_id, int)
            self.assertEqual(int_id, t)

        # test that method returns empty list for out-of-range intervals
        out_of_range = [tuple(('+', 'chr19', v)) for v in [0, -1, int(1e11)]]
        translated = [cgc.translate(*v) for v in out_of_range]
        self.assertTrue(all(isinstance(t, list) and not t for t in translated))

        # test that translate method raises for malformed input
        # this chromosome has no genes, thus won't be in dictionary; should return empty
        # list
        self.assertEqual([], cgc.translate('+', 'KI270751.1', 1000))


    def test_convert_gene_features_save_and_load(self):
        cgc = seqc.convert_features.ConvertGeneCoordinates.from_gtf(self.gtf)
        cgc.pickle(self.test_dir + 'pickled.p')
        cgc2 = seqc.convert_features.ConvertGeneCoordinates.from_pickle(
            self.test_dir + 'pickled.p')

        # check that they are the same type of object, holding the same data
        self.assertEqual(type(cgc), type(cgc2))
        self.assertEqual(cgc._data, cgc2._data)


class ThreeBitCellBarcodesTest(unittest.TestCase):
    test_dir = 'test_seqc/barcodes/'

    @classmethod
    def setUpClass(cls):

        if not os.path.isdir(cls.test_dir):
            os.makedirs(cls.test_dir)

    @params(*config.data_types)
    def test_create_cell_barcode_object(self, data_type):

        # create barcode files
        barcode_dir = config.barcode_prefix_pattern % data_type
        barcode_files = [barcode_dir + f for f in os.listdir(barcode_dir)
                         if f.endswith('_partial.txt')]
        cb = seqc.barcodes.CellBarcodes.from_files(*barcode_files)
        with open(self.test_dir + '%s.p' % data_type, 'wb') as f:
            pickle.dump(cb, f)

        # confirm loaded and dumped files contain same data
        with open(self.test_dir + '%s.p' % data_type, 'rb') as f:
            cb2 = pickle.load(f)
        self.assertEqual(cb.perfect_codes, cb2.perfect_codes)
        self.assertEqual(cb.error_codes, cb2.error_codes)

        # read all the barcodes
        barcode_data = []
        for f in barcode_files:
            with open(f) as fin:
                barcode_data.append([l.strip() for l in fin.readlines()])

        # call the appropriate processor
        tbp = seqc.encodings.ThreeBit.default_processors(data_type.replace('_', '-'))

        # get a random barcode from each file and an UMI of appropriate length
        barcodes = [random.choice(bcs) for bcs in barcode_data]
        umi = ''.join(np.random.choice(list('ACGT'), tbp._rmt_len))

        # convert to 3bit
        int_codes = [tbp.str2bin(bc) for bc in barcodes]
        cell = tbp.ints2int(int_codes)

        # check that this does the same thing as getting the cell fom _extract_cell
        # join barcodes by spacer if present
        if data_type == 'in_drop':
            str_code = tbp._w1.join(barcodes)
        else:
            str_code = ''.join(barcodes)

        seq = str_code + umi
        cell2, rmt, n_poly_t = tbp.process_forward_sequence(seq)
        self.assertEqual(cell, cell2)

        # check that this sequence is in the cell barcodes
        self.assertTrue(cb.perfect_match(cell), '%d is not a perfect match' % cell)
        self.assertTrue(cb.close_match(cell), '%d is not a close match' % cell)

    @classmethod
    def tearDownClass(cls):
        if os.path.isdir(cls.test_dir):
            shutil.rmtree(cls.test_dir)


class ThreeBitInDropTest(unittest.TestCase):
    def setUp(self):
        self.tb = seqc.encodings.ThreeBit.default_processors('in-drop')
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
        s = seqc.encodings.ThreeBit.str2bin(self.example_seqs[0])
        seq = seqc.encodings.ThreeBit.bin2str(s)
        self.assertTrue(seq == self.example_seqs[0])

    def test_bin2str_inverts_str2bin_with_ints2int(self):
        b1 = seqc.encodings.ThreeBit.str2bin(self.c1s[0])
        b2 = seqc.encodings.ThreeBit.str2bin(self.c1s[1])
        b12 = seqc.encodings.ThreeBit.ints2int([b1, b2])
        self.assertTrue(b12 == seqc.encodings.ThreeBit.str2bin(self.c1s[0] + self.c1s[1]))

    def test_str2bin_converts_forwards_not_backwards(self):
        s = seqc.encodings.ThreeBit.str2bin('ACGT')
        self.assertEqual(s, 0b100110101011)
        self.assertNotEqual(s, 0b011101110100)

    def test_ints2int_and_cells(self):
        s1 = 'AAAAAAAA'
        s2 = 'TTTTTTTT'
        seqs = [seqc.encodings.ThreeBit.str2bin(s) for s in [s1, s2]]
        c = seqc.encodings.ThreeBit.ints2int(seqs)
        cells = self.tb.split_cell(c)
        c1, c2 = [seqc.encodings.ThreeBit.bin2str(c) for c in cells]
        self.assertTrue(s1 == c1)
        self.assertTrue(s2 == c2)

    def test_extract_cell_spacer_extraction(self):

        results = [3446, 2478, 2869, 1894]
        for i, strseq in enumerate(self.example_seqs):
            # note this tracks code in extract cell; test_data will not reflect function
            # if extract cell is changed (not an ideal test_data)
            s = seqc.encodings.ThreeBit.str2bin(strseq)
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
            str_cell = seqc.encodings.ThreeBit.bin2str(cell)
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


class AlignSTARTest(unittest.TestCase):

    test_dir = 'test_seqc/alignment/%s'

    @classmethod
    def setUpClass(cls):
        for dtype in config.data_types:
            if not os.path.isdir(cls.test_dir % dtype):
                os.makedirs(cls.test_dir % dtype)
            check_merged_fastq(dtype)
        cls.n_threads = 7

    @params(*config.data_types)
    @unittest.skip('slow; genome loading takes time.')
    def test_alignment(self, dtype):

        # get merged file name and make sure it's a file
        merged = config.merged_pattern % dtype
        self.assertTrue(os.path.isfile(merged), '%s is not a file' % merged)

        samfile = seqc.align.STAR.align(merged, config.index, self.n_threads,
                                        self.test_dir % dtype)

        # check that file was created
        self.assertTrue(os.path.isfile(samfile))

        # create an iterator over the samefile
        samreader = seqc.sam.Reader(samfile)

        # check that all the reads aligned. There may be more than this number due to
        # multimapping
        with open(merged, 'rb') as f:
            n = int(sum(1 for line in f) / 4)

        aligned = 0
        for record in samreader:
            if not int(record.flag) & 4:
                aligned += 1
        self.assertGreaterEqual(aligned, n)

    def test_align_giberish(self):
        res = seqc.align.STAR.generate_unalignable_sequences(
            1000, 100, config.index, self.n_threads)
        # print(type(res))
        # print(len(res))
        # print(res[:5])

    @classmethod
    def tearDownClass(cls):
        if os.path.isdir('test_seqc/'):
            shutil.rmtree('test_seqc/')


class ParallelConstructH5Test(unittest.TestCase):

    test_dir = 'seqc_test/'
    h5_name = test_dir + 'test_seqc.h5'
    n_processes = 7

    @classmethod
    def setUpClass(cls):
        if not os.path.isdir(cls.test_dir):
            os.makedirs(cls.test_dir)
        for dtype in config.data_types:
            if not os.path.isfile(config.samfile_pattern % dtype):
                check_sam(dtype)
            assert os.path.isfile(config.samfile_pattern % dtype)

    @params(*config.data_types)
    def test_parallel_construct_h5(self, dtype):
        chunk_size = int(1e7)
        self.assertFalse(os.path.isfile(self.h5_name))
        samfile = config.samfile_pattern % dtype
        seqc.sam.to_h5(samfile, self.h5_name, self.n_processes, chunk_size,
                       config.gtf, fragment_length=1000)

        # check that the h5 file was successfully created.
        self.assertTrue(os.path.isfile(self.h5_name))

        # read the h5 file and check several characteristics:
        ra = seqc.arrays.ReadArray.from_h5(self.h5_name)

        # figure out the number of records in the input samfile
        fastq_reader = seqc.fastq.Reader(config.merged_pattern % dtype)
        n_records = len(fastq_reader)

        sam_reader = seqc.sam.Reader(samfile)
        n_multialignments = sum(1 for ma in sam_reader.iter_multialignments())
        n_alignments = len(sam_reader)

        # should have equal numbers of multialignments, records, and readarray records
        self.assertEqual(n_records, n_multialignments)
        self.assertEqual(n_multialignments, len(ra))

        # n_chunks = np.ceil(os.path.getsize(samfile) / chunk_size)

        # get expected minimum number of alignments if no reads multi-mapped
        # min_alignments = n_records - (n_chunks * 2) - 2

        # all reads should have a valid cell
        self.assertTrue(np.all(ra.data['valid_cell']))

        # all reads should have nonzero cell & rmt fields
        self.assertTrue(np.all(ra.data['rmt']))
        self.assertTrue(np.all(ra.data['cell']))

        # all reads should have average quality == 40 (due to how data was generated)
        self.assertTrue(np.all(ra.data['rev_quality'] == 40))
        self.assertTrue(np.all(ra.data['fwd_quality'] == 40))

        # todo:
        # all reads that are not ambiguous should have features; how many fail this?
        # we can check this by parsing the IntervalTree to see which overlaps have
        # multiple features. If they do, then they will not show up.
        n_with_features = sum(1 for f in ra.features if np.any(f))

        self.assertTrue(n_with_features > len(ra) * .98)

        # check data types
        self.assertEqual(ra.features.data.dtype, np.int64)

        # can parse the sam file to check that the right number of molecules are being
        # detected in the ReadArray.
        # count molecule
        # todo this is slow, refactor to single dictionary, expand out to df
        molecules = defaultdict(lambda: defaultdict(int))
        for ma in sam_reader.iter_multialignments():
            # record[0] contains the molecule/read information
            cell = ma[0].cell
            rmt = ma[0].rmt
            molecules[cell][rmt] += 1

        # get the same data from the ReadArray
        ra_molecules = defaultdict(lambda: defaultdict(int))
        data = np.hstack((ra.data['cell'][:, np.newaxis], ra.data['rmt'][:, np.newaxis]))
        for row in data:
            ra_molecules[row[0]][row[1]] += 1

        ra_df = pd.DataFrame(ra_molecules)
        sam_df = pd.DataFrame(molecules)

        # check that indices are equal
        self.assertEqual(sorted(ra_df.index), sorted(sam_df.index))
        self.assertEqual(sorted(ra_df.columns), sorted(sam_df.columns))

        # check total counts
        self.assertEqual(ra_df.sum().sum(), sam_df.sum().sum())

        # check column sums
        idx = sorted(ra_df.index)
        self.assertTrue(np.array_equal(ra_df.loc[idx].sum(axis=0),
                                       sam_df.loc[idx].sum(axis=0)))

        # check row sums
        clm = sorted(ra_df.columns)
        self.assertTrue(np.array_equal(ra_df[clm].sum(axis=1), sam_df[clm].sum(axis=1)))

    def tearDown(self):
        if os.path.isfile(self.h5_name):
            os.remove(self.h5_name)

    @classmethod
    def tearDownClass(cls):
        if os.path.isdir(cls.test_dir):
            shutil.rmtree(cls.test_dir)


class ConvertGeneCoordinatesTest(unittest.TestCase):

    test_dir = 'seqc_test/'
    h5_name = test_dir + 'test_seqc.h5'

    @classmethod
    def setUpClass(cls):
        """mimic the input of TestParallelConstructH5 exactly"""
        if not os.path.isdir(cls.test_dir):
            os.makedirs(cls.test_dir)
        for dtype in config.data_types:
            if not os.path.isfile(config.samfile_pattern % dtype):
                check_sam(dtype)
            assert os.path.isfile(config.samfile_pattern % dtype)


    @params(*config.data_types)
    def test_convert_gene_features(self, dtype):
        """
        the vast majority of the inputs are unique, therefore the vast majority should
        convert properly
        """

        # create the converter, allow the full gene to be used to create intervals!
        cgc = seqc.convert_features.ConvertGeneCoordinates.from_gtf(config.gtf, int(1e10))

        samfile = config.samfile_pattern % dtype
        rd = seqc.sam.Reader(samfile)
        fwd_res = []
        rev_res = []
        both = []
        total = 0
        res = []
        expected = []
        for records in rd.iter_multialignments():
            if len(records) == 1:
                if not int(records[0].flag) & 4:
                    strand = records[0].strand
                    chrom = records[0].rname
                    pos = int(records[0].pos)
                    try:
                        res.append(cgc.int2str_id(cgc.translate(strand, chrom, pos)[0]))
                    except IndexError:
                        res.append(None)
                    expected.append(records[0].qname.strip().split(':')[-1])

        total = len(res)
        n = 0
        for r, e in zip(res, expected):
            if r == e:
                n += 1
            elif r:
                pass
                # print(r, e)
        n_correct = sum(1 for r, e in zip(res, expected) if r == e)
        # print([(r, e) for (r, e) in zip(res[:15], expected[:15])])
        # print('percentage correct: %f' % (n_correct / len(res)))
        self.assertTrue(n_correct > 0.95, 'Of %d alignments, only %d converted' %
                        (total, n_correct))

    @classmethod
    def tearDownClass(cls):
        if os.path.isdir(cls.test_dir):
            shutil.rmtree(cls.test_dir)


class JaggedArrayTest(unittest.TestCase):

    # todo expand testing here.

    def generate_input_iterable(self, n):
        for i in range(n):
            yield [random.randint(0, 5) for _ in range(random.randint(0, 5))]

    def test_jagged_array(self):
        n = int(1e3)
        data = list(self.generate_input_iterable(n))
        data_size = 0
        for i in data:
            try:
                data_size += len(i)
            except TypeError:
                data_size += 1
        jarr = seqc.arrays.JaggedArray.from_iterable(data)
        self.assertTrue(jarr._data.dtype == np.int64)
        self.assertTrue(jarr._data.shape == (data_size,))


class UniqueArrayCreationTest(unittest.TestCase):
    """
    suspicion -> unique array creation is breaking something in the pipeline;
    could be an earlier method but lets check upstream so we know how each part is
    working.

    method: find the smallest real dataset that breaks it and use that to test_data.
    """

    test_dir = 'seqc_test/'
    h5_name = 'seqc_test/seqc_test.h5'

    @classmethod
    def setUpClass(cls):
        if not os.path.isdir(cls.test_dir):
            os.mkdir(cls.test_dir)

    @params(*config.data_types)
    # @unittest.skip('slow')
    def test_num_unique_samfile(self, dtype):
        h5 = seqc.sam.to_h5(config.samfile_pattern % dtype, self.h5_name, 7,
                            int(1e7), config.gtf, 1000)

        ra = seqc.arrays.ReadArray.from_h5(h5)

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
        # not yet tested explicitly, but implicitly working based on jagged.to_unique()
        fbool = ra.features.is_unique()
        data = ra._data[fbool]
        features = ra.features.to_unique(fbool)
        positions = ra.positions.to_unique(fbool)
        ua = seqc.arrays.UniqueReadArray(data, features, positions)
        self.assertTrue(ua.shape[0] == n_unique, 'Incorrect number of unique reads: '
                                                 '%d != %d' % (ua.shape[0], n_unique))

        # check if other filters are causing it to fail
        if dtype == 'drop_seq':
            ua2 = ra.to_unique(0)
        else:
            ua2 = ra.to_unique(3)
        self.assertTrue(ua2.shape == ua.shape, 'Incongruent number of unique reads: '
                                               '%d != %d' % (ua2.shape[0], ua.shape[0]))

        # figure out how many reads of each type should be in the sam file by parsing
        total_reads = 0
        unique_reads = 0
        unaligned_reads = 0
        multi_aligned_reads = 0
        rd = seqc.sam.Reader(config.samfile_pattern % dtype)
        for a in rd.iter_multialignments():
            if a.is_unmapped:
                unaligned_reads += 1
            elif a.is_uniquely_mapped:
                unique_reads += 1
            else:
                multi_aligned_reads += 1
            total_reads += 1
        self.assertEqual(unaligned_reads + unique_reads + multi_aligned_reads,
                         total_reads)

        # we lose a few genes (<5%) by requiring uniqueness
        self.assertEqual(len(ua), unique_reads)
        self.assertEqual(len(ra), total_reads)
        # self.assertTrue(len(ua) > len(ra) * .95, 'ua: %d << ra: %d' % (len(ua), len(ra)))

        # check that features and positions have the right type
        self.assertEqual(ua.features.dtype, np.int64)
        self.assertEqual(ua.positions.dtype, np.int64)

    @classmethod
    def tearDownClass(cls):
        if os.path.isdir(cls.test_dir):
            shutil.rmtree(cls.test_dir)


@unittest.skip('These tests are for to_experiment_new()')
class UniqueArrayToExperimentTest(unittest.TestCase):
    """Goal: understand if UniqueArray conversion is operating properly"""

    @classmethod
    def setUpClass(cls):
        cls.ua = {}
        for data_type in config.data_types:
            h5 = check_h5(data_type)
            ra = seqc.ReadArray.from_h5(h5)
            cls.ua[data_type] = ra.to_unique(0)
            cls.ua[data_type]._sort()

    @params(*config.data_types)
    def test_order_of_keys_does_not_matter_in_boundary_finding(self, data_type):

        order = self.ua[data_type]._sorted
        input_orders = list(permutations(['cell', 'rmt', 'features'], r=3))
        test_result = self.ua[data_type]._find_boundaries(order, *input_orders[0])
        for argument_ord in input_orders:
            res = self.ua[data_type]._find_boundaries(order, *argument_ord)
            self.assertTrue(np.array_equal(test_result, res))

    @params(*config.data_types)
    def test_find_boundaries(self, data_type):
        order = self.ua[data_type]._sorted
        molecule_boundaries = self.ua[data_type]._find_boundaries(
            order, 'cell', 'rmt', 'features')

        # lengths should be identical
        self.assertEqual(len(molecule_boundaries), len(order))

        # data_type should be a boolean
        self.assertEqual(molecule_boundaries.dtype, np.bool)

        # there should be as many ones as unique values
        boundaries = self.ua[data_type]._find_boundaries(order, 'cell')
        self.assertEqual(np.sum(boundaries.astype(int)),
                         len(np.unique(self.ua[data_type].data['cell'])))

        # make sure the wrong input raises a TypeError
        self.assertRaises(TypeError, self.ua[data_type]._find_boundaries, boundaries,
                          'cell')

    @params(*config.data_types)
    def test_get_counts(self, data_type):
        order = self.ua[data_type]._sorted
        molecule_boundaries = self.ua[data_type]._find_boundaries(
            order, 'cell', 'rmt', 'features')
        counts, index = self.ua[data_type]._get_counts(molecule_boundaries)

        # there should be fewer or equal numbers of counts than molecule boundaries,
        # but equal numbers of counts and indices
        self.assertLessEqual(len(counts), len(molecule_boundaries))
        self.assertEqual(len(counts), len(index))

        # there should be counts for as many as there are unique records
        # here's a simpler example
        boundaries = self.ua[data_type]._find_boundaries(order, 'cell')
        counts, index = self.ua[data_type]._get_counts(boundaries)
        self.assertEqual(len(counts), len(np.unique(self.ua[data_type].data['cell'])))

        # the sum of counts should be equal to the length of the input boundaries
        self.assertEqual(len(boundaries), np.sum(counts))

        # index should always be non-negative and less than the input length
        self.assertTrue(np.all((index >= 0) & (index < len(boundaries))))

    @params(*config.data_types)
    def test_filter_low_support_features(self, data_type):
        """this test is too weak; needs intermediate data with low support to test."""
        order = self.ua[data_type]._sorted
        molecule_boundaries = self.ua[data_type]._find_boundaries(
            order, 'cell', 'rmt', 'features')
        counts, index = self.ua[data_type]._get_counts(molecule_boundaries)

        supp_index, supp_counts = self.ua[data_type]._filter_low_support_features(
            index, counts, 1)
        self.assertEqual(len(supp_index), len(index))

        # when required_support is 1, all molecules should pass
        supported_idx, supported_counts = self.ua[data_type]._filter_low_support_features(
            index, counts, 1)
        self.assertTrue(np.array_equal(supported_idx, index))
        self.assertEqual(len(supported_counts), len(supported_idx))

        # when required_support is high, ValueError is raised if it causes all counts
        # to fail the filter
        self.assertRaises(ValueError, self.ua[data_type]._filter_low_support_features,
            index, counts, 10)

    @params(*config.data_types)
    def test_map_to_unique_index(self, data_type):

        # presents a set of indices which should multi-index cell or rmt, thus be
        # mappable down to a unique index
        order = self.ua[data_type]._sorted
        molecule_boundaries = self.ua[data_type]._find_boundaries(
            order, 'cell', 'rmt', 'features')
        counts, index = self.ua[data_type]._get_counts(molecule_boundaries)

        # test validity of input data; set of cells should be smaller than indexed cells
        cells = self.ua[data_type].data['cell']
        self.assertLess(len(set(cells)), len(cells))

        # test that np.unique generates the same list as set()
        id_set_np = np.unique(np.ravel(cells))
        id_set = set(cells)
        self.assertEqual(sorted(list(id_set_np)), sorted(list(id_set)))

        # since this is a map of ids to id_set, the returned new indices should be the
        # same shape as input ids, and the returned cells should be smaller than
        # the input cells (we required this above)
        index, cell_set = seqc.arrays.UniqueReadArray._map_to_unique_index(cells)
        self.assertEqual(len(index), len(cells))
        self.assertLess(len(cell_set), len(cells))

        # all new cells should be found in the original set, and all cells in the original
        # set should be found in the new cells.
        intersection = set(cells).intersection(cell_set)
        self.assertEqual(len(intersection), len(set(cells)))

        # returned index should be a vector of integers whose values are within the range
        # [0, len(cell_set))
        self.assertTrue(np.all((index >= 0) & (index < len(cell_set))))

    @params(*config.data_types)
    def test_reads_per_gene_per_cell_and_sparse_matrix_generation(self, data_type):
        """
        reads per cell is the simplest sparse matrix generation method, test together
        """
        order = self.ua[data_type]._sorted
        read_gene_cell_boundaries = self.ua[data_type]._find_boundaries(
            order, 'cell', 'features')
        read_counts, ua_gene_cell_positions = self.ua[data_type]._get_counts(
            read_gene_cell_boundaries)

        sm = self.ua[data_type]._generate_sparse_read_matrix(
            read_counts, ua_gene_cell_positions, order)

        # sm should be a SparseCounts object
        self.assertIsInstance(sm, seqc.SparseCounts)

        # sm should contain an equal number of counts as ua
        self.assertEqual(len(self.ua[data_type]), sm.counts.sum().sum())

        # the shape of sm should be equal to the set of all cells times all features
        all_cells = set(self.ua[data_type].data['cell'])
        all_features = set(self.ua[data_type].features)
        self.assertEqual(sm.counts.shape, (len(all_cells), len(all_features)))

        # when sorted, the indices should equal a sorted set of all cells and features in
        # ua
        self.assertEqual(sorted(all_cells), sorted(sm.index))
        self.assertEqual(sorted(all_features), sorted(sm.columns))

        # the row sums should equal the number of times a given cell was observed in ua
        cell_counts = pd.Series(np.zeros(len(all_cells)), index=all_cells)
        for c in self.ua[data_type].data['cell']:
            cell_counts[c] += 1

        # index of cell_counts and sm should be identical when sorted
        self.assertEqual(sorted(sm.index), sorted(cell_counts.index))

        # get cell counts from sm; should be equal shape to index
        sm_cell_counts = np.ravel(sm.counts.sum(axis=1))
        self.assertEqual(sm_cell_counts.shape, sm.index.shape)

        # check that cell sums are identical
        self.assertTrue(np.array_equal(sm_cell_counts, cell_counts[sm.index].values))

        # repeat the above process for genes
        gene_counts = pd.Series(np.zeros(len(all_features)), index=all_features)
        for f in self.ua[data_type].features:
            gene_counts[f] += 1

        # index of cell_counts and sm should be identical when sorted
        self.assertEqual(sorted(sm.columns), sorted(gene_counts.index))

        # get cell counts from sm; should be equal shape to index
        sm_gene_counts = np.ravel(sm.counts.sum(axis=0))
        self.assertEqual(sm_gene_counts.shape, sm.columns.shape)

        # check that cell sums are identical
        self.assertTrue(np.array_equal(sm_gene_counts, gene_counts[sm.columns].values))

    @params(*config.data_types)
    def test_molecules_per_cell_per_gene(self, data_type):

        ua = self.ua[data_type]
        order = ua._sorted

        # get number of records associated with each molecule
        molecule_boundaries = ua._find_boundaries(
            order, 'cell', 'rmt', 'features')
        molecule_counts, ua_molecule_positions = ua._get_counts(
            molecule_boundaries)

        # create molecules; allow filter to pass-through for now, and confirm that pass-
        # through is occuring
        res = ua._filter_low_support_features(
            ua_molecule_positions, molecule_counts, required_support=1)
        supported_molecule_positions, supported_molecule_counts = res
        self.assertTrue(np.array_equal(supported_molecule_positions,
                                       ua_molecule_positions))
        self.assertTrue(np.array_equal(supported_molecule_counts, molecule_counts))

        #################################################################################
        # This part doesn't matter for molecule generation, included in case bugs
        # get indices for the subset of unique molecules with counts > required_support

        # get reads per gene per cell by diffing on the original sort without considering
        # the RMT. This has the effect of counting the number of reads per gene, per cell.
        read_gene_cell_boundaries = ua._find_boundaries(
            order, 'cell', 'features')
        read_counts, ua_gene_cell_positions = ua._get_counts(
            read_gene_cell_boundaries)
        #################################################################################

        # if we call find_boundaries again with the SAME inputs (not the reduced one used
        # to find gene/cell boundaries), we should get a vector of all 1s. Note that we
        # need to pass the index created in find_counts() not find_boundaries() otherwise
        # the function gets the length wrong.
        ones = ua._find_boundaries(ua_molecule_positions, 'cell', 'rmt', 'features')
        self.assertEqual(sum(ones), len(ones))

        # get the gene/cell positions in molecule space
        molecule_gene_cell_boundaries = ua._find_boundaries(
            supported_molecule_positions, 'cell', 'features')

        # get the number of molecules per gene, per cell, and their indices
        molecule_gene_cell_counts, ua_molecule_gene_cell_positions = ua._get_counts(
            molecule_gene_cell_boundaries)

        # count the molecules in the original data; we should have the same number
        molecules = defaultdict(int)
        for i, row in enumerate(ua.data):
            key = (row[0], row[1], ua.features[i])
            molecules[key] += 1
        self.assertEqual(len(molecules), np.sum(molecule_gene_cell_counts))

        # check that the index is correctly identifying all cells
        cells = set(ua.data['cell'][ua_molecule_gene_cell_positions])
        all_cells = set(ua.data['cell'])
        self.assertEqual(cells, all_cells)

        # ... and all genes
        features = set(ua.features[ua_molecule_gene_cell_positions])
        all_features = set(ua.features)
        self.assertEqual(features, all_features)

        sm = ua._generate_sparse_molecule_matrix(
            molecule_gene_cell_counts, ua_molecule_gene_cell_positions)

    @params(*config.data_types)
    def test_molecules_sparse_matrix_generation(self, data_type):
        """
        test molecules per cell
        """

        # set a short-hand variable to the ua-appropriate data_type
        ua = self.ua[data_type]
        order = ua._sorted

        # get number of records associated with each molecule
        molecule_boundaries = ua._find_boundaries(
            order, 'cell', 'rmt', 'features')
        molecule_counts, ua_molecule_positions = ua._get_counts(
            molecule_boundaries)

        #################################################################################
        # This part doesn't matter for molecule generation, included in case bugs
        # get indices for the subset of unique molecules with counts > required_support

        # create molecules; allow filter to pass-through for now.
        res = ua._filter_low_support_features(
            ua_molecule_positions, molecule_counts, required_support=1)
        supported_molecule_positions, supported_molecule_counts = res

        # get reads per gene per cell by diffing on the original sort without considering
        # the RMT. This has the effect of counting the number of reads per gene, per cell.
        read_gene_cell_boundaries = ua._find_boundaries(
            order, 'cell', 'features')
        read_counts, ua_gene_cell_positions = ua._get_counts(
            read_gene_cell_boundaries)
        #################################################################################

        # molecules per cell can be calculated from reads per molecule by re-diffing
        # reads per molecule, without considering the rmt. This has the effect of
        # counting the number of RMTs per gene, per cell.
        molecule_gene_cell_boundaries = ua._find_boundaries(
            supported_molecule_positions, 'cell', 'features')
        molecule_gene_cell_counts, ua_molecule_gene_cell_positions = ua._get_counts(
            molecule_gene_cell_boundaries)

        sm = ua._generate_sparse_molecule_matrix(
            molecule_gene_cell_counts, ua_molecule_gene_cell_positions)

        #################################################################################
        # tests

        # sm should be a SparseCounts object
        self.assertIsInstance(sm, seqc.SparseCounts)

        # count molecules in ua; sm should contain this many molecules
        test_molecules = defaultdict(int)
        for i, row in enumerate(ua.data):
            key = (row[0], row[1], ua.features[i])
            test_molecules[key] += 1
        self.assertEqual(len(dict(test_molecules)), sm.counts.sum().sum())

        # the shape of the matrix should be equal to its indices and columns
        self.assertEqual(sm.counts.shape, (len(sm.index), len(sm.columns)))

        ################################ TEST GENES #####################################

        # the features in sm should be equal to the features in ua
        all_features = set(ua.features)
        self.assertEqual(sorted(sm.columns), sorted(all_features))

        # the row sums should equal the number of molecules observed for each cell
        gene_counts = pd.Series(np.zeros(len(all_features)), index=all_features)
        for cell, rmt, feature in test_molecules:
            gene_counts[feature] += 1

        # index of cell_counts and sm should be identical when sorted
        self.assertEqual(sorted(sm.columns), sorted(gene_counts.index))

        # convert to ndarray of correct dtype; check indices and types
        gene_counts = gene_counts[sm.columns].values.astype(np.int32)
        sm_gene_counts = np.ravel(sm.counts.sum(axis=0))
        self.assertEqual(sm_gene_counts.shape, sm.columns.shape)
        self.assertEqual(sm_gene_counts.dtype, gene_counts.dtype)

        # make sure the total number of molecules are all equal
        self.assertEqual(gene_counts.sum(), len(dict(test_molecules)))
        self.assertEqual(gene_counts.sum(), sm.counts.sum().sum())
        self.assertEqual(gene_counts.sum(), sm_gene_counts.sum())

        # print(np.sum(sm_gene_counts != gene_counts))
        res = sm_gene_counts - gene_counts
        # print(res[np.where(res)])

        # could the below problem be that the index isn't properly sorted? Nope!
        # self.assertTrue(np.array_equal(sorted(sm_gene_counts), sorted(gene_counts)))

        # todo gene_counts inconsistently 1 less than sm_gene_counts (per bin)
        self.assertTrue(np.array_equal(sm_gene_counts, gene_counts))

        ################################ TEST CELLS #####################################

        # the cells in sm should be equal to the cells in ua
        all_cells = set(ua.data['cell'])
        self.assertEqual(sorted(sm.index), sorted(all_cells))

        # the row sums should equal the number of molecules observed for each cell
        cell_counts = pd.Series(np.zeros(len(all_cells)), index=all_cells)
        for cell, rmt, feature in test_molecules:
            cell_counts[cell] += 1

        # index of cell_counts and sm should be identical when sorted
        self.assertEqual(sorted(sm.index), sorted(cell_counts.index))

        # convert to ndarray of correct dtype; check indices and types
        cell_counts = cell_counts[sm.index].values.astype(np.int32)
        sm_cell_counts = np.ravel(sm.counts.sum(axis=1))
        self.assertEqual(sm_cell_counts.shape, sm.index.shape)
        self.assertEqual(sm_cell_counts.dtype, cell_counts.dtype)

        # make sure the total number of molecules are all equal
        self.assertEqual(cell_counts.sum(), len(dict(test_molecules)))
        self.assertEqual(cell_counts.sum(), sm.counts.sum().sum())
        self.assertEqual(cell_counts.sum(), sm_cell_counts.sum())

        # print(np.sum(sm_cell_counts != cell_counts))
        res = sm_cell_counts - cell_counts
        # print(res[np.where(res)])

        # could the below problem be that the index isn't properly sorted? Nope!
        self.assertTrue(np.array_equal(sorted(sm_cell_counts), sorted(cell_counts)))

        # todo cell_counts inconsistently 1 less than sm_cell_counts (per bin)
        self.assertTrue(np.array_equal(sm_cell_counts, cell_counts))



        # repeat the above process for genes
        gene_counts = pd.Series(np.zeros(len(all_features)), index=all_features)
        for f in ua.features:
            gene_counts[f] += 1

        # index of gene_counts and sm should be identical when sorted
        self.assertEqual(sorted(sm.columns), sorted(gene_counts.index))

        # get cell counts from sm; should be equal shape to index
        sm_gene_counts = np.ravel(sm.counts.sum(axis=0))
        self.assertEqual(sm_gene_counts.shape, sm.columns.shape)

        # check that cell sums are identical
        self.assertTrue(np.array_equal(sm_gene_counts, gene_counts[sm.columns].values))

    # @params(*config.data_types)
    # def test_generate_experiment(self, data_type):
    #     pass

    # @params(*config.data_types)
    # @unittest.skip('commented out old code that was used to test this')
    # def test_to_experiment1_and_to_experiment2_give_same_result(self, data_type):
    #     """only works so long as old version of map_to_unique_index is used"""
    #     e1 = ua.to_experiment(2)
    #     e2 = ua.to_experiment2(2)
    #     # test that molecule indices are equal
    #     self.assertTrue(np.array_equal(e1.molecules.index, e2.molecules.index))
    #     self.assertTrue(np.array_equal(e1.molecules.columns, e2.molecules.columns))
    #
    #     # test that read indices are equal
    #     self.assertTrue(np.array_equal(e1.reads.index, e2.reads.index))
    #     self.assertTrue(np.array_equal(e1.reads.columns, e2.reads.columns))
    #
    #     # test that the same number of counts are being detected
    #     self.assertEqual(e1.molecules.counts.sum().sum(),
    #                      e2.molecules.counts.sum().sum())
    #     self.assertEqual(e1.reads.counts.sum().sum(),
    #                      e2.reads.counts.sum().sum())
    #
    #     # test that data is equal
    #     self.assertTrue(np.array_equal(e1.molecules.counts.todense(),
    #                                    e2.molecules.counts.todense()))
    #     self.assertTrue(np.array_equal(e1.reads.counts.todense(),
    #                                    e2.reads.counts.todense()))
    #
    # @params(*config.data_types)
    # @unittest.skip('commented out old code that was used to test this')
    # def test_map_to_unique_index2_and_manu_methods_give_same_result(self, data_type):
    #     """
    #     implements updated unique index mapping for both new (e4) and old (e3) methods
    #     """
    #     ids = np.random.randint(0, 10, (10000,))
    #     ac_indices, ac_unique = seqc.arrays.UniqueReadArray._map_to_unique_index2(ids)
    #     ms_indices, ms_unique = seqc.arrays.UniqueReadArray._map_to_unique_index_manu(ids)
    #
    #     self.assertTrue(np.array_equal(ac_indices, ms_indices))
    #     self.assertTrue(np.array_equal(ac_unique, ac_unique))
    #
    # @params(*config.data_types)
    # @unittest.skip('commented out old code that was used to test this')
    # def test_to_experiment_3_and_to_experiment4_give_same_result(self, data_type):
    #     """
    #     implements updated unique index mapping for both new (e4) and old (e3) methods
    #     """
    #     e1 = ua.to_experiment3(2)
    #     e2 = ua.to_experiment4(2)
    #     # test that molecule indices are equal
    #     self.assertTrue(np.array_equal(e1.molecules.index, e2.molecules.index))
    #     self.assertTrue(np.array_equal(e1.molecules.columns, e2.molecules.columns))
    #
    #     # test that read indices are equal
    #     self.assertTrue(np.array_equal(e1.reads.index, e2.reads.index))
    #     self.assertTrue(np.array_equal(e1.reads.columns, e2.reads.columns))
    #
    #     # test that the same number of counts are being detected
    #     self.assertEqual(e1.molecules.counts.sum().sum(),
    #                      e2.molecules.counts.sum().sum())
    #     self.assertEqual(e1.reads.counts.sum().sum(),
    #                      e2.reads.counts.sum().sum())
    #
    #     # test that data is equal
    #     self.assertTrue(np.array_equal(e1.molecules.counts.todense(),
    #                                    e2.molecules.counts.todense()))
    #     self.assertTrue(np.array_equal(e1.reads.counts.todense(),
    #                                    e2.reads.counts.todense()))


class CountingUniqueArrayTest(unittest.TestCase):

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
        self.simple_duplicate = create_unique_read_array(rmts * 2, cells * 2,
                                                         features * 2)

    def test_print_data(self):
        self.simple_unique._sort()
        data = append_fields(self.simple_unique.data, 'features',
                             self.simple_unique.features)

    def test_counting_filter_too_high(self):
        # seems to be sorting cell, feature, rmt
        self.assertRaises(ValueError, self.simple_unique.to_experiment, 1)

    def test_counting_simple(self):
        exp = self.simple_unique.to_experiment(0)
        rdata = np.array(exp.reads.counts.todense())
        mdata = np.array(exp.molecules.counts.todense())
        self.assertTrue(np.all(rdata == 3))
        self.assertTrue(np.all(mdata == 3))

    def test_counting_doublets(self):
        exp = self.simple_duplicate.to_experiment(0)
        rdata = np.array(exp.reads.counts.todense())
        mdata = np.array(exp.molecules.counts.todense())
        self.assertTrue(np.all(rdata == 6))
        self.assertTrue(np.all(mdata == 3))


@unittest.skip('These only work when create_experiment() -> create_experiment_old() and '
               'create_experiment_new() -> create_experiment()')
class ExperimentCreationTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        for dtype in config.data_types:
            check_h5(dtype)

    @params(*config.data_types)
    def test_create_experiment(self, dtype):
        h5 = config.h5_name_pattern % dtype
        ra = seqc.arrays.ReadArray.from_h5(h5)

        # check that dtypes are correct
        self.assertEqual(ra.features.data.dtype, np.int64)

        ua = ra.to_unique(0)

        # check that dtypes are correct
        self.assertEqual(ua.features.dtype, np.int64)

        exp = ua.to_experiment(0)
        self.assertEqual(len(ua), exp.reads.counts.sum().sum())

        # check that the types are correct
        self.assertEqual(exp.reads.columns.dtype, np.int64)
        self.assertEqual(exp.molecules.columns.dtype, np.int64)

        exp2 = ua.to_experiment_old(0)

        # create molecules from the uniquearray
        test_data = defaultdict(lambda: defaultdict(list))
        for i, row in enumerate(ua.data):
            cell, rmt, feature = row[0], row[1], ua.features[i]
            test_data[feature][cell].append(rmt)

        # convert to counts, reads
        reads = defaultdict(lambda: defaultdict(int))
        molecules = defaultdict(lambda: defaultdict(int))
        for feature in test_data:
            for cell, rmts in test_data[feature].items():
                reads[feature][cell] = len(rmts)
                molecules[feature][cell] = len(set(rmts))

        # make comparison DataFrames
        reads = pd.DataFrame(reads).fillna(0).astype(np.int32)
        molecules = pd.DataFrame(molecules).fillna(0).astype(np.int32)

        # index so they equal the exp2
        reads_matching_exp2 = reads.ix[exp2.reads.index, exp2.reads.columns].values
        molecules_matching_exp2 = molecules.ix[exp2.molecules.index,
                                               exp2.molecules.columns].values

        self.assertTrue(np.array_equal(reads_matching_exp2, exp2.reads.counts.todense()))
        self.assertTrue(np.array_equal(molecules_matching_exp2,
                                        exp2.molecules.counts.todense()))

        # todo these are to test refactored to_experiment()
        # self.assertTrue(np.array_equal(exp.reads.counts, exp2.reads.counts))
        # self.assertTrue(np.array_equal(exp.molecules.counts, exp2.molecules.counts))


class ThreeBitGeneralTest(unittest.TestCase):
    def test_3bit_mars_seq(self):
        self.assertRaises(NotImplementedError, seqc.encodings.ThreeBit.default_processors,
                          'mars-seq')

    def test_3bit_drop_seq(self):
        cell = 'ACGTACGTACGT'
        rmt = 'AACCGGTT'
        poly_t = ''
        sample_seq = cell + rmt
        tb = seqc.encodings.ThreeBit.default_processors('drop-seq')
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
        tb = seqc.encodings.ThreeBit.default_processors('cel-seq')
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
        tb = seqc.encodings.ThreeBit.default_processors('avo-seq')
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

        tb = seqc.encodings.ThreeBit.default_processors('avo-seq')

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
        bin_string = seqc.encodings.ThreeBit.str2bin(test_string)
        result = seqc.encodings.ThreeBit.gc_content(bin_string)
        self.assertEqual(result, expected_result)


class IndexGenerationTest(unittest.TestCase):
    # todo import these tests from scseq/seqdb
    pass


class ResolveAlignmentsTest(unittest.TestCase):
    """this may need more tests for more sophisticated input data."""
    # todo rewrite for standard locations

    def setUp(self):
        self.data_dir = '/'.join(seqc.__file__.split('/')[:-2]) + '/data/'
        self.expectations = self.data_dir + 'genome/mm38_chr19/p_coalignment.pckl'

    @unittest.skip('only run this if new ra disambiguation input is needed')
    def test_generate_disambiguation_read_array(self):
        cell_barcodes = self.data_dir + 'in_drop/barcodes/cb_3bit.p'
        save = self.data_dir + 'in_drop/disambiguation_ra_input.p'
        _ = seqc.generate_in_drop_read_array(self.expectations, cell_barcodes, 1000, 10,
                                             save=save)

    @unittest.skip('')
    def test_ra_disambiguate(self):
        expectations = self.data_dir + 'genome/mm38_chr19/p_coalignment_array.p'
        arr_pickle = self.data_dir + 'in_drop/disambiguation_ra_input.p'
        with open(arr_pickle, 'rb') as f:
            ra = pickle.load(f)
        ra.resolve_alignments(expectations)
        # 4 == complete disambiguation.
        vals, counts = np.unique(ra.data['disambiguation_results'], return_counts=True)
        self.assertTrue(np.array_equal(vals, np.array([2, 4])))


class GeneTableTest(unittest.TestCase):
    def setUp(self):
        self.genome_dir = '/'.join(seqc.__file__.split('/')[:-2]) + '/data/genome/'

    def test_gene_table(self):
        gt = seqc.convert_features.GeneTable(config.gtf)
        chromosome = 'chr19'
        start = 4007800
        end = 4007900
        strand = '-'
        genes = gt.coordinates_to_gene_ids(chromosome, start, end, strand)


# todo reimplement these tests
class SMARTSeqSamToCountTest(unittest.TestCase):
    def setUp(self):
        self.data_dir = '/'.join(seqc.__file__.split('/')[:-2]) + '/data/'

    @unittest.skip('not working atm.')
    def test_sam_to_count_single_file(self):
        samfile = self.data_dir + 'test_seqc_merge_drop_seq_fastq/Aligned.out.sam'
        gtf = self.data_dir + '/genome/mm38_chr19/annotations.gtf'
        coo, gene_index, ci = seqc.sam.to_count_single_file(samfile, gtf, 100)
        # print(repr(coo))
        # print(len(gene_index))

        samfile = self.data_dir + 'test_seqc_merge_in_drop_fastq/Aligned.out.sam'
        gtf = self.data_dir + '/genome/mm38_chr19/annotations.gtf'
        coo, gene_index, ci = seqc.sam.to_count_single_file(samfile, gtf, 100)
        # print(repr(coo))
        # print(len(gene_index))

    @unittest.skip('not working atm.')
    def test_sam_to_count_multiple_files(self):
        samfile = self.data_dir + 'test_seqc_merge_drop_seq_fastq/Aligned.out.sam'
        gtf = self.data_dir + '/genome/mm38_chr19/annotations.gtf'
        coo, gene_index, ci = seqc.sam.to_count_multiple_files([samfile, samfile], gtf,
                                                               100)
        # print(repr(coo))
        # print(len(gene_index))

        samfile = self.data_dir + 'test_seqc_merge_in_drop_fastq/Aligned.out.sam'
        gtf = self.data_dir + '/genome/mm38_chr19/annotations.gtf'
        coo, gene_index, ci = seqc.sam.to_count_multiple_files([samfile, samfile], gtf,
                                                               100)
        # print(repr(coo))
        # print(len(gene_index))


class GenerateSRATest(unittest.TestCase):
    def setUp(self):
        self.data_dir = '/'.join(seqc.__file__.split('/')[:-2]) + '/data/'
        self.sra_dir = self.data_dir + '.test_dump_sra/'
        if not os.path.isdir(self.sra_dir):
            os.makedirs(self.sra_dir)

    @unittest.skip("doesn't work")
    def test_md5sum(self):
        res = seqc.sra.SRAGenerator.md5sum(self.data_dir + 'in_drop/sample_data_r1.fastq')
        # print(res)

    @unittest.skip('')
    def test_generate_experiment_xml_file(self):
        seqc.sra.SRAGenerator.make_experiment_xml('mm38', 'DUMMY_SAMPLE',
                                                  fout_stem=self.sra_dir + 'EXP_XML')
        xml_data = xml.dom.minidom.parse(self.sra_dir + 'EXP_XML_experiment.xml')
        pretty = xml_data.toprettyxml()
        # print('\n', pretty)

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
        # print('\n', pretty)

    @unittest.skip('')
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


class GroupForErrorCorrectionTest(unittest.TestCase):

    def setUp(self):
        # create a dummy array with 10 records
        dtype = seqc.arrays.ReadArray._dtype
        dummy_data = np.zeros((20,), dtype=dtype)

        # we'll use these cell barcodes; no palindromes
        s2b = seqc.encodings.ThreeBit.str2bin
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
        i2i = seqc.encodings.ThreeBit.ints2int
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

    def test_ints2int(self):
        i2i = seqc.encodings.ThreeBit.ints2int
        for record in self.ra.data:
            i = i2i([int(record['cell']), int(record['rmt'])])
            self.assertTrue(i > 0)

    def test_group_for_error_correction(self):
        grp = self.ra.group_for_error_correction()
        # print(self.ra.data)
        # print(grp)
        # print(grp.keys())
        b2s = seqc.encodings.ThreeBit.bin2str
        for feature in grp.keys():
            for seq in grp[feature].keys():
                if seq:
                    b2s(seq)
                else:
                    pass  # not aligned


class ConvertSparseGeneIDsTest(unittest.TestCase):

    test_dir = 'test_seqc/'

    @classmethod
    def setUpClass(cls):
        for dtype in config.data_types:
            check_h5(dtype)
        if not os.path.isdir(cls.test_dir):
            os.makedirs(cls.test_dir)

    @params(*config.data_types)
    def test_convert_sparse_gene_ids(self, data_type):
        ra = seqc.ReadArray.from_h5(config.h5_name_pattern % data_type)
        gmap = config.gene_map_pattern % data_type
        if data_type == 'drop_seq':
            exp = seqc.Experiment.from_read_array(ra, 0, 2)
        else:
            exp = seqc.Experiment.from_read_array(ra, 3, 2)
        self.assertEqual(exp.reads.columns.dtype, np.int64)
        self.assertEqual(exp.molecules.columns.dtype, np.int64)
        exp.convert_gene_ids(gmap)
        self.assertIsInstance(exp.molecules.columns[0], str)
        self.assertIsInstance(exp.molecules.columns[0], str)

        # make sure the result is pickleable
        exp.to_npz(self.test_dir + 'test_seqc.npz')

    @classmethod
    def tearDownClass(cls):
        if os.path.isdir('test_seqc'):
            shutil.rmtree('test_seqc')


@unittest.skip('Not implemented')
class DownloadInputFilesTest(unittest.TestCase):

    test_dir = 'test_seqc/'

    @classmethod
    def setUpClass(cls):
        if not os.path.isdir(cls.test_dir):
            os.makedirs(cls.test_dir)

    def test_forward_fastq_basespace_pattern(self):
        forward_file = 'Day0-ligation-11-17_S1_L001_R1_001.fastq.gz'
        reverse_file = 'Day0-ligation-11-17_S1_L001_R2_001.fastq.gz'
        forward_pattern = r'_R1_.*?\.fastq\.gz'
        reverse_pattern = r'_R2_.*?\.fastq\.gz'
        self.assertTrue(re.search(forward_pattern, forward_file))
        self.assertFalse(re.search(reverse_pattern, forward_file))
        self.assertTrue(re.search(reverse_pattern, reverse_file))
        self.assertFalse(re.search(forward_pattern, reverse_file))

    @unittest.skip('extremely slow due to s3 download speeds')
    def test_download_input_files_incorrect_input_raises(self):
        data_type = 'in_drop'
        complete_kwargs = dict(
            output_dir=self.test_dir,
            forward_fastq=[config.s3_forward_fastq_pattern % data_type],
            reverse_fastq=[config.s3_reverse_fastq_pattern % data_type],
            merged=config.s3_merged_fastq_pattern % data_type,
            samfile=config.s3_sam_pattern % data_type
        )

        empty_kwargs = dict(
            output_dir=self.test_dir,
            forward_fastq=[],
            reverse_fastq=[],
            merged='',
            samfile='',
        )

        # check empty input
        self.assertRaises(ValueError, seqc.core.check_input_data, **empty_kwargs)

        # check wrong input
        wrong_input_type = complete_kwargs.copy()
        wrong_input_type['forward_fastq'] = 'input_type_should_be_list_not_str'
        self.assertRaises(TypeError, seqc.core.check_input_data, **wrong_input_type)

        # check passing multiple arguments raises
        self.assertRaises(ValueError, seqc.core.check_input_data, **complete_kwargs)

        # check passing single fastq input raises
        partial_fastq_input = empty_kwargs.copy()
        partial_fastq_input['forward_fastq'] = complete_kwargs['forward_fastq']
        self.assertRaises(ValueError, seqc.core.check_input_data, **partial_fastq_input)
        partial_fastq_input = empty_kwargs.copy()
        partial_fastq_input['reverse_fastq'] = complete_kwargs['reverse_fastq']
        self.assertRaises(ValueError, seqc.core.check_input_data, **partial_fastq_input)

        # check download fastq; directory exists -- raises
        download_fastq = empty_kwargs.copy()
        download_fastq['forward_fastq'] = complete_kwargs['forward_fastq']
        download_fastq['reverse_fastq'] = complete_kwargs['reverse_fastq']
        os.makedirs(self.test_dir + 'forward_fastq/')
        self.assertRaises(FileExistsError, seqc.core.check_input_data, **download_fastq)
        os.rmdir(self.test_dir + 'forward_fastq/')

        # set some directory names
        forward_dir = self.test_dir + 'forward_fastq/'
        reverse_dir = self.test_dir + 'reverse_fastq/'
        sam_dir = self.test_dir + 'sam/'
        merged_dir = self.test_dir + 'merged_fastq/'

        add_prefix = lambda prefix, list_: [prefix + l for l in list_]

        # test download fastq
        self.assertFalse(os.path.isdir(forward_dir))
        self.assertFalse(os.path.isdir(reverse_dir))
        fwd, rev, mgd, sam = seqc.core.check_input_data(**download_fastq)
        self.assertEqual(add_prefix(forward_dir, sorted(os.listdir(forward_dir))), fwd)
        self.assertEqual(add_prefix(reverse_dir, sorted(os.listdir(reverse_dir))), rev)

        self.assertFalse(os.path.isdir(sam_dir))
        self.assertFalse(os.path.isdir(merged_dir))

        # clean up the directories
        for d in (forward_dir, merged_dir, sam_dir, merged_dir):
            if os.path.isdir(d):
                shutil.rmtree(d)

        # test download merged fastq
        download_merged = empty_kwargs.copy()
        download_merged['merged'] = complete_kwargs['merged']
        self.assertTrue(download_merged['merged'])  # was showing mgd == ''
        self.assertFalse(os.path.isdir(merged_dir))
        fwd, rev, mgd, sam = seqc.core.check_input_data(**download_merged)
        self.assertEqual(add_prefix(merged_dir, sorted(os.listdir(merged_dir)))[0], mgd)

        # clean up the directories
        for d in (forward_dir, merged_dir, sam_dir, merged_dir):
            if os.path.isdir(d):
                shutil.rmtree(d)

        # test download merged fastq
        download_sam = empty_kwargs.copy()
        download_sam['samfile'] = complete_kwargs['samfile']
        self.assertFalse(os.path.isdir(sam_dir))
        fwd, rev, mgd, sam = seqc.core.check_input_data(**download_sam)
        self.assertEqual(add_prefix(sam_dir, sorted(os.listdir(sam_dir)))[0], sam)

        # clean up the directories
        for d in (forward_dir, merged_dir, sam_dir, merged_dir):
            if os.path.isdir(d):
                shutil.rmtree(d)

        # confirm that passing non-s3 links passes-through arguments as normal.
        complete_local_kwargs = dict(
            forward_fastq=['local_fastq_r1.fastq'],
            reverse_fastq=['local_fastq_r2.fastq'],
            merged='merged.fastq',
            samfile='alignments.sam'
        )

        # test each in sequence
        expected_results = ([], [], '', complete_local_kwargs['samfile'])
        local_sam_kwargs = empty_kwargs.copy()
        local_sam_kwargs['samfile'] = complete_local_kwargs['samfile']
        res = seqc.core.check_input_data(**local_sam_kwargs)
        self.assertEqual(expected_results, res)

        expected_results = ([], [], complete_local_kwargs['merged'], '')
        local_merged_kwargs = empty_kwargs.copy()
        local_merged_kwargs['merged'] = complete_local_kwargs['merged']
        res = seqc.core.check_input_data(**local_merged_kwargs)
        self.assertEqual(expected_results, res)

        expected_results = (complete_local_kwargs['forward_fastq'],
                            complete_local_kwargs['reverse_fastq'], '', '')
        local_fastq_kwargs = empty_kwargs.copy()
        local_fastq_kwargs['forward_fastq'] = complete_local_kwargs['forward_fastq']
        local_fastq_kwargs['reverse_fastq'] = complete_local_kwargs['reverse_fastq']
        res = seqc.core.check_input_data(**local_fastq_kwargs)
        self.assertEqual(expected_results, res)

    def test_download_base_space(self):
        """
        unittest to make sure that BaseSpace is downloading properly in the context
        of seqc.
        """
        raise NotImplementedError  # todo test

    @classmethod
    def tearDownClass(cls):
        if os.path.isdir(cls.test_dir):
            shutil.rmtree(cls.test_dir)


@unittest.skip('Not implemented')
class DownloadBaseSpaceTest(unittest.TestCase):
    """unittests to make sure BaseSpace is correctly functioning"""

    def test_download_base_space(self):
        raise NotImplementedError  # todo test


@unittest.skip('Not implemented')
class SSHToolsTest(unittest.TestCase):

    def test_wrong_ssh_key_raises(self):
        # is it possible to test this without creating an instance? Can one ssh into one's
        # own computer or something to test?
        raise NotImplementedError
        self.assertRaises(seqc.ssh_utils.SSHServer())


class ExperimentTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        for dtype in config.data_types:
            check_experiment(dtype)

    @params(*config.data_types)
    def test_equalize_features(self, data_type):
        e = seqc.Experiment.from_npz(config.experiment_pattern % data_type)
        new = e.equalize_features()
        self.assertTrue(np.array_equal(new.molecules.columns, new.reads.columns))

    @params(*config.data_types)
    def test_equalize_cells(self, data_type):
        e = seqc.Experiment.from_npz(config.experiment_pattern % data_type)
        new = e.equalize_cells()
        self.assertTrue(np.array_equal(new.molecules.index, new.reads.index))


class DenseCountsTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        for dtype in config.data_types:
            check_experiment(dtype)

    @params(*config.data_types)
    def test_downsample(self, data_type):
        e = seqc.Experiment.from_npz(config.experiment_pattern % data_type)
        dc = e.molecules.to_dense()  # no need for threshold, don't have many counts
        res = dc.downsample()


# @unittest.skip('too much printing at the moment')
class SparseCountsTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.data = {}
        for data_type in config.data_types:
            check_sam(data_type)
            check_experiment(data_type)

            # set poly-t requirements
            if data_type == 'drop_seq':
                n_poly_t_required = 0
            else:
                n_poly_t_required = 3

            # set required_support
            required_support = 0

            # molecules per gene per cell = molecules[cell][gene] = set(rmt)
            molecule_rmts = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))

            # reads per gene per cell = molecules[cell][gene] = #(int)
            read_counts = defaultdict(lambda: defaultdict(int))

            # process the samfile
            for r in seqc.sam.Reader(config.samfile_pattern % data_type):

                if r.optional_fields['NH'] != 1:
                    continue
                if r.is_unmapped:
                    continue
                if r.n_poly_t < n_poly_t_required:
                    continue

                gene = r.name.split(':')[-1]
                molecule_rmts[r.cell][gene][r.rmt] += 1
                read_counts[r.cell][gene] += 1

            # post-process the molecule counts
            molecule_counts = defaultdict(lambda: defaultdict(int))
            for cell in molecule_rmts.keys():
                for gene in molecule_rmts[cell].keys():
                    for rmt, count in molecule_rmts[cell][gene].items():
                        if count > required_support:  # filter based on required support
                            molecule_counts[cell][gene] += 1

            # compare to the count matrix, generated from the samfile
            h5 = seqc.sam.to_h5(config.samfile_pattern % data_type,
                                config.h5_name_pattern % data_type,
                                6, int(1e7), config.gtf)
            ra = seqc.ReadArray.from_h5(h5)

            # make sure the read array is the right size
            sam_reader = seqc.sam.Reader(config.samfile_pattern % data_type)
            n_records = sum(1 for _ in sam_reader.iter_multialignments())

            # todo 1-off error here; 39570 vs 39569
            assert len(ra) == n_records, '%d != %d' % (len(ra), n_records)
            # assert len(ra) == n_records * 0.95, '%d != %d' % (len(ra), n_records)

            # todo filters make a big difference
            ua = ra.to_unique(n_poly_t_required)

            # make sure the UniqueArray is the right size
            n_unique = 0
            for r in seqc.sam.Reader(config.samfile_pattern % data_type):
                if r.optional_fields['NH'] == 1:
                    n_unique += 1
            assert len(ua) >= (n_unique * .95), '%d !>= %d * .95' % (len(ua), n_unique)

            e = ua.to_experiment(required_support=required_support)

            # convert sparse counts to dataframe
            id_map = seqc.convert_features.ConvertGeneCoordinates.from_pickle(
                config.gene_map_pattern % data_type)
            e.molecules.convert_gene_ids(id_map)
            mol_df = e.molecules.to_dense().df
            e.reads.convert_gene_ids(id_map)
            read_df = e.reads.to_dense().df

            # convert dict of dicts to df
            test_mol_df = pd.DataFrame(molecule_counts).T
            test_read_df = pd.DataFrame(read_counts).T

            if len(test_mol_df.columns) != len(set(test_mol_df.columns)):
                print('%s: duplicate molecule features detected' % data_type)
            if len(test_read_df.columns) != len(set(test_read_df.columns)):
                print('%s: duplicate read features detected' % data_type)

            # remove duplicates
            test_mol_df = test_mol_df.T
            test_mol_df = test_mol_df.groupby(test_mol_df.index).sum().T

            test_read_df = test_read_df.T
            test_read_df = test_read_df.groupby(test_read_df.index).sum().T

            cls.data[data_type] = {
                'pipeline': {
                    'molecules': test_mol_df,
                    'reads': test_read_df
                }, 'sam': {
                    'molecules': mol_df,
                    'reads': read_df
                }
            }

    @params(*config.data_types)
    def test_plot_fraction_mitochondrial(self, data_type):
        self.e = seqc.Experiment.from_npz(config.experiment_pattern % data_type)
        self.assertRaises(TypeError, self.e.molecules.plot_fraction_mitochondrial,
                          molecules=True)

        # convert ids
        self.e.convert_gene_ids(config.gene_map_pattern % data_type)
        self.assertRaises(ValueError, self.e.molecules.plot_fraction_mitochondrial,
                          molecules=True)

        # todo add some mitochondrial RNA to the data somehow; otherwise cannot test.
        # these tests only test type-checking; very weak.

    @params(*config.data_types)
    def test_identifies_same_set_of_cells_for_reads(self, data_type):
        pipeline_cells = self.data[data_type]['pipeline']['reads'].index
        sam_cells = self.data[data_type]['sam']['reads'].index
        difference = set(pipeline_cells).symmetric_difference(sam_cells)
        self.assertEqual(len(difference), 0,
                         '%s: the following cells are found in either the processed '
                         'pipeline data or the sam file data but not both: %s'
                         % (data_type, repr(difference)))

    @params(*config.data_types)
    def test_identifies_same_set_of_cells_for_molecules(self, data_type):
        pipeline_cells = self.data[data_type]['pipeline']['molecules'].index
        sam_cells = self.data[data_type]['sam']['molecules'].index
        difference = set(pipeline_cells).symmetric_difference(sam_cells)
        self.assertEqual(len(difference), 0,
                         '%s: the following cells are found in either the processed '
                         'pipeline data or the sam file data but not both: %s'
                         % (data_type, repr(difference)))

    @params(*config.data_types)
    def test_identifies_same_set_of_genes_for_reads(self, data_type):
        pipeline_genes = self.data[data_type]['pipeline']['reads'].columns
        sam_genes = self.data[data_type]['sam']['reads'].columns
        difference = set(pipeline_genes).symmetric_difference(sam_genes)
        self.assertEqual(len(difference), 0,
                         '%s: the following genes are found in either the processed '
                         'pipeline data or the sam file data but not both: %s'
                         % (data_type, repr(difference)))

    @params(*config.data_types)
    def test_identifies_same_set_of_genes_for_molecules(self, data_type):
        pipeline_genes = self.data[data_type]['pipeline']['molecules'].columns
        sam_genes = self.data[data_type]['sam']['molecules'].columns
        difference = set(pipeline_genes).symmetric_difference(sam_genes)
        self.assertEqual(len(difference), 0,
                         '%s: the following genes are found in either the processed '
                         'pipeline data or the sam file data but not both: %s'
                         % (data_type, repr(difference)))

    @params(*config.data_types)
    def test_molecule_matrices_are_same_shape(self, data_type):
        pipeline_shape = self.data[data_type]['pipeline']['molecules'].shape
        sam_shape = self.data[data_type]['sam']['molecules'].shape
        self.assertEqual(pipeline_shape, sam_shape)

    @params(*config.data_types)
    def test_read_matrices_are_same_shape(self, data_type):
        pipeline_shape = self.data[data_type]['pipeline']['reads'].shape
        sam_shape = self.data[data_type]['sam']['reads'].shape
        self.assertEqual(pipeline_shape, sam_shape)

    @params(*config.data_types)
    def test_read_matrices_have_same_total_counts(self, data_type):
        pipeline_counts = self.data[data_type]['pipeline']['reads'].sum().sum()
        sam_counts = self.data[data_type]['sam']['reads'].sum().sum()
        self.assertEqual(pipeline_counts, sam_counts)

    @params(*config.data_types)
    def test_molecule_matrices_have_same_total_counts(self, data_type):
        pipeline_counts = self.data[data_type]['pipeline']['molecules'].sum().sum()
        sam_counts = self.data[data_type]['sam']['molecules'].sum().sum()
        self.assertEqual(pipeline_counts, sam_counts)


class SamReaderTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        for data_type in config.data_types:
            check_sam(data_type)

    # todo weak test; only tests running without errors
    @params(*config.data_types)
    def test_run(self, data_type):
        reader = seqc.sam.Reader(config.samfile_pattern % data_type)
        first_record = next(iter(reader))

        # check all the fields
        self.assertIsInstance(first_record.optional_fields['NH'], int)
        self.assertIsInstance(first_record.rmt, int)
        self.assertIsInstance(first_record.cell, int)
        self.assertIsInstance(first_record.n_poly_t, int)
        self.assertIsInstance(first_record.valid_cell, int)
        self.assertIsInstance(first_record.dust_score, int)
        self.assertIsInstance(first_record.fwd_quality, int)
        self.assertIsInstance(first_record.qname, str)
        self.assertIsInstance(first_record.flag, str)
        self.assertIsInstance(first_record.rname, str)
        self.assertIsInstance(first_record.pos, str)
        self.assertIsInstance(first_record.mapq, str)
        self.assertIsInstance(first_record.cigar, str)
        self.assertIsInstance(first_record.rnext, str)
        self.assertIsInstance(first_record.pnext, str)
        self.assertIsInstance(first_record.tlen, str)
        self.assertIsInstance(first_record.seq, str)
        self.assertIsInstance(first_record.qual, str)
        self.assertTrue(first_record.is_mapped)
        self.assertFalse(first_record.is_unmapped)


if __name__ == '__main__':
    suite = unittest.TestSuite()
    # suite.addTest(AlignSTARTest('test_align_giberish'))
    # suite.addTest(FastqGenerateTest('test_generate_reverse_complex'))
    suite.addTest(FastqGenerateTest('test_generate_typed_fastq'))
    runner = unittest.TextTestRunner()
    runner.run(suite)

    # nose2.main()
    # final cleanup
    if os.path.isdir('test_seqc'):
        shutil.rmtree('test_seqc')
    if os.path.isdir('seqc_test'):
        shutil.rmtree('seqc_test')
    # awful that we need to have this; can't find out how this is getting created
    if os.path.isdir('U'):
        shutil.rmtree('U')