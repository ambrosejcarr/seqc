__author__ = "Ambrose J. Carr"

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
from itertools import islice
from numpy.lib.recfunctions import append_fields
from io import StringIO

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
        gen_func(10000, prefix, config.fasta, config.gtf, barcodes=barcodes)

    return forward, reverse


def check_merged_fastq(data_type: str) -> str:
    """check if the required merged fastq files are present, else generate them"""
    forward, reverse = check_fastq(data_type)

    if not os.path.isfile(config.merged_pattern % data_type):
        forward, reverse = [forward], [reverse]
        exp_type = data_type.replace('_', '-')
        n_processors = 4
        output_dir = os.path.abspath('./test_data/%s/fastq/' % data_type)
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


def check_sam(data_type: str) -> str:
    """check if the required sam files are present, else generate them"""

    # replace any dashes with underscores
    data_type = data_type.replace('-', '_')

    # generation params
    n = 10000
    samfile = config.samfile_pattern % data_type
    prefix = config.seqc_dir + 'test_data/%s/fastq/seqc_test' % data_type
    barcodes = config.barcode_partial_serial_pattern % data_type

    expected_samfile = config.samfile_pattern % data_type
    sam_dir = '/'.join(expected_samfile.strip('/')[:-1])
    if not os.path.isdir(sam_dir):
        os.makedirs(sam_dir)

    if not os.path.isfile(samfile):
        gen_func = getattr(seqc.sam.GenerateSam, data_type)
        gen_func(n=n, filename=samfile, prefix=prefix, fasta=config.fasta, gtf=config.gtf,
                 index=config.index, barcodes=barcodes)

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
        return h5


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

    def test_reverse_three_prime(self):
        fasta = config.fasta
        gtf = config.gtf
        res = seqc.fastq.GenerateFastq._reverse_three_prime(
            5, read_length=30, fasta=fasta, gtf=gtf)
        self.assertTrue(True) # this test only tests that the data runs.

        # todo test that reads align
        # todo test that reads are converted properly


class FastqMergeTest(unittest.TestCase):
    test_dir = 'test_seqc_fastq_merge/'

    @classmethod
    def setUpClass(cls):
        for t in config.data_types:
            check_barcodes(t)
            check_fastq(t)
        if not os.path.isdir(cls.test_dir):
            os.makedirs(cls.test_dir)

    @params(*config.data_types)
    @unittest.skip('~6s; multiprocessing has some basic time-lag built in')
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

        # get the size of the original file:
        with open(forward[0]) as f:
            fsize = 2 * sum(1 for _ in f)

        # file object should have approximately same number of lines as input files
        # this number should be within 2 * 4 * n_million, where n_million is the number of
        # lines in the file
        margin = 2 * 4 * int(len(data) / (1e6 * 4))
        self.assertTrue(len(data) - margin <= len(data) <= fsize)

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

        # @classmethod
        # def tearDownClass(cls):
        #     if os.path.isdir(cls.test_dir):
        #         shutil.rmtree(cls.test_dir)
        #     pass

    @classmethod
    def tearDownClass(cls):
        if os.path.isdir(cls.test_dir):
            shutil.rmtree(cls.test_dir)


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

    output_dir = 'test_seqc/'

    def test_load_barcodes_wrong_input_type_raises(self):
        self.assertRaises(TypeError, seqc.core.check_and_load_barcodes,
                          10, self.output_dir)

    @params(*config.data_types)
    @unittest.skip('Downloads from s3, takes time, bandwidth')
    def test_load_barcodes_from_aws_link(self, dtype):
        barcode_link = config.barcode_serialized_link_pattern % dtype
        self.assertFalse(os.path.isdir('test_seqc/'))
        barcodes = seqc.core.check_and_load_barcodes(barcode_link, self.output_dir)
        self.assertIsInstance(barcodes, seqc.barcodes.CellBarcodes)

    @params(*config.data_types)
    @unittest.skip('This function takes a long time to load the serialized object. Think '
                   'about improving its run time.')
    def test_load_barcodes_local(self, dtype):
        if dtype == 'drop_seq':
            return  # None should be passed, not a link that does not target a file obj.

        barcode_file = config.barcode_serial_pattern % dtype
        barcodes = seqc.core.check_and_load_barcodes(barcode_file)
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
                          {}, {})

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
        # chr should be str
        self.assertRaises(TypeError, cgc.translate, '+', 19, 1000)
        # strand should be + or -
        self.assertRaises(ValueError, cgc.translate, 'plus', 'chr19', 1000)

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
        tbp = seqc.three_bit.ThreeBit.default_processors(data_type.replace('_', '-'))

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


class TestAlignSTAR(unittest.TestCase):

    test_dir = 'test_seqc/alignment/%s'

    @classmethod
    def setUpClass(cls):
        for dtype in config.data_types:
            if not os.path.isdir(cls.test_dir % dtype):
                os.makedirs(cls.test_dir % dtype)
            check_merged_fastq(dtype)

    @params(*config.data_types)
    @unittest.skip('slow; genome loading takes time.')
    def test_alignment(self, dtype):

        # get merged file name and make sure it's a file
        merged = config.merged_pattern % dtype
        self.assertTrue(os.path.isfile(merged), '%s is not a file' % merged)

        n_threads = 7
        samfile = seqc.align.STAR.align(merged, config.index, n_threads,
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

    @classmethod
    def tearDownClass(cls):
        if os.path.isdir('test_seqc/'):
            shutil.rmtree('test_seqc/')


class TestParallelConstructH5(unittest.TestCase):

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


    # @unittest.skip('slow; parallel has baked in waiting time')
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
        n_chunks = np.ceil(os.path.getsize(samfile) / chunk_size)

        # get expected minimum number of alignments if no reads multi-mapped
        min_alignments = n_records - (n_chunks * 2) - 2

        # check that at least the minimum number of reads are present
        self.assertGreaterEqual(len(ra), min_alignments)

        # check that at least the minimum number of reads align
        n_aligned = sum(ra.data['is_aligned'])
        # self.assertGreaterEqual(n_aligned, min_alignments)

        # all reads should have a valid cell
        self.assertTrue(np.all(ra.data['valid_cell']))

        # all reads should have nonzero cell & rmt fields
        self.assertTrue(np.all(ra.data['rmt']))
        self.assertTrue(np.all(ra.data['cell']))

        # all reads should have average quality == 40 (due to how data was generated)
        self.assertTrue(np.all(ra.data['rev_quality'] == 40))
        self.assertTrue(np.all(ra.data['fwd_quality'] == 40))

        # more than 98% of reads should have features.
        n_with_features = sum(1 for f in ra.features if np.any(f))
        self.assertTrue(n_with_features > len(ra) * .98)

    def tearDown(self):
        if os.path.isfile(self.h5_name):
            os.remove(self.h5_name)

    @classmethod
    def tearDownClass(cls):
        if os.path.isdir(cls.test_dir):
            shutil.rmtree(cls.test_dir)


class TestConvertGeneCoordinates(unittest.TestCase):

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
        for records in rd.iter_multialignments():
            if len(records) == 1:
                if not int(records[0].flag) & 4:
                    strand = records[0].strand
                    chrom = records[0].rname
                    pos = int(records[0].pos)
                    # forward = cgc.translate('+', chrom, pos)
                    # reverse = cgc.translate('-', chrom, pos)
                    # if forward and reverse:
                    #     both.append([1])
                    # elif forward:
                    #     fwd_res.append(forward)
                    # elif reverse:
                    #     rev_res.append(reverse)
                    # total += 1
                    res.append(cgc.translate(strand, chrom, pos))

        # print('plus: %d' % len(fwd_res))
        # print('minus: %d' % len(rev_res))
        # print('both: %d' % len(both))
        # print('total: %d' % total)
        # print('All features: %d' % (len(fwd_res) + len(rev_res) + len(both)))
        total = len(res)
        converted = sum(1 for r in res if r)
        self.assertTrue(total * .98 < converted, 'Of %d alignments, only %d converted' %
                        (total, converted))

    @classmethod
    def tearDownClass(cls):
        if os.path.isdir(cls.test_dir):
            shutil.rmtree(cls.test_dir)


class TestJaggedArray(unittest.TestCase):

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
        self.assertTrue(jarr._data.dtype == np.uint32)
        self.assertTrue(jarr._data.shape == (data_size,))


class TestUniqueArrayCreation(unittest.TestCase):
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

        self.assertTrue(len(ua) > len(ra) * .98)

    @classmethod
    def tearDownClass(cls):
        if os.path.isdir(cls.test_dir):
            shutil.rmtree(cls.test_dir)


class TestCountingUniqueArray(unittest.TestCase):

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


class TestExperimentCreation(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        for dtype in config.data_types:
            check_h5(dtype)

    @params(*config.data_types)
    def test_create_experiment(self, dtype):
        h5 = config.h5_name_pattern % dtype
        ra = seqc.arrays.ReadArray.from_h5(h5)
        ua = ra.to_unique(0)
        exp = ua.to_experiment(0)
        print('ReadCounts: %d' % exp.reads.counts.sum().sum())
        print('Molecule Counts: %d' % exp.molecules.counts.sum().sum())
        print('UniqueArray Length: %d' % len(ua))
        print('ReadArray length: %d' % len(ra))


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


class TestIndexGeneration(unittest.TestCase):
    # todo import these tests from scseq/seqdb
    pass


class TestResolveAlignments(unittest.TestCase):
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


# todo reimplement these tests
class TestSMARTSeqSamToCount(unittest.TestCase):
    def setUp(self):
        self.data_dir = '/'.join(seqc.__file__.split('/')[:-2]) + '/data/'

    @unittest.skip('not working atm.')
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

    @unittest.skip('not working atm.')
    def test_sam_to_count_multiple_files(self):
        samfile = self.data_dir + 'test_seqc_merge_drop_seq_fastq/Aligned.out.sam'
        gtf = self.data_dir + '/genome/mm38_chr19/annotations.gtf'
        coo, gene_index, ci = seqc.sam.to_count_multiple_files([samfile, samfile], gtf,
                                                               100)
        print(repr(coo))
        print(len(gene_index))

        samfile = self.data_dir + 'test_seqc_merge_in_drop_fastq/Aligned.out.sam'
        gtf = self.data_dir + '/genome/mm38_chr19/annotations.gtf'
        coo, gene_index, ci = seqc.sam.to_count_multiple_files([samfile, samfile], gtf,
                                                               100)
        print(repr(coo))
        print(len(gene_index))


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

    def test_ints2int(self):
        i2i = seqc.three_bit.ThreeBit.ints2int
        for record in self.ra.data:
            i = i2i([int(record['cell']), int(record['rmt'])])
            self.assertTrue(i > 0)

    def test_group_for_error_correction(self):
        grp = self.ra.group_for_error_correction()
        # print(self.ra.data)
        # print(grp)
        # print(grp.keys())
        b2s = seqc.three_bit.ThreeBit.bin2str
        for feature in grp.keys():
            for seq in grp[feature].keys():
                if seq:
                    b2s(seq)
                else:
                    pass  # not aligned


@unittest.skip('')
class TestDownloadInputFiles(unittest.TestCase):

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

@unittest.skip('')
class TestDownloadBaseSpace(unittest.TestCase):
    """unittests to make sure BaseSpace is correctly functioning"""

    def test_download_base_space(self):
        raise NotImplementedError  # todo test

if __name__ == '__main__':
    nose2.main()