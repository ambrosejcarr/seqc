__author__ = 'ambrose'

import hashlib
from memory_profiler import profile, memory_usage
from numpy.lib.recfunctions import append_fields
import unittest
import os
import pickle
import seqc
import tables as tb
from seqc import three_bit, fastq, align, sam, convert_features, barcodes, io_lib
from seqc import arrays
from io import StringIO
from itertools import product
from xml.etree import ElementTree as ET
from subprocess import Popen, PIPE
import cProfile
from pstats import Stats
import xml.dom.minidom
import logging
import socket
import re
import numpy as np
import pandas as pd
import random
import gzip
import ftplib
import shutil

################################ STATE CONFIGURATION ####################################

# this is the universal data dir for these tests
_seqc_dir = '/'.join(seqc.__file__.split('/')[:-3]) + '/'
_data_dir = '/'.join(seqc.__file__.split('/')[:-2]) + '/data/'

# set of dtypes
_data_types = ['in_drop', 'drop_seq']

# universal file patterns
_samfile_pattern = _seqc_dir + 'test/%s/.seqc_test/Aligned.out.sam'
_forward_pattern = _seqc_dir + 'test/%s/fastq/test_seqc_r1.fastq'
_reverse_pattern = _seqc_dir + 'test/%s/fastq/test_seqc_r2.fastq'
_barcode_pattern = _seqc_dir + 'test/%s/barcodes/barcodes.p'
_h5_name_pattern = _seqc_dir + 'test/%s/test_seqc.h5'

# universal index files
_gtf = _seqc_dir + 'test/genome/annotations.gtf'
_fasta = _seqc_dir + 'test/genome/mm38_chr19.fa'
_index = _seqc_dir + 'test/genome/'

# config parameters
_n_threads = 7

############################### SETUP FILE GENERATION ###################################


def check_fastq(data_type: str):
    """check if the required fastq files are present, else generate them"""

    # replace any dashes with underscores
    data_type = data_type.replace('-', '_')

    # get forward and reverse file ids
    forward = _forward_pattern % data_type
    reverse = _reverse_pattern % data_type

    if not all(os.path.isfile(f) for f in [forward, reverse]):
        gen_func = getattr(seqc.fastq.GenerateFastq, data_type)
        gen_func(10000, _seqc_dir + 'test/%s/fastq/test_seqc', _fasta, _gtf)

    return forward, reverse


def check_sam(data_type: str):
    """check if the required sam files are present, else generate them"""

    # replace any dashes with underscores
    data_type = data_type.replace('-', '_')

    # get sam pattern
    samfile = _samfile_pattern % data_type

    # generation params
    n = 10000
    prefix = _seqc_dir + 'test/%s/fastq/test_seqc'

    if not os.path.isfile(samfile):
        gen_func = getattr(seqc.sam.GenerateSam, data_type)
        gen_func(n=n, prefix=prefix, fasta=_fasta, gtf=_gtf, index=_index)

    return samfile


def check_index():
    """ensure that there is an index present. If not, download it."""
    gfiles = ['Genome', 'SA', 'SAindex', 'annotations.gtf']
    if not os.path.isdir(_index) or not all(os.path.isfile(_index + f) for f in gfiles):
        index_bucket = 'dplab-data'
        index_prefix = 'genomes/mm38_chr19/'
        seqc.io_lib.S3.download_files(
            bucket=index_bucket, key_prefix=index_prefix, output_prefix=_index,
            cut_dirs=2)


##################################### UNIT TESTS ########################################

class TestFastq




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
        jarr = arrays.JaggedArray.from_iterable(data_size, data)
        self.assertTrue(jarr._data.dtype == np.uint32)
        self.assertTrue(jarr._data.shape == (data_size,))
        print(jarr[10])


# todo keep
@unittest.skip('')
class TestThreeBitInDrop(unittest.TestCase):

    def setUp(self):
        self.tb = three_bit.ThreeBit.default_processors('in-drop')
        # note: don't test with palindromic sequences, these can produce unexpected
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
        s = three_bit.ThreeBit.str2bin(self.example_seqs[0])
        seq = three_bit.ThreeBit.bin2str(s)
        self.assertTrue(seq == self.example_seqs[0])

    def test_bin2str_inverts_str2bin_with_ints2int(self):
        b1 = three_bit.ThreeBit.str2bin(self.c1s[0])
        b2 = three_bit.ThreeBit.str2bin(self.c1s[1])
        b12 = three_bit.ThreeBit.ints2int([b1, b2])
        self.assertTrue(b12 == three_bit.ThreeBit.str2bin(self.c1s[0] + self.c1s[1]))

    def test_str2bin_converts_forwards_not_backwards(self):
        s = three_bit.ThreeBit.str2bin('ACGT')
        self.assertEqual(s, 0b100110101011)
        self.assertNotEqual(s, 0b011101110100)

    def test_ints2int_and_cells(self):
        s1 = 'AAAAAAAA'
        s2 = 'TTTTTTTT'
        seqs = [three_bit.ThreeBit.str2bin(s) for s in [s1, s2]]
        c = three_bit.ThreeBit.ints2int(seqs)
        cells = self.tb.split_cell(c)
        c1, c2 = [three_bit.ThreeBit.bin2str(c) for c in cells]
        self.assertTrue(s1 == c1)
        self.assertTrue(s2 == c2)

    def test_extract_cell_spacer_extraction(self):

        results = [3446, 2478, 2869, 1894]
        for i, strseq in enumerate(self.example_seqs):
            # note this tracks code in extract cell; test will not reflect function
            # if extract cell is changed (not an ideal test)
            s = three_bit.ThreeBit.str2bin(strseq)
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
            str_cell = three_bit.ThreeBit.bin2str(cell)
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
        self.assertRaises(NotImplementedError, three_bit.ThreeBit.default_processors,
                          'mars-seq')

    def test_3bit_drop_seq(self):
        cell = 'ACGTACGTACGT'
        rmt = 'AACCGGTT'
        poly_t = ''
        sample_seq = cell + rmt
        tb = three_bit.ThreeBit.default_processors('drop-seq')
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
        tb = three_bit.ThreeBit.default_processors('cel-seq')
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
        tb = three_bit.ThreeBit.default_processors('avo-seq')
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

        tb = three_bit.ThreeBit.default_processors('avo-seq')

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
        bin_string = three_bit.ThreeBit.str2bin(test_string)
        result = three_bit.ThreeBit.gc_content(bin_string)
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
        cb = barcodes.CellBarcodes(*barcode_files)
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
        tb = three_bit.ThreeBit.default_processors('in-drop')
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
        cb = barcodes.CellBarcodes(cb1, cb2)
        perfect_code = three_bit.ThreeBit.str2bin('AAACCC')
        error_code = three_bit.ThreeBit.str2bin('NTTCCC')
        self.assertTrue(cb.perfect_match(perfect_code))
        self.assertFalse(cb.perfect_match(error_code))
        self.assertTrue(cb.close_match(perfect_code))
        self.assertTrue(cb.close_match(error_code))

    def test_in_drop_error_code_error_type_identification(self):
        cb1 = StringIO('AAA\nTTT\n')
        cb2 = StringIO('CCC\nGGG\n')
        cb = barcodes.CellBarcodes(cb1, cb2)
        error_code = three_bit.ThreeBit.str2bin('NTTCCC')

        errors = cb.map_errors(error_code)
        self.assertEqual({'TN'}, set(errors))


# todo re-write
@unittest.skip('')
class TestFastq(unittest.TestCase):
    pass


# todo import these tests from scseq/seqdb
@unittest.skip('')
class TestIndexGeneration(unittest.TestCase):
    pass


# todo
# contains memory and time profiling examples for other unit tests; final is a
# relevant test
@unittest.skip('')
class TestSamToReadArray(unittest.TestCase):

    def setUp(self):
        self.data_dir = '/'.join(seqc.__file__.split('/')[:-2]) + '/data/'
        self.in_drop_samfile = (
            self.data_dir + 'test_seqc_merge_in_drop_fastq/Aligned.out.sam')
        self.drop_seq_samfile = (
            self.data_dir + 'test_seqc_merge_drop_seq_fastq/Aligned.out.sam')
        self.gtf = self.data_dir + 'genome/mm38_chr19/annotations.gtf'
        self.ftable, self.fpositions = sam.construct_feature_table(
            self.gtf, fragment_len=1000)
        self.h5name = self.data_dir + 'test_seqc_merge_in_drop_fastq/ra.h5'

    @unittest.skip('')
    def test_construct_read_array(self):
        a = arrays.ReadArray.from_samfile(
            self.in_drop_samfile, self.fpositions, self.ftable)
        array_lengths = (len(a.data), len(a.features), len(a.positions))
        self.assertTrue(len(set(array_lengths)) == 1)  # all arrays should have same #rows

    @unittest.skip('')
    def test_generate_indices_for_error_correction(self):
        a = arrays.ReadArray.from_samfile(
            self.in_drop_samfile, self.fpositions, self.ftable)
        _ = (len(a.data), len(a.features), len(a.positions))
        _ = a.group_for_error_correction(required_poly_t=4)

    @unittest.skip('')
    def test_save_and_load_readarray(self):
        ra = arrays.ReadArray.from_samfile(
            self.in_drop_samfile, self.fpositions, self.ftable)
        ra.save_h5(self.h5name)
        ra2 = arrays.ReadArray.from_h5(self.h5name)

        self.assertTrue(np.array_equiv(ra.data, ra2.data))
        self.assertTrue(np.array_equiv(ra._features._data, ra2._features._data))
        self.assertTrue(np.array_equiv(ra._features._index, ra2._features._index))
        self.assertTrue(np.array_equiv(ra._positions._data, ra2._positions._data))
        self.assertTrue(np.array_equiv(ra._positions._index, ra2._positions._index))

    def memory_usage(self, sizes, directory, index):

        # generate fastq files
        for s in sizes:
            fastq_filename = directory + str(s)
            generate_in_drop_fastq_data(s, fastq_filename)
        forward = [directory + str(s) + '_r1.fastq' for s in sizes]
        reverse = [directory + str(s) + '_r2.fastq' for s in sizes]

        # get cell barcodes
        with open(self.data_dir + 'in_drop/barcodes/cb_3bit.p', 'rb') as f:
            cb = pickle.load(f)

        # merge the data
        merged_files = []
        for f, r in zip(forward, reverse):
            fq, _ = fastq.merge_fastq([f], [r], 'in-drop', directory, cb)
            merged_files.append(fq)

        # align the data
        samfiles = align.STAR.align_multiple_files(
            merged_files, index, 7, directory
        )

        ft, fp = convert_features.construct_feature_table(self.gtf, 1000)

        logfile = open(directory + 'report.txt', 'w+')

        h5files = [directory + '%dm.h5' % s for s in sizes]

        # define and run functions
        def test_memory_usage(idx):
            memory_usage((arrays.ReadArray.from_samfile, (samfiles[idx], ft, fp)))
            # arr.save_h5(h5files[index])

        for i in range(len(samfiles)):
            test_memory_usage(i)

        logfile.close()

    @unittest.skip('')
    def test_ra_memory_usage_small(self):
        index = self.data_dir + 'genome/mm38_chr19/'
        working_directory = self.data_dir + 'test_ra_memory_usage/'
        if not os.path.isdir(working_directory):
            os.makedirs(working_directory)
        self.memory_usage([int(1e6)], working_directory, index)

    @unittest.skip('')
    def test_profile_mem_usage(self):
        samfile = ('/Users/ambrose/PycharmProjects/SEQC/src/data/test_ra_memory_usage/'
                   'merged_temp/Aligned.out.sam')
        ft, fp = convert_features.construct_feature_table(self.gtf, 1000)
        usage = memory_usage((arrays.ReadArray.from_samfile, (samfile, ft, fp)))
        print(np.array(usage))

    @unittest.skip('')
    def test_ra_memory_usage(self):
        data_dir = '/'.join(seqc.__file__.split('/')[:-2]) + '/data/'

        # generate fastq data for testing
        files = [data_dir + 'test_memory_usage/' + prefix for prefix in
                 ['1m', '2m', '4m', '8m', '16m']]
        index = data_dir + 'genome/mm38_chr19/'
        generate_in_drop_fastq_data(int(1e6), files[0])
        generate_in_drop_fastq_data(int(2e6), files[1])
        generate_in_drop_fastq_data(int(4e6), files[2])
        generate_in_drop_fastq_data(int(8e6), files[3])
        generate_in_drop_fastq_data(int(16e6), files[4])

        # align the data, generating sam files
        samfiles = align.STAR.align_multiple_files(
            files, index, 7, data_dir + 'test_memory_usage/')

        ft, fp = convert_features.construct_feature_table(self.gtf, 1000)

        logfile = open(data_dir + 'test_memory_usage/report.txt', 'w+')

        h5files = [f + '.h5' for f in files]

        @profile(stream=logfile)
        def test_1m():
            arr = arrays.ReadArray.from_samfile(samfiles[0], ft, fp)
            arr.save_h5(h5files[0])

        @profile(stream=logfile)
        def test_2m():
            arr = arrays.ReadArray.from_samfile(samfiles[1], ft, fp)
            arr.save_h5(h5files[1])

        @profile(stream=logfile)
        def test_4m():
            arr = arrays.ReadArray.from_samfile(samfiles[2], ft, fp)
            arr.save_h5(h5files[2])

        @profile(stream=logfile)
        def test_8m():
            arr = arrays.ReadArray.from_samfile(samfiles[3], ft, fp)
            arr.save_h5(h5files[3])

        @profile(stream=logfile)
        def test_16m():
            arr = arrays.ReadArray.from_samfile(samfiles[4], ft, fp)
            arr.save_h5(h5files[4])

        for f in [test_1m, test_2m, test_4m, test_8m, test_16m]:
            f()

        logfile.close()

    @unittest.skip('')
    def test_create_ra_object(self):
        samfile = ('/Users/ambrose/PycharmProjects/SEQC/src/data/test_ra_memory_usage/'
                   'merged_temp/Aligned.out.sam')
        fc = convert_features.ConvertFeatureCoordinates.from_gtf(self.gtf, 1000)
        res = arrays.ReadArray.from_samfile(samfile, fc)

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
        _ = generate_in_drop_read_array(self.expectations, cell_barcodes, 1000, 10,
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
        gt = convert_features.GeneTable(_gtf)
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
        io_lib.S3.upload_file(self.file1, bucket=bucket, key=key1)

        download_key = key1 + self.file1.split('/')[-1]
        fout = '.test_aws_s3/recovered_file1.txt'

        io_lib.S3.download_file(bucket=bucket, key=download_key, fout=fout)
        with open(fout, 'r') as f:
            self.assertEqual(f.read(), self.txt1)

        io_lib.S3.remove_file(bucket, download_key)

        os.remove('.test_aws_s3/recovered_file1.txt')

    def test_aws_s3_multiple_upload(self):
        # ResourceWarnings are generated by an interaction of io_lib.S3.list() and
        # the unittesting suite. These should not occur in normal usage.
        bucket = 'dplab-home'
        key_prefix = 'ajc2205/testing_aws/'
        file_prefix = '.test_aws_s3/*'

        # upload file1, file2, and file3
        io_lib.S3.upload_files(file_prefix, bucket, key_prefix)

        aws_dir = set(io_lib.S3.list(bucket, key_prefix))
        aws_files = {key_prefix + f for f in
                     ['test_aws_s3_file1.txt', 'test_aws_s3_file2.txt',
                      'test2/subfile.txt']}
        self.assertEqual(aws_files, aws_dir)

        # download the files again
        output_prefix = '.test_aws_s3/download_test/'
        io_lib.S3.download_files(bucket, key_prefix, output_prefix)
        all_downloaded_files = []
        for path, subdirs, files in os.walk(output_prefix):
            for name in files:
                all_downloaded_files.append(os.path.join(path, name))

        # remove the files that we uploaded
        io_lib.S3.remove_files(bucket, key_prefix)

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
        res = SRAGenerator.md5sum(self.data_dir + 'in_drop/sample_data_r1.fastq')
        print(res)

    @unittest.skip('')
    def test_generate_experiment_xml_file(self):

        SRAGenerator.make_experiment_xml('mm38', 'DUMMY_SAMPLE',
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
        SRAGenerator.make_run_xml(
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
        exp_xml = SRAGenerator.make_experiment_xml(
            reference, experiment_alias, platform_name, single_or_paired_end, instrument,
            fout_stem)
        run_xml = SRAGenerator.make_run_xml(
            forward_fastq, experiment_alias, forward_md5, fout_stem, data_block=None,
            reverse_fastq_files=reverse_fastq, reverse_checksum_results=reverse_md5)
        output_path = data_dir + 'test_generate_sra'
        SRAGenerator.fastq_load(file_directory, run_xml, exp_xml, output_path)


# @unittest.skip('')
class TestProcessSingleFileSCSEQExperiment(unittest.TestCase):

    def setUp(self):
        self.forward, self.reverse = check_fastq('in_drop')
        self.s3_bucket = 'dplab-home'
        self.s3_key = 'ajc2205/test_in_drop.npz'

    @unittest.skip('')
    def test_process_single_file_no_sra_download(self):

        # set some variables
        index_bucket = None
        index_key = None

        experiment_name = 'test_in_drop'
        s3_bucket = self.s3_bucket
        s3_key = self.s3_key
        cell_barcodes = ('/Users/ambrose/PycharmProjects/SEQC/src/data/in_drop/barcodes/'
                         'in_drop_barcodes.p')

        # set the index
        if not index:  # download the index
            index_dir = working_directory + 'index/'
            S3.download_files(bucket=index_bucket, key_prefix=index_key,
                              output_prefix=index_dir, no_cut_dirs=True)
            index = index_dir + index_key.lstrip('/')
        if not os.path.isdir(index):
            raise FileNotFoundError('Index does not lead to a directory')

        # merge fastq files
        merged_fastq, _ = fastq.merge_fastq(
            self.forward, self.reverse, 'in-drop', self.working_directory, cell_barcodes)

        # align the data
        sam_file = STAR.align(
            merged_fastq, index, n_threads, working_directory, reverse_fastq_file=None)

        # create the matrix
        gtf_file = index + 'annotations.gtf'
        coo, rowind, colind = sam_to_count_single_file(sam_file, gtf_file)

        numpy_archive = experiment_name + '.npz'
        with open(numpy_archive, 'wb') as f:
            np.savez(f, mat=coo, row=rowind, col=colind)

        # upload the matrix to amazon s3
        S3.upload_file(numpy_archive, s3_bucket, s3_key)

    def test_process_multiple_file_no_sra_download(self):
        # set some variables
        index = self.index
        working_directory = self.working_directory
        index_bucket = None
        index_key = None
        S3 = io_lib.S3
        STAR = align.STAR
        n_threads = 7
        sam_to_count_multiple_files = qc.sam_to_count_multiple_files
        experiment_name = 'test_in_drop'
        s3_bucket = self.s3_bucket
        s3_key = self.s3_key
        cell_barcodes = _barcode_pattern % dtype

        # potential issue: reverse should never map..
        forward = [self.forward[0]] * 3
        reverse = [self.reverse[0]] * 3

        # set the index
        if not index:  # download the index
            index_dir = working_directory + 'index/'
            S3.download_files(bucket=index_bucket, key_prefix=index_key,
                              output_prefix=index_dir, no_cut_dirs=True)
            index = index_dir + index_key.lstrip('/')
        if not os.path.isdir(index):
            raise FileNotFoundError('Index does not lead to a directory')

        # align the data
        sam_files = STAR.align_multiple_files(
            forward, index, n_threads, working_directory, reverse_fastq_files=reverse)

        # create the matrix
        gtf_file = index + 'annotations.gtf'
        coo, rowind, colind = sam_to_count_multiple_files(sam_files, gtf_file)

        numpy_archive = experiment_name + '.npz'
        with open(numpy_archive, 'wb') as f:
            np.savez(f, mat=coo, row=rowind, col=colind)

        # upload the matrix to amazon s3
        S3.upload_file(numpy_archive, s3_bucket, s3_key)


class TestGroupForErrorCorrection(unittest.TestCase):

    def setUp(self):
        # create a dummy array with 10 records
        dtype = seqc.arrays.ReadArray._dtype
        dummy_data = np.zeros((20,), dtype=dtype)

        # we'll use these cell barcodes; no palindromes
        s2b = three_bit.ThreeBit.str2bin
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
        i2i = three_bit.ThreeBit.ints2int
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

        positions = arrays.JaggedArray.from_iterable(positions)
        features = arrays.JaggedArray.from_iterable(features)

        self.ra = arrays.ReadArray(dummy_data, features, positions)

    @unittest.skip('')
    def test_ints2int(self):
        i2i = three_bit.ThreeBit.ints2int
        for record in self.ra.data:
            print(type(record['cell'].astype(int)))
            print(i2i([record['cell'].astype(int), record['rmt'].astype(int)]))

    # @unittest.skip('')
    def test_group_for_error_correction(self):
        grp = self.ra.group_for_error_correction()
        print(self.ra.data)
        print(grp)
        print(grp.keys())
        b2s = three_bit.ThreeBit.bin2str
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
        h5_name = 'test.h5'
        n_processes = 7
        chunk_size = 10000
        seqc.sam.to_h5(self.samfile, h5_name, n_processes, chunk_size,
                       self.gtf, fragment_length=1000)

    @unittest.skip('')
    def test_writing_non_parallel(self):
        h5_name = 'test.h5'
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
        fc = seqc.convert_features.ConvertFeatureCoordinates.from_gtf(self.gtf, 1000)

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
        h5_name = 'test.h5'
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
        fc = seqc.convert_features.ConvertFeatureCoordinates.from_gtf(self.gtf, 1000)

        # open all our tables
        h5f = seqc.sam.ReadArrayH5Writer(h5_name)
        h5f.create(expectedrows)
        itersam = seqc.sam._iterate_chunk(self.samfile, n=chunk_size)
        for chunk in itersam:
            processed = seqc.sam._process_chunk(chunk, fc)
            h5f.write(processed)
        h5f.close()

    def test_writing_parallel_writeobj(self):
        h5_name = 'test.h5'
        n_processes = 7
        chunk_size = 100000
        seqc.sam.to_h5(self.samfile, h5_name, n_processes, chunk_size,
                       self.gtf, fragment_length=1000)
        nlines = self.samfile


class TestCounting(unittest.TestCase):

    def setUp(self):

        # design a UniqueReadArray to test validity of method
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

        # test that sorting is working
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

    method: find the smallest real dataset that breaks it and use that to test.
    """

    def test_num_unique_samfile(self):
        fc = seqc.convert_features.ConvertFeatureCoordinates.from_gtf(_gtf, 1000)
        # todo fix this to use check_sam()
        ra = seqc.arrays.ReadArray.from_samfile(samfile, fc)

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
    import nose2
    nose2.main()
