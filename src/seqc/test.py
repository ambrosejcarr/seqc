__author__ = 'ambrose'

import unittest
import os
import pickle
import seqc
from seqc import three_bit, fastq, align, sam, qc, convert_features
from io import StringIO
from itertools import product
import re
import numpy as np
import pandas as pd
import random
import gzip
import ftplib


def generate_in_drop_fastq_data(n, prefix):
    data_dir = '/'.join(seqc.__file__.split('/')[:-2]) + '/data/'
    gfq = GenerateFastq()
    fwd_len = 50
    rev_len = 100
    barcodes = data_dir + 'in_drop/barcodes/in_drop_barcodes.p'
    umi_len = 6
    forward = gfq.generate_forward_fastq(n, fwd_len, barcodes, umi_len)
    forward = forward.read()  # consume the StringIO object
    reverse = gfq.generate_reverse_fastq(n, rev_len)
    reverse = reverse.read()  # consume the StringIO object
    with open(prefix + '_r1.fastq', 'w') as f:
        f.write(forward)
    with open(prefix + '_r2.fastq', 'w') as r:
        r.write(reverse)


def generate_in_drop_disambiguation_data(expectations, cell_barcodes, n, k, save=None):

    """generate N observations split between k ambiguous molecular models"""
    with open(expectations, 'rb') as f:
        expectations = pickle.load(f)

    # split the data between two cell barcodes and 2 rmts
    with open(cell_barcodes, 'rb') as f:
        cb = pickle.load(f)
    cells = np.random.choice(list(cb.perfect_codes), 2)

    # get two random rmts
    rmts = [''.join(np.random.choice(list('ACGT'), 6)) for _ in range(2)]
    rmts = [three_bit.ThreeBit.str2bin(r) for r in rmts]

    n_poly_t = 5
    valid_cell = 1
    trimmed_bases = 0
    fwd_quality = 40
    rev_quality = 40
    alignment_score = 50

    # get all the expectations that are not unique
    non_unique = {}
    for e, prob_dict in expectations.items():
        if len(prob_dict) > 1:
            non_unique[e] = prob_dict

    # get total number of models
    non_unique = pd.Series(non_unique)

    # select models
    models = np.random.choice(non_unique.index, size=k)

    # reads per model
    rpm = np.round(n / k)

    # create a container for the data
    arrays = []

    # create the data
    for m in models:
        cell = random.choice(cells)
        rmt = random.choice(rmts)
        features, probs = zip(*non_unique[m].items())
        counts = np.random.multinomial(rpm, probs)
        arr = sam.create_structured_array(sum(counts))

        i = 0
        for f, c in zip(features, counts):
            for _ in range(c):
                arr[i] = (cell, rmt, n_poly_t, valid_cell, trimmed_bases, rev_quality,
                          fwd_quality, f, tuple([0] * len(f)), alignment_score)
                i += 1
        arrays.append(arr)

    arrays = np.hstack(arrays)

    if isinstance(save, str):
        with open(save, 'wb') as f:
            pickle.dump(arrays, f)

    return arrays


class GenerateFastq(object):

    def __init__(self):
        data_dir = '/'.join(seqc.__file__.split('/')[:-2]) + '/data/'
        if os.path.isfile(data_dir + 'genome/mm38_chr19/chr19seq.p'):
            with open(data_dir + 'genome/mm38_chr19/chr19seq.p', 'rb') as f:
                self.sequence = pickle.load(f)
        else:
            self.prepare_chr19()
            with open(data_dir + 'genome/mm38_chr19/chr19seq.p', 'rb') as f:
                self.sequence = pickle.load(f)

    def generate_reverse_fastq(self, n, read_length):
        """generate an n-record in-drop fastq file with data from chromosome 19"""
        # TODO | incorporate error rates.
        names = range(n)
        name2 = '+'
        quality = 'I' * read_length
        seqs = self.generate_sequence_records(n, read_length)
        records = []
        for name, seq in zip(names, seqs):
            records.append('\n'.join(['@%d' % name, seq, name2, quality]))
        reverse_fastq = StringIO('\n'.join(records))
        return reverse_fastq

    @staticmethod
    def generate_forward_fastq(n, read_length, barcodes, umi_len):
        with open(barcodes, 'rb') as f:
            barcodes = list(pickle.load(f))
        names = range(n)
        name2 = '+'
        quality = 'I' * read_length
        records = []
        alphabet = np.array(['A', 'C', 'G', 'T'])
        for name in names:
            cb = random.choice(barcodes)
            umi = ''.join(np.random.choice(alphabet, umi_len))
            poly_a = (read_length - len(cb) - len(umi)) * 'T'
            records.append('\n'.join(['@%d' % name, cb + umi + poly_a, name2, quality]))
        forward_fastq = StringIO('\n'.join(records))
        return forward_fastq

    def generate_sequence_records(self, n, read_length):
        """generate a sequencing read from chr19"""
        idx = np.random.randint(len(self.sequence), size=n)
        sequences = []
        for i in idx:
            p = random.randint(0, max(len(self.sequence[i]) - read_length, 0))
            seq = self.sequence[i][p: p + read_length]
            # if the sequence isn't long enough, prepend 'A' nucleotides to fill length
            if len(seq) < read_length:
                seq = (read_length - len(seq)) * 'A' + seq
            sequences.append(seq)
        return sequences

    @classmethod
    def prepare_chr19(cls):
        """prepare a pickled chr19file for quick data synthesis

        Additionally, prepares a chr19 fasta, gtf, and cdna file for use with disambiguate
        """

        prefix = '/'.join(seqc.__file__.split('/')[:-2]) + '/data/genome/mm38_chr19/'
        if all(os.path.isfile(prefix + f) for f in ['mm38_chr19_cdna.fa', 'chr19seq.p',
                                                    'mm38_chr19.gtf']):
            return  # all necessary files already present

        # get cDNA file
        link = ('ftp://ftp.ensembl.org/pub/release-76/fasta/mus_musculus/cdna/'
                'Mus_musculus.GRCm38.cdna.all.fa.gz')


        # create prefix directory if it does not exist
        if not os.path.isdir(prefix):
            os.makedirs(prefix)

        # get the cdna file
        ip, *path, fa_name = link.split('/')[2:]  # [2:] -- eliminate leading 'ftp://'
        path = '/'.join(path)

        # check if file already exists
        if os.path.isfile(prefix + fa_name):
            pass  # file is already present, nothing to do.
        else:
            ftp = ftplib.FTP(ip)
            try:
                ftp.login()
                ftp.cwd(path)
                with open(prefix + fa_name, 'wb') as fout:
                    ftp.retrbinary('RETR %s' % fa_name, fout.write)
            finally:
                ftp.close()

        # get the gtf file
        link = ('ftp://ftp.ensembl.org/pub/release-76/gtf/mus_musculus/'
                'Mus_musculus.GRCm38.76.gtf.gz')
        ip, *path, gtf_name = link.split('/')[2:]  # [2:] -- eliminate leading 'ftp://'
        path = '/'.join(path)

        if os.path.isfile(prefix + gtf_name):
            pass  # file is already present, nothing to do.
        else:
            ftp = ftplib.FTP(ip)
            try:
                ftp.login()
                ftp.cwd(path)
                with open(prefix + gtf_name, 'wb') as fout:
                    ftp.retrbinary('RETR %s' % gtf_name, fout.write)
            finally:
                ftp.close()

        # get the genomic fasta data
        link = ('ftp://ftp.ensembl.org/pub/release-76/fasta/mus_musculus/dna/'
                'Mus_musculus.GRCm38.dna.primary_assembly.fa.gz')
        ip, *path, genome_name = link.split('/')[2:]  # [2:] -- eliminate leading 'ftp://'
        path = '/'.join(path)

        if os.path.isfile(prefix + genome_name):
            pass  # file is already present, nothing to do.
        else:
            ftp = ftplib.FTP(ip)
            try:
                ftp.login()
                ftp.cwd(path)
                with open(prefix + genome_name, 'wb') as fout:
                    ftp.retrbinary('RETR %s' % genome_name, fout.write)
            finally:
                ftp.close()

        # read the complete gtf
        with gzip.open(prefix + gtf_name, 'rt') as f:
            data = f.readlines()

        # get the chr19 specific data and write it to file for coalignment probabilities
        with open(prefix + 'mm38_chr19.gtf', 'w') as f:
            chr19 = [line for line in data if line.split('\t')[0] == '19']
            f.write(''.join(chr19))

        # get the chr19 transcript names to reduce the fasta file
        chr19_transcripts = []
        pattern = r'(.*?; transcript_id ")(.*?)("; .*)'
        for line in chr19:
            if line.split('\t')[2] == 'transcript':
                tx_id = re.match(pattern, line).group(2)
                chr19_transcripts.append(tx_id)
        chr19_transcripts = set(chr19_transcripts)

        # reduce the fasta file for coalignment probabilities
        with gzip.open(prefix + fa_name, 'rt') as fin:
            fa_data = fin.read().strip('>').split('>')

        # track relevant fasta records
        chr19_fa_records = []
        with open(prefix + 'mm38_chr19_cdna.fa', 'w') as fout:
            fout.write('>')
            for record in fa_data:
                if record.split()[0] in chr19_transcripts:
                    chr19_fa_records.append(''.join(record.split('\n')[1:]))
                    fout.write('>' + record)

        # finally, pickle a reduced fasta for synthetic data generation
        records = dict(zip(range(len(chr19_fa_records)), chr19_fa_records))

        # dump the processed file
        with open(prefix + 'chr19seq.p', 'wb') as f:
            pickle.dump(records, f)

    @staticmethod
    def prepare_barcodes(fout, barcode_files, spacer=None):
        barcodes = []
        for file_ in barcode_files:
            with open(file_, 'r') as f:
                barcodes.append([l.strip() for l in f.readlines()])

        if len(barcodes) == 1:
            barcodes = set(barcodes)
            with open(fout, 'wb') as f:
                pickle.dump(barcodes, f)
        else:
            try:
                barcodes = set(spacer.join(pair) for pair in product(*barcodes))
            except AttributeError:
                raise ValueError('must pass spacer argument for experiments with multiple'
                                 'barcodes')
            with open(fout, 'wb') as f:
                pickle.dump(barcodes, f)


@unittest.skip('')
class TestGenerateFastq(unittest.TestCase):

    def setUp(self):
        self.data_dir = '/'.join(seqc.__file__.split('/')[:-2]) + '/data/'

    def test_generate_seqs(self):
        gfq = GenerateFastq()
        sequences = gfq.generate_sequence_records(5, 100)
        self.assertTrue(all(len(s) == 100 for s in sequences))

    def test_generate_reverse_fastq_stringio(self):
        gfq = GenerateFastq()
        fq_records = gfq.generate_reverse_fastq(5, 100)
        data = fq_records.readlines()
        self.assertTrue(len(data) == 20)

    def test_generate_forward_fastq_stringio(self):
        gfq = GenerateFastq()
        fq_records = gfq.generate_forward_fastq(5, 100, self.data_dir +
                                                'in_drop/barcodes/in_drop_barcodes.p',
                                                6)
        data = fq_records.readlines()
        self.assertTrue(len(data) == 20)
        self.assertTrue(all(l.startswith('@') for l in data[::4]))
        self.assertTrue(all(l == '+\n' for l in data[2::4]))


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


@unittest.skip('')
class TestThreeBitCellBarcodes(unittest.TestCase):

    def setUp(self):
        self.barcode_dir = ('/'.join(seqc.__file__.split('/')[:-2]) +
                            '/data/in_drop/barcodes/')

    @unittest.skip('creation of the cb_3bit.p object is slow')
    def test_create_cb(self):

        # create barcode files
        barcode_files = [self.barcode_dir + 'cb1.txt', self.barcode_dir + 'cb2.txt']
        cb = three_bit.CellBarcodes(*barcode_files)
        with open(self.barcode_dir + 'cb_3bit.p', 'wb') as f:
            pickle.dump(cb, f)

    def test_created_cb(self):
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

    def test_error_code_construction_and_recognition(self):
        cb1 = StringIO('AAA\nTTT\n')
        cb2 = StringIO('CCC\nGGG\n')
        cb = three_bit.CellBarcodes(cb1, cb2)
        perfect_code = three_bit.ThreeBit.str2bin('AAACCC')
        error_code = three_bit.ThreeBit.str2bin('NTTCCC')
        self.assertTrue(cb.perfect_match(perfect_code))
        self.assertFalse(cb.perfect_match(error_code))
        self.assertTrue(cb.close_match(perfect_code))
        self.assertTrue(cb.close_match(error_code))

    # def test_error_code_different_barcode_set(self):
    #     cb1 = StringIO('AAA\nTTT\n')
    #     cb = three_bit.CellBarcodes('in-drop', cb1, cb2)

    def test_error_code_error_type_identification(self):
        cb1 = StringIO('AAA\nTTT\n')
        cb2 = StringIO('CCC\nGGG\n')
        cb = three_bit.CellBarcodes(cb1, cb2)
        error_code = three_bit.ThreeBit.str2bin('NTTCCC')

        errors = cb.map_errors(error_code)
        self.assertEqual({'TN'}, set(errors))


@unittest.skip('')
class TestFastq(unittest.TestCase):

    def setUp(self):
        data_dir = '/'.join(seqc.__file__.split('/')[:-2]) + '/data/'
        self.gfq = GenerateFastq()
        self.fwd_len = 50
        self.rev_len = 100
        self.barcodes = data_dir + 'in_drop/barcodes/in_drop_barcodes.p'
        self.umi_len = 6
        self.temp_dir = data_dir + 'test_seqc_merge_fastq'
        self.exp_type = 'in-drop'
        if not os.path.isdir(self.temp_dir):
            os.makedirs(self.temp_dir)
        self.cell_barcode_pickle = data_dir + 'in_drop/barcodes/cb_3bit.p'

    def test_remove_homopolymer(self):
        non_polymer = 'CGTACGATCGATAGCTAG'
        testseq = ('A' * 9) + non_polymer + ('T' * 12)
        testrecord = ['@name', testseq, '+name', 'I' * len(testseq)]
        r, trimmed_bases = fastq.remove_homopolymer(testrecord)
        self.assertEqual(non_polymer + '\n', r[1])
        self.assertTrue(len(non_polymer + '\n') == len(r[3]))
        self.assertEqual(trimmed_bases, 21)

    def test_process_record(self):
        n = 10000
        _ = self.gfq.generate_forward_fastq(
            n, self.fwd_len, self.barcodes, self.umi_len)
        _ = self.gfq.generate_reverse_fastq(n, self.rev_len)

    def test_merge_record(self):
        tbp = three_bit.ThreeBit.default_processors(self.exp_type)
        # cb = CellBarcodes(processor, *cell_barcode_files)
        with open(self.cell_barcode_pickle, 'rb') as f:
            cb = pickle.load(f)
        n = 1
        forward = self.gfq.generate_forward_fastq(
            n, self.fwd_len, self.barcodes, self.umi_len).readlines()
        reverse = self.gfq.generate_reverse_fastq(n, self.rev_len).readlines()
        fq = fastq.process_record(forward, reverse, tbp, cb)
        self.assertEqual(len(fq.split()), 4)

    @unittest.skip('Takes a long time')
    def test_merge_fastq(self):
        # takes approximately 6.35 minutes per million reads, of which 37% is
        # accorded to the estimation of sequence quality
        n = 10000
        self.forward = self.gfq.generate_forward_fastq(
            n, self.fwd_len, self.barcodes, self.umi_len)
        self.reverse = self.gfq.generate_reverse_fastq(n, self.rev_len)
        _ = fastq.merge_fastq(
            self.forward, self.reverse, self.exp_type, self.temp_dir,
            self.cell_barcode_pickle)

    @unittest.skip('this is not currently working.')
    def test_merge_fastq_threaded(self):
        n = 10000
        n_threads = 7
        self.forward = self.gfq.generate_forward_fastq(
            n, self.fwd_len, self.barcodes, self.umi_len)
        self.reverse = self.gfq.generate_reverse_fastq(n, self.rev_len)

        _ = fastq.merge_fastq_threaded(
            self.forward, self.reverse, n_threads, self.exp_type, self.temp_dir,
            self.cell_barcode_pickle)

    def tearDown(self):
        # eventually, should remove temp_dir
        pass


@unittest.skip('')
class TestTranslateFeature(unittest.TestCase):
    """
    I'm reasonably confident that this is working as intended, but the tests in this suite
    are not strong enough. Specifically, there is no test to determine if multialignments
    are working properly, and the data synthesis process is not working as it should.
    instead of looking at transcripts, we should look at exons for data generation, but
    need to constrain our search only to exons within fragment_length (1000bp) of the
    3' end of each transcript
    """

    def setUp(self):
        self.genome_dir = ('/'.join(seqc.__file__.split('/')[:-2]) +
                           '/data/genome/mm38_chr19/')
        self.gtf = self.genome_dir + 'annotations.gtf'
        self.fragment_len = 1000

    @unittest.skip('run only when new table and positions must be generated')
    def test_construct_tables(self):
        """saves speed for repeat testing"""
        ft, fp = convert_features.construct_feature_table(self.gtf, self.fragment_len)
        with open(self.genome_dir + 'feature_table_and_positions.pckl', 'wb') as f:
            pickle.dump((ft, fp), f)

    # @unittest.skip('')
    def test_translate_known_feature(self):
        with open(self.genome_dir + 'feature_table_and_positions.pckl', 'rb') as f:
            ft, fp = pickle.load(f)

        # pick a few random features
        with open(self.gtf, 'r') as f:
            gtf_data = [l.strip().split('\t') for l in f.readlines()]
            gtf_data = [l for l in gtf_data if l[2] == 'transcript']

        # generate some alignments to that feature
        pattern = re.compile(r'(^.*?scseq_id "SC)(.*?)(";.*?$)')
        random_features = [random.choice(gtf_data) for _ in range(1000)]
        fake_alignments = []
        for f in random_features:
            reference = f[0]
            strand = f[6]
            if strand == '+':
                end = int(f[4])
                # must be consistent with fragment length
                pos = random.randint(end - self.fragment_len, end)
            else:
                end = int(f[3])
                pos = np.random.randint(end, end + self.fragment_len)

            scid = int(re.match(pattern, f[-1]).group(2))
            fake_alignments.append(((reference, strand, pos, ft, fp), scid))

        # test if translate_feature recovers them
        n_successful = 0
        for a, scid in fake_alignments:
            try:
                pred_scid = sam.translate_feature(*a)
            except KeyError:
                pred_scid = 0
            if scid == pred_scid:
                n_successful += 1
        message = 'fails because not all positions in the transcript are found in exons'
        self.assertEqual(n_successful, 1000, message)


@unittest.skip('')
class TestAlign(unittest.TestCase):

    def setUp(self):
        data_dir = '/'.join(seqc.__file__.split('/')[:-2]) + '/data/'
        self.temp_dir = data_dir + 'test_seqc_merge_fastq/'
        self.fastq = data_dir + 'test_seqc_merge_fastq/merged_temp.fastq'
        self.index = data_dir + 'genome/mm38_chr19/'

    def test_align(self):
        star = align.STAR(self.temp_dir, 7, self.index, 'mm38')
        star.align(self.fastq)


@unittest.skip('')
class TestSEQC(unittest.TestCase):

    def test_set_up(self):
        exp_name = 'test_set_up'
        temp_dir = '.' + exp_name
        if not os.path.isdir(temp_dir):
            os.makedirs(temp_dir)
        if not temp_dir.endswith('/'):
            temp_dir += '/'

        self.assertEqual('.%s/' % exp_name, temp_dir)
        os.rmdir('.%s/' % exp_name)


@unittest.skip('')
class TestSamProcessing(unittest.TestCase):

    def setUp(self):
        data_dir = '/'.join(seqc.__file__.split('/')[:-2]) + '/data/'
        self.samfile = data_dir + 'test_seqc_merge_fastq/Aligned.out.sam'
        self.gtf = data_dir + 'genome/mm38_chr19/annotations.gtf'

    def test_process_alignments(self):
        n_threads = 4
        arr = sam.process_alignments(self.samfile, n_threads, self.gtf, fragment_len=1000)
        print(arr)


@unittest.skip('')
class TestPeekable(unittest.TestCase):

    def test_peekable_string(self):
        iterable = 'ACGTACGT'
        pk = sam.Peekable(iterable)
        self.assertEqual(pk.peek, 'A')
        first = next(pk)
        self.assertEqual(first, 'A')

    def test_peekable_fileobj(self):
        fileobj = StringIO('1\n2\n3\n4\n')
        pk = sam.Peekable(fileobj)
        first = next(pk)
        self.assertTrue(first == '1\n')
        self.assertTrue(pk.peek == '2\n')


@unittest.skip('')
class TestDisambiguation(unittest.TestCase):

    def setUp(self):
        self.data_dir = '/'.join(seqc.__file__.split('/')[:-2]) + '/data/'
        self.expectations = self.data_dir + 'genome/mm38_chr19/p_coalignment.pckl'

    @unittest.skip('only run this if new disambiguation input is needed')
    def test_create_test_data(self):
        cell_barcodes = self.data_dir + 'in_drop/barcodes/cb_3bit.p'
        save = self.data_dir + 'in_drop/disambiguation_input.p'
        self.data = generate_in_drop_disambiguation_data(
            self.expectations, cell_barcodes, 1000, 10, save=save)

    # todo this has a weak test for success. Fix to check specific disambiguation results
    def test_disambiguation(self):
        arr_pickle = self.data_dir + 'in_drop/disambiguation_input.p'
        with open(arr_pickle, 'rb') as f:
            arr = pickle.load(f)
        res, data = qc.disambiguate(arr, self.expectations)
        self.assertTrue(~np.all(data == 0))


@unittest.skip('')
class TestErrorCorrection(unittest.TestCase):

    def setUp(self):
        self.data_dir = '/'.join(seqc.__file__.split('/')[:-2]) + '/data/'
        self.expectations = self.data_dir + 'genome/mm38_chr19/p_coalignment.pckl'

    def test_create_cell_barcodes(self):
        cell_barcodes = self.data_dir + 'in_drop/barcodes/cb_3bit.p'
        save = self.data_dir + 'in_drop/disambiguation_input.p'
        data = generate_in_drop_disambiguation_data(
            self.expectations, cell_barcodes, 10000, 2, save=save)

        # mutate the first base to 'N' in 5% of cases
        for i in range(data.shape[0]):
            if np.random.uniform(0, 1) < 0.05:
                cell = data['cell'][i]
                data['cell'][i] = cell & 0b111

        res, err_rate = qc.correct_errors(data, cell_barcodes, 1, 0.2)
        self.assertTrue(~np.all(res == 0))
        print(err_rate)

    def test_error_correction(self):
        pass


if __name__ == '__main__':
    unittest.main()