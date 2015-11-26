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
import threading
from pyftpdlib.handlers import FTPHandler
from pyftpdlib.servers import FTPServer
from pyftpdlib.authorizers import DummyAuthorizer


# this is the universal data dir for these tests
_seqc_dir = '/'.join(seqc.__file__.split('/')[:-3]) + '/'
_data_dir = '/'.join(seqc.__file__.split('/')[:-2]) + '/data/'


class SRAGenerator:

    @classmethod
    def make_run_xml(cls, forward_fastq_files, experiment_alias, forward_checksum_results,
                     fout_stem=None, data_block=None, reverse_fastq_files=None,
                     reverse_checksum_results=None):
        """data_block should be a member name if the experiment is a pooled experiment,
        and this run is a demultiplexed member"""

        if not fout_stem:
            fout_stem = experiment_alias + '_run.xml'

        if reverse_fastq_files:
            input_files = [forward_fastq_files, reverse_fastq_files,
                           forward_checksum_results, reverse_checksum_results]
        else:
            input_files = [forward_fastq_files, forward_checksum_results]
        if not len(set(len(f) for f in input_files)) == 1:
            raise ValueError('Input files must be of equal length')
        n = len(input_files[0])

        run_set = ET.Element('RUN_SET')
        run_set.set('xmlns:xsi', 'http://www.w3.org/2001/XMLSchema-instance')
        run_set.set('xsi:noNamespaceSchemaLocation',
                    'ftp://ftp.sra.ebi.ac.uk/meta/xsd/sra_1_5/SRA.run.xsd')
        run = ET.SubElement(run_set, 'RUN', alias='RUN NAME',
                            center_name='Columbia University')
        experiment_ref = ET.SubElement(run, 'EXPERIMENT_REF')
        experiment_ref.set('alias', experiment_alias)

        if data_block:
            data_block_field = ET.SubElement(run, 'DATA_BLOCK')
            data_block_field.set('member_name', data_block)

        files = ET.SubElement(run, 'FILES')
        for i in range(n):
            forward_file = ET.SubElement(files, 'FILE')
            forward_file.set('filename', forward_fastq_files[i])
            forward_file.set('filetype', 'fastq')
            forward_file.set('checksum_method', 'MD5')
            forward_file.set('checksum', forward_checksum_results[i])
            if reverse_fastq_files:
                reverse_file = ET.SubElement(files, 'FILE')
                reverse_file.set('filename', reverse_fastq_files[i])
                reverse_file.set('filetype', 'fastq')
                reverse_file.set('checksum_method', 'MD5')
                reverse_file.set('checksum', reverse_checksum_results[i])
        tree = ET.ElementTree(run_set)
        tree.write(fout_stem + '_run.xml', method='xml')
        return fout_stem + '_run.xml'

    @classmethod
    def make_experiment_xml(
            cls, reference_alias, sample_alias, platform_name='ILLUMINA',
            single_or_paired_end='SINGLE', instrument_model='Illumina HiSeq 2500',
            fout_stem='SRA'):

        if not single_or_paired_end in ['SINGLE', 'PAIRED']:
            raise ValueError('single_or_paired_end must be one of "SINGLE" or "PAIRED"')

        valid_platform_names = ['LS454', 'ILLUMINA', 'COMPLETE_GENOMICS', 'PACBIO_SMRT',
                                'ION_TORRNET', 'OXFORD_NANOPORE', 'CAPILLARY']
        if not platform_name in valid_platform_names:
            raise ValueError('platform_name must be one of %s' %
                             repr(valid_platform_names))

        exp_set = ET.Element('EXPERIMENT_SET')
        exp_set.set('xmlns:xsi', 'http://www.w3.org/2001/XMLSchema-instance')
        exp_set.set('xsi:noNamespaceSchemaLocation',
                    'ftp://ftp.sra.ebi.ac.uk/meta/xsd/sra_1_5/SRA.experiment.xsd')
        exp = ET.SubElement(exp_set, 'EXPERIMENT')
        exp.set("alias", "DUMMY EXPERIMENT FOR TESTING")
        exp.set("center_name", "Columbia University")
        title = ET.SubElement(exp, 'TITLE')
        title.text = 'EXPERIMENT TITLE'
        study_ref = ET.SubElement(exp, 'STUDY_REF')
        study_ref.set('refname', '%s' % reference_alias)
        design = ET.SubElement(exp, 'DESIGN')
        design_description = ET.SubElement(design, 'DESIGN_DESCRIPTION')
        design_description.text = 'DETAILS ABOUT SETUP AND GOALS'
        sample_descriptor = ET.SubElement(design, 'SAMPLE_DESCRIPTOR')
        sample_descriptor.set('refname', '%s' % sample_alias)
        library_descriptor = ET.SubElement(design, 'LIBRARY_DESCRIPTOR')
        library_name = ET.SubElement(library_descriptor, 'LIBRARY_NAME')
        library_name.text = 'DUMMY NAME'
        library_strategy = ET.SubElement(library_descriptor, 'LIBRARY_STRATEGY')
        library_strategy.text = 'RNA-Seq'
        library_source = ET.SubElement(library_descriptor, 'LIBRARY_SOURCE')
        library_source.text = 'TRANSCRIPTOMIC'
        library_selection = ET.SubElement(library_descriptor, 'LIBRARY_SELECTION')
        library_selection.text = 'Oligo-dT'
        library_layout = ET.SubElement(library_descriptor, 'LIBRARY_LAYOUT')
        library_layout.text = single_or_paired_end
        platform = ET.SubElement(exp, 'PLATFORM')
        specific_platform = ET.SubElement(platform, platform_name)
        instrument = ET.SubElement(specific_platform, 'INSTRUMENT_MODEL')
        instrument.text = instrument_model
        tree = ET.ElementTree(exp_set)
        tree.write(fout_stem + '_experiment.xml', method='xml')
        return fout_stem + '_experiment.xml'

    # @classmethod
    # def create_xml_for_fastq_data(cls, forward_fastq, reverse_fastq=None):
    #     if reverse_fastq:
    #         single_or_paired_end = 'PAIRED'
    #     else:
    #         single_or_paired_end = 'SINGLE'

    @staticmethod
    def md5sum(fname):

        def hashfile(afile, hasher, blocksize=65536):
            buf = afile.read(blocksize)
            while len(buf) > 0:
                hasher.update(buf)
                buf = afile.read(blocksize)
            return hasher.digest()

        return hashfile(open(fname, 'rb'), hashlib.md5())

    @staticmethod
    def fastq_load(file_directory, run_xml, experiment_xml, output_path):
        cmd = ['fastq-load', '-r', run_xml, '-e', experiment_xml, '-o', output_path, '-i',
               file_directory]
        p = Popen(cmd, stderr=PIPE)
        _, err = p.communicate()
        if err:
            raise ChildProcessError(err)




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
    is_aligned = True

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
    arrays_ = []

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
                          fwd_quality, sam.ObfuscatedTuple(f),
                          sam.ObfuscatedTuple(tuple([0] * len(f))), is_aligned,
                          alignment_score)
                i += 1
        arrays_.append(arr)

    arrays_ = np.hstack(arrays_)

    if isinstance(save, str):
        with open(save, 'wb') as f:
            pickle.dump(arrays_, f)

    return arrays_


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
        reverse_fastq = StringIO('\n'.join(records) + '\n')
        return reverse_fastq

    @staticmethod
    def generate_forward_in_drop_fastq(n, read_length, barcodes_, umi_len):
        with open(barcodes_, 'rb') as f:
            barcodes_ = list(pickle.load(f))
        names = range(n)
        name2 = '+'
        quality = 'I' * read_length
        records = []
        alphabet = np.array(['A', 'C', 'G', 'T'])
        for name in names:
            cb = random.choice(barcodes_)
            umi = ''.join(np.random.choice(alphabet, umi_len))
            poly_a = (read_length - len(cb) - len(umi)) * 'T'
            records.append('\n'.join(['@%d' % name, cb + umi + poly_a, name2, quality]))
        forward_fastq = StringIO('\n'.join(records) + '\n')
        return forward_fastq

    @staticmethod
    def generate_forward_drop_seq_fastq(n):
        names = range(n)
        name2 = '+'
        quality = 'I' * 20
        records = []
        alphabet = np.array(['A', 'C', 'G', 'T'])
        for name in names:
            cb = ''.join(np.random.choice(alphabet, 12))
            umi = ''.join(np.random.choice(alphabet, 8))
            records.append('\n'.join(['@%d' % name, cb + umi, name2, quality]))
        forward_fastq = StringIO('\n'.join(records) + '\n')
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

        Additionally, prepares a chr19 fasta, gtf, and cdna file for use with resolve_alignments
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
        barcodes_ = []
        for file_ in barcode_files:
            with open(file_, 'r') as f:
                barcodes_.append([l.strip() for l in f.readlines()])

        if len(barcodes_) == 1:
            barcodes_ = set(barcodes_)
            with open(fout, 'wb') as f:
                pickle.dump(barcodes_, f)
        else:
            try:
                barcodes_ = set(spacer.join(pair) for pair in product(*barcodes_))
            except AttributeError:
                raise ValueError('must pass spacer argument for experiments with multiple'
                                 'barcodes')
            with open(fout, 'wb') as f:
                pickle.dump(barcodes_, f)


class DummyFTPClient(threading.Thread):
    """A threaded FTP server used for running tests.

    This is basically a modified version of the FTPServer class which
    wraps the polling loop into a thread.

    The instance returned can be used to start(), stop() and
    eventually re-start() the server.

    The instance can also launch a client using ftplib to navigate and download files.
    it will serve files from home.
    """
    handler = FTPHandler
    server_class = FTPServer

    def __init__(self, addr=None, home=None):

        try:
            host = socket.gethostbyname('localhost')
        except socket.error:
            host = 'localhost'

        threading.Thread.__init__(self)
        self.__serving = False
        self.__stopped = False
        self.__lock = threading.Lock()
        self.__flag = threading.Event()
        if addr is None:
            addr = (host, 0)

        if not home:
            home = os.getcwd()

        authorizer = DummyAuthorizer()
        authorizer.add_anonymous(home, perm='erl')
        # authorizer.add_anonymous(home, perm='elr')
        self.handler.authorizer = authorizer
        # lower buffer sizes = more "loops" while transfering data
        # = less false positives
        self.handler.dtp_handler.ac_in_buffer_size = 4096
        self.handler.dtp_handler.ac_out_buffer_size = 4096
        self.server = self.server_class(addr, self.handler)
        self.host, self.port = self.server.socket.getsockname()[:2]
        self.client = None

    def __repr__(self):
        status = [self.__class__.__module__ + "." + self.__class__.__name__]
        if self.__serving:
            status.append('active')
        else:
            status.append('inactive')
        status.append('%s:%s' % self.server.socket.getsockname()[:2])
        return '<%s at %#x>' % (' '.join(status), id(self))

    def generate_local_client(self):
        self.client = ftplib.FTP()
        self.client.connect(self.host, self.port)
        self.client.login()
        return self.client

    @property
    def running(self):
        return self.__serving

    def start(self, timeout=0.001):
        """Start serving until an explicit stop() request.
        Polls for shutdown every 'timeout' seconds.
        """
        if self.__serving:
            raise RuntimeError("Server already started")
        if self.__stopped:
            # ensure the server can be started again
            DummyFTPClient.__init__(self, self.server.socket.getsockname(), self.handler)
        self.__timeout = timeout
        threading.Thread.start(self)
        self.__flag.wait()

    def run(self):
        logging.basicConfig(filename='testing.log', level=logging.DEBUG)
        self.__serving = True
        self.__flag.set()
        while self.__serving:
            self.__lock.acquire()
            self.server.serve_forever(timeout=self.__timeout, blocking=False)
            self.__lock.release()
        self.server.close_all()

    def stop(self):
        """Stop serving (also disconnecting all currently connected
        clients) by telling the serve_forever() loop to stop and
        waits until it does.
        """
        if not self.__serving:
            raise RuntimeError("Server not started yet")
        self.__serving = False
        self.__stopped = True
        self.join()
        self.client.close()


@unittest.skip('')
class TestJaggedArray(unittest.TestCase):

    def generate_input_iterable(self, n):
        for i in range(n):
            yield [random.randint(0, 5) for _ in range(random.randint(0, 5))]

    def test_jagged_array(self):
        n = int(1e6)
        data = list(self.generate_input_iterable(n))
        data_size = sum(len(i) for i in data)
        jarr = arrays.JaggedArray.from_iterable(data_size, data)
        self.assertTrue(jarr._data.dtype == np.uint32)
        self.assertTrue(jarr._data.shape == (data_size,))
        print(jarr[10])


@unittest.skip('')
class TestGenerateFastq(unittest.TestCase):

    def setUp(self):
        self.data_dir = '/'.join(seqc.__file__.split('/')[:-2]) + '/data/'
        self.in_drop_barcodes = (
            self.data_dir + 'in_drop/barcodes/concatenated_string_in_drop_barcodes.p')
        if not os.path.isfile(self.in_drop_barcodes):
            cb1 = self.data_dir + 'in_drop/barcodes/cb1.txt'
            cb2 = self.data_dir + 'in_drop/barcodes/cb2.txt'
            GenerateFastq.prepare_barcodes(
                self.in_drop_barcodes, [cb1, cb2], spacer='GAGTGATTGCTTGTGACGCCTT')

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
        fq_records = gfq.generate_forward_in_drop_fastq(
            5, 100, self.in_drop_barcodes, 6)
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


@unittest.skip('')
class TestFastq(unittest.TestCase):

    def setUp(self):
        data_dir = '/'.join(seqc.__file__.split('/')[:-2]) + '/data/'
        self.gfq = GenerateFastq()
        self.fwd_len = 50
        self.rev_len = 100
        self.in_drop_string_barcodes = (
            data_dir + 'in_drop/barcodes/concatenated_string_in_drop_barcodes.p')
        self.umi_len = 6
        self.in_drop_temp_dir = data_dir + 'test_seqc_merge_in_drop_fastq'
        self.drop_seq_temp_dir = data_dir + 'test_seqc_merge_drop_seq_fastq'
        self.in_drop_processor = 'in-drop'
        self.drop_seq_processor = 'drop-seq'
        for d in [self.drop_seq_temp_dir, self.in_drop_temp_dir]:
            if not os.path.isdir(d):
                os.makedirs(d)
        self.in_drop_cell_barcode_pickle = data_dir + 'in_drop/barcodes/cb_3bit.p'
        self.drop_seq_cell_barcode_pickle = data_dir + 'drop_seq/barcodes/cb_3bit.p'

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
        _ = self.gfq.generate_forward_in_drop_fastq(
            n, self.fwd_len, self.in_drop_string_barcodes, self.umi_len)
        _ = self.gfq.generate_reverse_fastq(n, self.rev_len)

    def test_merge_in_drop_record(self):
        tbp = three_bit.ThreeBit.default_processors(self.in_drop_processor)
        with open(self.in_drop_cell_barcode_pickle, 'rb') as f:
            cb = pickle.load(f)
        n = 1
        forward = self.gfq.generate_forward_in_drop_fastq(
            n, self.fwd_len, self.in_drop_string_barcodes, self.umi_len).readlines()
        reverse = self.gfq.generate_reverse_fastq(n, self.rev_len).readlines()
        fq = fastq.process_record(forward, reverse, tbp, cb)
        self.assertEqual(len(fq.split()), 4)

    def test_merge_drop_seq_record(self):
        tbp = three_bit.ThreeBit.default_processors(self.drop_seq_processor)
        with open(self.drop_seq_cell_barcode_pickle, 'rb') as f:
            cb = pickle.load(f)
        n = 1
        forward = self.gfq.generate_forward_drop_seq_fastq(n).readlines()
        reverse = self.gfq.generate_reverse_fastq(n, self.rev_len).readlines()
        fq = fastq.process_record(forward, reverse, tbp, cb)
        self.assertEqual(len(fq.split()), 4)

    # @unittest.skip('Takes a long time')
    def test_merge_in_drop_fastq(self):
        # takes approximately 6.35 minutes per million reads, of which 37% is
        # accorded to the estimation of sequence quality
        n = 10000
        forward = self.gfq.generate_forward_in_drop_fastq(
            n, self.fwd_len, self.in_drop_string_barcodes, self.umi_len)
        reverse = self.gfq.generate_reverse_fastq(n, self.rev_len)
        _ = fastq.merge_fastq(
            [forward], [reverse], self.in_drop_processor, self.in_drop_temp_dir,
            self.in_drop_cell_barcode_pickle)

    # @unittest.skip('Takes a long time')
    def test_merge_drop_seq_fastq(self):
        # takes approximately 6.35 minutes per million reads, of which 37% is
        # accorded to the estimation of sequence quality
        n = 10000
        forward = self.gfq.generate_forward_drop_seq_fastq(n)
        reverse = self.gfq.generate_reverse_fastq(n, self.rev_len)
        _ = fastq.merge_fastq(
            [forward], [reverse], self.drop_seq_processor, self.drop_seq_temp_dir,
            self.drop_seq_cell_barcode_pickle)

    @unittest.skip('this is not currently working.')
    def test_merge_in_drop_fastq_threaded(self):
        n = 10000
        n_threads = 7
        forward = self.gfq.generate_forward_in_drop_fastq(
            n, self.fwd_len, self.in_drop_string_barcodes, self.umi_len)
        reverse = self.gfq.generate_reverse_fastq(n, self.rev_len)

        _ = fastq.merge_fastq_threaded(
            forward, reverse, n_threads, self.in_drop_processor, self.in_drop_temp_dir,
            self.in_drop_cell_barcode_pickle)


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
        self.in_drop_temp_dir = data_dir + 'test_seqc_merge_in_drop_fastq/'
        self.drop_seq_temp_dir = data_dir + 'test_seqc_merge_drop_seq_fastq/'
        self.in_drop_fastq = self.in_drop_temp_dir + 'merged_temp.fastq'
        self.drop_seq_fastq = self.drop_seq_temp_dir + 'merged_temp.fastq'
        self.index = data_dir + 'genome/mm38_chr19/'

    def test_align_in_drop(self):
        star = align.STAR(self.in_drop_temp_dir, 7, self.index, 'mm38')
        samfile = star.align(self.in_drop_fastq)
        self.assertEqual(samfile, self.in_drop_temp_dir + 'Aligned.out.sam')
ew
    def test_align_drop_seq(self):
        star = align.STAR(self.drop_seq_temp_dir, 7, self.index, 'mm38')
        samfile = star.align(self.drop_seq_fastq)
        self.assertEqual(samfile, self.drop_seq_temp_dir + 'Aligned.out.sam')

    def test_align_drop_seq_multiple_files(self):
        star = align.STAR(self.drop_seq_temp_dir, 7, self.index, 'mm38')
        samfiles = star.align_multiple_files([self.drop_seq_fastq, self.drop_seq_fastq,
                                              self.drop_seq_fastq])
        print(samfiles)


# todo import these tests from scseq/seqdb
@unittest.skip('')
class TestIndexGeneration(unittest.TestCase):
    pass


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
        self.in_drop_samfile = data_dir + 'test_seqc_merge_in_drop_fastq/Aligned.out.sam'
        self.drop_seq_samfile = (
            data_dir + 'test_seqc_merge_drop_seq_fastq/Aligned.out.sam')
        self.gtf = data_dir + 'genome/mm38_chr19/annotations.gtf'

    def test_process_in_drop_alignments(self):
        n_threads = 4
        arr = sam.process_alignments(self.in_drop_samfile, n_threads, self.gtf,
                                     fragment_len=1000)
        print(arr)

    def test_process_drop_seq_alignments(self):
        n_threads = 4
        arr = sam.process_alignments(self.drop_seq_samfile, n_threads, self.gtf,
                                     fragment_len=1000)
        print(arr)


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
        ft, fp = convert_features.construct_feature_table(self.gtf, 1000)
        res = arrays.ReadArray.from_samfile(samfile, ft, fp)
        res.save_h5(self.h5name)

    def test_profile_counts_creation(self):
        ra = arrays.ReadArray.from_h5(self.h5name)
        pr = cProfile.Profile()
        pr.enable()
        ra.to_sparse_counts(collapse_molecules=True, n_poly_t_required=0)
        pr.disable()
        p = Stats(pr)
        p.strip_dirs()
        p.sort_stats('cumtime')
        p.print_stats()

    # def tearDown(self):
    #     if os.path.isfile(self.h5name):
    #         os.remove(self.h5name)


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
class TestUnionFind(unittest.TestCase):

    def setUp(self):
        # get all combinations of single and multiple features with and without overlaps
        self.f1 = sam.ObfuscatedTuple((1, 2, 3))  # multiple features, overlaps
        self.f2 = sam.ObfuscatedTuple((1, 2))
        self.f3 = sam.ObfuscatedTuple((1,))  # one feature, overlaps
        self.f4 = sam.ObfuscatedTuple((3, 4))
        self.f5 = sam.ObfuscatedTuple((4,))
        self.f6 = sam.ObfuscatedTuple((6, 7))  # multiple features, no overlaps
        self.f7 = sam.ObfuscatedTuple((11,))  # one feature, no overlaps

    def test_union_find_on_single_feature_group(self):
        """This works with the new object because of its __iter__() method"""
        arr = np.array([self.f7], dtype=np.object)
        uf = qc.UnionFind()
        uf.union_all(arr)
        set_membership, sets = uf.find_all(arr)
        self.assertTrue(np.all(sets == np.array([0])))
        self.assertTrue(np.all(set_membership == np.array([0])))

    def test_union_find_on_multiple_feature_group(self):
        arr = np.array([self.f6], dtype=np.object)
        uf = qc.UnionFind()
        uf.union_all(arr)
        set_membership, sets = uf.find_all(arr)
        self.assertTrue(np.all(sets == np.array([0])))
        self.assertTrue(np.all(set_membership == np.array([0])))

    def test_union_find_on_multiple_feature_group_no_repitition_with_overlaps(self):
        arr = np.array([self.f1, self.f2, self.f3, self.f4, self.f5, self.f6, self.f7],
                       dtype=np.object)
        uf = qc.UnionFind()
        uf.union_all(arr)
        set_membership, sets = uf.find_all(arr)
        self.assertEqual(sorted(sets), [0, 1, 2])
        # there aren't guaranteed set associations -- all will have same number if they're
        # in the same component, but that number could be any of the sets 0, 1, 2
        s1 = set_membership[0]
        s2 = set_membership[-2]
        s3 = set_membership[-1]
        prediction = np.array([s1, s1, s1, s1, s1, s2, s3])
        # only 3 sets should be present, they should be 0, 1, and 2
        self.assertEqual(sorted(np.unique(set_membership)), [0, 1, 2])

        # they should find the right membership

        self.assertTrue(np.all(set_membership == prediction),
                        '%s != %s' % (repr(set_membership), repr(prediction)))

    def test_union_find_on_multiple_feature_group_with_repitition_with_overlaps(self):
        arr = np.array([self.f1, self.f2, self.f3, self.f3, self.f4, self.f4, self.f5,
                        self.f6, self.f7, self.f7], dtype=np.object)
        uf = qc.UnionFind()
        uf.union_all(arr)
        set_membership, sets = uf.find_all(arr)
        self.assertEqual(sorted(sets), [0, 1, 2])
        # there aren't guaranteed set associations -- all will have same number if they're
        # in the same component, but that number could be any of the sets 0, 1, 2
        s1 = set_membership[0]
        s2 = set_membership[-3]
        s3 = set_membership[-1]
        prediction = np.array([s1, s1, s1, s1, s1, s1, s1, s2, s3, s3])
        # only 3 sets should be present, they should be 0, 1, and 2
        self.assertEqual(sorted(np.unique(set_membership)), [0, 1, 2])

        # they should find the right membership

        self.assertTrue(np.all(set_membership == prediction),
                        '%s != %s' % (repr(set_membership), repr(prediction)))


@unittest.skip('')
class TestResolveAlignments(unittest.TestCase):
    """this may need more tests for more sophisticated input data."""

    def setUp(self):
        self.data_dir = '/'.join(seqc.__file__.split('/')[:-2]) + '/data/'
        self.expectations = self.data_dir + 'genome/mm38_chr19/p_coalignment.pckl'

    @unittest.skip('only run this if new disambiguation input is needed')
    def test_create_test_data(self):
        cell_barcodes = self.data_dir + 'in_drop/barcodes/cb_3bit.p'
        save = self.data_dir + 'in_drop/disambiguation_input.p'
        self.data = generate_in_drop_disambiguation_data(
            self.expectations, cell_barcodes, 1000, 10, save=save)

    @unittest.skip('only run this if new ra disambiguation input is needed')
    def test_generate_disambiguation_read_array(self):
        cell_barcodes = self.data_dir + 'in_drop/barcodes/cb_3bit.p'
        save = self.data_dir + 'in_drop/disambiguation_ra_input.p'
        _ = generate_in_drop_read_array(self.expectations, cell_barcodes, 1000, 10,
                                        save=save)

    @unittest.skip('')
    def test_disambiguation(self):
        arr_pickle = self.data_dir + 'in_drop/disambiguation_input.p'
        with open(arr_pickle, 'rb') as f:
            arr = pickle.load(f)
        res, data = qc.disambiguate(arr, self.expectations)
        self.assertTrue(np.all(res == 4))  # 4 == complete disambiguation.

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


@unittest.skip('')
class TestSaveCountsMatrices(unittest.TestCase):

    def setUp(self):
        self.data_dir = '/'.join(seqc.__file__.split('/')[:-2]) + '/data/'

    def test_save_counts_matrices(self):

        with open(self.data_dir + 'in_drop/disambiguation_input.p', 'rb') as f:
            arr = pickle.load(f)
        _ = qc.counts_matrix(arr, True)
        _ = qc.counts_matrix(arr, False)


@unittest.skip('')
class TestGeneTable(unittest.TestCase):

    def setUp(self):
        self.genome_dir = '/'.join(seqc.__file__.split('/')[:-2]) + '/data/genome/'

    def test_gene_table(self):
        gt = convert_features.GeneTable(self.genome_dir + 'mm38_chr19/annotations.gtf')
        chromosome = 'chr19'
        start = 4007800
        end = 4007900
        strand = '-'
        genes = gt.coordinates_to_gene_ids(chromosome, start, end, strand)
        print(genes)


@unittest.skip('')
class TestSamToCount(unittest.TestCase):

    def setUp(self):
        self.data_dir = '/'.join(seqc.__file__.split('/')[:-2]) + '/data/'

    def test_sam_to_count_single_file(self):
        samfile = self.data_dir + 'test_seqc_merge_drop_seq_fastq/Aligned.out.sam'
        gtf = self.data_dir + '/genome/mm38_chr19/annotations.gtf'
        coo, gene_index, ci = qc.sam_to_count_single_file(samfile, gtf, 100)
        print(repr(coo))
        print(len(gene_index))

        samfile = self.data_dir + 'test_seqc_merge_in_drop_fastq/Aligned.out.sam'
        gtf = self.data_dir + '/genome/mm38_chr19/annotations.gtf'
        coo, gene_index, ci = qc.sam_to_count_single_file(samfile, gtf, 100)
        print(repr(coo))
        print(len(gene_index))

    def test_sam_to_count_multiple_files(self):
        samfile = self.data_dir + 'test_seqc_merge_drop_seq_fastq/Aligned.out.sam'
        gtf = self.data_dir + '/genome/mm38_chr19/annotations.gtf'
        coo, gene_index, ci = qc.sam_to_count_multiple_files([samfile, samfile], gtf, 100)
        print(repr(coo))
        print(len(gene_index))

        samfile = self.data_dir + 'test_seqc_merge_in_drop_fastq/Aligned.out.sam'
        gtf = self.data_dir + '/genome/mm38_chr19/annotations.gtf'
        coo, gene_index, ci = qc.sam_to_count_multiple_files([samfile, samfile], gtf, 100)
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


@unittest.skip('')
class TestProcessSingleFileSCSEQExperiment(unittest.TestCase):
    """This is a longer, functional test which will run on 10,000 fastq records"""

    def setUp(self):
        # dummy up some fastq data
        data_dir = '/'.join(seqc.__file__.split('/')[:-2]) + '/data/'
        self.working_directory = data_dir + '.test_process_single_seqc_exp/'
        if not os.path.isdir(self.working_directory):
            os.mkdir(self.working_directory)
        # generate_in_drop_fastq_data(10000, self.cwd + 'testdata')

        # create an SRA file
        # todo this is a nightmare, so I'm skipping this step for now

        # unpack SRA file; will do when create SRA doesn't fail.

        # expect unpacked fastq files to be SRRxxxx_1.fastq and SRRxxxx_2.fastq if
        # paired-end or SRRxxxx.fastq if single-ended

        fastq_stem = '/Users/ambrose/PycharmProjects/SEQC/src/data/in_drop/'
        self.forward = [fastq_stem + 'sample_data_r1.fastq']
        self.reverse = [fastq_stem + 'sample_data_r2.fastq']

        self.index = data_dir + 'genome/mm38_chr19/'

        self.s3_bucket = 'dplab-home'
        self.s3_key = 'ajc2205/test_in_drop.npz'

    @unittest.skip('')
    def test_process_single_file_no_sra_download(self):

        # set some variables
        index = self.index
        working_directory = self.working_directory
        index_bucket = None
        index_key = None
        S3 = io_lib.S3
        STAR = align.STAR
        n_threads = 7
        sam_to_count_single_file = qc.sam_to_count_single_file
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

    # @unittest.skip('')
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
        _ = ('/Users/ambrose/PycharmProjects/SEQC/src/data/in_drop/barcodes/'
             'in_drop_barcodes.p')

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

    def tearDown(self):
        shutil.rmtree(self.working_directory)
        io_lib.S3.remove_file(self.s3_bucket, self.s3_key)
        os.remove('test_in_drop.npz')


@unittest.skip('')
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


@unittest.skip('')
class TestParallelConstructSam(unittest.TestCase):

    def setUp(self):
        data_dir = '/'.join(seqc.__file__.split('/')[:-2]) + '/data/'
        self.samfile = data_dir + 'in_drop/Aligned.out.sam'
        self.gtf = data_dir + 'genome/mm38_chr19/annotations.gtf'

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

@unittest.skip('')
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

    @unittest.skip('')
    def test_print_data(self):
        self.simple_unique._sort()
        data = append_fields(self.simple_unique.data, 'features', self.simple_unique.features)
        print(data[self.simple_unique._sorted][['cell', 'features', 'rmt']])

    @unittest.skip('')
    def test_counting_filter_too_high(self):

        # test that sorting is working
        # seems to be sorting cell, feature, rmt
        self.assertRaises(ValueError, self.simple_unique.to_experiment, 1)

    # @unittest.skip('')
    def test_counting_simple(self):
        exp = self.simple_unique.to_experiment(0)
        rdata = np.array(exp.reads.counts.todense())
        mdata = np.array(exp.molecules.counts.todense())
        print(pd.DataFrame(rdata, exp.reads.index, exp.reads.columns))
        print(pd.DataFrame(mdata, exp.molecules.index, exp.molecules.columns))

    def test_counting_doublets(self):
        exp = self.simple_duplicate.to_experiment(0)
        rdata = np.array(exp.reads.counts.todense())
        mdata = np.array(exp.molecules.counts.todense())
        print(pd.DataFrame(rdata, exp.reads.index, exp.reads.columns))
        print(pd.DataFrame(mdata, exp.molecules.index, exp.molecules.columns))


class TestUniqueArrayCreation(unittest.TestCase):
    """
    suspicion -> unique array creation is breaking something in the pipeline;
    could be an earlier method but lets check upstream so we know how each part is
    working.

    method: find the smallest real dataset that breaks it and use that to test.
    """

    def setUp(self):
        # verify that Aligned.out.sam doesn't produce correct results.
        self.samfile = _seqc_dir + '.test_download_index/Aligned.out.sam'
        self.gtf = _data_dir + 'genome/mm38_chr19/annotations.gtf'

    def test_num_unique_samfile(self):
        fc = seqc.convert_features.ConvertFeatureCoordinates.from_gtf(self.gtf, 1000)
        ra = seqc.arrays.ReadArray.from_samfile(self.samfile, fc)

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
    unittest.main(failfast=True, warnings='ignore')
