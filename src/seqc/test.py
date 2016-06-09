import random
import os
import nose2
import unittest
import seqc
from seqc import process_experiment


seqc_dir = '/'.join(seqc.__file__.split('/')[:-3]) + '/'


class TestProcessExperimentGeneral(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.human = 's3://dplab-data/genomes/hg38_phiX/'
        cls.mouse = 's3://dplab-data/genomes/mm38_phiX/'
        cls.email = 'testemail@nowhere.com'
        cls.output = 's3://dplab-data/seqc/test/{}/'
        cls.barcode_files = 's3://dplab-data/barcodes/{}/flat/'
        cls.barcode_fastq = 's3://dplab-data/seqc/test/{}/barcode/'
        cls.genomic_fastq = 's3://dplab-data/seqc/test/{}/genomic/'
        cls.log_name = seqc_dir + 'test/test_remote/seqc_{}.log'

    def test_recreate_command(self):
        platform = 'in_drop'
        args = [
            platform,
            '-o', self.output.format(platform),
            '-i', self.human,
            '--email-status', self.email,
            '-b', self.barcode_fastq.format(platform),
            '-g', self.genomic_fastq.format(platform),
            '--barcode-files', self.barcode_files.format(platform),
            '--log-name', self.log_name.format(platform),
        ]
        parsed_args = process_experiment.parse_args(args)
        print(process_experiment.recreate_command_line_arguments(parsed_args))


class TestRemoteProcessExperiment(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        os.makedirs(seqc_dir + 'test/test_remote', exist_ok=True)  # make sure test directory exists
        cls.human = 's3://dplab-data/genomes/hg38_phiX/'
        cls.mouse = 's3://dplab-data/genomes/mm38_phiX/'
        cls.email = input('provide an email address to receive test results: ')
        cls.output = 's3://dplab-data/seqc/test/{}/'
        cls.barcode_files = 's3://dplab-data/barcodes/{}/flat/'
        cls.barcode_fastq = 's3://dplab-data/seqc/test/{}/barcode/'
        cls.genomic_fastq = 's3://dplab-data/seqc/test/{}/genomic/'
        cls.log_name = seqc_dir + 'test/test_remote/seqc_{}.log'

    def test_in_drop(self):
        platform = 'in_drop'
        args = [
            platform,
            '-o', self.output.format(platform),
            '-i', self.human,
            '--email-status', self.email,
            '-b', self.barcode_fastq.format(platform),
            '-g', self.genomic_fastq.format(platform),
            '--barcode-files', self.barcode_files.format(platform),
            '--log-name', self.log_name.format(platform)
        ]
        try:
            process_experiment.main(args)
        except SystemExit:
            pass  # designed to exit when complete
        print("Initialization succeeded, wait for email to evaluate test results.")

    def test_in_drop_v2(self):
        platform = 'in_drop_v2'
        args = [
            platform,
            '-o', self.output.format(platform),
            '-i', self.human,
            '--email-status', self.email,
            '-b', self.barcode_fastq.format(platform),
            '-g', self.genomic_fastq.format(platform),
            '--barcode-files', self.barcode_files.format(platform),
        ]
        try:
            process_experiment.main(args)
        except SystemExit:
            pass  # designed to exit when complete
        print("Initialization succeeded, wait for email to evaluate test results.")

    def test_drop_seq(self):
        platform = 'drop_seq'
        args = [
            platform,
            '-o', self.output.format(platform),
            '-i', self.human,
            '--email-status', self.email,
            '-b', self.barcode_fastq.format(platform),
            '-g', self.genomic_fastq.format(platform),
        ]
        try:
            process_experiment.main(args)
        except SystemExit:
            pass  # designed to exit when complete
        print("Initialization succeeded, wait for email to evaluate test results.")

    def test_mars1_seq(self):
        platform = 'mars1_seq'
        args = [
            platform,
            '-o', self.output.format(platform),
            '-i', self.human,
            '--email-status', self.email,
            '-b', self.barcode_fastq.format(platform),
            '-g', self.genomic_fastq.format(platform),
            '--barcode-files', self.barcode_files.format(platform)
        ]
        try:
            process_experiment.main(args)
        except SystemExit:
            pass  # designed to exit when complete
        print("Initialization succeeded, wait for email to evaluate test results.")

    def test_mars2_seq(self):
        platform = 'mars2_seq'
        args = [
            platform,
            '-o', self.output.format(platform),
            '-i', self.human,
            '--email-status', self.email,
            '-b', self.barcode_fastq.format(platform),
            '-g', self.genomic_fastq.format(platform),
            '--barcode-files', self.barcode_files.format(platform)
        ]
        try:
            process_experiment.main(args)
        except SystemExit:
            pass  # designed to exit when complete
        print("Initialization succeeded, wait for email to evaluate test results.")


class TestLocalProcessExperiment(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.human = 's3://dplab-data/genomes/hg38_phiX'
        cls.mouse = 's3://dplab-data/genomes/mm38_phiX'
        cls.email = input('provide an email address to receive test results: ')
        cls.output = 's3://dplab-data/seqc/test/{platform}/'
        cls.barcode_files = 's3://dplab-data/barcodes/{platform}/flat/'
        cls.barcode_fastq = 's3://dplab-data/seqc/test/{platform}/barcode/'
        cls.genomic_fastq = 's3://dplab-data/seqc/test/{platform}/genomic/'

    def test_in_drop(self):
        platform = 'in_drop'
        args = [
            platform,
            '-o', self.output.format(platform),
            '-i', self.human,
            '--email-status', self.email,
            '-b', self.barcode_fastq.format(platform),
            '-g', self.genomic_fastq.format(platform),
            '--barcode-files', self.barcode_files.format(platform)
        ]
        process_experiment.main(args)
        print("Initialization succeeded, wait for email to evaluate test results.")

    def test_in_drop_v2(self):
        platform = 'in_drop_v2'
        args = [
            platform,
            '-o', self.output.format(platform),
            '-i', self.human,
            '--email-status', self.email,
            '-b', self.barcode_fastq.format(platform),
            '-g', self.genomic_fastq.format(platform),
            '--barcode-files', self.barcode_files.format(platform)
        ]
        process_experiment.main(args)
        print("Initialization succeeded, wait for email to evaluate test results.")

    def test_drop_seq(self):
        platform = 'drop_seq'
        args = [
            platform,
            '-o', self.output.format(platform),
            '-i', self.human,
            '--email-status', self.email,
            '-b', self.barcode_fastq.format(platform),
            '-g', self.genomic_fastq.format(platform),
        ]
        process_experiment.main(args)
        print("Initialization succeeded, wait for email to evaluate test results.")

    def test_mars1_seq(self):
        platform = 'mars1_seq'
        args = [
            platform,
            '-o', self.output.format(platform),
            '-i', self.human,
            '--email-status', self.email,
            '-b', self.barcode_fastq.format(platform),
            '-g', self.genomic_fastq.format(platform),
            '--barcode-files', self.barcode_files.format(platform)
        ]
        process_experiment.main(args)
        print("Initialization succeeded, wait for email to evaluate test results.")

    def test_mars2_seq(self):
        platform = 'mars2_seq'
        args = [
            platform,
            '-o', self.output.format(platform),
            '-i', self.human,
            '--email-status', self.email,
            '-b', self.barcode_fastq.format(platform),
            '-g', self.genomic_fastq.format(platform),
            '--barcode-files', self.barcode_files.format(platform)
        ]
        process_experiment.main(args)
        print("Initialization succeeded, wait for email to evaluate test results.")


class TestThreeBitEquivalence(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cb1_file = seqc_dir + 'data_test/in_drop/barcodes/cb1.txt'
        cb2_file = seqc_dir + 'data_test/in_drop/barcodes/cb2.txt'

        with open(cb1_file) as f1, open(cb2_file) as f2:
            cls.cb1_codes = [c.strip() for c in f1.readlines()]
            cls.cb2_codes = [c.strip() for c in f2.readlines()]

    def test_three_bit_version_1(self):
        """
        :return:
        """

        # test conversion of single cb1
        for _ in range(100):
            case = random.choice(self.cb1_codes).encode()
            coded = seqc.sequence.encodings.DNA3Bit.encode(case)
            decoded = seqc.sequence.encodings.DNA3Bit.decode(coded)
            self.assertEqual(case, decoded)

        # test conversion of single cb2
        for _ in range(100):
            case = random.choice(self.cb2_codes).encode()
            coded = seqc.sequence.encodings.DNA3Bit.encode(case)
            decoded = seqc.sequence.encodings.DNA3Bit.decode(coded)
            self.assertEqual(case, decoded)

        # test conversion of merged barcodes
        for _ in range(100):
            case1 = random.choice(self.cb1_codes).encode()
            case2 = random.choice(self.cb2_codes).encode()
            coded = seqc.sequence.encodings.DNA3Bit.encode(case1 + case2)
            decoded = seqc.sequence.encodings.DNA3Bit.decode(coded)
            self.assertEqual(case1 + case2, decoded)

    def test_three_bit_version_2(self):
        """
        :return:
        """

        # test conversion of single cb1
        for _ in range(100):
            case = random.choice(self.cb1_codes)
            coded = seqc.sequence.encodings.ThreeBit.str2bin(case)
            decoded = seqc.sequence.encodings.ThreeBit.bin2str(coded)
            self.assertEqual(case, decoded)

        # test conversion of single cb2
        for _ in range(100):
            case = random.choice(self.cb2_codes)
            coded = seqc.sequence.encodings.ThreeBit.str2bin(case)
            decoded = seqc.sequence.encodings.ThreeBit.bin2str(coded)
            self.assertEqual(case, decoded)

        # test conversion of merged barcodes
        for _ in range(100):
            case1 = random.choice(self.cb1_codes)
            case2 = random.choice(self.cb2_codes)
            coded = seqc.sequence.encodings.ThreeBit.str2bin(case1 + case2)
            decoded = seqc.sequence.encodings.ThreeBit.bin2str(coded)
            self.assertEqual(case1 + case2, decoded)

    def test_mixed_encoding(self):

        # test conversion of single cb1
        for _ in range(100):
            case = random.choice(self.cb1_codes).encode()
            coded = seqc.sequence.encodings.DNA3Bit.encode(case)
            decoded = seqc.sequence.encodings.ThreeBit.bin2str(coded).encode()
            self.assertEqual(case, decoded)

        # test conversion of single cb2
        for _ in range(100):
            case = random.choice(self.cb2_codes).encode()
            coded = seqc.sequence.encodings.DNA3Bit.encode(case)
            decoded = seqc.sequence.encodings.ThreeBit.bin2str(coded).encode()
            self.assertEqual(case, decoded)

        # test conversion of merged barcodes
        for _ in range(100):
            case1 = random.choice(self.cb1_codes).encode()
            case2 = random.choice(self.cb2_codes).encode()
            coded = seqc.sequence.encodings.DNA3Bit.encode(case1 + case2)
            decoded = seqc.sequence.encodings.ThreeBit.bin2str(coded).encode()
            self.assertEqual(case1 + case2, decoded)

        # test conversion and extraction of cb1, cb2
        for _ in range(100):  # todo fails
            case1 = random.choice(self.cb1_codes)
            case2 = random.choice(self.cb2_codes)
            coded = seqc.sequence.encodings.DNA3Bit.encode(
                    case1.encode() + case2.encode())
            case1_int = seqc.sequence.encodings.ThreeBit.c1_from_codes(coded)
            case2_int = seqc.sequence.encodings.ThreeBit.c2_from_codes(coded)
            case1_decoded = seqc.sequence.encodings.ThreeBit.bin2str(case1_int)
            case2_decoded = seqc.sequence.encodings.ThreeBit.bin2str(case2_int)
            self.assertEqual(case1, case1_decoded)
            self.assertEqual(case2, case2_decoded)


if __name__ == "__main__":
    nose2.main()
