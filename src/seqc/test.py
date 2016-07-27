import random
import os
import unittest
import seqc
import pandas as pd
import numpy as np
import paramiko
from seqc import process_experiment
import nose2

seqc_dir = '/'.join(seqc.__file__.split('/')[:-3]) + '/'


class TestProcessExperimentGeneral(unittest.TestCase):
    """
    Tests for portions of ProcessExperiment (the main package pipeline function).
    These tests should immediately return success or failure.
    """

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

    def format_args(self, platform):
        """Added for further input validation"""    # TODO implement said validation
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
        return args

    def test_recreate_command(self):
        parsed_args = process_experiment.parse_args(self.format_args('in_drop'))
        print(process_experiment.recreate_command_line_arguments(parsed_args))


class TestRemoteProcessExperiment(unittest.TestCase):
    """
    Complete tests for the remote running of process_experiment.py
    """

    @classmethod
    def setUpClass(cls):
        os.makedirs(seqc_dir + 'test/TestRemoteProcessExperiment', exist_ok=True)
        cls.human = 's3://dplab-data/genomes/hg38_phiX/'
        cls.mouse = 's3://dplab-data/genomes/mm38_phiX/'
        # cls.email = 'cyril-cros@hotmail.fr'
        cls.email = input('provide an email address to receive test results: ')
        cls.output = 's3://dplab-data/seqc/test/{}/'
        cls.barcode_files = 's3://dplab-data/barcodes/{}/flat/'
        cls.barcode_fastq = 's3://dplab-data/seqc/test/{}/barcode/'
        cls.genomic_fastq = 's3://dplab-data/seqc/test/{}/genomic/'
        cls.log_name = seqc_dir + 'test/TestRemoteProcessExperiment/seqc_{}.log'

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

    @unittest.skip("Not currently stable")
    def test_in_drop_v3(self):
        platform = 'in_drop_v3'
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
            print(args)
            process_experiment.main(args)
        except SystemExit:
            pass  # designed to exit when complete
        print("Initialization succeeded, wait for email to evaluate test results.")

    @unittest.skip("Not currently stable")
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

    @unittest.skip("Not currently stable")
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
    """
    Complete tests for local running of process_experiment.py
    """

    @classmethod
    def setUpClass(cls):
        os.makedirs(seqc_dir + 'test/TestLocalProcessExperiment', exist_ok=True)
        cls.human = 's3://dplab-data/genomes/hg38_chr19/'
        cls.mouse = 's3://dplab-data/genomes/mm38_chr19/'
        # cls.email = input('provide an email address to receive test results: ')
        cls.output = 'test/TestLocalProcessExperiment/{}/test'  # Note: not directory
        cls.barcode_files = 's3://dplab-data/barcodes/{}/flat/'
        cls.barcode_fastq = 's3://dplab-data/seqc/test/{}_chr19/barcode/'
        cls.genomic_fastq = 's3://dplab-data/seqc/test/{}_chr19/genomic/'

    def test_in_drop(self):
        platform = 'in_drop'
        os.makedirs(seqc_dir + 'test/TestLocalProcessExperiment/{}'.format(platform),
                    exist_ok=True)
        args = [
            platform,
            '-o', self.output.format(platform),
            '-i', self.human,
            # '--email-status', self.email,
            '-b', self.barcode_fastq.format(platform),
            '-g', self.genomic_fastq.format(platform),
            '--barcode-files', self.barcode_files.format(platform),
            '--local',
        ]
        process_experiment.main(args)
        print("Initialization succeeded, wait for email to evaluate test results.")

    def test_in_drop_v2(self):
        platform = 'in_drop_v2'
        os.makedirs(seqc_dir + 'test/TestLocalProcessExperiment/{}'.format(platform),
                    exist_ok=True)
        args = [
            platform,
            '-o', self.output.format(platform),
            '-i', self.human,
            # '--email-status', self.email,
            '-b', self.barcode_fastq.format(platform),
            '-g', self.genomic_fastq.format(platform),
            '--barcode-files', self.barcode_files.format(platform),
            '--local'
        ]
        process_experiment.main(args)
        print("Initialization succeeded, wait for email to evaluate test results.")

    @unittest.skip("Not currently stable")
    def test_in_drop_v3(self):
        platform = 'in_drop_v3'
        os.makedirs(seqc_dir + 'test/TestLocalProcessExperiment/{}'.format(platform),
                    exist_ok=True)
        args = [
            platform,
            '-o', self.output.format(platform),
            '-i', self.human,
            # '--email-status', self.email,
            '-b', self.barcode_fastq.format(platform),
            '-g', self.genomic_fastq.format(platform),
            '--barcode-files', self.barcode_files.format(platform),
            '--local',
        ]
        process_experiment.main(args)
        print("Initialization succeeded, wait for email to evaluate test results.")

    def test_drop_seq(self):
        platform = 'drop_seq'
        os.makedirs(seqc_dir + 'test/TestLocalProcessExperiment/{}'.format(platform),
                    exist_ok=True)
        args = [
            platform,
            '-o', self.output.format(platform),
            '-i', self.human,
            # '--email-status', self.email,
            '-b', self.barcode_fastq.format(platform),
            '-g', self.genomic_fastq.format(platform),
            '--local',
        ]
        process_experiment.main(args)
        print("Initialization succeeded, wait for email to evaluate test results.")

    @unittest.skip("Not currently stable")
    def test_mars1_seq(self):
        platform = 'mars1_seq'
        os.makedirs(seqc_dir + 'test/TestLocalProcessExperiment/{}'.format(platform),
                    exist_ok=True)
        args = [
            platform,
            '-o', self.output.format(platform),
            '-i', self.human,
            # '--email-status', self.email,
            '-b', self.barcode_fastq.format(platform),
            '-g', self.genomic_fastq.format(platform),
            '--barcode-files', self.barcode_files.format(platform),
            '--local',
        ]
        process_experiment.main(args)
        print("Initialization succeeded, wait for email to evaluate test results.")

    @unittest.skip("Not currently stable")
    def test_mars2_seq(self):
        platform = 'mars2_seq'
        os.makedirs(seqc_dir + 'test/TestLocalProcessExperiment/{}'.format(platform),
                    exist_ok=True)
        args = [
            platform,
            '-o', self.output.format(platform),
            '-i', self.human,
            # '--email-status', self.email,
            '-b', self.barcode_fastq.format(platform),
            '-g', self.genomic_fastq.format(platform),
            '--barcode-files', self.barcode_files.format(platform),
            '--local',
        ]
        process_experiment.main(args)
        print("Initialization succeeded, wait for email to evaluate test results.")


class TestRemoteSpeciesMixExperiment(unittest.TestCase):
    """
    Complete tests for local running of process_experiment.py
    """

    @classmethod
    def setUpClass(cls):
        os.makedirs(seqc_dir + 'test/TestRemoteSpeciesMixExperiment', exist_ok=True)
        cls.index = 's3://dplab-data/genomes/hg19_mm10/'
        cls.email = input('provide an email address to receive test results: ')
        cls.output = 's3://dplab-data/seqc/test/species_mix/'
        cls.barcode_files = 's3://dplab-data/barcodes/{}/flat/'
        cls.barcode_fastq = 's3://dplab-data/seqc/test/species_mix/barcode/'
        cls.genomic_fastq = 's3://dplab-data/seqc/test/species_mix/genomic/'

    def test_in_drop_v2(self):
        platform = 'in_drop_v2'
        os.makedirs(seqc_dir + 'test/TestRemoteSpeciesMixExperiment/{}'.format(platform),
                    exist_ok=True)
        args = [
            platform,
            '-o', self.output,
            '-i', self.index,
            '--email-status', self.email,
            '-b', self.barcode_fastq,
            '-g', self.genomic_fastq,
            '--barcode-files', self.barcode_files.format(platform),
            '--instance-type', 'r3',
            '--spot-bid', '1.0',
        ]
        process_experiment.main(args)
        print("Initialization succeeded, wait for email to evaluate test results.")


class TestingRemote(unittest.TestCase):
    """
    Testing EC2 management and SSH communications
    Should not be parallelized
    """  # TODO Spot bidding currently unsupported

    @classmethod
    def setUp(cls):
        cls.vol = 10   # GB
        cls.cluster = seqc.remote.ClusterServer()
        config_file = os.path.expanduser('~/.seqc/config')
        cls.cluster.configure_cluster(config_file, "r3")

    @classmethod
    def tearDown(cls):
        try:
            if cls.cluster.is_cluster_running():
                inst_id_cluster = cls.cluster.inst_id.id
                seqc.remote.terminate_cluster(inst_id_cluster)
            seqc.remote.cluster_cleanup()
        except Exception as e:
            print(e)

    def test010_bad_instance(self):
        bad_inst_type = "r2d2"
        self.cluster.inst_type = bad_inst_type
        self.cluster.spot_bid = 1
        with self.assertRaises(ValueError):
            self.cluster.create_spot_cluster(self.vol)
        self.cluster.spot_bid = None
        with self.assertRaises(ValueError):
            self.cluster.create_cluster()

    def test030_volume_creation(self):
        self.cluster.create_security_group()
        self.cluster.create_cluster()
        with self.assertRaises(ValueError):
            self.cluster.allocate_space(vol_size=-1)
        self.cluster.allocate_space(vol_size=self.vol)

    def test020_bad_ssh_config(self):
        self.cluster.create_security_group()
        self.cluster.create_cluster()
        with self.assertRaises(FileNotFoundError):
            self.cluster.keypath = os.path.expanduser("~/.seqc/config_foo")  # file not found
            self.cluster.connect_server()
        try:
            self.cluster.keypath = os.path.expanduser("~/.ssh/id_rsa_Github")   # another key
            self.cluster.connect_server()
        except (paramiko.AuthenticationException, paramiko.SSHException):
            pass
        except:
            self.fail("Bad authentication should have been caught")

    def test040_primed_for_remote_run(self):
        # This starts a remote job, but tearDown will kill it as soon as run_remote exits => != from above
        ready_made_cmd = ['drop_seq', '-o', 's3://dplab-data/seqc/test/drop_seq/', '-i',
                          's3://dplab-data/genomes/hg38_phiX/', '--email-status', 'cyril-cros@hotmail.fr', '-b',
                          's3://dplab-data/seqc/test/drop_seq/barcode/', '-g',
                          's3://dplab-data/seqc/test/drop_seq/genomic/']
        args = process_experiment.parse_args(ready_made_cmd)
        process_experiment.run_remote(args, volsize=self.vol)


class TestThreeBitEquivalence(unittest.TestCase):
    """
    Tests that show that the decode and encode functions in
    seqc.sequence.Encodings.DNA3bit properly regenerate the input value.
    """

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


class TestLogData(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        test_folder = 'test/TestLogData/'
        os.makedirs(test_folder, exist_ok=True)
        # cls.logfile = test_folder + 'seqc.log'
        cls.logfile = test_folder + 'JE9636_LYMPHNODE2.log'
        if not os.path.isfile(cls.logfile):
            seqc.io.S3.download_file(
                'dplab-data', 'seqc/test/in_drop_v2/seqc.log', cls.logfile)
        # note that this is a frozen instance of the summary that is printed in the log
        # at v0.1.8rc3
        cls.test_pattern = (
            '{divide}\nINPUT\n{divide}\n'
            'Total input reads:\t{n_fastq}\n'
            '{divide}\nALIGNMENT (% FROM INPUT)\n{divide}\n'
            'Total reads aligned:\t{n_sam} ({prop_al}%)\n'
            ' - Genomic alignments:\t{genomic} ({prop_gen}%)\n'
            ' - PhiX alignments:\t{phi_x} ({prop_phix}%)\n'
            ' - Transcriptome alignments:\t{trans} ({prop_trans}%)\n'
            '{divide}\nFILTERING (% FROM ALIGNMENT)\n{divide}\n'
            'Genomic alignments:\t{genomic} ({bad_gen}%)\n'
            'PhiX alignments:\t{phi_x} ({bad_phi}%)\n'
            'Incorrect barcodes:\t{wrong_cb} ({bad_cb}%)\n'
            'Missing cell barcodes/RMT:\t{no_cell} ({bad_cell}%)\n'
            'N present in RMT:\t{rmt_N} ({bad_rmtN}%)\n'
            'Insufficient poly(T):\t{poly_t} ({bad_polyt}%)\n'
            '{divide}\nCELL/MOLECULE COUNT DISTRIBUTION\n{divide}\n'
            'Total molecules:\t\t{tot_mc}\n'
            'Molecules lost:\t{mols_lost}\n'
            'Cells lost:\t{cells_lost}\n'
            'Cell description:\n{cell_desc}\n'
            '{divide}\nSUMMARY\n{divide}\n'
            'Total retained reads:\t{n_good} ({prop_good}%)\n'
            'Total reads unaligned:\t{lost_al} ({prop_un}%)\n'
            'Total reads filtered:\t{n_bad} ({prop_bad}%)\n'
            '{divide}\n')

    def test_string_to_regex(self):
        pattern2 = self.test_pattern.replace('{divide}', '-*?')
        pattern2 = pattern2.replace('(', '\(')
        pattern2 = pattern2.replace(')', '\)')
        pattern2 = pattern2.replace('{', '(?P<')
        pattern2 = pattern2.replace('}', '>.*?)')
        pattern = seqc.log.LogData.string_to_regex(self.test_pattern)
        assert pattern == pattern2, '{} != {}'.format(pattern, pattern2)

    def test_identify_duplicate_patterns(self):
        pattern = seqc.log.LogData.string_to_regex(self.test_pattern)
        d = seqc.log.LogData.identify_duplicate_patterns(pattern)
        assert dict(d) == {'genomic': 1, 'phi_x': 1}

    def test_replace_replicated_patterns(self):
        pattern = seqc.log.LogData.string_to_regex(self.test_pattern)
        duplicates = seqc.log.LogData.identify_duplicate_patterns(pattern)
        for k, v in duplicates.items():
            pattern = seqc.log.LogData.replace_replicated_patterns(pattern, k)
        # verify that no replicates remain
        residual = seqc.log.LogData.identify_duplicate_patterns(pattern)
        assert not residual

    def test_match_log(self):
        mo = seqc.log.LogData.match_log(self.logfile)
        assert mo

    def test_parse_special_fields(self):
        mo = seqc.log.LogData.match_log(self.logfile)
        matches = seqc.log.LogData.parse_special_fields(mo)
        assert all(matches)

    def test_dictionary_to_dataframe(self):
        mo = seqc.log.LogData.match_log(self.logfile)
        df = seqc.log.LogData.dictionary_to_dataframe(mo, 'test')
        assert isinstance(df, pd.DataFrame)

    def test_parse_log(self):
        df = seqc.log.LogData.parse_log(self.logfile)

    def test_parse_multiple(self):
        df = seqc.log.LogData.parse_multiple(os.path.expanduser(
            '~/google_drive/manuscripts/breast_cancer_immune/data/'),
            exclude='.*?seqc.log')
        assert not np.sum(np.isnan(df.values))

###########################################################################################

if __name__ == "__main__":
    # suite_Gen = unittest.makeSuite(TestProcessExperimentGeneral)
    # suite_Loc = unittest.makeSuite(TestLocalProcessExperiment)
    # suite_Rem = unittest.makeSuite(TestRemoteProcessExperiment)
    # suite_SSH = unittest.makeSuite(TestingRemote)
    # full_suite = unittest.TestSuite(suite_SSH)
    # runner_all_tests = unittest.TextTestRunner(descriptions=5, verbosity=5)
    # runner_all_tests.run(full_suite)
    nose2.main()

###########################################################################################



