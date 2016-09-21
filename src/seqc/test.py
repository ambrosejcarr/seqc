# import random
import os
import unittest
import gzip
from seqc import remote, log
from seqc.io import S3
from seqc.core import parser
import pandas as pd
import numpy as np
import paramiko
from seqc.core import process_experiment
from seqc.sequence import index, gtf
import ftplib
import nose2

import seqc
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
        """Added for further input validation"""
        args = [
            'run',
            platform,
            '-o', self.output.format(platform),
            '-i', self.human,
            '--email-status', self.email,
            '-b', self.barcode_fastq.format(platform),
            '-g', self.genomic_fastq.format(platform),
            '--barcode-files', self.barcode_files.format(platform),
            '--log-name', self.log_name.format(platform),
            '--no-terminate', 'on-success',
        ]
        return args

    def test_recreate_command(self):
        parsed_args = parser.parse_args(self.format_args('in_drop'))
        print(parsed_args)
        print(parser.generate_remote_cmdline_args(self.format_args('in_drop')))
        print(parser.recreate_cmdline_args(self.format_args('in_drop')))


class TestRemoteProcessExperiment(unittest.TestCase):
    """
    Complete tests for the remote running of process_experiment.py
    """

    @classmethod
    def setUpClass(cls):
        os.makedirs(seqc_dir + 'test/TestRemoteProcessExperiment', exist_ok=True)
        cls.human = 's3://dplab-data/genomes/hg38_phiX/'
        cls.mouse = 's3://dplab-data/genomes/mm38_phiX/'
        cls.email = input('provide an email address to receive test results: ')
        cls.output = 's3://dplab-data/seqc/test/{}/'
        cls.barcode_files = 's3://dplab-data/barcodes/{}/flat/'
        cls.barcode_fastq = 's3://dplab-data/seqc/test/{}/barcode/'
        cls.genomic_fastq = 's3://dplab-data/seqc/test/{}/genomic/'
        cls.log_name = seqc_dir + 'test/TestRemoteProcessExperiment/seqc_{}.log'

    def test_in_drop(self):
        platform = 'in_drop'
        argv = [
            'run',
            platform,
            '-o', self.output.format(platform),
            '-i', self.human,
            '--email-status', self.email,
            '-b', self.barcode_fastq.format(platform),
            '-g', self.genomic_fastq.format(platform),
            '--barcode-files', self.barcode_files.format(platform),
            '--log-name', self.log_name.format(platform),
            '--no-terminate', 'on-success'
        ]
        try:
            arguments = parser.parse_args(argv)
            process_experiment.run(argv, arguments)
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
            '--no-terminate', 'on-success'
        ]
        try:
            process_experiment.run(args)
        except SystemExit:
            pass  # designed to exit when complete
        print("Initialization succeeded, wait for email to evaluate test results.")

    def test_in_drop_v2_no_mt(self):
        platform = 'in_drop_v2'
        args = [
            platform,
            '-o', self.output.format(platform),
            '-i', self.human,
            '--email-status', self.email,
            '-b', self.barcode_fastq.format(platform),
            '-g', self.genomic_fastq.format(platform),
            '--barcode-files', self.barcode_files.format(platform),
            '--no-terminate', 'on-success',
            '--no-filter-mitochondrial-rna',
        ]
        try:
            process_experiment.run(args)
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
            '--no-terminate', 'on-success'
        ]
        try:
            process_experiment.run(args)
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
            '--no-terminate', 'on-success'
        ]
        try:
            print(args)
            process_experiment.run(args)
        except SystemExit:
            pass  # designed to exit when complete
        print("Initialization succeeded, wait for email to evaluate test results.")

    def test_ten_x(self):
        platform = 'ten_x'
        args = [
            platform,
            '-o', self.output.format(platform),
            '-i', self.mouse,
            '--email-status', self.email,
            '-b', self.barcode_fastq.format(platform),
            '-g', self.genomic_fastq.format(platform),
            '--barcode-files', self.barcode_files.format(platform),
            '--no-terminate', 'on-success'
        ]
        try:
            print(args)
            process_experiment.run(args)
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
            '--barcode-files', self.barcode_files.format(platform),
            '--no-terminate', 'on-success'
        ]
        try:
            process_experiment.run(args)
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
            '--barcode-files', self.barcode_files.format(platform),
            '--no-terminate', 'on-success'
        ]
        try:
            process_experiment.run(args)
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
        process_experiment.run(args)
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
        process_experiment.run(args)
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
        process_experiment.run(args)
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
        process_experiment.run(args)
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
        process_experiment.run(args)
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
        process_experiment.run(args)
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
        process_experiment.run(args)
        print("Initialization succeeded, wait for email to evaluate test results.")


class TestingRemote(unittest.TestCase):
    """
    Testing EC2 management and SSH communications
    Should not be parallelized
    """  # TODO Spot bidding currently unsupported
    cluster = None

    @classmethod
    def setUp(cls):
        cls.vol = 10   # GB
        cls.cluster = remote.ClusterServer()
        config_file = os.path.expanduser('~/.seqc/config')
        cls.cluster.configure_cluster(config_file, "r3")

    @classmethod
    def tearDown(cls):
        try:
            if cls.cluster.is_cluster_running():
                inst_id_cluster = cls.cluster.inst_id.id
                remote.terminate_cluster(inst_id_cluster)
            remote.cluster_cleanup()
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
            # file not found
            self.cluster.keypath = os.path.expanduser("~/.seqc/config_foo")
            self.cluster.connect_server()
        try:
            # another key
            self.cluster.keypath = os.path.expanduser("~/.ssh/id_rsa_Github")
            self.cluster.connect_server()
        except (paramiko.AuthenticationException, paramiko.SSHException):
            pass
        except:  # todo which exceptions should this catch?
            self.fail("Bad authentication should have been caught")

    def test040_primed_for_remote_run(self):
        # This starts a remote job, but tearDown will kill it as soon as run_remote
        # exits => != from above
        ready_made_cmd = ['drop_seq', '-o', 's3://dplab-data/seqc/test/drop_seq/', '-i',
                          's3://dplab-data/genomes/hg38_phiX/', '--email-status',
                          'cyril-cros@hotmail.fr', '-b',
                          's3://dplab-data/seqc/test/drop_seq/barcode/', '-g',
                          's3://dplab-data/seqc/test/drop_seq/genomic/']
        args = process_experiment.parse_args(ready_made_cmd)
        process_experiment.run_remote(args, volsize=self.vol)


class TestIndexCreation(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.outdir = os.environ['TMPDIR']
        # cls.ftp = ftplib.FTP('ftp.ensembl.org')

    def test_Index_raises_ValueError_when_organism_is_not_provided(self):
        self.assertRaises(ValueError, index.Index, organism='', additional_id_fields=[])

    def test_Index_raises_ValueError_when_organism_isnt_lower_case(self):
        self.assertRaises(ValueError, index.Index, organism='Homo_sapiens',
                          additional_id_fields=[])
        self.assertRaises(ValueError, index.Index, organism='Homo_Sapiens',
                          additional_id_fields=[])
        self.assertRaises(ValueError, index.Index, organism='hoMO_Sapiens',
                          additional_id_fields=[])

    def test_Index_raises_ValueError_when_organism_has_no_underscore(self):
        self.assertRaises(ValueError, index.Index, organism='homosapiens',
                          additional_id_fields=[])

    def test_Index_raises_TypeError_when_additional_id_fields_is_not_correct_type(self):
        self.assertRaises(TypeError, index.Index, organism='homo_sapiens',
                          additional_id_fields='not_an_array_tuple_or_list')
        self.assertRaises(TypeError, index.Index, organism='homo_sapiens',
                          additional_id_fields='')

    def test_False_evaluating_additional_id_fields_are_accepted_but_set_empty_list(self):
        idx = index.Index('homo_sapiens', [])
        self.assertEqual(idx.additional_id_fields, [])
        idx = index.Index('homo_sapiens', tuple())
        self.assertEqual(idx.additional_id_fields, [])
        idx = index.Index('homo_sapiens', np.array([]))
        self.assertEqual(idx.additional_id_fields, [])

    def test_converter_xml_contains_one_attribute_line_per_gene_list(self):
        idx = index.Index('homo_sapiens', ['hgnc_symbol', 'mgi_symbol'])
        self.assertEqual(idx._converter_xml.count('Attribute name'), 3)
        idx = index.Index('homo_sapiens', [])
        self.assertEqual(idx._converter_xml.count('Attribute name'), 1)

    def test_converter_xml_formats_genome_as_first_initial_plus_species(self):
        idx = index.Index('homo_sapiens', ['hgnc_symbol', 'mgi_symbol'])
        self.assertIn('hsapiens', idx._converter_xml)
        idx = index.Index('mus_musculus')
        self.assertIn('mmusculus', idx._converter_xml)

    def test_can_login_to_ftp_ensembl(self):
        with ftplib.FTP(host='ftp.ensembl.org') as ftp:
            ftp.login()

    def test_download_converter_gets_output_and_is_pandas_loadable(self):
        idx = index.Index('ciona_intestinalis', ['entrezgene'])
        filename = self.outdir + 'ci.csv'
        idx._download_conversion_file(filename)
        converter = pd.read_csv(filename, index_col=0)
        self.assertGreaterEqual(len(converter), 10)
        self.assertEqual(converter.shape[1], 1)
        os.remove(filename)  # cleanup

    def test_identify_newest_release_finds_a_release_which_is_gt_eq_85(self):
        idx = index.Index('ciona_intestinalis', ['entrezgene'])
        with ftplib.FTP(host='ftp.ensembl.org') as ftp:
            ftp.login()
            newest = idx._identify_newest_release(ftp)
        self.assertGreaterEqual(int(newest), 85)  # current=85, they only get bigger

    def test_identify_genome_file_finds_primary_assembly_when_present(self):
        idx = index.Index('homo_sapiens', ['entrezgene'])
        with ftplib.FTP(host='ftp.ensembl.org') as ftp:
            ftp.login()
            newest = idx._identify_newest_release(ftp)
            ftp.cwd('/pub/release-%d/fasta/%s/dna' % (newest, idx.organism))
            filename = idx._identify_genome_file(ftp.nlst())
        self.assertIn('primary_assembly', filename)

    def test_identify_genome_file_defaults_to_toplevel_when_no_primary_assembly(self):
        idx = index.Index('ciona_intestinalis', ['entrezgene'])
        with ftplib.FTP(host='ftp.ensembl.org') as ftp:
            ftp.login()
            newest = idx._identify_newest_release(ftp)
            ftp.cwd('/pub/release-%d/fasta/%s/dna' % (newest, idx.organism))
            filename = idx._identify_genome_file(ftp.nlst())
        self.assertIn('toplevel', filename)

    def test_download_fasta_file_gets_a_properly_formatted_file(self):
        idx = index.Index('ciona_intestinalis', ['entrezgene'])
        with ftplib.FTP(host='ftp.ensembl.org') as ftp:
            ftp.login()
            filename = self.outdir + 'ci.fa.gz'
            idx._download_fasta_file(ftp, filename)
        with gzip.open(filename, 'rt') as f:
            self.assertIs(f.readline()[0], '>')  # starting character for genome fa record
        os.remove(filename)

    def test_identify_annotation_file_finds_a_gtf_file(self):
        idx = index.Index('ciona_intestinalis', ['entrezgene'])
        with ftplib.FTP(host='ftp.ensembl.org') as ftp:
            ftp.login()
            newest = idx._identify_newest_release(ftp)
            ftp.cwd('/pub/release-%d/gtf/%s/' % (newest, idx.organism))
            filename = idx._identify_gtf_file(ftp.nlst(), newest)
        self.assertIsNotNone(filename)

    def test_download_gtf_file_gets_a_file_readable_by_seqc_gtf_reader(self):
        idx = index.Index('ciona_intestinalis', ['entrezgene'])
        with ftplib.FTP(host='ftp.ensembl.org') as ftp:
            ftp.login()
            filename = self.outdir + 'ci.gtf.gz'
            idx._download_gtf_file(ftp, filename)
        rd = gtf.Reader(filename)
        rc = next(rd.iter_genes())
        self.assertIsInstance(rc, gtf.Gene)
        os.remove(filename)

    def test_subset_genes_does_nothing_if_no_additional_fields_or_valid_biotypes(self):
        idx = index.Index('ciona_intestinalis')
        fasta_name = self.outdir + 'ci.fa.gz'
        gtf_name = self.outdir + 'ci.gtf.gz'
        conversion_name = self.outdir + 'ci_ids.csv'
        idx._download_ensembl_files(fasta_name, gtf_name, conversion_name)
        idx._subset_genes(conversion_name, gtf_name, self.outdir + 'test.csv',
                          valid_biotypes=None)
        self.assertFalse(os.path.isfile(self.outdir))

    def test_subset_genes_produces_a_reduced_annotation_file_when_passed_fields(self):
        organism = 'ciona_intestinalis'
        idx = index.Index(organism, ['entrezgene'])
        os.chdir(self.outdir)
        idx._download_ensembl_files()
        self.assertTrue(os.path.isfile('%s.fa.gz' % organism), 'fasta file not found')
        self.assertTrue(os.path.isfile('%s.gtf.gz' % organism), 'gtf file not found')
        self.assertTrue(os.path.isfile('%s_ids.csv' % organism), 'id file not found')

        idx._subset_genes()
        self.assertTrue(os.path.isfile('%s_multiconsortia.gtf' % organism))
        gr_subset = gtf.Reader('%s_multiconsortia.gtf' % organism)
        gr_complete = gtf.Reader('%s.gtf.gz' % organism)
        self.assertLess(
            len(gr_subset), len(gr_complete),
            'Subset annotation was not smaller than the complete annotation')

        # make sure only valid biotypes are returned
        complete_invalid = False
        valid_biotypes = {b'protein_coding', b'lincRNA'}
        for r in gr_complete.iter_genes():
            if r.attribute(b'gene_biotype') not in valid_biotypes:
                complete_invalid = True
                break
        self.assertTrue(complete_invalid)
        subset_invalid = False
        for r in gr_subset.iter_genes():
            if r.attribute(b'gene_biotype') not in valid_biotypes:
                subset_invalid = True
                break
        self.assertFalse(subset_invalid)
        self.assertGreater(len(gr_subset), 0)

    def test_create_star_index_produces_an_index(self):
        organism = 'ciona_intestinalis'
        idx = index.Index(organism, ['entrezgene'])
        os.chdir(self.outdir)
        idx._download_ensembl_files()
        idx._subset_genes()
        print(os.getcwd())
        idx._create_star_index()
        self.assertTrue(os.path.isfile('{outdir}/{organism}/SAindex'.format(
            outdir=self.outdir, organism=organism)))

    def test_upload_star_index_correctly_places_index_on_s3(self):
        os.chdir(self.outdir)
        organism = 'ciona_intestinalis'
        idx = index.Index(organism, ['entrezgene'])
        index_directory = organism + '/'
        idx._download_ensembl_files()
        idx._subset_genes()
        idx._create_star_index()
        idx._upload_index(index_directory, 's3://dplab-data/genomes/ciona_intestinalis/')
        # watch output for success; will list uploads

    def test_create_index_produces_and_uploads_an_index(self):
        os.chdir(self.outdir)
        organism = 'ciona_intestinalis'
        idx = index.Index(organism, ['entrezgene'])
        idx.create_index(s3_location='s3://dplab-data/genomes/%s/' % idx.organism)


    @classmethod
    def tearDownClass(cls):
        pass


class TestLogData(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        test_folder = 'test/TestLogData/'
        os.makedirs(test_folder, exist_ok=True)
        # cls.logfile = test_folder + 'log'
        cls.logfile = test_folder + 'JE9636_LYMPHNODE2.log'
        if not os.path.isfile(cls.logfile):
            S3.download_file('dplab-data', 'seqc/test/in_drop_v2/log', cls.logfile)
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
        pattern = log.LogData.string_to_regex(self.test_pattern)
        assert pattern == pattern2, '{} != {}'.format(pattern, pattern2)

    def test_identify_duplicate_patterns(self):
        pattern = log.LogData.string_to_regex(self.test_pattern)
        d = log.LogData.identify_duplicate_patterns(pattern)
        assert dict(d) == {'genomic': 1, 'phi_x': 1}

    def test_replace_replicated_patterns(self):
        pattern = log.LogData.string_to_regex(self.test_pattern)
        duplicates = log.LogData.identify_duplicate_patterns(pattern)
        for k, v in duplicates.items():
            pattern = log.LogData.replace_replicated_patterns(pattern, k)
        # verify that no replicates remain
        residual = log.LogData.identify_duplicate_patterns(pattern)
        assert not residual

    def test_match_log(self):
        mo = log.LogData.match_log(self.logfile)
        assert mo

    def test_parse_special_fields(self):
        mo = log.LogData.match_log(self.logfile)
        matches = log.LogData.parse_special_fields(mo)
        assert all(matches)

    def test_dictionary_to_dataframe(self):
        mo = log.LogData.match_log(self.logfile)
        df = log.LogData.dictionary_to_dataframe(mo, 'test')
        assert isinstance(df, pd.DataFrame)

    def test_parse_log(self):
        log.LogData.parse_log(self.logfile)  # called for side-effect

    def test_parse_multiple(self):
        df = log.LogData.parse_multiple(os.path.expanduser(
            '~/google_drive/manuscripts/breast_cancer_immune/data/'),
            exclude='.*?log')
        assert not np.sum(np.isnan(df.values))

#########################################################################################

if __name__ == "__main__":
    # suite_Gen = unittest.makeSuite(TestProcessExperimentGeneral)
    # suite_Loc = unittest.makeSuite(TestLocalProcessExperiment)
    # suite_Rem = unittest.makeSuite(TestRemoteProcessExperiment)
    # suite_SSH = unittest.makeSuite(TestingRemote)
    # full_suite = unittest.TestSuite(suite_SSH)
    # runner_all_tests = unittest.TextTestRunner(descriptions=5, verbosity=5)
    # runner_all_tests.run(full_suite)
    nose2.main()

#########################################################################################
