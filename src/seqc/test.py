import os
import unittest
import gzip
import pandas as pd
import numpy as np
import ftplib
import nose2
from nose2.tools import params
from seqc.sequence import index, gtf
from seqc.core import main
from seqc.read_array import ReadArray
import seqc
import logging

logging.basicConfig()

seqc_dir = '/'.join(seqc.__file__.split('/')[:-3]) + '/'

# fill and uncomment these variables to avoid having to provide input to tests
TEST_BUCKET = "seqc-public"  # None
EMAIL = None
# RSA_KEY = None

# define some constants for testing
BARCODE_FASTQ = 's3://seqc-public/test/%s/barcode/'  # platform
GENOMIC_FASTQ = 's3://seqc-public/test/%s/genomic/'  # platform
MERGED = 's3://seqc-public/test/%s/%s_merged.fastq.gz'  # platform, platform
SAMFILE = 's3://seqc-public/test/%s/Aligned.out.bam'  # platform
INDEX = 's3://seqc-public/genomes/hg38_chr19/'
LOCAL_OUTPUT = os.environ['TMPDIR'] + 'seqc/%s/test'  # test_name
REMOTE_OUTPUT = './test'
UPLOAD = 's3://%s/seqc_test/%s/'  # bucket_name, test_folder
PLATFORM_BARCODES = 's3://seqc-public/barcodes/%s/flat/'  # platform


def makedirs(output):
    """make directories based on OUTPUT"""
    stem, prefix = os.path.split(output)
    os.makedirs(stem, exist_ok=True)


class TestSEQC(unittest.TestCase):

    email = globals()['EMAIL'] if 'EMAIL' in globals() else None
    bucket = globals()['TEST_BUCKET'] if 'TEST_BUCKET' in globals() else None
    rsa_key = globals()['RSA_KEY'] if 'RSA_KEY' in globals() else None
    if rsa_key is None:
        try:
            rsa_key = os.environ['AWS_RSA_KEY']
        except KeyError:
            pass

    def check_parameters(self):
        if self.email is None:
            self.email = input('please provide an email address for SEQC to mail '
                               'results: ')
        if self.bucket is None:
            self.bucket = input('please provide an amazon s3 bucket to upload test '
                                'results: ')
            self.bucket.rstrip('/')
        if self.rsa_key is None:
            self.rsa_key = input('please provide an RSA key with permission to create '
                                 'aws instances: ')

    # @params('in_drop', 'in_drop_v2', 'drop_seq', 'ten_x', 'mars_seq')
    def test_local(self, platform='in_drop_v2'):
        """test seqc after pre-downloading all files"""
        with open('seqc_log.txt', 'w') as f:
            f.write('Dummy log. nose2 captures input, so no log is produced. This causes '
                    'pipeline errors.\n')
        test_name = 'test_no_aws_%s' % platform
        makedirs(LOCAL_OUTPUT % test_name)
        if self.email is None:
            self.email = input('please provide an email address for SEQC to mail results: ')
        argv = [
            'run',
            platform,
            '-o', test_name,
            '-i', INDEX,
            '-b', BARCODE_FASTQ % platform,
            '-g', GENOMIC_FASTQ % platform,
            '--barcode-files', PLATFORM_BARCODES % platform,
            '-e', self.email,
            '--local']
        main.main(argv)
        os.remove('./seqc_log.txt')  # clean up the dummy log we made.

    # @params('in_drop', 'in_drop_v2', 'drop_seq', 'ten_x', 'mars_seq')
    def test_remote_from_raw_fastq(self, platform='ten_x_v2'):
        test_name = 'test_remote_%s' % platform
        self.check_parameters()
        argv = [
            'run',
            platform,
            '-o', REMOTE_OUTPUT,
            '-u', UPLOAD % (self.bucket, test_name),
            '-i', INDEX,
            '-e', self.email,
            '-b', BARCODE_FASTQ % platform,
            '-g', GENOMIC_FASTQ % platform,
            '--instance-type', 'c4.large',
            '--spot-bid', '1.0',
            '-k', self.rsa_key, '--debug']
        if platform != 'drop_seq':
            argv += ['--barcode-files', PLATFORM_BARCODES % platform]
        main.main(argv)

    # @params('in_drop', 'in_drop_v2', 'drop_seq', 'ten_x', 'mars_seq')
    def test_remote_from_merged(self, platform='in_drop_v2'):
        test_name = 'test_remote_%s' % platform
        self.check_parameters()
        argv = [
            'run',
            platform,
            '-o', REMOTE_OUTPUT,
            '-u', UPLOAD % (self.bucket, test_name),
            '-i', INDEX,
            '-e', self.email,
            '-m', MERGED % (platform, platform),
            '-k', self.rsa_key,
            '--instance-type', 'c4.large',
            # '--spot-bid', '1.0'
        ]
        if platform != 'drop_seq':
            argv += ['--barcode-files', PLATFORM_BARCODES % platform]
        main.main(argv)

    # @params('in_drop', 'in_drop_v2', 'drop_seq', 'ten_x', 'mars_seq')
    def test_remote_from_samfile(self, platform='in_drop_v2'):
        test_name = 'test_remote_%s' % platform
        self.check_parameters()
        argv = [
            'run',
            platform,
            '-o', REMOTE_OUTPUT,
            '-u', UPLOAD % (self.bucket, test_name),
            '-i', INDEX,
            '-e', self.email,
            '-a', SAMFILE % platform,
            '-k', self.rsa_key,
            '--instance-type', 'c4.8xlarge',
            '--debug',
            # '--spot-bid', '1.0'
        ]
        if platform != 'drop_seq':
            argv += ['--barcode-files', PLATFORM_BARCODES % platform]
        main.main(argv)


class TestIndex(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.outdir = os.environ['TMPDIR']

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
        self.assertEqual(idx.additional_id_types, [])
        idx = index.Index('homo_sapiens', tuple())
        self.assertEqual(idx.additional_id_types, [])
        idx = index.Index('homo_sapiens', np.array([]))
        self.assertEqual(idx.additional_id_types, [])

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
        if 'TEST_BUCKET' in globals():
            bucket = globals()['TEST_BUCKET']
        else:
            bucket = input('please provide an amazon s3 bucket to upload test results: ')
        organism = 'ciona_intestinalis'
        idx = index.Index(organism, ['entrezgene'])
        index_directory = organism + '/'
        idx._download_ensembl_files()
        idx._subset_genes()
        idx._create_star_index()
        idx._upload_index(index_directory, 's3://%s/genomes/ciona_intestinalis/' % bucket)

    def test_create_index_produces_and_uploads_an_index(self):
        if 'TEST_BUCKET' in globals():
            bucket = globals()['TEST_BUCKET']
        else:
            bucket = input('please provide an amazon s3 bucket to upload test results: ')
        organism = 'ciona_intestinalis'
        idx = index.Index(organism, ['entrezgene'], self.outdir)
        idx.create_index(s3_location='s3://%s/genomes/%s/' % (bucket, idx.organism))


class TestReadArrayCreation(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        platform = 'in_drop_v2'
        # cls.bamfile = LOCAL_OUTPUT % platform + '_bamfile.bam'
        # cls.annotation = LOCAL_OUTPUT % platform + '_annotations.gtf'
        # S3.download(SAMFILE % platform, cls.bamfile, recursive=False)
        # S3.download(INDEX + 'annotations.gtf', cls.annotation, recursive=False)

        cls.bamfile = os.path.expanduser('~/Downloads/mm_test_short.bam')
        cls.annotation = os.path.expanduser('~/Downloads/annotations.gtf')
        cls.summary = os.path.expanduser('~/Downloads/mm_test_summary.txt')
        cls.total_input_reads = 12242659
        cls.translator = gtf.GeneIntervals(cls.annotation, 10000)

    def test_read_array_creation(self):
        ra = ReadArray.from_alignment_file(self.bamfile, self.translator,
                                           required_poly_t=0)
        print(repr(ra.data.shape[0]))
        print(repr(ra.genes))
        ra.save(os.path.expanduser('~/Downloads/test.ra'))


class TestTranslator(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.annotation = os.path.expanduser('~/Downloads/annotations.gtf')

    def test_construct_translator(self):
        translator = gtf.GeneIntervals(self.annotation)
        print(len(translator._chromosomes_to_genes))

    def get_length_of_gtf(self):
        rd = gtf.Reader(self.annotation)
        # print(len(rd))
        print(sum(1 for _ in rd.iter_transcripts()))


#########################################################################################

if __name__ == "__main__":
    nose2.main()

#########################################################################################
