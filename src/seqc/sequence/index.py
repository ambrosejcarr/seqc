from abc import ABCMeta, abstractmethod
from ftplib import FTP
from subprocess import check_call
import pandas as pd
import numpy as np
from seqc.sequence import gtf
from seqc.alignment import star
from seqc.io import S3


class AbstractIndex:

    __metaclass__ = ABCMeta

    @abstractmethod
    @property
    def organism(self) -> str:
        """Must be defined by subclasses. Should be the lower case genus and species,
        joined by an underscore. For example, human is homo_sapiens and mouse is
        mus_musculus.

        For implementation examples, see Mouse.organism or Human.organism, defined below.

        :return str: organism
        """
        pass

    @abstractmethod
    @property
    def additional_id_fields(self) -> list:
        """Must be defined by subclasses. Should be a list of id fields that whose
        presence declares the associated ENSEMBL gene id a valid identifier. If multiple
        such fields are passed, a gene can be defined in any of the provided fields, and
        it will be declared valid.

        The effect of these fields is to limit the ENSEMBL genes to a subset of genes
        which are also defined by other consortia. The rationale for requiring multiple
        definitions is that ENSEMBL has very relaxed standards, with many of its genes
        being defined based on predicted locations without any biological evidence, or,
        more importantly, any associated biological information. These such genes are
        often uninformative as a result, and are better excluded from the index.

        For example definitions, see Human.additional_id_fields or
        Mouse.additional_id_fields, defined below

        :return [str]: list of string id fields
        """
        pass

    @property
    def converter_xml(self) -> str:
        """Generate The xml query to download an ENSEMBL BioMART file mapping
        ENSEMBL gene ids to any identifiers implemented in self.additional_id_fields
        """
        attributes = ''.join(
            '<Attribute name = "%s" />' % f for f in self.additional_id_fields)
        genus, species = self.organism.split('_')
        genome_name = genus[0] + species
        xml = (
            '<?xml version = "1.0" encoding = "UTF-8"?>'
            '<!DOCTYPE Query>'
            '<Query virtualSchemaName = "default" formatter = "CSV" header = "0" '
            'uniqueRows = "0" count = "" datasetConfigVersion = "0.6">'
            '<Dataset name = "{genome}_gene_ensembl" interface = "default">'
            '<Attribute name = "ensembl_gene_id" />'
            '{attr}'
            '</Dataset>'
            '</Query>'.format(genome=genome_name, attr=attributes))
        return xml

    @staticmethod
    def identify_genome_file(files: [str]) -> str:
        """Identify and return the soft-masked genome file from a list of files"""
        for f in files:
            if '.dna_sm.primary_assembly' in f:
                return f

    @staticmethod
    def identify_gtf_file(files: [str], newest: int) -> str:
        """Identify and return the basic gtf file from a list of annotation files"""
        for f in files:
            if f.endswith('.%d.gtf.gz' % newest):
                return f

    @staticmethod
    def identify_newest_release(open_ftp: FTP) -> int:
        """Identify the most recent genome release given an open link to ftp.ensembl.org

        :param FTP open_ftp: open FTP link to ftp.ensembl.org
        """
        open_ftp.cwd('/pub')
        releases = [f for f in open_ftp.nlst() if 'release' in f]
        newest = max(int(r[r.find('-'):]) for r in releases)
        return newest

    @classmethod
    def download_fasta_file(cls, ftp: FTP, download_name: str) -> None:
        """download the fasta file for cls.organism from ftp, an open Ensembl FTP server

        :param FTP ftp: open FTP link to ENSEMBL
        :param str download_name: filename for downloaded fasta file
        """
        newest = cls.identify_newest_release(ftp)
        ftp.cwd('/pub/release-%d/fasta/%s/dna' % (newest, cls.organism))
        ensembl_fasta_filename = cls.identify_genome_file(ftp.nlst())
        with open(download_name, 'wb') as f:
            ftp.retrbinary('RETR %s' % ensembl_fasta_filename, f.write)

    @classmethod
    def download_gtf_file(cls, ftp, download_name) -> None:
        """download the gtf file for cls.organism from ftp, an open Ensembl FTP server

        :param FTP ftp: open FTP link to ENSEMBL
        :param str download_name: filename for downloaded gtf file
        """
        newest = cls.identify_newest_release(ftp)
        ftp.cwd('/pub/release-%d/gtf/%s/' % (newest, cls.organism))
        ensembl_gtf_filename = cls.identify_gtf_file(ftp.nlst(), newest)
        with open(download_name, 'wb') as f:
            ftp.retrbinary('RETR %s' % ensembl_gtf_filename, f.write)

    @classmethod  # todo remove wget dependency
    def download_conversion_file(cls, download_name: str) -> None:
        """download a conversion file from BioMART that maps ENSEMBL ids to additional
        ids defined in cls.additional_id_fields

        :param download_name: name for the downloaded file
        """
        cmd = 'wget -O %s http://www.ensembl.org/biomart/martservice?query=%s' % (
            download_name, cls.converter_xml
        )
        err = check_call(cmd, shell=True)
        if err:
            raise ChildProcessError('conversion file download failed: %s' % err)

    @classmethod
    def download_ensembl_files(cls, fasta_name: str=None, gtf_name: str=None,
                               conversion_name: str=None) -> None:
        """download the fasta, gtf, and id_mapping file for the organism defined in
        cls.organism

        :param fasta_name: name for the downloaded fasta file
        :param gtf_name: name for the downloaded gtf file
        :param conversion_name: name for the downloaded conversion file
        """

        if fasta_name is None:
            fasta_name = './%s.fa.gz' % cls.organism
        if gtf_name is None:
            gtf_name = './%s.fa.gz' % cls.organism
        if conversion_name is None:
            conversion_name = './%s_ids.csv' % cls.organism

        with FTP(host='ftp.ensembl.org') as ftp:
            cls.download_fasta_file(ftp, fasta_name)
            cls.download_gtf_file(ftp, gtf_name)

        cls.download_conversion_file(conversion_name)

    @classmethod
    def subset_genes(
            cls,
            conversion_file: str=None,
            gtf_file: str=None,
            truncated_annotation: str=None):
        """
        Remove any annotation from the annotation_file that is not also defined by at
        least one additional identifer present in conversion file.

        The effect of these fields is to limit the ENSEMBL genes to a subset of genes
        which are also defined by other consortia. The rationale for requiring multiple
        definitions is that ENSEMBL has very relaxed standards, with many of its genes
        being defined based on predicted locations without any biological evidence, or,
        more importantly, any associated biological information. These such genes are
        often uninformative as a result, and are better excluded from the index.

        :param conversion_file: file location of the conversion file
        :param gtf_file: file location of the annotation file
        :param truncated_annotation: name for the generated output file
        """
        if gtf_file is None:
            gtf_file = './%s.fa.gz' % cls.organism
        if conversion_file is None:
            conversion_file = './%s_ids.csv' % cls.organism
        if truncated_annotation is None:
            truncated_annotation = './%s_multiconsortia.gtf' % cls.organism

        # extract valid ensembl ids from the conversion file
        c = pd.read_csv(conversion_file, index_col=[0])
        valid_ensembl_ids = set(c[np.any(~c.isnull().values, axis=1)].index)

        # remove any invalid ids from the annotation file
        gr = gtf.Reader(gtf_file)
        with open(truncated_annotation, 'wb') as f:
            for line_fields in gr:
                record = gtf.Record(line_fields)
                if record.attribute(b'gene_id').decode() in valid_ensembl_ids:
                    f.write(bytes(record))

    @classmethod
    def create_star_index(
            cls,
            fasta_file: str=None,
            gtf_file: str=None,
            genome_dir: str=None,
            read_length: int=75) -> None:
        """Create a new STAR index for the associated genome

        :param fasta_file:
        :param gtf_file:
        :param genome_dir:
        :param read_length:
        :return:
        """
        if fasta_file is None:
            fasta_file = './%s.fa.gz' % cls.organism
        if gtf_file is None:
            gtf_file = './%s.gtf.gz' % cls.organism
        if genome_dir is None:
            genome_dir = cls.organism
        star.create_index(fasta_file, gtf_file, genome_dir, read_length)

    @classmethod
    def upload_index(cls, index_directory: str, s3_upload_location: str) -> None:
        """Upload the newly constructed index to s3 at s3_upload_location

        :param index_directory: folder containing index
        :param s3_upload_location: location to upload index on s3
        """
        if not index_directory.endswith('/'):
            index_directory += '/'
        if not s3_upload_location.endswith('/'):
            s3_upload_location += '/'
        bucket, *dirs = s3_upload_location.replace('s3://', '').split('/')
        key_prefix = '/'.join(dirs)
        S3.upload_files(file_prefix=index_directory, bucket=bucket, key_prefix=key_prefix)


class Human(AbstractIndex):

    @property
    def organism(self):
        return 'homo_sapiens'

    @property
    def additional_id_fields(self):
        return ['hgnc_symbol', 'entrezgene']


class Mouse(AbstractIndex):

    @property
    def organism(self):
        return 'mus_musculus'

    @property
    def additional_id_fields(self):
        return ['mgi_symbol', 'entrezgene']
