import os
from ftplib import FTP
from subprocess import check_call
import pandas as pd
import numpy as np
from seqc.sequence import gtf
from seqc.alignment import star
from seqc.io import S3


class Index:

    def __init__(self, organism, additional_id_types=None, index_folder_name='.'):
        """Create an Index object for organism, requiring that a valid annotation have
        both an ENSEMBL id and at least one additional id provided by an
        additional_id_field (if provided)

        :param organism: name of the organism. Must be lower-case genus_species. For
          example, human=homo_sapiens, mouse=mus_musculus, etc.
        :param additional_id_types: (default: None) Names of additional ID types whose
          presence declares the associated ENSEMBL gene id a valid identifier. If multiple
          such fields are passed, a gene can be defined in any of the provided fields, and
          it will be declared valid. If no fields are passed, no ID filtering will be
          performed, and all genes with ENSEMBL ids will be considered valid.

          The effect of these fields is to limit the ENSEMBL genes to a subset of genes
          which are also defined by other consortia. The rationale for requiring multiple
          definitions is that ENSEMBL has very relaxed standards, with many of its genes
          being defined based on predicted locations without any biological evidence, or,
          more importantly, any associated biological information. These such genes are
          often uninformative as a result, and are better excluded from the index.

          example fields for human are ['hgnc_symbol', 'entrezgene'], example fields for
          mouse are ['mgi_symbol', 'entrezgene']

          to find out the specific spelling and case of different id fields you may want
          to use, you will need to generate sample XML queries from ENSEMBL BioMART.
          This README describes how to accomplish this:
          http://useast.ensembl.org/info/data/biomart/biomart_restful.html
        """

        # check organism input
        if not organism:
            raise ValueError(
                'organism must be formatted as genus_species in all lower case')
        elif not isinstance(organism, str):
            raise TypeError('organism must be a string')
        elif any([('_' not in organism) or (organism.lower() != organism)]):
            raise ValueError(
                'organism must be formatted as genus_species in all lower case')
        self._organism = organism

        # check additional_id_fields argument
        if not (isinstance(additional_id_types, (list, tuple, np.ndarray)) or
                additional_id_types is None):
            raise TypeError(
                'if provided, additional id fields must be a list, tuple, or numpy '
                'array')
        if additional_id_types:
            self._additional_id_types = additional_id_types
        else:
            self._additional_id_types = []

        # todo type checks
        self.index_folder_name = index_folder_name

    @property
    def organism(self) -> str:
        return self._organism

    @property
    def additional_id_types(self) -> list:
        return self._additional_id_types

    @property
    def _converter_xml(self) -> str:
        """Generate The xml query to download an ENSEMBL BioMART file mapping
        ENSEMBL gene ids to any identifiers implemented in self.additional_id_fields
        """
        attributes = ''.join(
            '<Attribute name = "%s" />' % f for f in self.additional_id_types)
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
            '</Query>\''.format(genome=genome_name, attr=attributes))
        return xml

    @staticmethod
    def _identify_genome_file(files: [str]) -> str:
        """Identify and return the soft-masked primary assembly file from a list of fasta
        files. If the primary assembly is not present, default to the top-level file,
        which should always be present.

        :param files: list of fasta files obtained from the ENSEMBL ftp server
        :return str: name of the correct genome file"""
        for f in files:
            if '.dna_sm.primary_assembly' in f:
                return f
        for f in files:
            if f.endswith('.dna_sm.toplevel.fa.gz'):
                return f
        raise FileNotFoundError('could not find the correct fasta file in %r' % files)

    @staticmethod
    def _identify_gtf_file(files: [str], newest: int) -> str:
        """Identify and return the basic gtf file from a list of annotation files"""
        for f in files:
            if f.endswith('.%d.gtf.gz' % newest):
                return f

    @staticmethod
    def _identify_newest_release(open_ftp: FTP) -> int:
        """Identify the most recent genome release given an open link to ftp.ensembl.org

        This returns the second-to-last release, which is typically the last
        reliable/stable ENSEMBL release.

        :param FTP open_ftp: open FTP link to ftp.ensembl.org
        """
        open_ftp.cwd('/pub')
        releases = [f for f in open_ftp.nlst() if 'release' in f]
        newest = max(int(r[r.find('-') + 1:]) for r in releases)
        return newest - 1

    def _download_fasta_file(self, ftp: FTP, download_name: str) -> None:
        """download the fasta file for cls.organism from ftp, an open Ensembl FTP server

        :param FTP ftp: open FTP link to ENSEMBL
        :param str download_name: filename for downloaded fasta file
        """
        newest = self._identify_newest_release(ftp)
        ftp.cwd('/pub/release-%d/fasta/%s/dna' % (newest, self.organism))
        ensembl_fasta_filename = self._identify_genome_file(ftp.nlst())
        with open(download_name, 'wb') as f:
            ftp.retrbinary('RETR %s' % ensembl_fasta_filename, f.write)

    def _download_gtf_file(self, ftp, download_name) -> None:
        """download the gtf file for cls.organism from ftp, an open Ensembl FTP server

        :param FTP ftp: open FTP link to ENSEMBL
        :param str download_name: filename for downloaded gtf file
        """
        newest = self._identify_newest_release(ftp)
        ftp.cwd('/pub/release-%d/gtf/%s/' % (newest, self.organism))
        ensembl_gtf_filename = self._identify_gtf_file(ftp.nlst(), newest)
        with open(download_name, 'wb') as f:
            ftp.retrbinary('RETR %s' % ensembl_gtf_filename, f.write)

    # todo remove wget dependency
    def _download_conversion_file(self, download_name: str) -> None:
        """download a conversion file from BioMART that maps ENSEMBL ids to additional
        ids defined in cls.additional_id_fields

        :param download_name: name for the downloaded file
        """
        cmd = ('wget -O %s \'http://www.ensembl.org/biomart/martservice?query=%s > '
               '/dev/null 2>&1' % (download_name, self._converter_xml))

        err = check_call(cmd, shell=True)
        if err:
            raise ChildProcessError('conversion file download failed: %s' % err)

    def _download_ensembl_files(
            self, fasta_name: str=None, gtf_name: str=None,
            conversion_name: str=None) -> None:
        """download the fasta, gtf, and id_mapping file for the organism defined in
        cls.organism

        :param fasta_name: name for the downloaded fasta file
        :param gtf_name: name for the downloaded gtf file
        :param conversion_name: name for the downloaded conversion file
        """

        if fasta_name is None:
            fasta_name = '%s/%s.fa.gz' % (self.index_folder_name, self.organism)
        if gtf_name is None:
            gtf_name = '%s/%s.gtf.gz' % (self.index_folder_name, self.organism)
        if conversion_name is None:
            conversion_name = '%s/%s_ids.csv' % (self.index_folder_name, self.organism)

        with FTP(host='ftp.ensembl.org') as ftp:
            ftp.login()
            self._download_fasta_file(ftp, fasta_name)
            self._download_gtf_file(ftp, gtf_name)

        self._download_conversion_file(conversion_name)

    def _subset_genes(
            self,
            conversion_file: str=None,
            gtf_file: str=None,
            truncated_annotation: str=None,
            valid_biotypes=(b'protein_coding', b'lincRNA')):
        """
        Remove any annotation from the annotation_file that is not also defined by at
        least one additional identifer present in conversion file.

        The effect of these fields is to limit the ENSEMBL genes to a subset of genes
        which are also defined by other consortia. The rationale for requiring multiple
        definitions is that ENSEMBL has very relaxed standards, with many of its genes
        being defined based on predicted locations without any biological evidence, or,
        more importantly, any associated biological information. These such genes are
        often uninformative as a result, and are better excluded from the index.

        valid_biotypes removes genes that are of biotypes that single-cell sequencing
        is unlikely to detect. For example, miRNA are rarely poly-adenylated, and are
        of a size that they are often removed with primers. In our experience, the only
        biotypes that are worth considering are protein coding genes and lincRNA, the
        defaults for this function.

        :param conversion_file: file location of the conversion file
        :param gtf_file: file location of the annotation file
        :param truncated_annotation: name for the generated output file
        :param list(bytes) valid_biotypes: only accept genes of this biotype.
        """
        if not (self.additional_id_types or valid_biotypes):  # nothing to be done
            return

        # change to set for efficiency
        if all(isinstance(t, str) for t in valid_biotypes):
            valid_biotypes = set((t.encode() for t in valid_biotypes))
        elif all(isinstance(t, bytes) for t in valid_biotypes):
            valid_biotypes = set(valid_biotypes)
        else:
            raise TypeError('mixed-type biotypes detected. Please pass valid_biotypes '
                            'as strings or bytes objects (but not both).')

        if gtf_file is None:
            gtf_file = '%s/%s.gtf.gz' % (self.index_folder_name, self.organism)
        if conversion_file is None:
            conversion_file = '%s/%s_ids.csv' % (self.index_folder_name, self.organism)
        if truncated_annotation is None:
            truncated_annotation = '%s/%s_multiconsortia.gtf' % (
                self.index_folder_name, self.organism)

        # extract valid ensembl ids from the conversion file
        c = pd.read_csv(conversion_file, index_col=[0])
        valid_ensembl_ids = set(c[np.any(~c.isnull().values, axis=1)].index)

        # remove any invalid ids from the annotation file
        gr = gtf.Reader(gtf_file)
        with open(truncated_annotation, 'wb') as f:
            for line_fields in gr:
                record = gtf.Record(line_fields)
                if (record.attribute(b'gene_id').decode() in valid_ensembl_ids and
                        record.attribute(b'gene_biotype') in valid_biotypes):
                    f.write(bytes(record))

    def _create_star_index(
            self,
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
            fasta_file = '%s/%s.fa.gz' % (self.index_folder_name, self.organism)
        if gtf_file is None:
            if os.path.isfile('%s/%s_multiconsortia.gtf' % (
                    self.index_folder_name, self.organism)):
                gtf_file = '%s/%s_multiconsortia.gtf' % (
                    self.index_folder_name, self.organism)
            else:
                gtf_file = '%s/%s.gtf.gz' % (self.index_folder_name, self.organism)
        if genome_dir is None:
            genome_dir = '%s/%s' % (self.index_folder_name, self.organism)
        star.create_index(fasta_file, gtf_file, genome_dir, read_length)

    @staticmethod
    def _upload_index(index_directory: str, s3_upload_location: str) -> None:
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

    def create_index(
            self, valid_biotypes=('protein_coding', 'lincRNA'), s3_location: str=None):
        """create an optionally upload an index

        :param valid_biotypes: gene biotypes that do not match values in this list will
          be discarded from the annotation and will not appear in final count matrices
        :param s3_location: optional, s3 location to upload the index to.
        :return:
        """
        if self.index_folder_name is not '.':
            os.makedirs(self.index_folder_name, exist_ok=True)
        self._download_ensembl_files()
        self._subset_genes(valid_biotypes=valid_biotypes)
        self._create_star_index()
        if s3_location:
            self._upload_index('%s/%s' % (self.index_folder_name, self.organism),
                               s3_location)
