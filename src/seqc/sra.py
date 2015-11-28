__author__ = 'ambrose'

from xml.etree import ElementTree as ET
from subprocess import Popen, PIPE
import hashlib


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

