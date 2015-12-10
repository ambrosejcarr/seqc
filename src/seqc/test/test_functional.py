__author__ = 'ambrose'

import seqc
import nose2
import os
from copy import copy
import json
import sys
import logging
from seqc.test.test_unit import config

logging.getLogger('nose2').setLevel(logging.CRITICAL)


def test_seqc_raw_fastq_input():

    # create test directory if not present

    test_dir = config.seqc_dir + 'test_data/'
    if not os.path.isdir(test_dir):
        os.mkdir(test_dir)

    for data_type in config.data_types:

        seqc.log.info('SEQC functest for %s' % data_type)

        # define some variables
        output_prefix = config.seqc_dir + 'test_data/%s/functest/test_seqc' % data_type
        output_dir = config.seqc_dir + 'test_data/%s/functest/' % data_type
        # output_dir = testing_dir + data_type + '/'

        n_threads = config.n_threads

        # check if output_directory exists
        if not os.path.isdir(output_dir):
            seqc.log.info('SEQC test_data: creating testing directory')
            os.makedirs(output_dir)

        # check for necessary genome index files
        if not os.path.isdir(config.index):
            seqc.log.info('SEQC test_data: downloading mouse chr19 genome index files')
            bucket, key_prefix = seqc.io.S3.split_link(config.index_link)
            prefix = config.index
            seqc.io.S3.download_files(bucket=bucket, key_prefix=key_prefix,
                                      output_prefix=prefix, cut_dirs=2)

        # define some file names
        forward = config.forward_pattern % data_type
        reverse = config.reverse_pattern % data_type
        barcodes = config.barcode_partial_serial_pattern % data_type

        ######################## correct special barcode cases ##########################

        if data_type == 'drop_seq':
            barcodes = ''

        # get barcode data
        if barcodes:  # skip for types lacking barcode files
            # get directory
            *path, filename = barcodes.strip('\n').split('/')
            barcode_dir = '/'.join(path) + '/'

            # create directory if missing
            seqc.log.info('Creating %s barcode directory.' % data_type)
            if not os.path.isdir(barcode_dir):
                os.makedirs(barcode_dir)

            # download file if missing
            if not os.path.isfile(barcodes):
                seqc.log.info('Downloading %s barcode files' % data_type)
                bucket, key = seqc.io.S3.split_link(
                    config.barcode_partial_serial_link_pattern % data_type)
                seqc.io.S3.download_file(bucket=bucket, key=key, fout=barcodes)

        # define some genome files

        # get fastq data
        fastq_files = [forward, reverse]

        # check if there are any fastq files for the given experiment; if not skip
        # this section
        if any(fastq_files):

            # check if files must be created, or if they already exist in the directory
            if not all([True if os.path.isfile(f) or f is None else False
                        for f in fastq_files]):

                # check that fastq directory exists; get the first filename that exists
                # and use it to check/create necessary path.
                for fastq in fastq_files:
                    if fastq:
                        break
                if not fastq:
                    raise ValueError('no fastq file detected')
                *path, _ = fastq.strip('\n').split('/')
                fastq_dir = '/'.join(path) + '/'
                if not os.path.isdir(fastq_dir):
                    os.makedirs(fastq_dir)

                # get generator for correct fastq data type and create fastq files
                generate_fastq = getattr(seqc.fastq.GenerateFastq, data_type)
                prefix = fastq_dir + 'test_seqc'
                generate_fastq(10000, prefix, config.fasta, config.gtf, barcodes=barcodes,
                               tag_type='gene_id')

        args = [
            data_type.replace('_', '-'),
            '-i', config.index,
            '-n', str(n_threads),
            '-o', output_prefix,
        ]
        if barcodes:
            args.extend(['-b', barcodes])
        if forward:
            args.extend(['-f', forward,])
        if reverse:
            args.extend(['-r', reverse])

        # yield a test_data
        yield seqc_raw_fastq_input, (args,)


def seqc_raw_fastq_input(args):
    """
    This function should be an exact mimic of the __main__ section in the SEQC script
    EXCEPT you must pass args to seqc.core.parse_args(parser, args) because there are
    no sys.argv values to extract.
    """

    def initialize_logging():
        arg_copy = copy(kwargs)
        del arg_copy['func']  # function is not serializable

        seqc.log.info('SEQC version: %s' % seqc.__version__)
        seqc.log.info('SEQC working directory: %s' % os.getcwd())
        seqc.log.info('Passed command line arguments: %s' %
                      json.dumps(arg_copy, separators=(',', ': '), indent=4))

    parser = seqc.core.create_parser()
    kwargs = seqc.core.parse_args(parser, args)

    if kwargs['remote']:
        try:
            # log command line arguments for debugging
            initialize_logging()
            seqc.log.info('Starting REMOTE run')

            # run remotely
            del kwargs['remote']
            seqc.core.run_remote(kwargs)
        except:
            seqc.log.exception()
            raise

    elif kwargs['email_status']:
        try:
            # log command line arguments for debugging
            arg_copy = copy(kwargs)
            initialize_logging()

            # run locally
            func = kwargs['func']
            func(**kwargs)
            seqc.cluster_utils.upload_results(
                kwargs['output_prefix'], kwargs['email_status'], kwargs['aws_upload_key'])
        except:
            seqc.log.exception()
            email_body = b'Process interrupted -- see attached error message'
            seqc.cluster_utils.email_user(attachment='seqc.log', email_body=email_body,
                                          email_address=kwargs['email_status'])
            sys.exit(1)
        finally:
            # todo add cluster cleanup here
            # Note that this finally loop is executed regardless
            # of whether cluster SETUP has succeeded; there may not be any cluster.
            # Program accordingly.
            pass

    else:  # not email_status and not remote
        try:
            # log command line arguments for debugging
            arg_copy = copy(kwargs)
            initialize_logging()

            # run locally
            func = kwargs['func']
            func(**kwargs)
        except:
            seqc.log.exception()
            raise

    # ADDITIONAL TESTING OF FUNCTION OUTPUTS
    # get number of lines in fastq file
    n_records = 0
    for f in kwargs['forward']:
        with open(f) as fastq_in:
            n_records += sum(1 for line in fastq_in.readlines()) / 4

    # get number of lines in merged file
    prefix = kwargs['output_prefix'] + '/'
    merged = prefix + 'merged.fastq'
    with open(merged) as f:
        n_merged = sum(1 for line in f) / 4

    # (1) test that n_merged is within expected range
    margin = 2 * 4 * int(n_merged / (1e6 * 4))
    assert n_records - margin <= n_merged <= n_records

    # get number of alignments in samfile
    samfile_name = ('/'.join(kwargs['output_prefix'].split('/')[:-1]) +
                    '/test_seqc/Aligned.out.sam')
    rd = seqc.sam.Reader(samfile_name)
    n_alignments = sum(1 for ma in rd.iter_multialignments())

    # (2) check that the number of alignments is the same as the number of input records
    assert n_alignments >= n_merged * .99, '%d != %d' % (n_alignments, n_records)

    # (3) open the counts matrix and count up the reads that are considered valid
    npz_name = kwargs['output_prefix'] + '_sp_counts.npz'
    exp = seqc.analyze.Experiment.from_npz(npz_name)


if __name__ == "__main__":
    seqc.log.setup_logger('seqc_functest.log')
    try:
        nose2.main(exit=False, module=__name__, verbosity=5)
    except:
        seqc.log.exception()
        raise
