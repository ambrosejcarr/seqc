__author__ = 'ambrose'

import seqc
import nose2
import os
from copy import copy
import json
import sys
import logging

logging.getLogger('nose2').setLevel(logging.CRITICAL)


def test_seqc_raw_fastq_input():

    # create test directory if not present
    package_dir = '/'.join(seqc.__file__.split('/')[:-3]) + '/'
    test_dir = package_dir + 'test_data/'
    if not os.path.isdir(test_dir):
        os.mkdir(test_dir)

    for data_type in ['drop_seq', 'in_drop']:

        seqc.log.info('SEQC functest for %s' % data_type)

        # define some variables
        output_prefix = package_dir + 'test_data/%s/seqc_test' % data_type
        output_dir = package_dir + 'test_data/%s/' % data_type
        # output_dir = testing_dir + data_type + '/'
        index = package_dir + 'test_data/genome/'
        n_threads = 7

        # check if output_directory exists
        if not os.path.isdir(output_dir):
            seqc.log.info('SEQC test_data: creating testing directory')
            os.mkdir(output_dir)

        # check for necessary genome index files
        if not os.path.isdir(index):
            seqc.log.info('SEQC test_data: downloading mouse chr19 genome index files')
            bucket = 'dplab-data'
            key_prefix = 'genomes/mm38_chr19/'
            prefix = index
            seqc.io.S3.download_files(bucket=bucket, key_prefix=key_prefix,
                                      output_prefix=prefix, cut_dirs=2)

        # define some file names
        forward = output_dir + 'fastq/seqc_test_r1.fastq'
        reverse = output_dir + 'fastq/seqc_test_r2.fastq'
        barcodes = output_dir + 'barcodes/barcodes.p'

        ########################### correct special barcode cases ###########################

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
                bucket = 'dplab-data'
                key = 'barcodes/%s/barcodes.p' % data_type
                seqc.io.S3.download_file(bucket=bucket, key=key, fout=barcodes)

        # define some genome files
        gtf = index + 'annotations.gtf'
        fasta = index + 'mm38_chr19.fa'

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
                prefix = fastq_dir + 'seqc_test'
                generate_fastq(10000, prefix, fasta, gtf, barcodes=barcodes,
                               tag_type='gene_id')

        args = [
            data_type.replace('_', '-'),
            '-i', index,
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
    seqc.log.setup_logger()

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


if __name__ == "__main__":
    seqc.log.setup_logger()
    try:
        nose2.main(exit=False, module=__name__, verbosity=5)
    except:
        seqc.log.exception()
        raise
