#!/usr/local/bin/python3

__author__ = 'ambrose'

from seqc.io_lib import S3, GEO
from seqc.align import STAR
from seqc.qc import sam_to_count_single_file
from seqc.fastq import merge_fastq
import os
import json
import numpy as np
import argparse
import logging
from datetime import datetime


def setup_logger():
    """create a simple log file in the cwd to track progress and any errors"""
    logging.basicConfig(filename='seqc.log', level=logging.DEBUG)


def log_info(message):
    """print a timestamped update for the user"""
    logging.info(datetime.now().strftime("%Y-%m-%d %H:%M:%S") + ':' + message)


# todo think about adding paired-end support
def parse_args():

    p = argparse.ArgumentParser()
    p.add_argument('-s', '--srp', help='FTP link to SRP experiment to download from GEO',
                   metavar='S', required=True, type=str)
    p.add_argument('-n', '--n-threads', help='number of threads to use', metavar='N',
                   required=True, type=int)
    p.add_argument('--s3-bucket', help='s3 bucket to upload data matrix', metavar='B',
                   required=True, type=str)
    p.add_argument('--s3-key', help='s3 key to upload data matrix', metavar='K',
                   required=True, type=str)
    p.add_argument('-e', '--experiment-name', help='stem for output data matrix',
                   metavar='E', required=True, type=str)
    p.add_argument('-i', '--index', help='location of directory containing star index',
                   metavar='I', type=str)
    p.add_argument('--index-bucket', help='s3 bucket for star index', metavar='IB',
                   type=str)
    p.add_argument('--index-key', help='s3 key for star index', metavar='IK', type=str)
    p.add_argument('-w', '--working-directory', metavar='W', type=str,
                   help='temporary working directory for script', required=True)
    p.add_argument('-t', '--experiment-type', help='pre-processor name', metavar='T',
                   type=str)
    args = vars(p.parse_args())
    return args


def main(srp, n_threads, s3_bucket, s3_key, cell_barcodes,
         experiment_type, experiment_name, index=None, index_bucket=None, index_key=None,
         working_directory='./'):

    if not working_directory.endswith('/'):
        working_directory += '/'

    # set the index
    if not index:  # download the index
        log_info('Downloading Index from S3')
        index_dir = working_directory + 'index/'
        S3.download_files(bucket=index_bucket, key_prefix=index_key,
                          output_prefix=index_dir, cut_dirs=False)
        index = index_dir + index_key.lstrip('/')
    if not os.path.isdir(index):
        log_info('Using local index')
        raise FileNotFoundError('Index does not lead to a directory')

    # download the data
    log_info('Downloading SRA Data')
    files = GEO.download_srp(srp, working_directory, min(n_threads, 10), verbose=False,
                             clobber=False)

    # unpack the .sra files into forward and reverse fastq files
    log_info('Unpacking SRA to fastq')
    forward, reverse = GEO.extract_fastq(files, n_threads, working_directory,
                                         verbose=False)

    # merge fastq files
    log_info('Extracting cell information and merging fastq files')
    merged_fastq, _ = merge_fastq(
        forward, reverse, experiment_type, working_directory, cell_barcodes)

    # align the data
    log_info('Aligning fastq records')
    sam_file = STAR.align(merged_fastq, index, n_threads, working_directory,
                          paired_end=False)

    # create the matrix
    log_info('Creating counts matrix')
    gtf_file = index + 'annotations.gtf'
    coo, rowind, colind = sam_to_count_single_file(sam_file, gtf_file)

    log_info('Saving counts matrix')
    numpy_archive = experiment_name + '.npz'
    with open(numpy_archive, 'wb') as f:
        np.savez(f, mat=coo, row=rowind, col=colind)

    log_info('Uploading counts matrix to S3')
    # upload the matrix to amazon s3
    S3.upload_file(numpy_archive, s3_bucket, s3_key)

    log_info('Run completed')

if __name__ == "__main__":
    kwargs = parse_args()
    setup_logger()
    log_info('Passed command line arguments: %s' %
             json.dumps(kwargs, separators=(',', ': '), indent=4))
    try:
        main(**kwargs)
    except:
        logging.exception(datetime.now().strftime("%Y-%m-%d %H:%M:%S") + ':main:')
        raise

