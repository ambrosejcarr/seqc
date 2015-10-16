#!/usr/local/bin/python3

__author__ = 'ambrose'

from seqc.io_lib import S3, GEO
from seqc.align import STAR
from seqc.qc import sam_to_count_multiple_files
import os
import numpy as np
import argparse


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
    args = vars(p.parse_args())
    return args


def main(srp, n_threads, s3_bucket, s3_key, experiment_name, index_key=None,
         index_bucket=None, index=None, working_directory='./'):

    # set the index
    if not index:  # download the index
        index_dir = working_directory + 'index/'
        S3.download_files(bucket=index_bucket, key_prefix=index_key,
                          output_prefix=index_dir, cut_dirs=False)
        index = index_dir + index_key.lstrip('/')
    if not os.path.isdir(index):
        raise FileNotFoundError('Index does not lead to a directory')

    # download the data
    files = GEO.download_srp(srp, working_directory, min(n_threads, 10), verbose=False,
                             clobber=False)

    # unpack the .sra files into forward and reverse fastq files
    forward, reverse = GEO.extract_fastq(files, n_threads)

    # align the data
    sam_files = STAR.align_multiple_files(
        forward, index, n_threads, working_directory, reverse_fastq_files=reverse)

    # create the matrix
    gtf_file = index_dir + 'annotations.gtf'
    coo, rowind, colind = sam_to_count_multiple_files(sam_files, gtf_file)

    numpy_archive = experiment_name + '.npz'
    with open(numpy_archive, 'wb') as f:
        np.savez(f, mat=coo, row=rowind, col=colind)

    # upload the matrix to amazon s3
    S3.upload_file(numpy_archive, s3_bucket, s3_key)


if __name__ == "__main__":
    kwargs = parse_args()
    main(**kwargs)

