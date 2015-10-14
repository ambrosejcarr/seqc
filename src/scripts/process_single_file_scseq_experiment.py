#!/usr/local/bin/python3

__author__ = 'ambrose'

from seqc.io_lib import S3, GEO
from seqc.align import STAR
from seqc.qc import sam_to_count_single_file
from seqc.fastq import merge_fastq
from glob import glob
import os
import numpy as np


def main(srp, n_threads, output_directory, index_key, index_bucket, s3_bucket, s3_key,
         cell_barcodes, experiment_type, experiment_name):
    # download the index
    index_dir = output_directory + 'index/'
    S3.download_files(bucket=index_bucket, key_prefix=index_key, output_prefix=index_dir)

    # todo index_dir will not be correct unless we cut the proper number of directories
    # are these cut by default based on how I wrote the program?

    # download the data
    GEO.download_srp(srp, output_directory, min(n_threads, 10), verbose=False,
                     clobber=False)
    # get the downloaded .sra files
    files = [f for f in os.listdir(output_directory) if os.path.isfile(f) and
             f.endswith('.sra')]

    # unpack the .sra files into forward and reverse fastq files
    # todo this needs an output prefix!
    GEO.extract_fastq(files, n_threads)

    # get the fastq files
    forward = glob('*_1.fastq')
    assert(len(forward) == 1)
    reverse = glob('*_2.fastq')
    assert(len(reverse) == 1)

    # merge fastq files
    merged_fastq = merge_fastq(forward, reverse, experiment_type, output_directory,
                               cell_barcodes)

    # align the data
    # todo make STAR arguments into static/class methods where appropriate.
    star = STAR(output_directory, n_threads, index_dir)
    # todo aligner needs to support paired-end reads! create new functions for this
    # todo these functions should output .sam filenames
    sam_file = star.align(merged_fastq)

    # todo estimate read length
    read_length = 100

    # create the matrix
    gtf_file = index_dir + 'annotations.gtf'
    coo, rowind, colind = sam_to_count_single_file(sam_file, gtf_file, read_length)

    numpy_archive = experiment_name + '.npz'
    with open(numpy_archive, 'wb') as f:
        np.savez(f, mat=coo, row=rowind, col=colind)

    # upload the matrix to amazon s3
    S3.upload_file(numpy_archive, s3_bucket, s3_key)

if __name__ == "__main__":
    main()
