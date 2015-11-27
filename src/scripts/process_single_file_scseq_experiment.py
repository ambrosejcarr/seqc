#!/usr/local/bin/python3

__author__ = 'ambrose'

import seqc
from scipy.sparse import coo_matrix
from collections import defaultdict
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

def sam_to_count_single_file(sam_file, gtf_file):
    """cannot separate file due to operating system limitations. Instead, implement
    a mimic of htseq-count that uses the default 'union' approach to counting, given the
    same gtf file"""

    # get conversion table, all possible genes for the count matrix
    gt = seqc.convert_features.GeneTable(gtf_file)
    all_genes = gt.all_genes()

    # map genes to ids
    n_genes = len(all_genes)
    gene_to_int_id = dict(zip(sorted(all_genes), range(n_genes)))
    cell_to_int_id = {'no_cell': 0}
    cell_number = 1
    read_count = defaultdict(int)

    # add metadata fields to mimic htseq output; remember to remove these in the final
    # analysis
    gene_to_int_id['ambiguous'] = n_genes
    gene_to_int_id['no_feature'] = n_genes + 1
    gene_to_int_id['not_aligned'] = n_genes + 2

    # estimate the average read length
    with open(sam_file, 'r') as f:
        sequences = []
        line = f.readline()
        while line.startswith('@'):
            line = f.readline()
        while len(sequences) < 100:
            sequences.append(f.readline().strip().split('\t')[9])
        read_length = round(np.mean([len(s) for s in sequences]))

    # pile up counts
    with open(sam_file) as f:
        for record in f:

            # discard headers
            if record.startswith('@'):
                continue
            record = record.strip().split('\t')

            # get cell id, discard if no cell, add new cell if not found.
            cell = record[0].split(':')[0]
            if cell == 0:
                int_cell_id = 0
            else:
                try:
                    int_cell_id = cell_to_int_id[cell]
                except KeyError:
                    cell_to_int_id[cell] = cell_number
                    int_cell_id = cell_number
                    cell_number += 1

            # get start, end, chrom, strand
            flag = int(record[1])
            if flag & 4:
                int_gene_id = n_genes + 2  # not aligned
            else:
                chromosome = record[2]
                if flag & 16:
                    strand = '-'
                    end = int(record[3])
                    start = end - read_length
                else:
                    strand = '+'
                    start = int(record[3])
                    end = start + read_length

                try:
                    genes = gt.coordinates_to_gene_ids(chromosome, start, end, strand)
                except KeyError:
                    continue  # todo include these non-chromosome scaffolds instead of
                    # discarding
                if len(genes) == 1:
                    int_gene_id = gene_to_int_id[genes[0]]
                if len(genes) == 0:
                    int_gene_id = n_genes + 1
                if len(genes) > 1:
                    int_gene_id = n_genes
            read_count[(int_cell_id, int_gene_id)] += 1

    # create sparse matrix
    cell_row, gene_col = zip(*read_count.keys())
    data = list(read_count.values())
    m = cell_number
    n = n_genes + 3

    coo = coo_matrix((data, (cell_row, gene_col)), shape=(m, n), dtype=np.int32)
    gene_index = np.array(sorted(all_genes) + ['ambiguous', 'no_feature', 'not_aligned'],
                          dtype=object)
    cell_index = np.array(['no_cell'] + list(range(1, cell_number)), dtype=object)

    return coo, gene_index, cell_index


def main(srp, n_threads, s3_bucket, s3_key, cell_barcodes,
         experiment_type, experiment_name, index=None, index_bucket=None, index_key=None,
         working_directory='./'):

    if not working_directory.endswith('/'):
        working_directory += '/'

    # set the index
    if not index:  # download the index
        log_info('Downloading Index from S3')
        index_dir = working_directory + 'index/'
        seqc.io_lib.S3.download_files(bucket=index_bucket, key_prefix=index_key,
                          output_prefix=index_dir, cut_dirs=False)
        index = index_dir + index_key.lstrip('/')
    if not os.path.isdir(index):
        log_info('Using local index')
        raise FileNotFoundError('Index does not lead to a directory')

    # download the data
    log_info('Downloading SRA Data')
    files = seqc.io_lib.GEO.download_srp(
        srp, working_directory, min(n_threads, 10), verbose=False,
                             clobber=False)

    # unpack the .sra files into forward and reverse fastq files
    log_info('Unpacking SRA to fastq')
    forward, reverse = seqc.io_lib.GEO.extract_fastq(files, n_threads, working_directory,
                                                     verbose=False)

    # merge fastq files
    log_info('Extracting cell information and merging fastq files')
    merged_fastq, _ = merge_fastq(
        forward, reverse, experiment_type, working_directory, cell_barcodes)

    # align the data
    log_info('Aligning fastq records')
    sam_file = seqc.align.STAR.align(merged_fastq, index, n_threads, working_directory,
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
    seqc.io_lib.S3.upload_file(numpy_archive, s3_bucket, s3_key)

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

