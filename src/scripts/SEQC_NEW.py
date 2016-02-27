#!/usr/local/bin/python3

import os
import sys
import argparse
import seqc
import multiprocessing
import pickle


def parse_args():
    p = argparse.ArgumentParser(description='Process Single-Cell RNA Sequencing Data')
    p.add_argument('platform',
                   choices=['in_drop', 'drop_seq', 'mars1_seq',
                            'mars2_seq', 'in_drop_v2'],
                   help='which platform are you merging annotations from?')
    p.add_argument('-g', '--genomic-fastq', nargs='+', metavar='G',
                   help='fastq file(s) containing genomic information')
    p.add_argument('-b', '--barcode-fastq', nargs='?', metavar='B',
                   help='fastq file(s) containing barcode information')
    p.add_argument('-o', '--output-stem', metavar='O', help='file stem for output files '
                   'e.g. ./seqc_output/tumor_run5')
    p.add_argument('-i', '--index', metavar='I', help='Folder containing index files for '
                   'alignment and resolution of ambiguous reads.')

    p.add_argument('-v', '--version', action='version',
                   version='{} {}'.format(p.prog, seqc.__version__))

    f = p.add_argument_group('filter arguments')
    f.add_argument('--max-insert-size', metavar='F', help='maximum paired-end insert '
                   'size that is considered a valid record', default=1000)
    f.add_argument('--min-poly-t', metavar='T', help='minimum size of poly-T tail that '
                   'is required for a barcode to be considered a valid record', default=3)
    return p.parse_args()


def main():

    # do a bit of argument checking
    args = parse_args()
    if args.output_stem.endswith('/'):
        print('-o/--output-stem should not be a directory')
        sys.exit(2)

    # set-up
    seqc.log.setup_logger()
    seqc.log.args(args)
    n_threads = multiprocessing.cpu_count() - 1

    seqc.log.info('Merging genomic reads and barcode annotations.')
    merge_function = getattr(seqc.sequence.merge_functions, args.platform)
    merged = seqc.sequence.fastq.merge_paired(
            merge_function=merge_function,
            fout=args.output_stem + '_merged.fastq',
            genomic=args.genomic_fastq,
            barcode=args.barcode_fastq)

    seqc.log.info('Aligning merged fastq records.')
    *base_directory, stem = args.output_stem.split('/')
    alignment_directory = '/'.join(base_directory) + '/alignments/'
    os.makedirs(alignment_directory, exist_ok=True)
    samfile = seqc.align.STAR.align(merged, args.index, n_threads, alignment_directory)

    seqc.log.info('Filtering alignments and constructing record database.')
    fc = seqc.convert_features.ConvertFeatureCoordinates.from_gtf(
            args.index + 'annotations.gtf',
            args.max_insert_size)
    h5db = seqc.arrays.ReadArray.from_samfile(samfile, fc)

    seqc.log.info('Resolving ambiguous alignments.')
    h5db.resolve_alignments_old(index=args.index, required_poly_t=args.min_poly_t)

    seqc.log.info('Storing corrected record database.')
    h5db.save_h5(args.output_stem + '.h5')

    seqc.log.info('Generating molecule x cell and read x cell matrices.')
    mols, mrow, mcol = h5db.unique_features_to_sparse_counts(
        collapse_molecules=True,
        n_poly_t_required=args.min_poly_t)
    reads, rrow, rcol = h5db.unique_features_to_sparse_counts(
        collapse_molecules=False,
        n_poly_t_required=args.min_poly_t)
    matrices = {'molecules': {'matrix': mols, 'row_ids': mrow, 'col_ids': mcol},
                'reads': {'matrix': reads, 'row_ids': rrow, 'col_ids': rcol}}
    with open(args.output_stem + '_read_and_mol_matrices.p', 'wb') as f:
        pickle.dump(matrices, f)

    seqc.log.info('Run complete.')


if __name__ == '__main__':
    main()
