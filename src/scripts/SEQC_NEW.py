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
                   choices=['in_drop_patched', 'drop_seq', 'mars1_seq',
                            'mars2_seq', 'in_drop_v2_patched'],
                   help='which platform are you merging annotations from?')
    p.add_argument('-g', '--genomic-fastq', nargs='*', metavar='G', default=[],
                   help='fastq file(s) containing genomic information')
    p.add_argument('-b', '--barcode-fastq', nargs='*', metavar='B', default=[],
                   help='fastq file(s) containing barcode information')
    p.add_argument('-m', '--merged-fastq', nargs='?', metavar='M', default='',
                   help='fastq file containing genomic information annotated with '
                        'barcode data')
    p.add_argument('-s', '--samfile', nargs='?', metavar='S', default='',
                   help='sam file containing aligned, merged fastq records.')
    p.add_argument('--basespace', metavar='BS', help='BaseSpace sample ID. '
                   'Identifies a sequencing run to download and process.')
    p.add_argument('--basespace-token', metavar='BT', help='BaseSpace access '
                   'token')

    p.add_argument('-o', '--output-stem', metavar='O', help='file stem for output files '
                   'e.g. ./seqc_output/tumor_run5')
    p.add_argument('-i', '--index', metavar='I', help='Folder or s3 link to folder '
                   'containing index files for alignment and resolution of ambiguous '
                   'reads.')

    p.add_argument('-v', '--version', action='version',
                   version='{} {}'.format(p.prog, seqc.__version__))

    f = p.add_argument_group('filter arguments')
    f.add_argument('--max-insert-size', metavar='F', help='maximum paired-end insert '
                   'size that is considered a valid record', default=1000)
    f.add_argument('--min-poly-t', metavar='T', help='minimum size of poly-T tail that '
                   'is required for a barcode to be considered a valid record', default=3)
    return p.parse_args()


def main():

    seqc.log.setup_logger()
    try:
        args = parse_args()
        seqc.log.args(args)

        # do a bit of argument checking
        if args.output_stem.endswith('/'):
            print('-o/--output-stem should not be a directory')
            sys.exit(2)

        # split output_stem into path and prefix
        output_dir, output_prefix = os.path.split(args.output_stem)

        if args.basespace:
            seqc.log.info('BaseSpace link provided for fastq argument. Downloading input '
                          'data.')
            if not args.basespace_token:
                raise ValueError(
                    'If the --basespace argument is used, the --basespace-token argument '
                    'must also be provided in order to gain access to the basespace '
                    'repository')
            args.barcode_fastq, args.genomic_fastq = seqc.io.BaseSpace.download_sample(
                args.platform, args.basespace, output_dir, args.basespace_token)

        # check if the index must be downloaded
        if not args.index.startswith('s3://'):
            if not os.path.isdir(args.index):
                raise ValueError('provided index: "%s" is neither an s3 link or a valid '
                                 'filepath' % args.index)
        else:
            try:
                seqc.log.info('AWS s3 link provided for index. Downloading index.')
                bucket, prefix = seqc.io.S3.split_link(args.index)
                args.index = output_dir + 'index/'  # set index  based on s3 download
                cut_dirs = prefix.count('/')
                seqc.io.S3.download_files(bucket, prefix, args.index, cut_dirs)
            except FileNotFoundError:
                raise FileNotFoundError('No index file or folder was identified at the '
                                        'specified s3 index location: %s' % args.index)
            except FileExistsError:
                pass  # file is already present.

        n_processes = multiprocessing.cpu_count() - 1  # get number of processors

        # determine where the script should start:
        merge = True
        align = True
        if args.samfile:
            merge = False
            align = False
        if args.merged_fastq:
            merge = False

        if merge:
            seqc.log.info('Merging genomic reads and barcode annotations.')
            merge_function = getattr(seqc.sequence.merge_functions, args.platform)
            args.merged_fastq = seqc.sequence.fastq.merge_paired(
                    merge_function=merge_function,
                    fout=args.output_stem + '_merged.fastq',
                    genomic=args.genomic_fastq,
                    barcode=args.barcode_fastq)

        if align:
            seqc.log.info('Aligning merged fastq records.')
            *base_directory, stem = args.output_stem.split('/')
            alignment_directory = '/'.join(base_directory) + '/alignments/'
            os.makedirs(alignment_directory, exist_ok=True)
            args.samfile = seqc.align.STAR.align(
                    args.merged_fastq, args.index, n_processes, alignment_directory)

        seqc.log.info('Filtering alignments and constructing record database.')
        # fc = seqc.convert_features.ConvertFeatureCoordinates.from_gtf(
        #         args.index + 'annotations.gtf',
        #         args.max_insert_size)
        fc = seqc.gtf.SCID_set(args.index + 'annotations.gtf')
        fc.slice_SCIDs(start=-args.max_insert_size)  # limit to 3'
        fc.create_interval_tree_scid()
        h5db = seqc.arrays.SimpleReadArray.from_samfile(args.samfile, fc)
        h5db.save_h5(args.output_stem + '.h5')

        seqc.log.info('Resolving ambiguous alignments.')
        # h5db.resolve_alignments_old(expectations=args.index + 'p_coalignment_array.p',
        #                             required_poly_t=args.min_poly_t)
        h5db.resolve_alignments(index=args.index, required_poly_t=args.min_poly_t)

        seqc.log.info('Storing corrected record database.')
        h5db.save_h5(args.output_stem + '.h5')

        seqc.log.info('Generating molecule x cell and read x cell matrices.')
        # mols, mrow, mcol = h5db.unique_features_to_sparse_counts(
        # mols, mrow, mcol = h5db.unique_features_to_sparse_counts_cmp(
        #     collapse_molecules=True,
        #     n_poly_t_required=args.min_poly_t)
        # # reads, rrow, rcol = h5db.unique_features_to_sparse_counts(
        # reads, rrow, rcol = h5db.unique_features_to_sparse_counts_cmp(
        #     collapse_molecules=False,
        #     n_poly_t_required=args.min_poly_t)
        # matrices = {'molecules': {'matrix': mols, 'row_ids': mrow, 'col_ids': mcol},
        #             'reads': {'matrix': reads, 'row_ids': rrow, 'col_ids': rcol}}
        # with open(args.output_stem + '_read_and_mol_matrices.p', 'wb') as f:
        #     pickle.dump(matrices, f)
        matrices = h5db.to_sparse_counts_diff(
                required_support=2, min_poly_t=args.min_poly_t)
        with open(args.output_stem + '_read_and_mol_matrices_diff.p', 'wb') as f:
            pickle.dump(matrices, f)

        seqc.log.info('Run complete.')
    except:
        seqc.log.exception()
        raise


if __name__ == '__main__':
    main()
