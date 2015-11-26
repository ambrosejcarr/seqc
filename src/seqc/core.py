__author__ = 'ambrose'

import argparse
import pickle
import os
import shutil
import seqc
import sys
import json


def create_parser():

    ################################# DEFINE ARGS ######################################

    parser = argparse.ArgumentParser()

    # add subparsers for each library construction method
    subparsers = parser.add_subparsers(help='library construction method types',
                                       dest='subparser_name')
    parse_in_drop = subparsers.add_parser('in-drop', help='in-drop help')
    parse_drop_seq = subparsers.add_parser('drop-seq', help='drop-seq help')
    parse_mars_seq = subparsers.add_parser('mars-seq', help='mars-seq help')
    parse_cel_seq = subparsers.add_parser('cel-seq', help='cel-seq help')
    parse_avo_seq = subparsers.add_parser('avo-seq', help='avo-seq help')
    parse_strt_seq = subparsers.add_parser('strt-seq', help='strt-seq help')

    # get a list of parsers, set-up pipeline function for each parser
    subparser_list = [parse_in_drop, parse_mars_seq, parse_cel_seq,
                      parse_avo_seq, parse_strt_seq, parse_drop_seq]
    default_functions = [in_drop, mars_seq, cel_seq, avo_seq, strt_seq, drop_seq]
    for p, f in zip(subparser_list, default_functions):
        p.set_defaults(func=f)

    # set barcode default for drop_seq:
    parse_drop_seq.set_defaults(barcodes=[])

    # set required arguments for all parsers
    for i, p in enumerate(subparser_list):
        r = p.add_argument_group('Required Arguments')
        # todo needs to take input from user on what organism it should be
        r.add_argument('-i', '--index', metavar='I', help='local or s3 location of star '
                       'alignment index folder. This folder will be created if it does '
                       'not exist', default=None)
        r.add_argument('-n', '--n-threads', help='number of threads to run', metavar='N',
                       type=int, default=None)
        r.add_argument('-o', '--output-file', metavar='O', default=None,
                       help='stem of filename in which to store output')

        # for all experiments except drop-seq, barcodes are a required input argument
        if i < 5:
            r.add_argument('-b', '--barcodes', metavar='B', default=None,
                           help='local or s3 location of serialized barcode object.')

        i = p.add_argument_group(
            title='Input Files',
            description='pass one input file type: sam (-s), raw fastq (-f, [-r]), or '
                        'processed fastq (-m)')

        i.add_argument('-f', '--forward', help='forward fastq file(s)', metavar='F',
                       nargs='*', default=None)
        i.add_argument('-r', '--reverse', help='reverse fastq file(s)', metavar='R',
                       nargs='*', default=None)
        i.add_argument('-s', '--samfile', metavar='S', nargs='?', default=None,
                       help='sam file(s) containing aligned, pre-processed reads')
        i.add_argument('-m', '--merged-fastq', metavar='M', default=None,
                       help='fastq file containing merged, pre-processed records')

        # disambiguation arguments
        d = p.add_argument_group('Optional arguments for disambiguation')
        d.add_argument('-l', '--frag-len', metavar='L', type=int, default=1000,
                       help='the number of bases from the 3 prime end to '
                       'consider when determining trancript overlaps')

        # alignment arguments
        a = p.add_argument_group('Optional arguments for STAR aligner')
        a.add_argument('--star-args', metavar='SA', nargs='+', default={},
                       help='additional arguments for STAR. Pass as arg=value without '
                            'leading "--". e.g. runMode=alignReads')
        a.add_argument('--list-default-star-args', default=False, action='store_true',
                       help='list SEQDB default args for the STAR aligner')

    # print help on call of "SEQC" with missing -h
    if len(sys.argv) == 1 and sys.argv[0].endswith('SEQC'):
        parser.print_help()
        sys.exit(2)

    # add a sub-parser for building the index
    pindex = subparsers.add_parser('index', help='SEQC index functions')
    pindex.add_argument('-b', '--build', action='store_true', default=False,
                        help='build a SEQC index')
    pindex.add_argument('-t', '--test', action='store_true', default=False,
                        help='test a SEQC index')
    pindex.add_argument('-o', '--organism', required=True, nargs='+', metavar='O',
                        help='build index for these organism(s)')
    pindex.add_argument('-i', '--index', help='name of folder where index should be '
                        'built or containing the index to be verified',
                        required=True, metavar='I')
    pindex.add_argument('-n', '--n-threads', type=int, default=4, help='number of threads'
                        ' to use when building index', metavar='N')
    pindex.add_argument('--phix', help='add phiX to the genome index and GTF file.',
                        action='store_true', default=False)

    # allow user to check version
    parser.add_argument('-v', '--version', help='print version and exit',
                        action='store_true', default=False)

    return parser


def parse_args(parser, args=None):

    arguments = parser.parse_args(args)

    if arguments.version:
        print('SEQC version: %s' % seqc.__version__)
        exit(2)

    if arguments.subparser_name == 'index':
        if arguments.build and not arguments.test:
            arguments.func = seqc.align.STAR.build_index
        elif arguments.test and not arguments.build:
            arguments.func = seqc.align.STAR.test_index
        else:
            print('SEQC index: error: one but not both of the following arguments must '
                  'be provided: -b/--build, -t/--test')
            sys.exit(2)
    else:

        # list star args if requested, then exit
        if arguments.list_default_star_args:
            printable_args = json.dumps(
                seqc.align.STAR.default_alignment_args('$FASTQ', '$N_THREADS', '$INDEX',
                                                       './$EXP_NAME/'),
                separators=(',', ': '), indent=4)
            print(printable_args)
            sys.exit(2)

        # check that at least one input argument was passed:
        check = [arguments.forward, arguments.reverse, arguments.samfile,
                 arguments.merged_fastq]
        if not any(check):
            print('SEQC %s: error: one or more of the following arguments must be '
                  'provided: -f/--forward, -r/--reverse, -m/--merged-fastq, -s/--sam' %
                  arguments.subparser_name)
            sys.exit(2)
        required = [arguments.output_file, arguments.index, arguments.n_threads]
        if not arguments.subparser_name == 'drop-seq':
            if not all(required + [arguments.barcodes]):
                print('SEQC %s: error: the following arguments are required: -i/--index, '
                      '-n/--n-threads, -o/--output-file, -b/--barcodes')
                sys.exit(2)
        else:
            if not all(required):
                print('SEQC %s: error: the following arguments are required: -i/--index, '
                      '-n/--n-threads, -o/--output-file')

    return vars(arguments)


def set_up(output_file, index, barcodes):
    """
    create temporary directory for run, find gtf in index, and create or load a
    serialized barcode object
    """

    # temporary directory should be made in the same directory as the output prefix
    *stem, final_dir = output_file.split('/')
    temp_dir = '/'.join(stem + ['.' + final_dir])

    # create temporary directory based on experiment name
    if not os.path.isdir(temp_dir):
        os.makedirs(temp_dir)
    if not temp_dir.endswith('/'):
        temp_dir += '/'

    # check that index exists. If index is an aws link, download the index
    if index.startswith('s3://'):
        seqc.log.info('AWS s3 link provided for index. Downloading index.')
        bucket, prefix = seqc.io_lib.S3.split_link(index)
        index = temp_dir + 'index/'
        cut_dirs = prefix.count('/')
        seqc.io_lib.S3.download_files(bucket, prefix, index, cut_dirs)

    # obtain gtf file from index argument
    gtf = index + 'annotations.gtf'
    if not os.path.isfile(gtf):
        raise FileNotFoundError('no file named "annotations.gtf" found in index: %s' % index)

    # get cell barcode files
    if not barcodes:
        cb = seqc.barcodes.DropSeqCellBarcodes()
    else:
        if barcodes.startswith('s3://'):
            seqc.log.info('AWS s3 link provided for barcodes. Downloading barcodes')
            bucket, key = seqc.io_lib.S3.split_link(barcodes)
            output_file = temp_dir + 'barcodes.p'
            try:
                barcodes = seqc.io_lib.S3.download_file(bucket, key, output_file)
            except FileExistsError:
                barcodes = output_file
                pass  # already have the file in this location from a previous run

        # now, we should definitely have the binary file. Load it.
        with open(barcodes, 'rb') as f:
            cb = pickle.load(f)

    return temp_dir, gtf, cb, index

# todo commented out until dependencies are installed automatically.
# todo should be coded to pep8 standards
# def htqc(output_dir, forward, reverse, r, t):
#     seqc.log.info('getting QC')
#     def run_htqc():
#         call(["ht-sample","-P","-q","-o",output_dir + "sample","-r",str(r),"-z","-i"] + [i for i in forward] + [i for i in reverse])
#         call(["ht-stat","-P","-o",output_dir,"-t",str(t),"-z","-i",output_dir + "sample_1.fastq.gz",output_dir + "sample_2.fastq.gz"])
#         call(["ht-stat-draw.pl", "--dir", output_dir])
#         os.remove(output_dir + "sample_1.fastq.gz")
#         os.remove(output_dir + "sample_2.fastq.gz")
#     p = Process(target=run_htqc, args=(), daemon=True)
#     p.start()


def merge(forward, reverse, samfile, merged_fastq, processor, temp_dir, cb, n_threads):
    if (forward or reverse) and not (samfile or merged_fastq):

        # pre-process reads
        seqc.log.info('Merging fastq files.')
        merged_fastq = seqc.fastq.merge_fastq(forward, reverse, processor, temp_dir, cb,
                                              n_threads)
    return merged_fastq


def align(merged_fastq, samfile, star_args, temp_dir, n_threads, index, output_file):
    if merged_fastq and not samfile:
        # process any additional arguments for star passed from the command line
        if star_args:
            kwargs = {}
            for arg in star_args:
                k, v = arg.split('=')
                kwargs['--' + k] = v

        # align fastq files
        seqc.log.info('Aligning merged fastq file')
        seqc.align.STAR.align(merged_fastq, index, n_threads, temp_dir, **kwargs)
        samfile = temp_dir + 'Aligned.out.sam'

        # copy alignment summary
        shutil.copyfile(temp_dir + 'Log.final.out', output_file + '_alignment_summary.txt')
    return samfile


def process_samfile(samfile, output_file, n_threads, gtf, frag_len):
    seqc.log.info('Post-processing alignments')
    h5_name = output_file + '.h5'
    chunk_size = int(2e8)  # ~ 10GB of data in memory at any time
    h5name = seqc.sam.to_h5(samfile, h5_name, n_threads, chunk_size, gtf, frag_len)
    read_array = seqc.arrays.ReadArray.from_h5(h5name)
    return read_array


def correct_errors():
    # try:
    #     info('Correcting errors')
    #     if not barcode_dir.endswith('/'):
    #         barcode_dir += '/'
    #     barcode_files = glob(barcode_dir + '*')
    #     corrected, nerr, err_rate = correct_errors(
    #         unambiguous, barcode_files, processor_name, meta, p_val=0.1)
    #
    #     # store metadata from these runs in post-processing summary
    #     # post_processing_summary(h5db, exp_name, dis_meta, nerr)
    #
    #     # store error rates
    #     # store_error_rate(h5db, exp_name, err_rate)
    #
    # except:  # fail gracefully on error_correction; not all methods can do this!
    #     logging.exception(datetime.now().strftime("%Y-%m-%d %H:%M:%S") + ':main:')
    #     corrected = unambiguous
    raise NotImplementedError


def resolve_alignments(index, arr, n, output_file):
    seqc.log.info('Resolving ambiguous alignments.')
    try:
        # todo
        # I should store this in a jagged array as well. Indexing is 1/2 as fast
        # but it would still be very quick, and VERY memory efficient.
        expectations = index + 'p_coalignment_array.p'
        arr.resolve_alignments(expectations, required_poly_t=n)
    except:
        seqc.log.info('Caught error in resolve_alignments(), saving data')
        arr.save_h5(output_file + '.h5')
        raise


def save_counts_matrices(output_file, arr, n):
    seqc.log.info('Generating Gene x Cell SparseMatrices for reads and molecules.')
    uniq = arr.to_unique(n_poly_t_required=n)
    del arr
    experiment = uniq.to_experiment()
    experiment.to_npz(output_file + '_sp_counts.npz')


def select_cells():
    ################################### SELECT CELLS ####################################
    # # filter cells with < threshold # reads. todo | add gc, length, maybe size biases
    # mols, meta = filter_cells_mols(mols, meta, 10)
    # reads, meta = filter_cells_reads(reads, meta, 100)
    #
    # # store files
    # if location:
    #     location += '/'
    # with open(location + exp_name + '_read_and_mol_frames.p', 'wb') as f:
    #     pickle.dump(((reads, rrow, rcol), (mols, mrow, mcol)), f)
    raise NotImplementedError


def store_results(output_file, arr):
    if not output_file.endswith('.h5'):
        output_file += '.h5'
    seqc.log.info('Storing processed data in %s' % output_file)
    arr.save_h5(output_file)


def run_complete():
    seqc.log.info('Run complete.')


def clean_up(temp_dir):
    seqc.log.info('Run succeeded. Cleaning up and terminating.')
    if os.path.isdir(temp_dir):
        shutil.rmtree(temp_dir)
    if os.path.isdir('_STARtmp'):
        shutil.rmtree('_STARtmp')


def in_drop(output_file, forward, reverse, samfile, merged_fastq, subparser_name, index,
            n_threads, frag_len, star_args, barcodes, **kwargs):

    temp_dir, gtf, cb, index = set_up(output_file, index, barcodes)

    # htqc(temp_dir, forward, reverse, 1000, 1)

    merged_fastq = merge(forward, reverse, samfile, merged_fastq, subparser_name,
                         temp_dir, cb, n_threads)

    samfile = align(merged_fastq, samfile, star_args, temp_dir, n_threads, index,
                    output_file)

    arr = process_samfile(samfile, output_file, n_threads, gtf, frag_len)

    resolve_alignments(index, arr, n=0, output_file=output_file)

    store_results(output_file, arr)

    save_counts_matrices(output_file, arr, n=3)

    run_complete()


def drop_seq(output_file, forward, reverse, samfile, merged_fastq, subparser_name, index,
             n_threads, frag_len, star_args, barcodes, **kwargs):

    temp_dir, gtf, cb, index = set_up(output_file, index, barcodes)

    # htqc(temp_dir, forward, reverse, 1000, 1)

    merged_fastq = merge(forward, reverse, samfile, merged_fastq, subparser_name,
                         temp_dir, cb, n_threads)

    samfile = align(merged_fastq, samfile, star_args, temp_dir, n_threads, index,
                    output_file)

    arr = process_samfile(samfile, output_file, n_threads, gtf, frag_len)

    resolve_alignments(index, arr, n=0, output_file=output_file)

    arr.save_h5(output_file + '.h5')

    save_counts_matrices(output_file, arr, n=0)

    run_complete()


def mars_seq():
    raise NotImplementedError


def cel_seq():
    raise NotImplementedError


def avo_seq():
    raise NotImplementedError


def strt_seq():
    """
    This method assumes STRT-seq datasets are being retrieved from GEO, where Sten
    Linnarsson has demultiplexed the data into multiple, individual fastq files."""
    raise NotImplementedError
