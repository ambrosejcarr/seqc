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
    parse_drop_seq.set_defaults(barcodes='')

    # set required arguments for all parsers
    for i, p in enumerate(subparser_list):
        r = p.add_argument_group('Required Arguments')
        # todo needs to take input from user on what organism it should be
        r.add_argument('-i', '--index', metavar='I', help='local or s3 location of star '
                       'alignment index folder. This folder will be created if it does '
                       'not exist', default=None)
        r.add_argument('-n', '--n-threads', help='number of threads to run', metavar='N',
                       type=int, default=None)
        r.add_argument('-o', '--output-prefix', metavar='O', default=None,
                       help='stem of filename in which to store output')
        r.add_argument('-c', '--cluster-name', metavar='C', default=None,
                       help='optional name for aws cluster')
        r.add_argument('-k', '--aws-upload-key', metavar='K', default=None,
                       help='upload processed results to this AWS folder. Required if '
                            '--remote is passed')
        r.add_argument('--remote', default=False, action='store_true',
                       help='run the requested SEQC command remotely')
        r.add_argument('--email-status', default='', metavar='E',
                       help='email results to this address')

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
    pindex.add_argument('-c', '--cluster-name', metavar='C', default=None,
                        help='optional name for aws cluster')
    pindex.add_argument('-k', '--aws-upload-key', metavar='K', default=None,
                        help='upload constructed index to this AWS folder. Required if '
                             '--remote is passed')
    pindex.add_argument('--remote', default=False, action='store_true',
                        help='run the requested SEQC command remotely')
    pindex.add_argument('--email-status', default='', metavar='E',
                        help='email results to this address')


    # allow user to check version
    parser.add_argument('-v', '--version', help='print version and exit',
                        action='store_true', default=False)
    # allow indication of remote run.


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
        required = [arguments.output_prefix, arguments.index, arguments.n_threads]
        if not arguments.subparser_name == 'drop-seq':
            if not all(required + [arguments.barcodes]):
                print('SEQC %s: error: the following arguments are required: -i/--index, '
                      '-n/--n-threads, -o/--output-prefix, -b/--barcodes')
                sys.exit(2)
        else:
            if not all(required):
                print('SEQC %s: error: the following arguments are required: -i/--index, '
                      '-n/--n-threads, -o/--output-prefix')

    if arguments.remote:
        if not arguments.aws_upload_key:
            print('SEQC: %s: error: if requesting a remote run with --remote, '
                  '-k/--aws-upload-key must be specified. Otherwise results are '
                  'discarded when the cluster is terminated.' % arguments.subparser_name)
            sys.exit(2)
        if not arguments.email_status:
            print('SEQC: %s: error: if requesting a remote run with --remote, '
                  '--email-status must specify the email address that updates or errors '
                  'should be sent to.')
            sys.exit(2)

    return vars(arguments)


def run_remote(kwargs: dict) -> None:
    """
    todo document me!

    args:
    -----

    returns:
    --------
    None
    """
    cmd = 'SEQC '

    # get the positional argument; doesn't need a '--' prefix
    positional = kwargs['subparser_name']
    clustname = kwargs['cluster_name']

    cmd += positional + ' '
    del kwargs['subparser_name']
    del kwargs['func']
    del kwargs['cluster_name']

    for k, v in kwargs.items():
        if isinstance(v, list):
            v = ' '.join(v)  # lists of input files should be merged with whitespace
        if v:
            edit_k = '-'.join(k.split('_'))
            cmd += '--%s %s ' % (edit_k, v)
    print('cmd: %s' %cmd)

    # set up remote cluster here, finishes all the way through gitpull
    cluster = seqc.cluster_utils.ClusterServer()
    cluster.cluster_setup(clustname)
    cluster.serv.connect()
    seqc.log.info('Remote server set-up complete.')

    # todo
    # figure out how to download files from basespace here; hard-coded inputs right now
    cluster.serv.put_file('/Users/kristyc/PycharmProjects/seqc/src/scripts/'
                          'short_f1.fastq','short_f1.fastq')
    cluster.serv.put_file('/Users/kristyc/PycharmProjects/seqc/src/scripts/'
                          'short_r1.fastq','short_r1.fastq')
    # cluster.serv.put_file('/Users/kristyc/PycharmProjects/seqc/src/scripts/notify.py',
    #                       'notify.py')
    # cluster.serv.exec_command('mv notify.py /data/software/notify.py')
    cluster.serv.exec_command('mv short_f1.fastq /data/software/short_f1.fastq')
    cluster.serv.exec_command('mv short_r1.fastq /data/software/short_r1.fastq')
    seqc.log.info('Remote file download complete.')

    #running SEQC on the cluster
    # cmdstring = "SEQC in-drop --forward short_f1.fastq --index
    # s3://dplab-data/genomes/mm38/ --frag-len 1000 --reverse short_r1.fastq --n-threads
    # 30 --barcodes s3://dplab-data/sc-seq/allon/barcodes/in_drop_barcodes.p
    # --output-prefix /data/software/qq"
    seqc.log.info('Beginning remote run.')
    print('cmd: %s' %cmd)
    cluster.serv.exec_command(cmd)
    #cluster.serv.exec_command('nohup %s > dev/null 2>&1 &' % cmd)
    # cluster.serv.exec_command('nohup %s > /data/software/nohup.txt' % cmd)
    # seqc.log.info('Terminating local client. Email will be sent when remote run '
    #               'completes')

    #shut down cluster after finished
    # todo | have local program exit while process still runs remotely
    # print('shutting down cluster...')
    # if out:
    #     cluster.terminate_cluster()
    #     sys.exit("process complete -- program exiting")
    # elif err:
    #     cluster.terminate_cluster()
    #     sys.exit("process complete -- program exiting")


def fix_output_paths(output_prefix: str) -> (str, str):
    """
    Returns an output prefix and output directory with absolute paths

    args:
    -----
    output_prefix: str; prefix for all seqc output files; there will also be a
      folder created with the same name to house some output files with non-unique names

    returns:
    --------
    absolute_output_prefix, absolute_output_directory
    """

    seqc.util.check_type(output_prefix, str, 'output_prefix')

    if output_prefix.endswith('/'):
        raise ValueError('Invalid output_prefix: "%s". output prefix must be a prefix,'
                         'not a directory name' % output_prefix)

    output_prefix = os.path.expanduser(output_prefix)
    output_prefix = os.path.abspath(output_prefix)

    # create the directory tree if it doesn't already exist
    output_directory = output_prefix + '/'
    if not os.path.isdir(output_directory):
        os.makedirs(output_directory)

    return output_prefix, output_directory


def check_index(index: str, output_dir: str='') -> (str, str):
    """
    Checks the provided index parameter.

    If index resembles a file path, makes sure it is present. If it resembles an amazon
    s3 link, downloads the index to output_dir. Checks the provided index for critical
    files that are necessary to run SEQC.

    args:
    -----
    index: str; either file path or s3 link pointing to the seqc index directory
    output_dir: str or None; directory where downloaded index should be placed.

    returns:
    --------
    index, gtf
    """

    seqc.util.check_type(index, str, 'index')
    seqc.util.check_type(output_dir, str, 'output_dir')

    if not index.endswith('/'):
        index += '/'

    critical_index_files = ['SA', 'SAindex', 'Genome', 'annotations.gtf',
                            'p_coalignment.pckl']

    if not index.startswith('s3://'):  # index is a file path
        if not os.path.isdir(index):
            raise ValueError('provided index: "%s" is neither an s3 link or a valid '
                             'filepath' % index)
        else:
            pass  # index points to a valid folder

    else:  # index is an aws link
        try:
            seqc.log.info('AWS s3 link provided for index. Downloading index.')
            bucket, prefix = seqc.io.S3.split_link(index)
            index = output_dir + 'index/'  # set index directory based on s3 download
            cut_dirs = prefix.count('/')
            seqc.io.S3.download_files(bucket, prefix, index, cut_dirs)
        except FileNotFoundError:  # index does not exist in the specified location
            raise FileNotFoundError('No index file or folder was identified at the '
                                    'specified s3 index location: %s' % index)
        except FileExistsError:
            pass  # file is already present.

    # check that the index contains the necessary files to run SEQC
    for f in critical_index_files:
        if not os.path.isfile(index + f):
            raise FileNotFoundError('Index is missing critical file "%s". Please '
                                    'regenerate the index.')

    # obtain gtf file from index argument
    gtf = index + 'annotations.gtf'

    return index, gtf


def check_and_load_barcodes(
        barcodes: str='', output_dir='') -> seqc.barcodes.CellBarcodes:
    """
    check if barcodes points to a valid s3 link or file object. If it does, load and
    return the CellBarcodes object.
    """

    seqc.util.check_type(barcodes, str, 'barcodes')

    # get cell barcode files
    if not barcodes:
        return seqc.barcodes.DropSeqCellBarcodes()
    elif barcodes.startswith('s3://'):
        seqc.log.info('AWS s3 link provided for barcodes. Downloading barcodes')
        bucket, key = seqc.io.S3.split_link(barcodes)
        output_prefix = output_dir + 'barcodes.p'
        try:
            barcodes = seqc.io.S3.download_file(bucket, key, output_prefix)
        except FileExistsError:
            barcodes = output_prefix
            pass  # already have the file in this location from a previous run
        except FileNotFoundError:
            raise FileNotFoundError('No barcode file was identified at the '
                                    'specified s3 barcodes location: %s' % barcodes)
    elif not os.path.isfile(barcodes):
        raise FileNotFoundError('No barcode file was found at %s' % barcodes)

    with open(barcodes, 'rb') as f:
        cb = pickle.load(f)

    if not isinstance(cb, seqc.barcodes.CellBarcodes):
        raise TypeError('specified barcodes file %s did not contain a CellBarcodes '
                        'object' % barcodes)

    return cb


def set_up(output_prefix, index, barcodes):
    """
    create temporary directory for run, find gtf in index, and create or load a
    serialized barcode object
    """

    # temporary directory should be made in the same directory as the output prefix
    *stem, final_prefix = output_prefix.split('/')
    output_dir = '/'.join(stem + ['.' + final_prefix])

    # create temporary directory based on experiment name
    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)
    if not output_dir.endswith('/'):
        output_dir += '/'

    # check that index exists. If index is an aws link, download the index
    if index.startswith('s3://'):
        seqc.log.info('AWS s3 link provided for index. Downloading index.')
        bucket, prefix = seqc.io.S3.split_link(index)
        index = output_dir + 'index/'
        cut_dirs = prefix.count('/')
        seqc.io.S3.download_files(bucket, prefix, index, cut_dirs)

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
            bucket, key = seqc.io.S3.split_link(barcodes)
            output_prefix = output_dir + 'barcodes.p'
            try:
                barcodes = seqc.io.S3.download_file(bucket, key, output_prefix)
            except FileExistsError:
                barcodes = output_prefix
                pass  # already have the file in this location from a previous run

        # now, we should definitely have the binary file. Load it.
        with open(barcodes, 'rb') as f:
            cb = pickle.load(f)

    return output_dir, gtf, cb, index

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


def align(merged_fastq, samfile, star_args, temp_dir, n_threads, index, output_prefix):
    if merged_fastq and not samfile:
        # process any additional arguments for star passed from the command line

        kwargs = {}
        if star_args:
            for arg in star_args:
                k, v = arg.split('=')
                kwargs['--' + k] = v

        # align fastq files
        seqc.log.info('Aligning merged fastq file')
        seqc.align.STAR.align(merged_fastq, index, n_threads, temp_dir, **kwargs)
        samfile = temp_dir + 'Aligned.out.sam'

        # copy alignment summary
        shutil.copyfile(temp_dir + 'Log.final.out', output_prefix + '_alignment_summary.txt')
    return samfile


def process_samfile(samfile, output_prefix, n_threads, gtf, frag_len):
    seqc.log.info('Post-processing alignments')
    h5_name = output_prefix + '.h5'
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


def resolve_alignments(index, arr, n, output_prefix):
    seqc.log.info('Resolving ambiguous alignments.')
    try:
        # todo
        # I should store this in a jagged array as well. Indexing is 1/2 as fast
        # but it would still be very quick, and VERY memory efficient.
        expectations = index + 'p_coalignment_array.p'
        arr.resolve_alignments(expectations, required_poly_t=n)
    except:
        seqc.log.info('Caught error in resolve_alignments(), saving data')
        arr.save_h5(output_prefix + '.h5')
        raise


def save_counts_matrices(output_prefix, arr, n):
    seqc.log.info('Generating Gene x Cell SparseMatrices for reads and molecules.')
    uniq = arr.to_unique(n_poly_t_required=n)
    del arr
    experiment = uniq.to_experiment()
    experiment.to_npz(output_prefix + '_sp_counts.npz')


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


def store_results(output_prefix, arr):
    if not output_prefix.endswith('.h5'):
        output_prefix += '.h5'
    seqc.log.info('Storing processed data in %s' % output_prefix)
    arr.save_h5(output_prefix)


def run_complete():
    seqc.log.info('Run complete.')


def clean_up(temp_dir):
    seqc.log.info('Run succeeded. Cleaning up and terminating.')
    if os.path.isdir(temp_dir):
        shutil.rmtree(temp_dir)
    if os.path.isdir('_STARtmp'):
        shutil.rmtree('_STARtmp')


def in_drop(output_prefix, forward, reverse, samfile, merged_fastq, subparser_name, index,
            n_threads, frag_len, star_args, barcodes, **kwargs):

    output_prefix, output_dir = fix_output_paths(output_prefix)

    index, gtf = check_index(index, output_dir)

    cb = check_and_load_barcodes(barcodes, output_dir)

    # htqc(temp_dir, forward, reverse, 1000, 1)

    merged_fastq = merge(forward, reverse, samfile, merged_fastq, subparser_name,
                         output_dir, cb, n_threads)

    samfile = align(merged_fastq, samfile, star_args, output_dir, n_threads, index,
                    output_prefix)

    arr = process_samfile(samfile, output_prefix, n_threads, gtf, frag_len)

    resolve_alignments(index, arr, n=0, output_prefix=output_prefix)

    store_results(output_prefix, arr)

    save_counts_matrices(output_prefix, arr, n=3)

    run_complete()


def drop_seq(output_prefix, forward, reverse, samfile, merged_fastq, subparser_name, index,
             n_threads, frag_len, star_args, barcodes, **kwargs):

    output_prefix, output_dir = fix_output_paths(output_prefix)

    index, gtf = check_index(index, output_dir)

    cb = check_and_load_barcodes(barcodes, output_dir)

    # htqc(temp_dir, forward, reverse, 1000, 1)

    merged_fastq = merge(forward, reverse, samfile, merged_fastq, subparser_name,
                         output_dir, cb, n_threads)

    samfile = align(merged_fastq, samfile, star_args, output_dir, n_threads, index,
                    output_prefix)

    arr = process_samfile(samfile, output_prefix, n_threads, gtf, frag_len)

    resolve_alignments(index, arr, n=0, output_prefix=output_prefix)

    arr.save_h5(output_prefix + '.h5')

    save_counts_matrices(output_prefix, arr, n=0)

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
