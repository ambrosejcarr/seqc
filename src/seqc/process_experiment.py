#!/usr/local/bin/python3

import argparse
import configparser
import multiprocessing
import os
import shutil
import pickle
import sys
from subprocess import Popen
from copy import copy

import numpy as np
import pandas as pd

import seqc
from seqc.exceptions import ConfigurationError, ArgumentParserError
from seqc.io import check_s3links, obtain_size


class NewArgumentParser(argparse.ArgumentParser):
    """
    Simple wrapper for ArgumentParser that allows flags to be caught before the presence
    of 'required' arguments are tested. Allows us to add flags, for example, for
    checking the progress of existing SEQC runs.
    """
    def error(self, message):
        # checks to see whether user wants to check remote experiment status
        if '--check-progress' in sys.argv[1:]:
            seqc.remote.check_progress()
        else:
            print(message)
        sys.exit(0)


def parse_args(args):
    """
    command line argument parser for process_experiment

    :param args: list of command-line arguments (either passed as a list, or retrieved
      from sys.argv.
    :returns: args, namespace object, output of ArgumentParser.parse_args()
    """
    p = NewArgumentParser(description='Process Single-Cell RNA Sequencing Data')
    p.add_argument('platform',
                   choices=['in_drop', 'drop_seq', 'mars1_seq',
                            'mars2_seq', 'in_drop_v2', 'in_drop_v3'],
                   help='which platform are you merging annotations from?')

    a = p.add_argument_group('required arguments')
    a.add_argument('-o', '--output-stem', metavar='O', required=True,
                   help='If remote=True (default), an s3 link to the directory where all '
                        'output data should be uploaded when the run completes (e.g. '
                        's3://my-bucket/seqc_output_run_x/). If remote=False, the '
                        'complete path and prefix for the output data '
                        '(e.g. /Users/me/data/seqc/run_x_output). Cannot have a terminal '
                        '/.')
    a.add_argument('-i', '--index', metavar='I', required=True,
                   help='Local folder or s3 link to a directory containing the STAR '
                        'index used for alignment.')
    a.add_argument('--barcode-files', nargs='*', metavar='BF', default=[],
                   help='Either (a) an s3 link to a folder containing only barcode '
                        'files, or (b) the full file path of each file on the local '
                        'machine.')
    a.add_argument('--email-status', metavar='E', default=None,
                   help='Email address to receive run summary or errors when running '
                        'remotely. Optional if running locally.')

    i = p.add_argument_group('input arguments')
    i.add_argument('-g', '--genomic-fastq', nargs='*', metavar='G', default=[],
                   help='List of fastq file(s) containing genomic information, or an s3 '
                        'link to a directory containing only genomic fastq file(s).')
    i.add_argument('-b', '--barcode-fastq', nargs='*', metavar='B', default=[],
                   help='List of fastq file(s) containing barcode information, or an s3 '
                        'link to a directory containing only barcode fastq file(s).')
    i.add_argument('-m', '--merged-fastq', nargs='?', metavar='M', default='',
                   help='Filename or s3 link to a fastq file containing genomic '
                        'information annotated with barcode data.')
    i.add_argument('-s', '--samfile', nargs='?', metavar='S', default='',
                   help='Filename or s3 link to a .sam file containing aligned, merged '
                        'sequence records.')
    i.add_argument('-r', '--read-array', nargs='?', metavar='RA', default='',
                   help='Filename or s3 link to a ReadArray (.h5) archive containing '
                        'processed sam records.')
    i.add_argument('--basespace', metavar='BS',
                   help='BaseSpace sample ID. The string of numbers indicating the id '
                        'of the BaseSpace sample. (e.g. if the link to the sample is '
                        'https://basespace.illumina.com/sample/34000253/0309, '
                        'then --basespace would be 34000253.')

    f = p.add_argument_group('filter arguments')
    f.add_argument('--max-insert-size', metavar='F',
                   help='maximum paired-end insert size that is considered a valid '
                        'record. For multialignment correction. Not currently used.',
                   default=1000)
    f.add_argument('--min-poly-t', metavar='T',
                   help='minimum size of poly-T tail that is required for a barcode to '
                        'be considered a valid record (default=None, automatically '
                        'estimates the parameter from the sequence length)',
                   default=None, type=int)
    f.add_argument('--max-dust-score', metavar='D',
                   help='maximum complexity score for a read to be considered valid. '
                        '(default=10, higher scores indicate lower complexity.)')
    f.add_argument('--max-ed', metavar='ED',
                   help='Maximum hamming distance for correcting barcode errors. Should '
                        'be set to 1 when running in_drop legacy barcodes. (Default=2).',
                   default=2, type=int)
    f.add_argument('--singleton-weight', metavar='SW',
                   help='Weight to apply to singletons in the count matrix. Float '
                        'between 0 and 1, default=1 (all molecules get full weight)',
                   default=1.0, type=float)

    s = p.add_argument_group('alignment arguments')
    s.add_argument('--star-args', default=None, nargs='*',
                   help='additional arguments that should be passed to the STAR '
                        'aligner. For example, to set the maximum allowable times for a '
                        'read to align to 20, one would set '
                        '--star-args outFilterMultimapNmax=20. Additional arguments can '
                        'be provided as a white-space separated list.')

    r = p.add_argument_group('remote run arguments')
    r.set_defaults(remote=True)
    r.set_defaults(check=False)
    r.add_argument('--local', dest="remote", action="store_false",
                   help='Run SEQC locally instead of initiating on AWS EC2 servers.')
    r.add_argument('--aws', default=False, action='store_true',
                   help='Automatic flag; no need for user specification.')
    r.add_argument('--no-terminate', default='False',
                   help='Do not terminate the EC2 instance after program completes. If '
                        '"on-success" is specified, the EC2 instance does not terminate '
                        'in case the SEQC run throws an error.')
    r.add_argument('--check-progress', dest="check", action="store_true",
                   help='Check progress of all currently running remote SEQC runs.')
    r.add_argument('--instance-type', default='c4',
                   help='AWS instance (c3, c4, r3) to run SEQC remotely. Default=c4.')
    r.add_argument('--spot-bid', type=float, default=None,
                   help='float, Amount to bid for a spot instance. Default=None (will '
                        'reserve a non-spot instance). WARNING: using spot instances '
                        'will cause your instance to terminate if instance prices exceed '
                        'your spot bid during runtime.')
    r.add_argument('--log-name', type=str, default='seqc.log',
                   help='Output log name (default=seqc.log)')

    p.add_argument('-v', '--version', action='version',
                   version='{} {}'.format(p.prog, seqc.__version__))

    try:
        return p.parse_args(args)
    except ArgumentParserError:
        raise


def recreate_command_line_arguments(args):
    """
    Helper function to recreate command line arguments for passing to the remote run that
    is often called by the local process. This function is aware that 'platform' is
    a positional argument, and '--remote' is a flag; otherwise arguments are
    reconstructed from the dictionary.

    This function may require adjustment when new flags or positional arguments are added
    to process_experiment.py
    """
    args = copy(args)  # don't break original object
    platform = args.platform
    del args.platform  # we explicitly create a positional parameter
    del args.remote  # this is set by --local (or default)
    cmd = 'process_experiment.py {} '.format(platform)
    for key, value in vars(args).items():
        if value:
            if isinstance(value, list):
                value = ' '.join(value)
            key = key.replace('_', '-')
            cmd += '--{} {} '.format(key, value)
    cmd += '--local --aws'
    return cmd


def run_remote(args, volsize):
    """
    Mirror the local arguments to an AWS server and execute the run there. When complete,
    terminates the local process.

    :param args: simple namespace object; output of parse_args()
    :param volsize: estimated volume needed for run
    """
    seqc.log.notify('Beginning remote SEQC run...')

    # recreate remote command, but instruct it to run locally on the server.
    cmd = recreate_command_line_arguments(args)

    # set up remote cluster
    cluster = seqc.remote.ClusterServer()
    volsize = int(np.ceil(volsize/1e9))

    try:  # if anything goes wrong during cluster setup, clean up the instance
        cluster.cluster_setup(volsize, args.instance_type, spot_bid=args.spot_bid)
        cluster.serv.connect()

        seqc.log.notify('Beginning remote run.')
        if args.output_stem.endswith('/'):
            args.output_stem = args.output_stem[:-1]
        # writing name of instance in local machine to keep track of instance
        with open(os.path.expanduser('~/.seqc/instance.txt'), 'a') as f:
            _, run_name = os.path.split(args.output_stem)
            f.write('%s:%s\n' % (cluster.inst_id.instance_id, run_name))

        # writing name of instance in /data/instance.txt for clean up
        inst_path = '/data/instance.txt'
        cluster.serv.exec_command(
            'echo {instance_id} > {inst_path}'.format(
                inst_path=inst_path, instance_id=str(cluster.inst_id.instance_id)))
        cluster.serv.exec_command('sudo chown -R ubuntu /home/ubuntu/.seqc/')
        cluster.serv.put_file(os.path.expanduser('~/.seqc/config'),
                              '/home/ubuntu/.seqc/config')
        cluster.serv.exec_command('cd /data; nohup {cmd} > /dev/null 2>&1 &'
                                  ''.format(cmd=cmd))

        # check that process_experiment.py is actually running on the cluster
        out, err = cluster.serv.exec_command('ps aux | grep process_experiment.py')
        res = ' '.join(out)
        if '/usr/local/bin/process_experiment.py' not in res:
            raise ConfigurationError('Error executing SEQC on the cluster!')
        seqc.log.notify('Terminating local client. Email will be sent when remote run '
                        'completes. Please use "process_experiment.py --check-progress" '
                        'to monitor the status of the remote run.')
    except Exception as e:
        seqc.log.notify('Error {e} occurred during cluster setup!'.format(e=e))
        if cluster.cluster_is_running():
            seqc.log.notify('Cleaning up instance {id} before exiting...'.format(
                id=cluster.inst_id.instance_id))
            seqc.remote.terminate_cluster(cluster.inst_id.instance_id)
        raise  # re-raise original exception


def check_executables() -> (bool, bool):
    """
    checks whether pigz and mutt are installed on the machine of the
    current seqc run. returns True/False for both

    :returns pigz, email: (bool, bool) Booleans indicating if pigz and mutt are installed.
    """
    pigz = False
    email = False
    if shutil.which('pigz'):
        pigz = True
    if shutil.which('mutt'):
        email = True
    return pigz, email


def check_arguments(args, basespace_token: str) -> float:
    """
    checks data input through the command line arguments and throws
    an error if the provided input data is invalid.

    additionally, this function obtains a rough estimate of how much
    volume storage is needed for the overall SEQC run.

    :param args: Namespace object, output from ArgumentParser.parse_args()
    :param basespace_token: str, OAuth token for BaseSpace.
    :returns total: float, estimated Kb of Volume space needed to run SEQC remotely.
    """

    # make sure only one filetype has been passed
    multi_input_error_message = ('Only one input type (-s, -m, -b/-g, or --basespace) '
                                 'should be passed to SEQC.')
    unpaired_fastq_error_message = ('If either genomic or barcode fastq files are '
                                    'provided, both must be provided.')

    # simplify variables and parse out the ones needed
    barcodes = args.barcode_files
    index_dir = args.index
    barcode_fastq = args.barcode_fastq
    genomic_fastq = args.genomic_fastq
    merged = args.merged_fastq
    samfile = args.samfile
    read_array = args.read_array
    basespace = args.basespace

    if args.spot_bid is not None:
        if args.spot_bid < 0:
            raise ValueError('"{bid}" must be a non-negative float! Exiting.'.format(
                bid=args.spot_bid))

    # check to make sure that --email-status is passed with remote run
    if args.remote and not args.email_status:
        raise ValueError('Please supply the --email-status flag for a remote SEQC run.')
    if args.instance_type not in ['c3', 'c4', 'r3']:
        raise ValueError('All AWS instance types must be either c3, c4, or r3.')
    if args.no_terminate not in ['True', 'true', 'False', 'false', 'on-success']:
        raise ValueError('the --no-terminate flag must be either True, False, '
                         'or on-success.')

    # make sure at least one input has been passed
    if not any([barcode_fastq, genomic_fastq, merged, samfile, basespace, read_array]):
        raise ValueError('At least one input argument must be passed to SEQC.')
    if not barcodes:
        if args.platform != 'drop_seq':
            raise ValueError('Barcode files are required for this SEQC run.')

    # keep track of which files need to be checked
    seqc_input = barcodes + [index_dir]

    # keep track of how much space is needed given input
    # using worst-case estimates to make sure we don't run out of space
    cushion = 5e10
    if args.instance_type == 'r3':
        cushion = 9e10
    total = 0

    if barcode_fastq or genomic_fastq:
        if not all((barcode_fastq, genomic_fastq)):
            raise ValueError(unpaired_fastq_error_message)
        if any((merged, samfile, basespace, read_array)):
            raise ValueError(multi_input_error_message)
        seqc_input = seqc_input + barcode_fastq + genomic_fastq
        check_s3links(seqc_input)

        # checking size of input file
        input_fastq = barcode_fastq + genomic_fastq
        for item in input_fastq:
            total += obtain_size(item)
        total += (total * 14) + cushion
    if samfile:
        if any((merged, barcode_fastq, genomic_fastq, basespace, read_array)):
            raise ValueError(multi_input_error_message)
        seqc_input += [samfile]
        check_s3links(seqc_input)

        # checking size of input file
        total += obtain_size(samfile)
        total += (total * 2) + 2e10
    if merged:
        if any((samfile, barcode_fastq, genomic_fastq, basespace, read_array)):
            raise ValueError(multi_input_error_message)
        seqc_input += [merged]
        check_s3links(seqc_input)

        # checking size of input file
        total += obtain_size(merged)
        total += (total * 13) + cushion
    if basespace:
        if any((samfile, merged, barcode_fastq, genomic_fastq, read_array)):
            raise ValueError(multi_input_error_message)
        if not basespace_token or basespace_token == 'None':
            raise ValueError(
                'If the --basespace argument is used, the BaseSpace token must be '
                'specified in the seqc config file.')
        seqc_input += [basespace]
        check_s3links(seqc_input)
    if read_array:
        if any((samfile, merged, barcode_fastq, genomic_fastq, basespace)):
            raise ValueError(multi_input_error_message)
        seqc_input += [read_array]
        seqc_input.remove(index_dir)
        check_s3links(seqc_input)

        # checking size of input
        for item in seqc_input:
            total += obtain_size(item)
        total += 1e10
    if basespace:
        seqc.io.BaseSpace.check_sample(basespace, basespace_token)
        # checking size of input file
        total = seqc.io.BaseSpace.check_size(basespace, basespace_token)
        total += (total * 14) + cushion

    # return total size needed for EBS volume
    return total


def read_config(config_file: str) -> dict:
    """
    :param config_file: str, location of seqc configuration file
    :returns config: dict, dictionary of arguments and values in the configuration file
    """
    # todo add additional checks to make sure correct parameters are filled in!
    config_file = os.path.expanduser('~/.seqc/config')
    config = configparser.ConfigParser()
    if not config.read(config_file):
        raise ConfigurationError('Please run ./configure (found in the seqc '
                                 'directory) before attempting to run '
                                 'process_experiment.py.')
    return config


def get_basespace_data(args, output_dir: str, basespace_token: str) -> (str, str):
    """
    If --basespace argument is passed, download BaseSpace data

    :param args: Namespace object, output of ArgumentParser.parse_args()
    :param output_dir: str, output directory
    :param basespace_token: str, OAuth token for BaseSpace authentication
    :return barcode_fastq, genomic_fastq:
    """
    seqc.log.info('BaseSpace link provided for fastq argument. Downloading '
                  'input data.')

    # making extra directories for BaseSpace download, changing permissions
    bspace_dir = output_dir + '/Data/Intensities/BaseCalls/'
    bf = Popen(['sudo', 'mkdir', '-p', bspace_dir])
    bf.communicate()
    if args.aws:  # changing permissions is unnecessary if local run
        bf2 = Popen(['sudo', 'chown', '-c', 'ubuntu', bspace_dir])
        bf2.communicate()
    barcode_fastq, genomic_fastq = seqc.io.BaseSpace.download(
        args.platform, args.basespace, output_dir, basespace_token)
    return barcode_fastq, genomic_fastq


def update_directories_for_aws(output_stem: str, output_prefix: str) -> (
        str, str, str, str):
    """
    if the --aws argument is passed, this function returns updated directories and
     prefixes for output files

    :param output_stem: str, stem for output files
    :param output_prefix: str, prefix for output files
    :returns output_stem, output_dir, output_prefix: (str, str, str), stem, directory and
     prefix for output files generated by an aws SEQC run.
    """
    aws_upload_key = output_stem
    if not aws_upload_key.endswith('/'):
        aws_upload_key += '/'
    output_stem = '/data/' + output_prefix
    output_dir, output_prefix = os.path.split(output_stem)
    return aws_upload_key, output_stem, output_dir, output_prefix


def get_s3_fastq(fastq_file: list, output_dir: str) -> list:
    """
    Checks if -g/--genomic-fastq was passed. If it was, downloads the necessary files if
    any passed arguments were s3 links.

    :param fastq_file: list, a list of fastq files or s3 links to fastq files
    :param output_dir: directory in which to download fastq files
    :returns fastq_file: list, filename(s) of local genomic_fastq files.
    """
    if fastq_file:
        if not fastq_file[0].startswith('s3://'):
            for gf in fastq_file:
                if not os.path.isfile(gf):
                    raise ValueError('Provided genomic fastq files: "[%s]" is '
                                     'neither an s3 link or a valid filepath' %
                                     ', '.join(map(str, fastq_file)))
        else:
            try:
                seqc.log.info('Downloading genomic fastq files from Amazon s3 link.')
                if fastq_file[0].endswith('/'):
                    # s3 directory specified, download all files recursively
                    bucket, prefix = seqc.io.S3.split_link(fastq_file[0])
                    cut_dirs = prefix.count('/')
                    fastq_file = seqc.io.S3.download_files(
                        bucket=bucket, key_prefix=prefix, output_prefix=output_dir,
                        cut_dirs=cut_dirs)
                else:
                    # individual s3 links provided, download each fastq file
                    downloaded_files = []
                    for s3link in fastq_file:
                        bucket, prefix = seqc.io.S3.split_link(s3link)
                        _, fname = os.path.split(prefix)
                        fname = output_dir + '/' + fname
                        seqc.io.S3.download_file(bucket, prefix, fname)
                        downloaded_files.append(fname)
                    fastq_file = sorted(downloaded_files)
                # note that this is printed regardless of whether a file is downloaded
                seqc.log.info('Genomic fastq files [%s] successfully installed.' %
                              ', '.join(map(str, fastq_file)))
            except FileExistsError:
                pass  # file is already present.
    return fastq_file


def get_index(output_dir: str, index: str, read_array: str, samfile: str) -> str:
    """
    Checks the index parameter, downloading files from s3 if an s3 link was passed

    :param output_dir: str, directory that the index will be downloaded to if necessary
    :param index: value of -i/--index parameter
    :param read_array: str, value of -r/--read-array parameter
    :param samfile: str, value of -s/--samfile parameter
    :return index: str, path to downloaded or local index
    """
    if not index.startswith('s3://'):
        if not os.path.isdir(index):
            raise ValueError('Provided index: "%s" is neither an s3 link or a valid '
                             'filepath' % index)
    else:
        if not read_array:  # does not require index
            seqc.log.info('AWS s3 link provided for index. Downloading index.')
            bucket, prefix = seqc.io.S3.split_link(index)
            index = output_dir + '/index/'  # set index  based on s3 download
            cut_dirs = prefix.count('/')
            # install whole index
            if not samfile:
                try:
                    seqc.io.S3.download_files(bucket=bucket, key_prefix=prefix,
                                              output_prefix=index,
                                              cut_dirs=cut_dirs)
                except FileExistsError:
                    pass  # file is already present
            else:  # samfile provided, only download annotations file
                try:
                    annotations_file = index + 'annotations.gtf'
                    seqc.io.S3.download_file(bucket, prefix + 'annotations.gtf',
                                             annotations_file)
                except FileExistsError:
                    pass  # file is already present
    return index


def determine_start_point(args) -> (bool, bool, bool):
    """
    determine where seqc should start based on which parameters were passed.

    :param args: Namespace object, result of ArgumentParser.parse_args()
    :returns merge, align, process_samfile: (bool, bool, bool) indicates whether merging,
      alignment, or processing samfiles should be carried out.
    """
    merge = True
    align = True
    process_samfile = True
    if args.read_array:
        merge = False
        align = False
        process_samfile = False
    if args.samfile:
        merge = False
        align = False
    if args.merged_fastq:
        merge = False
    return merge, align, process_samfile


def merge_fastq_files(
        platform: str, barcode_fastq: str, merged_fastq: str, output_stem: str,
        genomic_fastq: str, pigz: bool) -> (str, int):
    """
    annotates genomic fastq with barcode information; merging the two files.

    :param platform:
    :param barcode_fastq:
    :param merged_fastq:
    :param output_stem:
    :param genomic_fastq:
    :param pigz:
    :returns merged_fastq, fastq_records: (str, int) name of merged fastq file and the
      number of fastq records that were processed.
    """

    seqc.log.info('Merging genomic reads and barcode annotations.')
    merge_function = getattr(seqc.sequence.merge_functions, platform)
    merged_fastq, fastq_records = seqc.sequence.fastq.merge_paired(
        merge_function=merge_function,
        fout=output_stem + '_merged.fastq',
        genomic=genomic_fastq,
        barcode=barcode_fastq)

    # delete genomic/barcode fastq files after merged.fastq creation
    seqc.log.info('Removing original fastq file for memory management.')
    delete_fastq = ' '.join(['rm'] + genomic_fastq + barcode_fastq)
    seqc.io.ProcessManager(delete_fastq).run_all()

    # zip merged fastq file
    seqc.log.info('Gzipping merged fastq file.')
    if pigz:
        pigz_zip = "pigz --best -k {fname}".format(fname=merged_fastq)
    else:
        pigz_zip = "gzip -kf {fname}".format(fname=merged_fastq)
    pigz_proc = seqc.io.ProcessManager(pigz_zip)
    pigz_proc.run_all()
    pigz_proc.wait_until_complete()  # prevents slowing down STAR alignment
    return merged_fastq, fastq_records


def align_fastq_records(
        merged_fastq, output_dir, output_stem, star_args, index, n_processes,
        aws_upload_key, pigz) -> (str, str, seqc.io.ProcessManager):
    """
    Align fastq records.

    :param merged_fastq: str, path to merged .fastq file
    :param output_dir: str, directory for output files
    :param output_stem: str, stem for output files
    :param star_args: dict, extra keyword arguments for STAR
    :param index: str, file path to directory containing STAR index
    :param n_processes: int, number of STAR processes to initiate
    :param aws_upload_key: str, location to upload files
    :param pigz: bool, True if pigz is installed, else False
    :return samfile, input_data, manage_merged: (str, str, seqc.io.ProcessManager)
      name of .sam file containing aligned reads, indicator of which data was used as
      input, and a ProcessManager for merged fastq files
    """
    input_data = 'start'
    seqc.log.info('Aligning merged fastq records.')
    if merged_fastq.startswith('s3://'):
        input_data = 'merged'
        bucket, prefix = seqc.io.S3.split_link(merged_fastq)
        _, fname = os.path.split(prefix)
        merged_fastq = seqc.io.S3.download_file(bucket, prefix, output_dir + '/' + fname)
        if merged_fastq.endswith('.gz'):
            if pigz:
                cmd = 'pigz -d {fname}'.format(fname=merged_fastq)
            else:
                cmd = 'gunzip {fname}'.format(fname=merged_fastq)
            gunzip_proc = seqc.io.ProcessManager(cmd)
            gunzip_proc.run_all()
            gunzip_proc.wait_until_complete()  # finish before STAR alignment
            merged_fastq = merged_fastq.strip('.gz')
        seqc.log.info('Merged fastq file %s successfully installed from S3.' %
                      merged_fastq)
    *base_directory, stem = output_stem.split('/')
    alignment_directory = '/'.join(base_directory) + '/alignments/'
    os.makedirs(alignment_directory, exist_ok=True)
    if star_args is not None:
        star_kwargs = dict(a.strip().split('=') for a in star_args)
    else:
        star_kwargs = {}
    samfile = seqc.alignment.star.align(
        merged_fastq, index, n_processes, alignment_directory,
        **star_kwargs)

    # Removing or uploading the merged fastq file depending on start point
    if input_data == 'merged':
        seqc.log.info('Removing merged.fastq file for memory management.')
        rm_cmd = 'rm {merged_file}'.format(merged_file=merged_fastq)
        seqc.io.ProcessManager(rm_cmd).run_all()
        manage_merged = None
    else:
        if aws_upload_key:
            seqc.log.info('Uploading gzipped merged fastq file to S3.')
            merge_upload = 'aws s3 mv {fname} {s3link}'.format(
                fname=merged_fastq + '.gz', s3link=aws_upload_key)
            manage_merged = seqc.io.ProcessManager(merge_upload)
            manage_merged.run_all()
        else:
            manage_merged = None
    return samfile, input_data, manage_merged


def create_or_download_read_array(
        process_samfile, samfile, output_dir, index, read_array, output_stem,
        aws_upload_key, input_data) -> (str, str, seqc.io.ProcessManager, int):
    """
    Create or download a ReadArray object.

    :param process_samfile: bool, True if the samfile should be processed
    :param samfile: str, filename of .sam file
    :param output_dir: str, directory for output files
    :param index: str, directory containing index files
    :param read_array: str, filename, s3 link to ReadArray, or None (if generating from
      .sam file)
    :param output_stem: str, stem for output files
    :param aws_upload_key: str, key where aws files should be uploaded
    :param input_data: str, indicator variable that contains the identity of the input
      data
    :returns ra, input_data, manage_samfile: (str, str, seqc.io.ProcessManager, int)
      ReadArray object, indicator object, samfile ProcessManager, and number of processed
      sam records
    """
    manage_samfile = None
    sam_records = None
    if process_samfile:
        seqc.log.info('Filtering aligned records and constructing record database.')
        if samfile.startswith('s3://'):
            input_data = 'samfile'
            bucket, prefix = seqc.io.S3.split_link(samfile)
            _, fname = os.path.split(prefix)
            samfile = seqc.io.S3.download_file(bucket, prefix, output_dir + '/' + fname)
            seqc.log.info('Samfile %s successfully installed from S3.' % samfile)
        ra, sam_records = seqc.core.ReadArray.from_samfile(
            samfile, index + 'annotations.gtf')
        read_array = output_stem + '.h5'
        ra.save(read_array)

        # converting sam to bam and uploading to S3, else removing samfile
        if input_data == 'samfile':
            seqc.log.info('Removing samfile for memory management.')
            rm_samfile = 'rm {fname}'.format(fname=samfile)
            seqc.io.ProcessManager(rm_samfile).run_all()
        else:
            if aws_upload_key:
                seqc.log.info('Converting samfile to bamfile and uploading to S3.')
                bamfile = output_dir + '/alignments/Aligned.out.bam'
                convert_sam = 'samtools view -bS -o {bamfile} {samfile}' \
                    .format(bamfile=bamfile, samfile=samfile)
                upload_bam = 'aws s3 mv {fname} {s3link}'.format(
                    fname=bamfile, s3link=aws_upload_key)
                manage_samfile = seqc.io.ProcessManager(convert_sam, upload_bam)
                manage_samfile.run_all()
    else:
        if read_array.startswith('s3://'):
            input_data = 'readarray'
            bucket, prefix = seqc.io.S3.split_link(read_array)
            _, fname = os.path.split(prefix)
            read_array = seqc.io.S3.download_file(
                bucket, prefix, output_dir + '/' + fname)
            seqc.log.info('Read array %s successfully installed from S3.' %
                          read_array)
        ra = seqc.core.ReadArray.load(read_array)
    return ra, input_data, manage_samfile, sam_records, read_array


def download_barcodes(platform: str, barcode_files: list, output_dir: str) -> list:
    """
    if an s3 link was passed and the platform was not drop-seq, download barcode files

    :param platform: str, type of experiment
    :param barcode_files: list, files containing barcode information or s3 links
    :param output_dir: str, directory for output files
    :returns barcode_files: list, files containing barcodes
    """
    if platform != 'drop_seq':
        if not barcode_files[0].startswith('s3://'):
            for cb in barcode_files:
                if not os.path.isfile(cb):
                    raise ValueError('provided barcode files: "[%s]" is neither '
                                     'an s3 link or a valid filepath' %
                                     ', '.join(map(str, barcode_files)))
        else:
            try:
                seqc.log.info('AWS s3 link provided for barcodes. Downloading files.')
                if not barcode_files[0].endswith('/'):
                    barcode_files[0] += '/'
                bucket, prefix = seqc.io.S3.split_link(barcode_files[0])
                cut_dirs = prefix.count('/')
                barcode_files = seqc.io.S3.download_files(
                    bucket=bucket, key_prefix=prefix, output_prefix=output_dir,
                    cut_dirs=cut_dirs)
                seqc.log.info('Barcode files [%s] successfully installed.' %
                              ', '.join(map(str, barcode_files)))
            except FileExistsError:
                pass  # file is already present.
    return barcode_files


def upload_data_and_notify_user(
        email_status: str, aws_upload_key: str, align: bool, input_data: str,
        manage_merged: seqc.io.ProcessManager, process_samfile: bool,
        merged_fastq: str, samfile: str, read_array: str,
        sparse_proc: seqc.io.ProcessManager, sparse_csv: str,
        manage_ra: seqc.io.ProcessManager, summary: dict, fastq_records: int,
        sam_records: int, total_molecules: int, mols_lost: int, cells_lost: int,
        manage_samfile: seqc.io.ProcessManager, cell_description: pd.Series,
        output_stem: str, log_name: str, email: bool) -> None:
    """
    Uploads data and notifies the user of the termination of a run

    :param email_status: str, email address to email results
    :param aws_upload_key: str, location to upload files on aws
    :param align: bool, whether or not alignment occurred
    :param input_data: str, indicator, stores the type of input data
    :param manage_merged: seqc.io.ProcessManager, an upload manager for merged_fastq
    :param process_samfile: bool, whether or not samfile was processed
    :param merged_fastq: str, name of merged fastq file
    :param samfile: str, name of sam file
    :param read_array: str, name of stored ReadArray
    :param sparse_proc: seqc.io.ProcessManager, an upload manager for sparse counts file
    :param sparse_csv: str, name of sparse_csv file
    :param manage_ra: seqc.io.ProcessManager, an upload manager for the ReadArray
    :param summary: dict, object containing summary statistics
    :param fastq_records: int, number of processed fastq records
    :param sam_records: int, number of processed sam records
    :param total_molecules: int, number of molecules
    :param mols_lost: int, number of molecules lost
    :param cells_lost: int, number of cells lost
    :param manage_samfile: seqc.io.ProcessManager, an upload manager for the sam file
    :param cell_description: pd.Series, a summary of cells and cell counts
    :param output_stem: str, stem for file output
    :param log_name: str, name of SEQC.log file
    :param email: bool, whether email should be sent
    :return None:
    """

    seqc.log.info('Starting file upload onto %s.' % aws_upload_key)

    if email_status and aws_upload_key is not None:
        # make sure that all other files are uploaded before termination
        if align and input_data != 'merged':
            manage_merged.wait_until_complete()
            seqc.log.info('Successfully uploaded %s to the specified S3 '
                          'location "%s"' % (merged_fastq, aws_upload_key))
        if process_samfile:
            if input_data != 'samfile':
                manage_samfile.wait_until_complete()
                seqc.log.info('Successfully uploaded %s to the specified S3 '
                              'location "%s"' % (samfile, aws_upload_key))
            manage_ra.wait_until_complete()
            seqc.log.info('Successfully uploaded %s to the specified S3 '
                          'location "%s"' % (read_array, aws_upload_key))
            sparse_proc.wait_until_complete()
            seqc.log.info('Successfully uploaded %s to the specified S3 '
                          'location "%s"' % (sparse_csv + '.gz', aws_upload_key))

        # upload count matrix and alignment summary at the very end
        if summary:
            if fastq_records:
                summary['n_fastq'] = fastq_records
            else:
                summary['n_fastq'] = 'NA'
            if sam_records:
                summary['n_sam'] = sam_records
            else:
                summary['n_sam'] = 'NA'
            summary['total_mc'] = total_molecules
            summary['mols_lost'] = mols_lost
            summary['cells_lost'] = cells_lost
            summary['cell_desc'] = cell_description

        # todo: run summary will not be reported if n_fastq or n_sam = NA
        seqc.remote.upload_results(
            output_stem, email_status, aws_upload_key, input_data,
            summary, log_name, email)


def main(args: list=None) -> None:
    """
    Execute a SEQC run.

    :param args: list of arguments. Normally either extracted from sys.argv or passed as
      a list from a testing function.
    :return: None
    """

    args = parse_args(args)
    if args.aws:
        args.log_name = args.log_name.split('/')[-1]
        seqc.log.setup_logger('/data/' + args.log_name)
    else:
        seqc.log.setup_logger(args.log_name)
    try:
        err_status = False
        pigz, email = check_executables()
        if not email and not args.remote:
            seqc.log.notify('mutt was not found on this machine; an email will not '
                            'be sent to the user upon termination of SEQC run.')
        seqc.log.args(args)

        # read in config file, make sure it exists
        config = read_config('~/.seqc/config')

        # extract basespace token and make sure args.index is a directory
        basespace_token = config['BaseSpaceToken']['base_space_token']
        if not args.index.endswith('/'):
            args.index += '/'

        # check arguments if aws flag is not set (arguments would already be checked)
        if not args.aws:
            total_size = check_arguments(args, basespace_token)
        else:
            total_size = None

        # check whether output is an s3 or local directory, split output_stem
        if args.output_stem.endswith('/'):
            if args.output_stem.startswith('s3://'):
                output_dir, output_prefix = os.path.split(args.output_stem[:-1])
            else:
                raise ValueError('-o/--output-stem should not be a directory for local '
                                 'SEQC runs.')
        else:
            output_dir, output_prefix = os.path.split(args.output_stem)

        # check output directory if remote run
        if args.remote:
            if not args.output_stem.startswith('s3://'):
                raise ValueError('-o/--output-stem must be an s3 link for remote SEQC '
                                 'runs.')
            seqc.remote.cluster_cleanup()
            run_remote(args, total_size)
            sys.exit(0)

        if args.aws:
            aws_upload_key, args.output_stem, output_dir, output_prefix = \
                update_directories_for_aws(args.output_stem, output_prefix)
        else:
            aws_upload_key = None

        # download data if necessary
        if args.basespace:
            args.barcode_fastq, args.genomic_fastq = get_basespace_data(
                args, output_dir, basespace_token)

        # check for remote fastq file links
        args.genomic_fastq = get_s3_fastq(args.genomic_fastq, output_dir)
        args.barcode_fastq = get_s3_fastq(args.barcode_fastq, output_dir)

        # check if the index must be downloaded
        args.index = get_index(output_dir, args.index, args.read_array, args.samfile)

        n_processes = multiprocessing.cpu_count() - 1  # get number of processors

        input_data = 'start'  # which data was used as input?

        # determine where the script should start:
        merge, align, process_samfile = determine_start_point(args)

        if merge:
            if args.min_poly_t is None:  # estimate min_poly_t if it was not provided
                args.min_poly_t = seqc.filter.estimate_min_poly_t(
                    args.barcode_fastq, args.platform)
                seqc.log.notify('Estimated min_poly_t={!s}'.format(args.min_poly_t))

            args.merged_fastq, fastq_records = merge_fastq_files(
                args.platform, args.barcode_fastq, args.merged_fastq, args.output_stem,
                args.genomic_fastq, pigz)
        else:
            fastq_records = None

        if align:
            args.samfile, input_data, manage_merged = align_fastq_records(
                args.merged_fastq, output_dir, args.output_stem, args.star_args,
                args.index, n_processes, aws_upload_key, pigz)

        ra, input_data, manage_samfile, sam_records, h5_file = \
            create_or_download_read_array(
            process_samfile, args.samfile, output_dir, args.index, args.read_array,
            args.output_stem, aws_upload_key, input_data)
        args.read_array = h5_file

        args.barcode_files = download_barcodes(
            args.platform, args.barcode_files, output_dir)

        # SEQC was started from input other than fastq files
        if args.min_poly_t is None:
            args.min_poly_t = 0
            seqc.log.notify('Warning: SEQC started from step other than unmerged fastq '
                            'with empty --min-poly-t parameter. Continuing with '
                            '--min-poly-t=0.')

        # correct errors
        seqc.log.info('Correcting cell barcode and RMT errors')
        correct_errors = getattr(seqc.correct_errors, args.platform)
        # for in-drop and mars-seq, summary is a dict. for drop-seq, it may be None
        cell_counts, _, summary = correct_errors(
            ra, args.barcode_files, reverse_complement=False,
            required_poly_t=args.min_poly_t, max_ed=args.max_ed,
            singleton_weight=args.singleton_weight)

        # uploading read array to S3 if created, else removing read array
        if input_data == 'readarray':
            seqc.log.info('Removing .h5 file for memory management.')
            rm_ra = 'rm {fname}'.format(fname=args.read_array)
            seqc.io.ProcessManager(rm_ra).run_all()
        else:
            if aws_upload_key:
                seqc.log.info('Uploading read array to S3.')
                upload_ra = 'aws s3 mv {fname} {s3link}'.format(fname=args.read_array,
                                                                s3link=aws_upload_key)
                manage_ra = seqc.io.ProcessManager(upload_ra)
                manage_ra.run_all()

        # generate count matrix
        seqc.log.info('Creating sparse count matrices')
        matrices = seqc.correct_errors.convert_to_matrix(cell_counts)
        with open(args.output_stem + '_read_and_count_matrices.p', 'wb') as f:
            pickle.dump(matrices, f)
        seqc.log.info('Successfully generated sparse count matrix.')

        seqc.log.info('filtering, summarizing cell molecule count distribution, and '
                      'generating dense count matrix')
        e = seqc.core.Experiment.from_count_matrices(
            args.output_stem + '_read_and_count_matrices.p')

        # todo @ajc speed up this function; is slow
        gene_id_map = seqc.core.Experiment.create_gene_id_to_official_gene_symbol_map(
            args.index + 'annotations.gtf')
        e = e.ensembl_gene_id_to_official_gene_symbol(gene_id_map=gene_id_map)
        dense, total_molecules, mols_lost, cells_lost, cell_description = (
            e.create_filtered_dense_count_matrix())
        df = dense.molecules
        df[df == 0] = np.nan  # for sparse_csv saving
        sparse_csv = args.output_stem + '_counts.csv'
        df.to_csv(sparse_csv)

        # gzipping sparse_csv and uploading to S3
        if pigz:
            sparse_zip = "pigz --best {fname}".format(fname=sparse_csv)
        else:
            sparse_zip = "gzip {fname}".format(fname=sparse_csv)
        sparse_upload = 'aws s3 mv {fname} {s3link}'.format(fname=sparse_csv+'.gz',
                                                            s3link=aws_upload_key)
        if aws_upload_key:
            sparse_proc = seqc.io.ProcessManager(sparse_zip, sparse_upload)
        else:
            sparse_proc = seqc.io.ProcessManager(sparse_zip)
        sparse_proc.run_all()

        # in this version, local runs won't be able to upload to S3
        # and also won't get an e-mail notification.
        if args.aws:
            upload_data_and_notify_user(
                args.email_status, aws_upload_key, align, input_data, manage_merged,
                process_samfile, args.merged_fastq, args.samfile, args.read_array,
                sparse_proc, sparse_csv, manage_ra, summary, fastq_records, sam_records,
                total_molecules, mols_lost, cells_lost, manage_samfile, cell_description,
                args.output_stem, args.log_name, email)

    except BaseException as e:
        if not isinstance(e, SystemExit):
            seqc.log.exception()
        err_status = True
        if args.email_status and not args.remote:
            email_body = 'Process interrupted -- see attached error message'
            if email:
                if args.aws:
                    attachment = '/data/' + args.log_name
                else:
                    attachment = args.log_name
                seqc.remote.email_user(attachment=attachment, email_body=email_body,
                                       email_address=args.email_status)
        raise  # re-raise exception

    finally:
        if not args.remote:  # Is local
            if args.no_terminate == 'on-success':
                if err_status:
                    args.no_terminate = 'True'
                else:
                    args.no_terminate = 'False'
            if args.no_terminate in ['False', 'false']:  # terminate = True
                fpath = '/data/instance.txt'
                if os.path.isfile(fpath):
                    with open(fpath, 'r') as f:
                        inst_id = f.readline().strip('\n')
                    seqc.remote.terminate_cluster(inst_id)
                else:
                    seqc.log.info('File containing instance id is unavailable!')
            else:
                seqc.log.info('Not terminating cluster -- user responsible for cleanup')


if __name__ == '__main__':
    main(sys.argv[1:])
