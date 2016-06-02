#!/usr/local/bin/python3

import argparse
import multiprocessing
import os
import sys
import pickle
import seqc
import configparser
import boto3
from subprocess import Popen, check_output
import numpy as np


class ConfigurationError(Exception):
    pass


class ArgumentParserError(Exception):
    pass


class NewArgumentParser(argparse.ArgumentParser):
    def error(self, message):
        # checks to see whether user wants to check remote experiment status
        if '--check-progress' in sys.argv[1:]:
            seqc.remote.check_progress()
        else:
            print(message)
        sys.exit(0)


def parse_args(args):
    p = NewArgumentParser(description='Process Single-Cell RNA Sequencing Data')
    p.add_argument('platform',
                   choices=['in_drop', 'drop_seq', 'mars1_seq',
                            'mars2_seq', 'in_drop_v2'],
                   help='which platform are you merging annotations from?')

    a = p.add_argument_group('required arguments')
    a.add_argument('--barcode-files', nargs='*', metavar='BF', default=list(),
                   help='Either (a) an s3 link to a folder containing only barcode '
                        'files, or (b) the full file path of each file on the local '
                        'machine.')
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
                        'of the BaseSpace sample. (e.g. if the link to the sample is'
                        'https://basespace.illumina.com/sample/34000253/0309, '
                        '--basespace=34000253.')

    f = p.add_argument_group('filter arguments')
    f.add_argument('--max-insert-size', metavar='F',
                   help='maximum paired-end insert size that is considered a valid '
                        'record. For multialignment correction. Not currently used.',
                   default=1000)
    f.add_argument('--min-poly-t', metavar='T',
                   help='minimum size of poly-T tail that is required for a barcode to '
                        'be considered a valid record (default=3)',
                   default=3, type=int)
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

    r = p.add_argument_group('run options')
    r.set_defaults(remote=True)
    r.set_defaults(check=False)
    r.add_argument('--local', dest="remote", action="store_false",
                   help='Run SEQC locally instead of initiating on AWS EC2 servers.')
    r.add_argument('--aws', default=False, action='store_true',
                   help='Automatic flag; no need for user specification.')
    r.add_argument('--email-status', metavar='E', default=None,
                   help='Email address to receive run summary or errors when running '
                        'remotely.')
    r.add_argument('--no-terminate', default=False, action='store_true',
                   help='Do not terminate the EC2 instance after program completes.')
    r.add_argument('--check-progress', dest="check", action="store_true",
                   help='Check progress of all currently running remote SEQC runs.')
    r.add_argument('--star-args', default=None, nargs='*',
                   help='additional arguments that should be passed to the STAR '
                        'aligner. For example, to set the maximum allowable times for a '
                        'read to align to 20, one would set '
                        '--star-args outFilterMultimapNmax=20. Additional arguments can '
                        'be provided as a white-space separated list.')
    r.add_argument('--instance-type', default='c4',
                   help='AWS instance (c3, c4, r3) to run SEQC remotely. Default=c4.')
    r.add_argument('--spot-bid', type=float, default=None,
                   help='Amount to bid for a spot instance. Default=None (will reserve a '
                        'non-spot instance). WARNING: using spot instances will cause '
                        'your instance to terminate if instance prices exceed your spot '
                        'bid during runtime.')

    p.add_argument('-v', '--version', action='version',
                   version='{} {}'.format(p.prog, seqc.__version__))

    try:
        return p.parse_args(args)
    except ArgumentParserError:
        raise


def run_remote(stem: str, volsize: int, aws_instance:str, spot_bid=None) -> None:
    """
    :param stem: output_prefix from main()
    :param volsize: estimated volume needed for run
    :param aws_instance: instance type from user
    """
    seqc.log.notify('Beginning remote SEQC run...')

    # recreate remote command, but instruct it to run locally on the server.
    cmd = 'process_experiment.py ' + ' '.join(sys.argv[1:]) + ' --local --aws'

    # set up remote cluster
    cluster = seqc.remote.ClusterServer()
    volsize = int(np.ceil(volsize/(1e9)))

    # if anything goes wrong during cluster setup, clean up the instance
    try:
        cluster.cluster_setup(volsize, aws_instance)
        cluster.serv.connect()

        seqc.log.notify('Beginning remote run.')
        if stem.endswith('/'):
            stem = stem[:-1]
        # writing name of instance in local machine to keep track of instance
        with open(os.path.expanduser('~/.seqc/instance.txt'), 'a') as f:
            _, run_name = os.path.split(stem)
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
        sys.exit(2)


def cluster_cleanup():
    """checks for all security groups that are unused and deletes
    them prior to each SEQC run. Updates ~/.seqc/instance.txt as well"""

    cmd = 'aws ec2 describe-instances --output text'
    cmd2 = 'aws ec2 describe-security-groups --output text'
    in_use = set([x for x in check_output(cmd.split()).split() if x.startswith(b'sg-')])
    all_sgs = set([x for x in check_output(cmd2.split()).split() if x.startswith(b'sg-')])
    to_delete = list(all_sgs - in_use)

    # set up ec2 resource from boto3
    ec2 = boto3.resource('ec2')

    # iteratively remove unused SEQC security groups
    for sg in to_delete:
        sg_id = sg.decode()
        sg_name = ec2.SecurityGroup(sg_id).group_name
        if 'SEQC-' in sg_name:
            seqc.remote.remove_sg(sg_id)

    # check which instances in instances.txt are still running
    inst_file = os.path.expanduser('~/.seqc/instance.txt')

    if os.path.isfile(inst_file):
        with open(inst_file, 'r') as f:
            seqc_list = [line.strip('\n') for line in f]
        if seqc_list:
            with open(inst_file, 'w') as f:
                for i in range(len(seqc_list)):
                    try:
                        entry = seqc_list[i]
                        inst_id = entry.split(':')[0]
                        instance = ec2.Instance(inst_id)
                        if instance.state['Name'] == 'running':
                            f.write('%s\n' % entry)
                    except:
                        continue
    else:
        pass  # instances.txt file has not yet been created


def check_s3links(input_args: list):
    """determine if valid arguments were passed before initiating run,
    specifically whether s3 links exist
    :param input_args: list of files that should be checked
    """

    s3 = boto3.resource('s3')
    for infile in input_args:
        try:
            if infile.startswith('s3://'):
                if not infile.endswith('/'):  # check that s3 link for file exists
                    bucket, key = seqc.io.S3.split_link(infile)
                    s3.meta.client.head_object(Bucket=bucket, Key=key)
                else:
                    cmd = 'aws s3 ls ' + infile  # directory specified in s3 link
                    res = check_output(cmd.split())
                    if b'PRE ' in res:  # subdirectories present
                        raise ValueError
        except:
            seqc.log.notify('Error: Provided s3 link "%s" does not contain the proper '
                            'input files to SEQC.' % infile)
            sys.exit(2)


def check_arguments(args, basespace_token: str):
    """
    checks data input through the command line arguments and throws
    an error if the provided input data is invalid.

    additionally, this function obtains a rough estimate of how much
    volume storage is needed for the overall SEQC run.
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

    # check to make sure that --email-status is passed with remote run
    if args.remote and not args.email_status:
        raise ValueError('Please supply the --email-status flag for a remote SEQC run.')
    if args.instance_type not in ['c3', 'c4', 'r3']:
        raise ValueError('All AWS instance types must be either c3, c4, or r3.')

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


def s3bucket_download(s3link: str, outdir: str):
    """
    recursively downloads files from given S3 link
    :param s3link: link to s3 bucket that holds files to download
    :param outdir: output directory where files will be downloaded
    returns sorted list of downloaded files
    """

    bucket, prefix = seqc.io.S3.split_link(s3link)
    cut_dirs = prefix.count('/')
    downloaded_files = seqc.io.S3.download_files(bucket, prefix, outdir, cut_dirs)
    return sorted(downloaded_files)


def s3files_download(s3links: list, outdir: str):
    """
    downloads each file in list of s3 links
    :param s3links: list of s3 links to be downloaded
    :param outdir: output directory where files will be downloaded
    returns sorted list of downloaded files
    """

    fnames = []
    for s3link in s3links:
        bucket, prefix = seqc.io.S3.split_link(s3link)
        _, fname = os.path.split(prefix)
        fname = outdir + '/' + fname
        seqc.io.S3.download_file(bucket, prefix, fname)
        fnames.append(fname)
    return sorted(fnames)


def obtain_size(item):
    cmd = 'aws s3 ls --summarize --recursive ' + item + ' | grep "Total Size"'
    obj_size = int(check_output(cmd, shell=True).decode().split()[-1])
    return obj_size


def main(args: list = None):
    seqc.log.setup_logger()
    args = parse_args(args)
    try:
        seqc.log.args(args)

        # read in config file, make sure it exists
        config_file = os.path.expanduser('~/.seqc/config')
        config = configparser.ConfigParser()
        if not config.read(config_file):
            raise ConfigurationError('Please run ./configure (found in the seqc '
                                     'directory) before attempting to run '
                                     'process_experiment.py.')

        if args.spot_bid is not None:
            if args.spot_bid < 0:
                seqc.log.notify('"{bid}" must be a non-negative float! Exiting.'.format(
                    spot=args.spot_bid))
                sys.exit(2)

        # extract basespace token and make sure args.index is a directory
        basespace_token = config['BaseSpaceToken']['base_space_token']
        if not args.index.endswith('/'):
            args.index += '/'

        if not args.aws:  # if args.aws, arguments would already have passed checks.
            total_size = check_arguments(args, basespace_token)
        else:
            total_size = None

        # check whether output is an s3 or local directory, split output_stem
        if args.output_stem.endswith('/'):
            if args.output_stem.startswith('s3://'):
                output_dir, output_prefix = os.path.split(args.output_stem[:-1])
            else:
                seqc.log.notify('-o/--output-stem should not be a directory for local '
                                'SEQC runs.')
                sys.exit(2)
        else:
            output_dir, output_prefix = os.path.split(args.output_stem)

        # check output directory if remote run
        if args.remote:
            if not args.output_stem.startswith('s3://'):
                raise ValueError('-o/--output-stem must be an s3 link for remote SEQC '
                                 'runs.')
            cluster_cleanup()
            run_remote(args.output_stem, total_size, args.instance_type,
                       spot_bid=args.spot_bid)
            sys.exit()

        if args.aws:
            aws_upload_key = args.output_stem
            if not aws_upload_key.endswith('/'):
                aws_upload_key += '/'
            args.output_stem = '/data/' + output_prefix
            output_dir, output_prefix = os.path.split(args.output_stem)
        else:
            # todo: need to fix usage of aws_upload_key for local runs
            aws_upload_key = None

        # download data if necessary
        if args.basespace:
            seqc.log.info('BaseSpace link provided for fastq argument. Downloading '
                          'input data.')
            # making extra directories for BaseSpace download, changing permissions
            bspace_dir = output_dir + '/Data/Intensities/BaseCalls/'
            bf = Popen(['sudo', 'mkdir', '-p', bspace_dir])
            bf.communicate()
            if args.aws:  # changing permissions is unnecessary if local run
                bf2 = Popen(['sudo', 'chown', '-c', 'ubuntu', bspace_dir])
                bf2.communicate()
            args.barcode_fastq, args.genomic_fastq = seqc.io.BaseSpace.download(
                args.platform, args.basespace, output_dir, basespace_token)

        # check for genomic fastq files
        if args.genomic_fastq:
            if not args.genomic_fastq[0].startswith('s3://'):
                for gf in args.genomic_fastq:
                    if not os.path.isfile(gf):
                        raise ValueError('Provided genomic fastq files: "[%s]" is '
                                         'neither an s3 link or a valid filepath' %
                                         ', '.join(map(str, args.genomic_fastq)))
            else:
                try:
                    seqc.log.info('Downloading genomic fastq files from Amazon s3 link.')
                    if args.genomic_fastq[0].endswith('/'):
                        # s3 directory specified, download all files recursively
                        args.genomic_fastq = s3bucket_download(args.genomic_fastq[0],
                                                           output_dir)
                    else:
                        # individual s3 links provided, download each fastq file
                        args.genomic_fastq = s3files_download(args.genomic_fastq,
                                                              output_dir)
                    seqc.log.info('Genomic fastq files [%s] successfully installed.' %
                                  ', '.join(map(str, args.genomic_fastq)))
                except FileExistsError:
                    pass  # file is already present.

        # todo: could make this less redundant. same code block but for barcode fastq
        if args.barcode_fastq:
            if not args.barcode_fastq[0].startswith('s3://'):
                for bf in args.barcode_fastq:
                    if not os.path.isfile(bf):
                        raise ValueError('provided genomic fastq files: "[%s]" is '
                                         'neither an s3 link or a valid filepath' %
                                         ', '.join(map(str, args.barcode_fastq)))
            else:
                try:
                    seqc.log.info('Downloading barcode fastq files from Amazon s3 link.')
                    if args.barcode_fastq[0].endswith('/'):
                        args.barcode_fastq = s3bucket_download(args.barcode_fastq[0],
                                                               output_dir)
                    else:
                        args.barcode_fastq = s3files_download(args.barcode_fastq,
                                                              output_dir)
                    seqc.log.info('Barcode fastq files [%s] successfully installed.' %
                                  ', '.join(map(str, args.barcode_fastq)))
                except FileExistsError:
                    pass  # file is already present.

        # check if the index must be downloaded
        if not args.index.startswith('s3://'):
            if not os.path.isdir(args.index):
                raise ValueError('Provided index: "%s" is neither an s3 link or a valid '
                                 'filepath' % args.index)
        else:
            if not args.read_array:
                seqc.log.info('AWS s3 link provided for index. Downloading index.')
                bucket, prefix = seqc.io.S3.split_link(args.index)
                args.index = output_dir + '/index/'  # set index  based on s3 download
                cut_dirs = prefix.count('/')
                # install whole index
                if not args.samfile:
                    try:
                        seqc.io.S3.download_files(bucket, prefix, args.index, cut_dirs)
                    except FileExistsError:
                        pass  # file is already present
                else:  # samfile provided, only download annotations file
                    try:
                        annotations_file = args.index + 'annotations.gtf'
                        seqc.io.S3.download_file(bucket, prefix + 'annotations.gtf',
                                                 annotations_file)
                    except FileExistsError:
                        pass  # file is already present

        n_processes = multiprocessing.cpu_count() - 1  # get number of processors

        # determine where the script should start:
        input_data = 'start'
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

        if merge:
            seqc.log.info('Merging genomic reads and barcode annotations.')
            merge_function = getattr(seqc.sequence.merge_functions, args.platform)
            args.merged_fastq = seqc.sequence.fastq.merge_paired(
                merge_function=merge_function,
                fout=args.output_stem + '_merged.fastq',
                genomic=args.genomic_fastq,
                barcode=args.barcode_fastq)

            # deleting genomic/barcode fastq files after merged.fastq creation
            seqc.log.info('Removing original fastq file for memory management.')
            delete_fastq = ' '.join(['rm'] + args.genomic_fastq + args.barcode_fastq)
            seqc.io.ProcessManager(delete_fastq).run_all()

            # wait until merged fastq file is zipped before alignment
            seqc.log.info('Gzipping merged fastq file.')
            pigz_zip = "pigz --best -k {fname}".format(fname=args.merged_fastq)
            pigz_proc = seqc.io.ProcessManager(pigz_zip)
            pigz_proc.run_all()
            pigz_proc.wait_until_complete()  # prevents slowing down STAR alignment

        if align:
            seqc.log.info('Aligning merged fastq records.')
            if args.merged_fastq.startswith('s3://'):
                input_data = 'merged'
                args.merged_fastq = s3files_download([args.merged_fastq], output_dir)[0]
                if args.merged_fastq.endswith('.gz'):
                    cmd = 'pigz -d {fname}'.format(fname=args.merged_fastq)
                    gunzip_proc = seqc.io.ProcessManager(cmd)
                    gunzip_proc.run_all()
                    gunzip_proc.wait_until_complete()  # finish before STAR alignment
                    args.merged_fastq = args.merged_fastq.strip('.gz')
                seqc.log.info('Merged fastq file %s successfully installed from S3.' %
                              args.merged_fastq)
            *base_directory, stem = args.output_stem.split('/')
            alignment_directory = '/'.join(base_directory) + '/alignments/'
            os.makedirs(alignment_directory, exist_ok=True)
            if args.star_args is not None:
                star_kwargs = dict(a.strip().split('=') for a in args.star_args)
            else:
                star_kwargs = {}
            args.samfile = seqc.alignment.star.align(
                args.merged_fastq, args.index, n_processes, alignment_directory,
                **star_kwargs)

            # Removing or uploading the merged fastq file depending on start point
            if input_data == 'merged':
                seqc.log.info('Removing merged.fastq file for memory management.')
                rm_cmd = 'rm {merged_file}'.format(merged_file=args.merged_fastq)
                seqc.io.ProcessManager(rm_cmd).run_all()
            else:
                seqc.log.info('Uploading gzipped merged fastq file to S3.')
                merge_upload = 'aws s3 mv {fname} {s3link}'.format(
                    fname=args.merged_fastq+'.gz', s3link=aws_upload_key)
                manage_merged = seqc.io.ProcessManager(merge_upload)
                manage_merged.run_all()

        if process_samfile:
            seqc.log.info('Filtering aligned records and constructing record database.')
            if args.samfile.startswith('s3://'):
                input_data = 'samfile'
                args.samfile = s3files_download([args.samfile], output_dir)[0]
                seqc.log.info('Samfile %s successfully installed from S3.' % args.samfile)
            ra = seqc.core.ReadArray.from_samfile(
                args.samfile, args.index + 'annotations.gtf')
            args.read_array = args.output_stem + '.h5'
            ra.save(args.read_array)

            # converting sam to bam and uploading to S3, else removing samfile
            if input_data == 'samfile':
                seqc.log.info('Removing samfile for memory management.')
                rm_samfile = 'rm {fname}'.format(fname=args.samfile)
                seqc.io.ProcessManager(rm_samfile).run_all()
            else:
                seqc.log.info('Converting samfile to bamfile and uploading to S3.')
                bamfile = output_dir + '/alignments/Aligned.out.bam'
                convert_sam = 'samtools view -bS -o {bamfile} {samfile}'.\
                    format(bamfile=bamfile, samfile=args.samfile)
                upload_bam = 'aws s3 mv {fname} {s3link}'.format(fname=bamfile,
                                                                 s3link=aws_upload_key)
                manage_samfile = seqc.io.ProcessManager(convert_sam, upload_bam)
                manage_samfile.run_all()
        else:
            if args.read_array.startswith('s3://'):
                input_data = 'readarray'
                args.read_array = s3files_download([args.read_array], output_dir)[0]
                seqc.log.info('Read array %s successfully installed from S3.' %
                              args.read_array)
            ra = seqc.core.ReadArray.load(args.read_array)

        # check if barcode files need to be downloaded
        if args.platform != 'drop_seq':
            if not args.barcode_files[0].startswith('s3://'):
                for cb in args.barcode_files:
                    if not os.path.isfile(cb):
                        raise ValueError('provided barcode files: "[%s]" is neither '
                                         'an s3 link or a valid filepath' %
                                         ', '.join(map(str, args.barcode_files)))
            else:
                try:
                    seqc.log.info('AWS s3 link provided for barcodes. Downloading files.')
                    if not args.barcode_files[0].endswith('/'):
                        args.barcode_files[0] += '/'
                    args.barcode_files = s3bucket_download(args.barcode_files[0], output_dir)
                    seqc.log.info('Barcode files [%s] successfully installed.' %
                                  ', '.join(map(str, args.barcode_files)))
                except FileExistsError:
                    pass  # file is already present.

        # correct errors
        seqc.log.info('Correcting cell barcode and RMT errors')
        correct_errors = getattr(seqc.correct_errors, args.platform)
        cell_counts, _ = correct_errors(
            ra, args.barcode_files, reverse_complement=False,
            required_poly_t=args.min_poly_t, max_ed=args.max_ed,
            singleton_weight=args.singleton_weight)

        # uploading read array to S3 if created, else removing read array
        if input_data == 'readarray':
            seqc.log.info('Removing .h5 file for memory management.')
            rm_ra = 'rm {fname}'.format(fname=args.read_array)
            seqc.io.ProcessManager(rm_ra).run_all()
        else:
            seqc.log.info('Uploading read array to S3.')
            upload_ra = 'aws s3 mv {fname} {s3link}'.format(fname=args.read_array,
                                                            s3link=aws_upload_key)
            manage_ra = seqc.io.ProcessManager(upload_ra)
            manage_ra.run_all()

        seqc.log.info('Creating count matrices')
        matrices = seqc.correct_errors.convert_to_matrix(cell_counts)
        with open(args.output_stem + '_read_and_count_matrices.p', 'wb') as f:
            pickle.dump(matrices, f)
        seqc.log.info('Successfully generated count matrix.')

        # in this version, local runs won't be able to upload to S3
        # and also won't get an e-mail notification.
        if args.aws:
            seqc.log.info('Starting file upload onto %s.' % aws_upload_key)

            if args.email_status:
                # make sure that all other files are uploaded before termination
                if align and input_data != 'merged':
                    manage_merged.wait_until_complete()
                    seqc.log.info('Successfully uploaded %s to the specified S3 '
                                  'location "%s"' % (args.merged_fastq, aws_upload_key))
                if process_samfile:
                    if input_data != 'samfile':
                        manage_samfile.wait_until_complete()
                        seqc.log.info('Successfully uploaded %s to the specified S3 '
                                      'location "%s"' % (args.samfile, aws_upload_key))
                    manage_ra.wait_until_complete()
                    seqc.log.info('Successfully uploaded %s to the specified S3 '
                                  'location "%s"' % (args.read_array, aws_upload_key))

                # upload count matrix and alignment summary at the very end
                seqc.remote.upload_results(
                    args.output_stem, args.email_status, aws_upload_key, input_data)

    except:
        seqc.log.exception()
        if args.email_status and not args.remote:
            email_body = 'Process interrupted -- see attached error message'
            seqc.remote.email_user(attachment='seqc.log', email_body=email_body,
                                   email_address=args.email_status)

        raise

    finally:
        if not args.remote:  # Is local
            if not args.no_terminate:  # terminate = True
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
