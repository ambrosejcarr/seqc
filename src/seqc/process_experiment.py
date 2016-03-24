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


class ConfigurationError(Exception):
    pass


def parse_args(args):
    p = argparse.ArgumentParser(description='Process Single-Cell RNA Sequencing Data')
    p.add_argument('platform',
                   choices=['in_drop', 'drop_seq', 'mars1_seq',
                            'mars2_seq', 'in_drop_v2'],
                   help='which platform are you merging annotations from?')

    a = p.add_argument_group('required arguments')
    a.add_argument('--barcode-files', nargs='*', metavar='BF', default=list(),
                   help='text file(s) containing valid cell barcodes (one barcode per '
                        'line)')
    a.add_argument('-o', '--output-stem', metavar='O', required=True,
                   help='file stem for output files e.g. ./seqc_output/tumor_run5')
    a.add_argument('-i', '--index', metavar='I', required=True,
                   help='Folder or s3 link to folder containing index files for '
                        'alignment and resolution of ambiguous reads.')

    i = p.add_argument_group('input arguments')
    i.add_argument('-g', '--genomic-fastq', nargs='*', metavar='G', default=[],
                   help='fastq file(s) containing genomic information')
    i.add_argument('-b', '--barcode-fastq', nargs='*', metavar='B', default=[],
                   help='fastq file(s) containing barcode information')
    i.add_argument('-m', '--merged-fastq', nargs='?', metavar='M', default='',
                   help='fastq file containing genomic information annotated with '
                        'barcode data')
    i.add_argument('-s', '--samfile', nargs='?', metavar='S', default='',
                   help='sam file containing aligned, merged fastq records.')
    i.add_argument('-r', '--read-array', nargs='?', metavar='RA', default='',
                   help='ReadArray archive containing processed sam records')
    i.add_argument('--basespace', metavar='BS',
                   help='BaseSpace sample ID. Identifies a sequencing run to download '
                        'and process.')

    f = p.add_argument_group('filter arguments')
    f.add_argument('--max-insert-size', metavar='F',
                   help='maximum paired-end insert size that is considered a valid '
                        'record',
                   default=1000)
    f.add_argument('--min-poly-t', metavar='T',
                   help='minimum size of poly-T tail that is required for a barcode to '
                        'be considered a valid record',
                   default=3)
    f.add_argument('--max-dust-score', metavar='D', help='maximum complexity score for a '
                                                         'read to be considered valid')

    r = p.add_argument_group('run options')
    r.set_defaults(remote=True)
    r.add_argument('--local', dest="remote", action="store_false",
                   help='run SEQC locally instead of initiating on AWS EC2 servers')
    r.add_argument('--aws', default=False, action='store_true',
                   help='automatic flag; no need for user specification')
    r.add_argument('--email-status', metavar='E', default=None,
                   help='email address to receive run summary when running remotely')
    r.add_argument('--cluster-name', default=None, metavar='C',
                   help='optional name for EC2 instance')
    r.add_argument('--no-terminate', default=False, action='store_true',
                   help='do not terminate the EC2 instance after program completes')
    r.add_argument('--aws-upload-key', default=None, metavar='A',
                   help='location to upload results')
    r.add_argument('--reverse-complement', default=False, action='store_true',
                   help='indicates that provided barcode files contain reverse '
                        'complements of what will be found in the sequencing data')

    p.add_argument('-v', '--version', action='version',
                   version='{} {}'.format(p.prog, seqc.__version__))

    return p.parse_args(args)


def check_arguments():
    """determine if valid arguments were passed before initiating run

    :return:
    """
    # check: at least one input was provided
    # check: if s3 link, key exist
    # check: if basespace id, id exists
    raise NotImplementedError


def run_remote(name: str, stem: str) -> None:
    """
    :param name: cluster name if provided by user, otherwise None
    """
    seqc.log.notify('Beginning remote SEQC run...')

    # recreate remote command, but instruct it to run locally on the server.
    cmd = 'process_experiment.py ' + ' '.join(sys.argv[1:]) + ' --local --aws'

    # set up remote cluster
    cluster = seqc.remote.ClusterServer()
    cluster.cluster_setup(name)
    cluster.serv.connect()

    seqc.log.notify('Beginning remote run.')
    # writing name of instance in local machine to keep track of instance
    with open(os.path.expanduser('~/.seqc/instance.txt'), 'a') as f:
        _, run_name = os.path.split(stem)
        f.write('%s:%s' % (cluster.inst_id.instance_id, run_name))

    # writing name of instance in /path/to/output_dir/instance.txt for clean up
    inst_path = '/data/instance.txt'
    cluster.serv.exec_command(
        'echo {instance_id} > {inst_path}'.format(inst_path=inst_path, instance_id=str(
            cluster.inst_id.instance_id)))
    cluster.serv.exec_command('mkdir /home/ubuntu/.seqc')
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
                    'completes.')


def cluster_cleanup():
    """checks for all security groups that are unused and deletes
    them prior to each SEQC run. Updates ~/.seqc/instance.txt as well"""
    cmd = 'aws ec2 describe-instances --output text'
    cmd2 = 'aws ec2 describe-security-groups --output text'
    in_use = set([x for x in check_output(cmd.split()).split() if x.startswith(b'sg-')])
    all_sgs = set([x for x in check_output(cmd2.split()).split() if x.startswith(b'sg-')])
    to_delete = list(all_sgs - in_use)
    for sg in to_delete:
        seqc.remote.remove_sg(sg.decode())

    # check which instances in instances.txt are still running
    ec2 = boto3.resource('ec2')
    inst_file = os.path.expanduser('~/.seqc/instance.txt')

    if os.path.isfile(inst_file):
        with open(inst_file, 'r') as f:
            seqc_list = [line.strip('\n') for line in f]
        if seqc_list:
            with open(inst_file, 'w') as f:
                for i in range(len(seqc_list)):
                    entry = seqc_list[i]
                    inst_id = entry.split(':')[0]
                    instance = ec2.Instance(inst_id)
                    # todo: could instance potentially be None?
                    if instance.state['Name'] == 'running':
                        f.write('%s\n' % entry)
    else:
        pass  # instances.txt file has not yet been created


def check_input_data(barcode_fastq: list, genomic_fastq: list, merged: str,
                     samfile: str, basespace: str):
    """
    checks the validity of data inputs.
    args:
    -----
    barcode_fastq, genomic_fastq: Lists of fastq files OR single aws s3 links
     specifying folders holding barcode and genomic fastq files
    samfile: an alignment file or an aws s3 link specifying the location of a samfile
    merged: a merged fastq file or an aws s3 link specifying the location of a merged
     fastq file
    basespace: the sample id necessary to download files from BaseSpace
    """

    # make sure at least one input has been passed
    if not any([barcode_fastq, genomic_fastq, merged, samfile, basespace]):
        raise ValueError('At least one input argument must be passed to SEQC.')

    # make sure only one filetype has been passed
    multi_input_error_message = ('Only one input type (-s, -m, -r/-f, or --basespace) '
                                 'should be passed to SEQC.')
    unpaired_fastq_error_message = ('If either genomic or barcode fastq files are '
                                    'provided, both must be provided.')

    if barcode_fastq or genomic_fastq:
        if not all((barcode_fastq, genomic_fastq)):
            raise ValueError(unpaired_fastq_error_message)
        if any((merged, samfile)):
            raise ValueError(multi_input_error_message)
    if samfile:
        if any((merged, barcode_fastq, genomic_fastq)):
            raise ValueError(multi_input_error_message)
    if merged:
        if any((samfile, barcode_fastq, genomic_fastq)):
            raise ValueError(multi_input_error_message)


def s3bucket_download(s3link: str, outdir: str):
    """recursively downloads files from given S3 link"""
    bucket, prefix = seqc.io.S3.split_link(s3link)
    cut_dirs = prefix.count('/')
    downloaded_files = seqc.io.S3.download_files(bucket, prefix, outdir, cut_dirs)
    return sorted(downloaded_files)


def s3files_download(s3links: list, outdir: str):
    """downloads each file in list of s3 links"""
    fnames = []
    for s3link in s3links:
        bucket, prefix = seqc.io.S3.split_link(s3link)
        _, fname = os.path.split(prefix)
        fname = outdir + '/' + fname
        seqc.io.S3.download_file(bucket, prefix, fname)
        fnames.append(fnames)
    return sorted(fnames)


def main(args: list = None):
    seqc.log.setup_logger()
    args = parse_args(args)
    try:
        seqc.log.args(args)
        check_input_data(args.barcode_fastq, args.genomic_fastq, args.merged_fastq,
                         args.samfile, args.basespace)

        # split output_stem into path and prefix, read in config file
        output_dir, output_prefix = os.path.split(args.output_stem)
        config_file = os.path.expanduser('~/.seqc/config')
        config = configparser.ConfigParser()
        if not config.read(config_file):
            raise ConfigurationError('Please run ./configure (found in the seqc '
                                     'directory) before attempting to run '
                                     'process_experiment.py.')

        if args.remote:
            cluster_cleanup()
            run_remote(args.cluster_name, args.output_stem)
            sys.exit()

        if args.aws:
            args.output_stem = '/data/' + output_prefix
            output_dir, output_prefix = os.path.split(args.output_stem)

        # do a bit of argument checking
        if args.output_stem.endswith('/'):
            seqc.log.notify('-o/--output-stem should not be a directory.')
            sys.exit(2)

        # download data if necessary
        if args.basespace:
            seqc.log.info('BaseSpace link provided for fastq argument. Downloading '
                          'input data.')
            basespace_token = config['BaseSpaceToken']['base_space_token']
            if basespace_token == 'None':
                raise ValueError(
                    'If the --basespace argument is used, the BaseSpace token must be '
                    'specified in the seqc config file.')

            # making extra directories for BaseSpace download, changing permissions
            # todo: this is unnecessary if being run locally
            bspace_dir = output_dir + '/Data/Intensities/BaseCalls/'
            bf = Popen(['sudo', 'mkdir', '-p', bspace_dir])
            bf.communicate()
            bf2 = Popen(['sudo', 'chown', '-c', 'ubuntu', bspace_dir])
            bf2.communicate()
            args.barcode_fastq, args.genomic_fastq = seqc.io.BaseSpace.download(
                args.platform, args.basespace, output_dir, basespace_token)

        # check for genomic fastq files
        if not args.genomic_fastq[0].startswith('s3://'):
            for gf in args.genomic_fastq:
                if not os.path.isfile(gf):
                    raise ValueError('provided genomic fastq files: "[%s]" is neither '
                                     'an s3 link or a valid filepath' %
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
                    args.genomic_fastq = s3files_download(args.genomic_fastq, output_dir)
                    # todo: need to log
            except FileNotFoundError:
                raise FileNotFoundError('No fastq files were found at the specified '
                                        's3 location: [%s]' %
                                        ', '.join(map(str, args.genomic_fastq)))
            except FileExistsError:
                pass  # file is already present.

        # same check for barcode fastq files --> todo: could make this less redundant
        if not args.barcode_fastq[0].startswith('s3://'):
            for bf in args.barcode_fastq:
                if not os.path.isfile(bf):
                    raise ValueError('provided genomic fastq files: "[%s]" is neither '
                                     'an s3 link or a valid filepath' %
                                     ', '.join(map(str, args.barcode_fastq)))
        else:
            try:
                seqc.log.info('Downloading barcode fastq files from Amazon s3 link.')
                if args.barcode_fastq[0].endswith('/'):
                    args.barcode_fastq = s3bucket_download(args.barcode_fastq[0],
                                                           output_dir)
                else:
                    args.barcode_fastq = s3files_download(args.barcode_fastq, output_dir)
                    # todo: need to log
            except FileNotFoundError:
                raise FileNotFoundError('No fastq files were found at the specified '
                                        's3 location: [%s]' %
                                        ', '.join(map(str, args.barcode_fastq)))
            except FileExistsError:
                pass  # file is already present.

        # check if the index must be downloaded
        if not args.index.startswith('s3://'):
            if not os.path.isdir(args.index):
                raise ValueError('provided index: "%s" is neither an s3 link or a valid '
                                 'filepath' % args.index)
        else:
            try:
                seqc.log.info('AWS s3 link provided for index. Downloading index.')
                bucket, prefix = seqc.io.S3.split_link(args.index)
                args.index = output_dir + '/index/'  # set index  based on s3 download
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

        if align:
            # todo: need to check for S3 installation of args.merged_fastq
            seqc.log.info('Aligning merged fastq records.')
            *base_directory, stem = args.output_stem.split('/')
            alignment_directory = '/'.join(base_directory) + '/alignments/'
            os.makedirs(alignment_directory, exist_ok=True)
            args.samfile = seqc.alignment.star.align(
                args.merged_fastq, args.index, n_processes, alignment_directory)

        if process_samfile:
            # todo: need to check for S3 installation of args.samfile
            seqc.log.info('Filtering aligned records and constructing record database')
            ra = seqc.arrays.ReadArray.from_samfile(
                args.samfile, args.index + 'annotations.gtf')
            ra.save(args.output_stem + '.h5')
        else:
            # todo: need to check for S3 installation of args.read_array
            ra = seqc.arrays.ReadArray.load(args.read_array)

        seqc.log.info('Correcting cell barcode and RMT errors')
        # check if barcode files need to be downloaded
        if not args.barcode_files[0].startswith('s3://'):
            for cb in args.barcode_files:
                if not os.path.isfile(cb):
                    raise ValueError('provided barcode files: "[%s]" is neither '
                                     'an s3 link or a valid filepath' %
                                     ', '.join(map(str, args.barcode_files)))
        else:
            try:
                seqc.log.info('AWS s3 link provided for barcodes. Downloading files.')
                args.barcode_files = s3bucket_download(args.barcode_files[0], output_dir)
            except FileNotFoundError:
                raise FileNotFoundError('No barcode files were found at the specified '
                                        's3 location: %s' % args.barcode_files[0])
            except FileExistsError:
                pass  # file is already present.
        cell_counts, _ = seqc.correct_errors.correct_errors(
            ra, args.barcode_files, reverse_complement=args.reverse_complement)

        seqc.log.info('Creating count matrices')
        matrices = seqc.correct_errors.convert_to_matrix(cell_counts)
        with open(args.output_stem + '_read_and_count_matrices.p', 'wb') as f:
            pickle.dump(matrices, f)

        seqc.log.info('Run complete.')

        if args.email_status:
            seqc.remote.upload_results(
                args.output_stem, args.email_status, args.aws_upload_key)

    except:
        pass
        seqc.log.exception()
        if args.email_status and not args.remote:
            email_body = 'Process interrupted -- see attached error message'
            seqc.remote.email_user(attachment='seqc.log', email_body=email_body,
                                   email_address=args.email_status)
        raise

    finally:
        if not args.remote:  # Is local
            if not args.no_terminate:  # terminate = True
                fpath = output_dir + '/instance.txt'
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
