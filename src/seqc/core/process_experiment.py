#!/usr/local/bin/python3

import multiprocessing
import os
import pickle
import sys
import numpy as np
from seqc import remote, log, io, filter, correct_errors, platforms
from seqc.sequence import fastq
from seqc.alignment import star
from seqc.read_array import ReadArray
from seqc.exceptions import ConfigurationError, retry_boto_call
from seqc.sequence.gtf import (ensembl_gene_id_to_official_gene_symbol,
                               create_gene_id_to_official_gene_symbol_map)
from seqc.filter import create_filtered_dense_count_matrix
from seqc.sparse_frame import SparseFrame
from seqc.core import parser, verify, config, download, upload, execution_control


def run_remote(args, argv, volsize) -> None:
    """
    Mirror the local arguments from a seqc.core.process_experiment call to an AWS server
    and execute the run there. When complete, terminates the local process unless
    otherwise specified by the args.no_terminate argument.

    :param argv: original argument list as received from sys.argv or passed to a main()
      object
    :param args: simple namespace object; output of parse_args()
    :param volsize: estimated volume needed for run
    """
    log.notify('Beginning remote SEQC run...')

    # recreate remote command, but instruct it to run locally on the server.
    cmd = parser.generate_remote_cmdline_args(argv)
    log.print_exact_command_line(cmd)
    # set up remote cluster
    cluster = remote.ClusterServer()
    volsize = int(np.ceil(volsize/1e9))

    try:  # if anything goes wrong during cluster setup, clean up the instance
        cluster.setup_cluster(volsize, args.instance_type, spot_bid=args.spot_bid)
        cluster.serv.connect()

        log.notify('Beginning remote run.')
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
        cluster.serv.exec_command('cd /data; nohup {cmd} > debug.log 2>&1 &'
                                  ''.format(cmd=cmd))

        # check that process_experiment.py is actually running on the cluster
        out, err = cluster.serv.exec_command('ps aux | grep process_experiment.py')
        res = ' '.join(out)
        if '/usr/local/bin/process_experiment.py' not in res:
            raise ConfigurationError('Error executing SEQC on the cluster!')
        log.notify('Terminating local client. Email will be sent when remote run '
                   'completes. Please use "process_experiment.py --check-progress" '
                   'to monitor the status of the remote run.')
    except Exception as e:
        log.notify('Error {e} occurred during cluster setup!'.format(e=e))
        if cluster.is_cluster_running():
            log.notify('Cleaning up instance {id} before exiting...'.format(
                id=cluster.inst_id.instance_id))
            remote.terminate_cluster(cluster.inst_id.instance_id)
        raise  # re-raise original exception


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
        platform, barcode_fastq: [str], output_stem: str,
        genomic_fastq: [str], pigz: bool) -> (str, int):
    """
    annotates genomic fastq with barcode information; merging the two files.

    :param platform: class from platforms.py that defines the characteristics of the
      data being processed
    :param barcode_fastq: list of str names of fastq files containing barcode information
    :param output_stem: str, stem for output files
    :param genomic_fastq: list of str names of fastq files containing genomic information
    :param pigz: bool, indicates if pigz is installed
    :returns merged_fastq, fastq_records: (str, int) name of merged fastq file and the
      number of fastq records that were processed.
    """

    log.info('Merging genomic reads and barcode annotations.')
    merged_fastq, fastq_records = fastq.merge_paired(
        merge_function=platform.merge_function,
        fout=output_stem + '_merged.fastq',
        genomic=genomic_fastq,
        barcode=barcode_fastq)

    # delete genomic/barcode fastq files after merged.fastq creation
    log.info('Removing original fastq file for memory management.')
    delete_fastq = ' '.join(['rm'] + genomic_fastq + barcode_fastq)
    io.ProcessManager(delete_fastq).run_all()

    # zip merged fastq file
    log.info('Gzipping merged fastq file.')
    if pigz:
        pigz_zip = "pigz --best -k {fname}".format(fname=merged_fastq)
    else:
        pigz_zip = "gzip -kf {fname}".format(fname=merged_fastq)
    pigz_proc = io.ProcessManager(pigz_zip)
    pigz_proc.run_all()
    pigz_proc.wait_until_complete()  # prevents slowing down STAR alignment
    return merged_fastq, fastq_records


def align_fastq_records(
        merged_fastq, output_dir, output_stem, star_args, index, n_processes,
        aws_upload_key, pigz) -> (str, str, io.ProcessManager):
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
    :return samfile, input_data, manage_merged: (str, str, io.ProcessManager)
      name of .sam file containing aligned reads, indicator of which data was used as
      input, and a ProcessManager for merged fastq files
    """
    input_data = 'start'
    log.info('Aligning merged fastq records.')
    if merged_fastq.startswith('s3://'):
        input_data = 'merged'
        bucket, prefix = io.S3.split_link(merged_fastq)
        _, fname = os.path.split(prefix)
        merged_fastq = io.S3.download_file(bucket, prefix, output_dir + '/' + fname)
        if merged_fastq.endswith('.gz'):
            if pigz:
                cmd = 'pigz -d {fname}'.format(fname=merged_fastq)
            else:
                cmd = 'gunzip {fname}'.format(fname=merged_fastq)
            gunzip_proc = io.ProcessManager(cmd)
            gunzip_proc.run_all()
            gunzip_proc.wait_until_complete()  # finish before STAR alignment
            merged_fastq = merged_fastq.strip('.gz')
        log.info('Merged fastq file %s successfully installed from S3.' % merged_fastq)
    *base_directory, stem = output_stem.split('/')
    alignment_directory = '/'.join(base_directory) + '/alignments/'
    os.makedirs(alignment_directory, exist_ok=True)
    if star_args is not None:
        star_kwargs = dict(a.strip().split('=') for a in star_args)
    else:
        star_kwargs = {}
    samfile = star.align(
        merged_fastq, index, n_processes, alignment_directory,
        **star_kwargs)

    # Removing or uploading the merged fastq file depending on start point
    if input_data == 'merged':
        log.info('Removing merged.fastq file for memory management.')
        rm_cmd = 'rm {merged_file}'.format(merged_file=merged_fastq)
        io.ProcessManager(rm_cmd).run_all()
        manage_merged = None
    else:
        if aws_upload_key:
            log.info('Uploading gzipped merged fastq file to S3.')
            merge_upload = 'aws s3 mv {fname} {s3link}'.format(
                fname=merged_fastq + '.gz', s3link=aws_upload_key)
            manage_merged = io.ProcessManager(merge_upload)
            manage_merged.run_all()
        else:
            manage_merged = None
    return samfile, input_data, manage_merged


def create_or_download_read_array(
        process_samfile, samfile, output_dir, index, read_array, output_stem,
        aws_upload_key, input_data) -> (str, str, io.ProcessManager, int):
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
    :returns ra, input_data, manage_samfile: (str, str, io.ProcessManager, int)
      ReadArray object, indicator object, samfile ProcessManager, and number of processed
      sam records
    """
    manage_samfile = None
    sam_records = None
    if process_samfile:
        log.info('Constructing record database.')
        if samfile.startswith('s3://'):
            input_data = 'samfile'
            bucket, prefix = io.S3.split_link(samfile)
            _, fname = os.path.split(prefix)
            samfile = io.S3.download_file(bucket, prefix, output_dir + '/' + fname)
            log.info('Samfile %s successfully installed from S3.' % samfile)
        ra, sam_records = ReadArray.from_samfile(samfile, index + 'annotations.gtf')
        read_array = output_stem + '.h5'
        #ra.save(read_array)

        # converting sam to bam and uploading to S3, else removing samfile
        if input_data == 'samfile':
            log.info('Removing samfile for memory management.')
            rm_samfile = 'rm {fname}'.format(fname=samfile)
            io.ProcessManager(rm_samfile).run_all()
        else:
            if aws_upload_key:
                log.info('Converting samfile to bamfile and uploading to S3.')
                bamfile = output_dir + '/alignments/Aligned.out.bam'
                convert_sam = 'samtools view -bS -o {bamfile} {samfile}' \
                    .format(bamfile=bamfile, samfile=samfile)
                upload_bam = 'aws s3 mv {fname} {s3link}'.format(
                    fname=bamfile, s3link=aws_upload_key)
                manage_samfile = io.ProcessManager(convert_sam, upload_bam)
                manage_samfile.run_all()
    else:
        raise ValueError('Must provide sam file')
        return
#        if read_array.startswith('s3://'):
#            input_data = 'readarray'
#            bucket, prefix = io.S3.split_link(read_array)
#            _, fname = os.path.split(prefix)
#            read_array = io.S3.download_file(
#                bucket, prefix, output_dir + '/' + fname)
#            log.info('Read array %s successfully installed from S3.' % read_array)
#        ra = ReadArray.load(read_array)
    return ra, input_data, manage_samfile, sam_records, read_array
    

def generate_count_matrices(args, cell_counts, pigz, aws_upload_key):
    """generate sparse count matrix, filter, and generate dense count matrix
    :param args:
    :param cell_counts:
    :param pigz:
    :param aws_upload_key:
    :return:
    """

    log.info('Creating sparse count matrices')
    matrices = correct_errors.convert_to_matrix(cell_counts)
    molecules = SparseFrame(
        matrices['molecules']['matrix'],
        matrices['molecules']['row_ids'],
        matrices['molecules']['col_ids'])
    reads = SparseFrame(
        matrices['reads']['matrix'],
        matrices['reads']['row_ids'],
        matrices['reads']['col_ids'])

    with open(args.output_stem + '_read_and_count_matrices.p', 'wb') as f:
        pickle.dump(matrices, f)
    log.info('Successfully generated sparse count matrix.')

    log.info('filtering, summarizing cell molecule count distribution, and '
             'generating dense count matrix')

    # todo speed up this function
    gene_id_map = create_gene_id_to_official_gene_symbol_map(
        args.index + 'annotations.gtf')

    # translate symbols to gene names
    reads = ensembl_gene_id_to_official_gene_symbol(reads, gene_id_map=gene_id_map)
    molecules = ensembl_gene_id_to_official_gene_symbol(
        molecules, gene_id_map=gene_id_map)

    # get dense molecules
    dense, total_molecules, mols_lost, cells_lost, cell_description = (
        create_filtered_dense_count_matrix(
            molecules, reads, plot=True, figname=args.output_stem + '_filters.png',
            filter_mitochondrial_rna=args.filter_mitochondrial_rna))

    # save sparse csv
    dense[dense == 0] = np.nan  # for sparse_csv saving
    sparse_csv = args.output_stem + '_counts.csv'
    dense.to_csv(sparse_csv)

    # gzipping sparse_csv and uploading to S3
    if pigz:
        sparse_zip = "pigz --best {fname}".format(fname=sparse_csv)
    else:
        sparse_zip = "gzip -f {fname}".format(fname=sparse_csv)
    sparse_upload = 'aws s3 mv {fname} {s3link}'.format(fname=sparse_csv + '.gz',
                                                        s3link=aws_upload_key)
    if aws_upload_key:
        sparse_proc = io.ProcessManager(sparse_zip, sparse_upload)
    else:
        sparse_proc = io.ProcessManager(sparse_zip)
    sparse_proc.run_all()

    return (sparse_proc, sparse_csv, total_molecules, mols_lost, cells_lost,
            cell_description)


def main(argv: list=None) -> None:

    args = parser.parse_args(argv)
    if args.aws:
        args.log_name = args.log_name.split('/')[-1]
        log.setup_logger('/data/' + args.log_name)
    else:
        log.setup_logger(args.log_name)

    print('beginning execution_control')
    with execution_control.cleanup(args):
        print('in_execution_control')
        log.print_exact_command_line('process_experiment.py ' + " ".join(argv))
        pigz, email = verify.executables('pigz', 'mutt')
        if not email and not args.remote:
            log.notify('mutt was not found on this machine; an email will not be sent to '
                       'the user upon termination of SEQC run.')
        log.args(args)

        # read in config file, make sure it exists
        configuration = config.read_config('~/.seqc/config')

        # extract basespace token and make sure args.index is a directory
        basespace_token = configuration['BaseSpaceToken']['base_space_token']
        if not args.index.endswith('/'):
            args.index += '/'

        # check arguments if aws flag is not set (arguments would already be checked)
        if not args.aws:
            total_size = verify.arguments(args, basespace_token)
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
            remote.cluster_cleanup()
            run_remote(args, argv, total_size)
            sys.exit(0)

        if args.aws:
            aws_upload_key, args.output_stem, output_dir, output_prefix = (
                upload.update_directories_for_aws(args.output_stem, output_prefix))
        else:
            aws_upload_key = None

        # download data if necessary
        if args.basespace:
            args.barcode_fastq, args.genomic_fastq = download.basespace(
                args, output_dir, basespace_token)

        # check for remote fastq file links
        args.genomic_fastq = download.s3_fastq(
            args.genomic_fastq, output_dir, 'genomic')
        args.barcode_fastq = download.s3_fastq(
            args.barcode_fastq, output_dir, 'barcode')

        # check if the index must be downloaded
        args.index = download.index(
            output_dir, args.index, args.read_array, args.samfile)

        n_processes = multiprocessing.cpu_count() - 1  # get number of processors

        input_data = 'start'  # which data was used as input?

        # determine where the script should start:
        merge, align, process_samfile = determine_start_point(args)

        # check if the platform name provided is supported by seqc
        platform_name = verify.platform_name(args.platform)
        platform = platforms.AbstractPlatform.factory(platform_name)  # returns platform class

        if merge:
            if args.min_poly_t is None:  # estimate min_poly_t if it was not provided
                args.min_poly_t = filter.estimate_min_poly_t(
                    args.barcode_fastq, platform)
                log.notify('Estimated min_poly_t={!s}'.format(args.min_poly_t))

            args.merged_fastq, fastq_records = merge_fastq_files(
                platform, args.barcode_fastq, args.output_stem,
                args.genomic_fastq, pigz)
        else:
            fastq_records = None

        if align:
            args.samfile, input_data, manage_merged = align_fastq_records(
                args.merged_fastq, output_dir, args.output_stem, args.star_args,
                args.index, n_processes, aws_upload_key, pigz)
        else:
            manage_merged = None

        ra, input_data, manage_samfile, sam_records, h5_file = (
            create_or_download_read_array(
                process_samfile, args.samfile, output_dir, args.index, args.read_array,
                args.output_stem, aws_upload_key, input_data))
        args.read_array = h5_file

        args.barcode_files = download.barcodes(
            platform_name, args.barcode_files, output_dir)

        # SEQC was started from input other than fastq files
        if args.min_poly_t is None:
            args.min_poly_t = 0
            log.notify('Warning: SEQC started from step other than unmerged fastq with '
                       'empty --min-poly-t parameter. Continuing with --min-poly-t=0.')
        if args.max_dust_score is None:
            args.max_dust_score = 10
            log.notify('Warning: --max-dust-score parameter was not supplied, continuing '
                       'with --max-dust-score=10.')

        log.info('Read array after reading sam file: {}'.format(ra))
        # Apply filters
        ra.apply_filters(required_poly_t=args.min_poly_t, max_dust_score = args.max_dust_score)
        log.info('Read array after filtering: {}'.format(ra))
        # Correct barcodes
        if platform.check_barcodes:
            error_rate = ra.apply_barcode_correction(platform, args.barcode_files, reverse_complement=True, max_ed=args.max_ed)
            log.info('Read array after barcode correction: {}'.format(ra))
        else:
            error_rate = None
            log.info('Skipping barcode correction')
        
        # Resolve multimapping
        ra.resolve_alignments(args.index)
        log.info('Read array after multialignment resolution: {}'.format(ra))
        # correct errors
        log.info('Filterring errors')
        # for in-drop and mars-seq, summary is a dict. for drop-seq, it may be None
        cell_counts = platform.correct_errors(ra, error_rate, singleton_weight=args.singleton_weight)
        log.info('Read array after error correction: {}'.format(ra))
        
        files = ra.to_count_matrix(args.output_stem)
##############
        bucket, key = io.S3.split_link(aws_upload_key)
        for item in files:
            try:
                retry_boto_call(io.S3.upload_file)(item, bucket, key)
                item_name = item.split('/')[-1]
                log.info('Successfully uploaded %s to the specified S3 location '
                         '"%s%s".' % (item, aws_upload_key, item_name))
            except FileNotFoundError:
                log.notify('Item %s was not found! Continuing with upload...' % item)

        #raise StandardError('All new code finished running')
        
###################
        #TODO: I need to create summary but from the ra and not throuygh error corrcetion like before
        summary = None

        # uploading read array to S3 if created, else removing read array
        if input_data == 'readarray':
            log.info('Removing .h5 file for memory management.')
            rm_ra = 'rm {fname}'.format(fname=args.read_array)
            io.ProcessManager(rm_ra).run_all()
            manage_ra = None
        else:
            if aws_upload_key:
                log.info('Uploading read array to S3.. Link is %s' % aws_upload_key)
                upload_ra = 'aws s3 mv {fname} {s3link}'.format(fname=args.read_array,
                                                                s3link=aws_upload_key)
                manage_ra = io.ProcessManager(upload_ra)
                manage_ra.run_all()
                

        #sparse_proc, sparse_csv, total_mols, mols_lost, cells_lost, cell_descr = (
        #    generate_count_matrices(args, cell_counts, pigz, aws_upload_key))
 #       return
        
        # in this version, local runs won't be able to upload to S3
        # and also won't get an e-mail notification.
        if args.aws:
            upload.data_and_notify(
                args.email_status, aws_upload_key, align, input_data, manage_merged,
                process_samfile, args.merged_fastq, args.samfile, args.read_array,
                sparse_proc, sparse_csv, manage_ra, summary, fastq_records, sam_records,
                total_mols, mols_lost, cells_lost, manage_samfile, cell_descr,
                args.output_stem, args.log_name, email)


if __name__ == '__main__':
    main(sys.argv[1:])
