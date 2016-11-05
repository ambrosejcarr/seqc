#!/usr/local/bin/python3

import sys
from seqc.core import parser, verify
from seqc import ec2


def run(args) -> None:
    """Run SEQC on the files provided in args, given specifications provided on the
    command line

    :param args: parsed argv, produced by seqc.parser(). This function is only called
      when args.subprocess_name is "run".
    """

    # functions to be run remotely must import all required arguments
    import os
    import pickle
    import multiprocessing
    import numpy as np
    from seqc import log, ec2, platforms, filter, io, correct_errors
    from seqc.sequence import fastq
    from seqc.alignment import star
    from seqc.filter import create_filtered_dense_count_matrix
    from seqc.read_array import ReadArray
    from seqc.sparse_frame import SparseFrame
    from seqc.core import verify, download, upload

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
            pigz_zip = "pigz --best -k -f {fname}".format(fname=merged_fastq)
        else:
            pigz_zip = "gzip -kf {fname}".format(fname=merged_fastq)
        pigz_proc = io.ProcessManager(pigz_zip)
        pigz_proc.run_all()
        pigz_proc.wait_until_complete()  # prevents slowing down STAR alignment
        return merged_fastq, fastq_records

    def align_fastq_records(
            merged_fastq, output_dir, star_args, index, n_processes,
            aws_upload_key, pigz) -> (str, str, io.ProcessManager):
        """
        Align fastq records.

        :param merged_fastq: str, path to merged .fastq file
        :param output_dir: str, stem for output files
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
            # bucket, prefix = io.S3.split_link(merged_fastq)
            # _, fname = os.path.split(prefix)
            merged_fastq = io.S3.download(merged_fastq, output_dir + '/')
            if merged_fastq.endswith('.gz'):
                if pigz:
                    cmd = 'pigz -d -f {fname}'.format(fname=merged_fastq)
                else:
                    cmd = 'gunzip -f {fname}'.format(fname=merged_fastq)
                gunzip_proc = io.ProcessManager(cmd)
                gunzip_proc.run_all()
                gunzip_proc.wait_until_complete()  # finish before STAR alignment
                merged_fastq = merged_fastq.replace('.gz', '')
            log.info(
                'Merged fastq file %s successfully installed from S3.' % merged_fastq)
        alignment_directory = output_dir + '/alignments/'
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
            rm_cmd = 'rm %s' % merged_fastq
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
            log.info('Filtering aligned records and constructing record database.')
            if samfile.startswith('s3://'):
                samfile = io.S3.download(samfile, output_dir)
                log.info('Samfile %s successfully installed from S3.' % samfile)
            ra, sam_records = ReadArray.from_samfile(
                samfile, index + 'annotations.gtf')
            read_array = output_stem + '.h5'
            ra.save(read_array)

            # converting sam to bam and uploading to S3, else removing samfile
            if input_data == 'samfile':
                log.info('Removing samfile for memory management.')
                rm_samfile = 'rm %s' % samfile
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
            if read_array.startswith('s3://'):
                input_data = 'readarray'
                read_array = io.S3.download(read_array, output_dir)
                log.info('Read array %s successfully installed from S3.' % read_array)
            ra = ReadArray.load(read_array)
        return ra, input_data, manage_samfile, sam_records, read_array

    def generate_count_matrices(args, cell_counts, pigz, aws_upload_key):
        """generate sparse count matrix, filter, and generate dense count matrix
        :param args:  parsed command-line arguments
        :param dict cell_counts: error correction results
        :param bool pigz: True if pigz is installed, else False
        :param aws_upload_key: location on aws to upload results
        :return tuple:
          ProcessManager sparse_proc: manager for uploading SparseMatrix results
          str sparse_csv: filtered count matrix, saved in sparse csv format
          int total_molecules: number of molecules recovered pre-filtering
          int mols_lost: number of molecules removed by filters
          int cells_lost: number of cells removed by filters
          pd.Series cell_description: descriptive statistics of the retained data
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

        with open(args.output_prefix + '_read_and_count_matrices.p', 'wb') as f:
            pickle.dump(matrices, f)
        log.info('Successfully generated sparse count matrix.')

        log.info('filtering, summarizing cell molecule count distribution, and '
                 'generating dense count matrix')

        # get dense molecules
        dense, total_molecules, mols_lost, cells_lost, cell_description = (
            create_filtered_dense_count_matrix(
                molecules, reads, plot=True, figname=args.output_prefix + '_filters.png',
                filter_mitochondrial_rna=False))  # todo fix hardcoding False for filter

        # save sparse csv
        dense[dense == 0] = np.nan  # for sparse_csv saving
        sparse_csv = args.output_prefix + '_counts.csv'
        dense.to_csv(sparse_csv)

        # gzipping sparse_csv and uploading to S3
        if pigz:
            sparse_zip = "pigz --best -f {fname}".format(fname=sparse_csv)
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

    ######################### MAIN FUNCTION BEGINS HERE ################################

    log.setup_logger(args.log_name)
    with ec2.instance_clean_up(
            email=args.email, upload=args.upload_prefix, log_name=args.log_name):
        pigz, mutt = verify.executables('pigz', 'mutt')
        if mutt:
            log.notify('mutt executable identified, email will be sent when run '
                       'terminates. ')
        else:
            log.notify('mutt was not found on this machine; an email will not be sent to '
                       'the user upon termination of SEQC run.')
        log.args(args)

        output_dir, output_prefix = os.path.split(args.output_prefix)
        if not output_dir:
            output_dir = '.'

        # download basespace data if necessary
        if args.basespace:
            args.barcode_fastq, args.genomic_fastq = download.basespace(
                args, output_dir, args.basespace_token)

        # check for remote fastq file links
        args.genomic_fastq = download.s3_data(
            args.genomic_fastq, output_dir + '/genomic_fastq/')
        args.barcode_fastq = download.s3_data(
            args.barcode_fastq, output_dir + '/barcode_fastq/')

        # check if the index must be downloaded
        print(args.index)
        if any((args.samfile, args.read_array)):
            index_link = args.index + 'annotations.gtf'
        else:
            index_link = args.index
        download.s3_data([index_link], output_dir + '/index/')
        args.index = output_dir + '/index/'
        print(args.index)

        # check if barcode files must be downloaded
        args.barcode_files = download.s3_data(
            args.barcode_files, output_dir + '/barcodes/')

        n_processes = multiprocessing.cpu_count() - 1  # get number of processors

        # determine where the script should start:
        input_data = 'start'  # which data was used as input?
        merge, align, process_samfile = determine_start_point(args)

        # check if the platform name provided is supported by seqc
        platform_name = verify.platform_name(args.platform)
        platform = platforms.AbstractPlatform.factory(platform_name)  # returns platform

        if merge:
            if args.min_poly_t is None:  # estimate min_poly_t if it was not provided
                args.min_poly_t = filter.estimate_min_poly_t(
                    args.barcode_fastq, platform)
                log.notify('Estimated min_poly_t={!s}'.format(args.min_poly_t))

            args.merged_fastq, fastq_records = merge_fastq_files(
                platform, args.barcode_fastq, args.output_prefix,
                args.genomic_fastq, pigz)
        else:
            fastq_records = None

        if align:
            args.samfile, input_data, manage_merged = align_fastq_records(
                args.merged_fastq, output_dir, args.star_args,
                args.index, n_processes, args.upload_prefix, pigz)
        else:
            manage_merged = None

        ra, input_data, manage_samfile, sam_records, h5_file = (
            create_or_download_read_array(
                process_samfile, args.samfile, output_dir, args.index, args.read_array,
                args.output_prefix, args.upload_prefix, input_data))
        args.read_array = h5_file

        # SEQC was started from input other than fastq files
        if args.min_poly_t is None:
            args.min_poly_t = 0
            log.notify('Warning: SEQC started from step other than unmerged fastq with '
                       'empty --min-poly-t parameter. Continuing with --min-poly-t=0.')
        if args.max_dust_score is None:
            args.max_dust_score = 10
            log.notify('Warning: --max-dust-score parameter was not supplied, continuing '
                       'with --max-dust-score=10.')

        files = []
        files += ra.to_count_matrix(args.output_prefix + '_phase1_')
        log.info('Read array after reading sam file: {}'.format(ra))
        # Apply filters
        ra.apply_filters(
            required_poly_t=args.min_poly_t, max_dust_score=args.max_dust_score)
        files += ra.to_count_matrix(args.output_prefix + '_phase2_')
        log.info('Read array after filtering: {}'.format(ra))
        # Correct barcodes
        if platform.check_barcodes:
            error_rate = ra.apply_barcode_correction(
                platform, args.barcode_files, reverse_complement=False,
                max_ed=args.max_ed)
            files += ra.to_count_matrix(args.output_prefix+'_phase3_')
            log.info('Read array after barcode correction: {}'.format(ra))
        else:
            error_rate = None
            log.info('Skipping barcode correction')
        
        # Resolve multimapping
        ra.resolve_alignments(args.index)
        files += ra.to_count_matrix(args.output_prefix + '_phase4_')
        log.info('Read array after multialignment resolution: {}'.format(ra))
        # correct errors
        log.info('Filterring errors')
        # for in-drop and mars-seq, summary is a dict. for drop-seq, it may be None
        platform.correct_errors(ra, error_rate, singleton_weight=args.singleton_weight)
        files += ra.to_count_matrix(args.output_prefix + '_phase5_')
        log.info('Read array after error correction: {}'.format(ra))
        
        if args.upload_prefix:
            # Upload count matrices files
            bucket, key = io.S3.split_link(args.upload_prefix)
            for item in files:
                try:
                    ec2.Retry(retries=5)(io.S3.upload_file)(item, bucket, key)
                    item_name = item.split('/')[-1]
                    log.info('Successfully uploaded %s to the specified S3 location '
                             '"%s%s".' % (item, args.upload_prefix, item_name))
                except FileNotFoundError:
                    log.notify('Item %s was not found! Continuing with upload...' % item)

            # uploading read array to S3 if created, else removing read array
            if input_data == 'readarray':
                log.info('Removing .h5 file for memory management.')
                rm_ra = 'rm {fname}'.format(fname=args.read_array)
                io.ProcessManager(rm_ra).run_all()

# This is commented out until we have a better way to save the ReadArray
#        else:
#            if aws_upload_key:
#                log.info('Uploading read array to S3.. Link is %s' % aws_upload_key)
#                upload_ra = 'aws s3 mv {fname} {s3link}'.format(fname=args.read_array,
#                                                                s3link=aws_upload_key)
#                manage_ra = io.ProcessManager(upload_ra)
#                manage_ra.run_all()


        #sparse_proc, sparse_csv, total_mols, mols_lost, cells_lost, cell_descr = (
        #    generate_count_matrices(args, cell_counts, pigz, aws_upload_key))
 #       return

        #TODO: I need to create summary but from the ra and not throuygh error corrcetion like before
        #TODO: generate a printout of these for the summary
        sparse_proc = None
        sparse_csv = 0
        total_mols = 0
        mols_lost = 0
        cells_lost = 0
        cell_descr = None
        summary = None
        manage_ra = None

        if args.upload_prefix:
            upload.data_and_notify(
                args.email, args.upload_prefix, align, input_data, manage_merged,
                process_samfile, args.merged_fastq, args.samfile, args.read_array,
                sparse_proc, sparse_csv, manage_ra, summary, fastq_records, sam_records,
                total_mols, mols_lost, cells_lost, manage_samfile, cell_descr,
                args.output_prefix, args.log_name, mutt)


def index(args):  # todo needs a volume attached, how to specify amount?
    """create an index for SEQC.

    :param args: parsed arguments. This function is only called if subprocess_name is
      'index'
    """

    # functions to be pickled and run remotely must import all their own modules
    from seqc import ec2, log
    from seqc.sequence.index import Index

    log.setup_logger(args.log_name)
    with ec2.instance_clean_up(args.email, args.upload, log_name=args.log_name):
        idx = Index(args.organism, args.additional_id_types)
        idx.create_index(args.upload_location)


def progress(args):
    """check the progress of a seqc experiment

    :param args:
    :return:
    """
    raise NotImplementedError


def main(argv):
    """Check arguments, then call the appropriate sub-module

    Created to allow the main pipeline to be tested from the earliest entry point
    (command-line arguments).

    :param argv: output of sys.argv[1:]
    """
    arguments = parser.parse_args(argv)
    this_module = sys.modules[__name__]
    func = getattr(this_module, arguments.subparser_name)
    assert func is not None
    verification_func = getattr(verify, arguments.subparser_name)
    verified_args = verification_func(arguments)
    if arguments.remote:
        remote_args = {
            k: getattr(verified_args, k) for k in
            ('rsa_key', 'instance_type', 'spot_bid', 'volume_size') if
            getattr(verified_args, k)}
        ec2.AWSInstance(synchronous=False, **remote_args)(func)(verified_args)
    else:
        func(verified_args)


if __name__ == '__main__':
    main(sys.argv[1:])
