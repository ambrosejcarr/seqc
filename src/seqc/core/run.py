
def run(args) -> None:
    """Run SEQC on the files provided in args, given specifications provided on the
    command line

    :param args: parsed argv, produced by seqc.parser(). This function is only called
      when args.subprocess_name is "run".
    """

    # import inside module for pickle functionality
    # top 2 only needed for post-filtering
    # from seqc.filter import create_filtered_dense_count_matrix  # todo can remove this?
    # from seqc.sparse_frame import SparseFrame  # todo can remove this?

    import os
    import multiprocessing
    from seqc import log, ec2, platforms, filter, io
    from seqc.sequence import fastq
    from seqc.alignment import star
    from seqc.email import email_user
    from seqc.read_array import ReadArray
    from seqc.core import verify, download

    def determine_start_point(arguments) -> (bool, bool, bool):
        """
        determine where seqc should start based on which parameters were passed.

        :param arguments: Namespace object, result of ArgumentParser.parse_args()
        :returns merge, align, process_samfile: indicates whether merging, alignment, and
          processing samfiles should be executed.
        """
        if arguments.read_array:
            return False, False, False
        if arguments.samfile:
            return False, False, True
        if arguments.merged_fastq:
            return False, True, True
        else:
            return True, True, True

    def download_input(dir_, arguments):
        """parse input arguments and download any necessary data

        :param str dir_: directory to download data to
        :param arguments: namespace object from argparse
        :return args: updated namespace object reflecting local file paths of downloaded
          files
        """
        # download basespace data if necessary
        if arguments.basespace:
            arguments.barcode_fastq, arguments.genomic_fastq = io.BaseSpace.download(
                arguments.platform, arguments.basespace, dir_, arguments.basespace_token)

        # check for remote fastq file links
        arguments.genomic_fastq = download.s3_data(
            arguments.genomic_fastq, dir_ + '/genomic_fastq/')
        arguments.barcode_fastq = download.s3_data(
            arguments.barcode_fastq, dir_ + '/barcode_fastq/')

        # get merged fastq file, unzip if necessary
        arguments.merged_fastq = (
            download.s3_data([arguments.merged_fastq], dir_ + '/')[0] if
            arguments.merged_fastq is not None else None)
        if arguments.merged_fastq and arguments.merged_fastq.endswith('.gz'):
            if pigz:
                unzip_command = 'pigz -d -f %s' % arguments.merged_fastq
            else:
                unzip_command = 'gunzip -f %s' % arguments.merged_fastq
            gunzip_proc = io.ProcessManager(unzip_command)
            gunzip_proc.run_all()
            gunzip_proc.wait_until_complete()
            arguments.merged_fastq = arguments.merged_fastq.replace('.gz', '')
            log.info(
                'Merged fastq file %s successfully installed from S3 and unzipped.' %
                arguments.merged_fastq)

        # check if the index must be downloaded
        if any((arguments.samfile, arguments.read_array)):
            index_link = arguments.index + 'annotations.gtf'
        else:
            index_link = arguments.index
        download.s3_data([index_link], dir_ + '/index/')
        arguments.index = dir_ + '/index/'

        # check if barcode files must be downloaded
        arguments.barcode_files = download.s3_data(
            arguments.barcode_files, dir_ + '/barcodes/')

        # check if samfile needs downloading
        if arguments.samfile:
            arguments.samfile = download.s3_data([arguments.samfile], dir_ + '/')[0]

        # check if readarray needs downloading
        if arguments.read_array:
            arguments.read_array = download.s3_data([arguments.read_array], dir_ + '/')[0]

        return arguments

    def merge_fastq_files(
            technology_platform, barcode_fastq: [str], output_stem: str,
            genomic_fastq: [str]) -> (str, int):
        """annotates genomic fastq with barcode information; merging the two files.

        :param technology_platform: class from platforms.py that defines the
          characteristics of the data being processed
        :param barcode_fastq: list of str names of fastq files containing barcode
          information
        :param output_stem: str, stem for output files
        :param genomic_fastq: list of str names of fastq files containing genomic
          information
        :returns merged_fastq, fastq_records: (str, int) name of merged fastq file and the
          number of fastq records that were processed.
        """

        log.info('Merging genomic reads and barcode annotations.')
        merged_fastq, fastq_records = fastq.merge_paired(
            merge_function=technology_platform.merge_function,
            fout=output_stem + '_merged.fastq',
            genomic=genomic_fastq,
            barcode=barcode_fastq)

        # delete genomic/barcode fastq files after merged.fastq creation
        log.info('Removing original fastq file for memory management.')
        delete_fastq = ' '.join(['rm'] + genomic_fastq + barcode_fastq)
        io.ProcessManager(delete_fastq).run_all()

        return merged_fastq, fastq_records

    def align_fastq_records(
            merged_fastq, dir_, star_args, star_index, n_proc,
            aws_upload_key) -> (str, str, io.ProcessManager):
        """
        Align fastq records.

        :param merged_fastq: str, path to merged .fastq file
        :param dir_: str, stem for output files
        :param star_args: dict, extra keyword arguments for STAR
        :param star_index: str, file path to directory containing STAR index
        :param n_proc: int, number of STAR processes to initiate
        :param aws_upload_key: str, location to upload files, or None if seqc was
          initiated from a merged fastq file.
        :return samfile, input_data, upload_manager: (str, str, io.ProcessManager)
          name of .sam file containing aligned reads, indicator of which data was used as
          input, and a ProcessManager for merged fastq files
        """
        log.info('Aligning merged fastq records.')
        alignment_directory = dir_ + '/alignments/'
        os.makedirs(alignment_directory, exist_ok=True)
        if star_args is not None:
            star_kwargs = dict(a.strip().split('=') for a in star_args)
        else:
            star_kwargs = {}
        samfile = star.align(
            merged_fastq, star_index, n_proc, alignment_directory,
            **star_kwargs)

        if aws_upload_key:
            log.info('Gzipping merged fastq file.')
            if pigz:
                pigz_zip = "pigz --best -k -f {fname}".format(fname=merged_fastq)
            else:
                pigz_zip = "gzip -kf {fname}".format(fname=merged_fastq)
            pigz_proc = io.ProcessManager(pigz_zip)
            pigz_proc.run_all()
            pigz_proc.wait_until_complete()  # prevents slowing down STAR alignment
            merged_fastq += '.gz'  # reflect gzipped nature of file

            log.info('Uploading gzipped merged fastq file to S3.')
            merge_upload = 'aws s3 mv {fname} {s3link}'.format(
                fname=merged_fastq, s3link=aws_upload_key)
            upload_manager = io.ProcessManager(merge_upload)
            upload_manager.run_all()
        else:
            log.info('Removing merged fastq file for memory management.')
            rm_merged = 'rm %s' % merged_fastq
            io.ProcessManager(rm_merged).run_all()

            upload_manager = None
        return samfile, upload_manager

    def create_read_array(
            samfile, dir_, index, aws_upload_key) -> (str, str, io.ProcessManager, int):
        """Create or download a ReadArray object.

        :param samfile: str, filename of .sam file
        :param dir_: str, directory for output files
        :param index: str, directory containing index files
        :param aws_upload_key: str, key where aws files should be uploaded
        :returns read_array, input_data, upload_manager: ReadArray object,
          samfile ProcessManager, and number of processed sam records
        """
        log.info('Filtering aligned records and constructing record database.')
        read_array, sam_record_count = ReadArray.from_samfile(
            samfile, index + 'annotations.gtf')

        # TODO save the created ReadArray object

        # converting sam to bam and uploading to S3, else removing samfile
        if aws_upload_key:
            log.info('Converting samfile to bamfile and uploading to S3.')
            bamfile = dir_ + '/alignments/Aligned.out.bam'
            convert_sam = 'samtools view -bS -o {bamfile} {samfile}' \
                .format(bamfile=bamfile, samfile=samfile)
            upload_bam = 'aws s3 mv {fname} {s3link}'.format(
                fname=bamfile, s3link=aws_upload_key)
            upload_manager = io.ProcessManager(convert_sam, upload_bam)
            upload_manager.run_all()
        else:
            log.info('Removing samfile for memory management.')
            rm_samfile = 'rm %s' % samfile
            io.ProcessManager(rm_samfile).run_all()
            upload_manager = None
        return read_array, upload_manager, sam_record_count

    # ######################## MAIN FUNCTION BEGINS HERE ################################

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

        # check if the platform name provided is supported by seqc
        # todo move into verify for run
        platform_name = verify.platform_name(args.platform)
        platform = platforms.AbstractPlatform.factory(platform_name)  # returns platform

        n_processes = multiprocessing.cpu_count() - 1  # get number of processors

        merge, align, process_samfile = determine_start_point(args)

        args = download_input(output_dir, args)

        if merge:
            if args.min_poly_t is None:  # estimate min_poly_t if it was not provided
                args.min_poly_t = filter.estimate_min_poly_t(
                    args.barcode_fastq, platform)
                log.notify('Estimated min_poly_t={!s}'.format(args.min_poly_t))

            args.merged_fastq, input_fastq_records = merge_fastq_files(
                platform, args.barcode_fastq, args.output_prefix,
                args.genomic_fastq)
        else:
            input_fastq_records = None

        if align:
            upload_merged = args.upload_prefix if merge else None
            args.samfile, manage_merged = align_fastq_records(
                args.merged_fastq, output_dir, args.star_args,
                args.index, n_processes, upload_merged)
        else:
            manage_merged = None

        if process_samfile:
            upload_samfile = args.upload_prefix if align else None
            ra, manage_samfile, sam_records = create_read_array(
                args.samfile, output_dir, args.index, upload_samfile)
        else:
            manage_samfile, sam_records = None, None
            ra = ReadArray.load(args.read_array)

        # SEQC was started from input other than fastq files
        if args.min_poly_t is None:
            args.min_poly_t = 0
            log.notify('Warning: SEQC started from step other than unmerged fastq with '
                       'empty --min-poly-t parameter. Continuing with --min-poly-t=0.')

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
            error_rate = ra.apply_barcode_correction(platform, args.barcode_files,
                                                     max_ed=args.max_ed)
            files += ra.to_count_matrix(args. output_prefix + '_phase3_')
            log.info('Read array after barcode correction: {}'.format(ra))
        else:
            error_rate = None
            log.info('Skipping barcode correction')

        # Resolve multimapping
        ra.resolve_alignments(args.index)
        files += ra.to_count_matrix(args.output_prefix + '_phase4_')
        log.info('Read array after multialignment resolution: {}'.format(ra))
        # correct errors
        log.info('Filtering errors')
        # for in-drop and mars-seq, summary is a dict. for drop-seq, it may be None
        platform.correct_errors(ra, error_rate, singleton_weight=args.singleton_weight)
        files += ra.to_count_matrix(args.output_prefix + '_phase5_')
        log.info('Read array after error correction: {}'.format(ra))

        # filter non-cells
        sp_reads, sp_mols = ra.to_count_matrix(
            sparse_frame=True, genes_to_symbols=args.index + 'annotations.gtf')

        ra.save(args.output_prefix + '.h5')
        log.notify('ReadArray saved.')

        # todo fold all this output into a summary page
        cell_filter_figure = 'cell_filters.png'
        sp_csv, total_molecules, molecules_lost, cells_lost, cell_description = (
            filter.create_filtered_dense_count_matrix(
                sp_mols, sp_reads, plot=True, figname=cell_filter_figure))
        dense_csv = args.output_prefix + '_dense.csv'
        sp_csv.to_csv(dense_csv)

        # get alignment summary
        if os.path.isfile(output_dir + '/alignments/Log.final.out'):
            os.rename(output_dir + '/alignments/Log.final.out',
                      output_dir + '/alignment_summary.txt')
            files += [output_dir + '/alignment_summary.txt']

        files += [dense_csv, cell_filter_figure, args.output_prefix + '.h5']

        if args.upload_prefix:
            # Upload count matrices files, logs, and return
            bucket, key = io.S3.split_link(args.upload_prefix)
            for item in files:
                try:
                    ec2.Retry(retries=5)(io.S3.upload_file)(item, bucket, key)
                    item_name = item.split('/')[-1]
                    log.info('Successfully uploaded %s to the specified S3 location '
                             '"%s%s".' % (item, args.upload_prefix, item_name))
                except FileNotFoundError:
                    log.notify('Item %s was not found! Continuing with upload...' % item)

            # todo read-array does not exist, no need to remove. waste of time
            # # uploading read array to S3 if created, else removing read array
            # if not process_samfile:
            #     log.info('Removing .h5 file for memory management.')
            #     rm_ra = 'rm {fname}'.format(fname=args.read_array)
            #     io.ProcessManager(rm_ra).run_all()

        if manage_merged:
            manage_merged.wait_until_complete()
            log.info('Successfully uploaded %s to the specified S3 location "%s"' %
                     (args.merged_fastq, args.upload_prefix))
        if manage_samfile:
            manage_samfile.wait_until_complete()
            log.info('Successfully uploaded %s to the specified S3 location "%s"'
                     % (args.samfile, args.upload_prefix))

        # todo we have no read_array now
        # if manage_ra:
        #     manage_ra.wait_until_complete()
        #     log.info('Successfully uploaded %s to the specified S3 location "%s"' %
        #              (args.read_array, args.aws_upload_key))

        log.info('SEQC run complete. Cluster will be terminated')

        # upload logs
        if args.upload_prefix:
            # Upload count matrices files, logs, and return
            bucket, key = io.S3.split_link(args.upload_prefix)
            for item in [args.log_name, './nohup.log']:
                try:
                    ec2.Retry(retries=5)(io.S3.upload_file)(item, bucket, key)
                    item_name = item.split('/')[-1]
                    log.info('Successfully uploaded %s to the specified S3 location '
                             '"%s%s".' % (item, args.upload_prefix, item_name))
                except FileNotFoundError:
                    log.notify('Item %s was not found! Continuing with upload...' % item)

        # TODO: I need to create summary but from the read_array instead of via error
        # correction
        # TODO: generate a printout of these for the summary
        if mutt:
            email_body = (
                '<font face="Courier New, Courier, monospace">'
                'SEQC RUN COMPLETE.\n\n'
                'The run log has been attached to this email and '
                'results are now available in the S3 location you specified: '
                '"%s"\n\n' % args.upload_prefix)
            email_body = email_body.replace('\n', '<br>').replace('\t', '&emsp;')
            email_user(args.log_name, email_body, args.email)
