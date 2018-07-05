

def run(args) -> None:
    """Run SEQC on the files provided in args, given specifications provided on the
    command line

    :param args: parsed argv, produced by seqc.parser(). This function is only called
      when args.subprocess_name is "run".
    """

    # import inside module for pickle functionality
    # top 2 only needed for post-filtering

    import os
    import multiprocessing
    from seqc import log, ec2, platforms, io
    from seqc.sequence import fastq
    from seqc.alignment import star
    from seqc.email_ import email_user
    from seqc.read_array import ReadArray
    from seqc.core import verify, download
    from seqc import filter
    from seqc.sequence.gtf import GeneIntervals
    from seqc.summary.summary import Section, Summary
    import numpy as np
    import scipy.io
    from shutil import copyfile
    from seqc.summary.summary import MiniSummary
    from seqc.stats.mast import run_mast
    import logging
    logger = logging.getLogger('weasyprint')
    logger.handlers = []  # Remove the default stderr handler
    logger.setLevel(100)
    logger.addHandler(logging.FileHandler('weasyprint.log'))

    def determine_start_point(arguments) -> (bool, bool, bool):
        """
        determine where seqc should start based on which parameters were passed.

        :param arguments: Namespace object, result of ArgumentParser.parse_args()
        :returns merge, align, process_bamfile: indicates whether merging, alignment, and
          processing bamfiles should be executed.
        """
        if arguments.read_array:
            return False, False, False
        if arguments.alignment_file:
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

        # check if the index must be downloaded
        if any((arguments.alignment_file, arguments.read_array)):
            index_link = arguments.index + 'annotations.gtf'
        else:
            index_link = arguments.index
        download.s3_data([index_link], dir_ + '/index/')
        arguments.index = dir_ + '/index/'

        # check if barcode files must be downloaded
        arguments.barcode_files = download.s3_data(
            arguments.barcode_files, dir_ + '/barcodes/')

        # check if alignment_file needs downloading
        if arguments.alignment_file:
            arguments.alignment_file = download.s3_data(
                [arguments.alignment_file], dir_ + '/')[0]

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
        :returns str merged_fastq: name of merged fastq file
        """

        log.info('Merging genomic reads and barcode annotations.')
        merged_fastq = fastq.merge_paired(
            merge_function=technology_platform.merge_function,
            fout=output_stem + '_merged.fastq',
            genomic=genomic_fastq,
            barcode=barcode_fastq)

        # delete genomic/barcode fastq files after merged.fastq creation
        log.info('Removing original fastq file for memory management.')
        delete_fastq = ' '.join(['rm'] + genomic_fastq + barcode_fastq)
        io.ProcessManager(delete_fastq).run_all()

        return merged_fastq

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
        :return bamfile, input_data, upload_manager: (str, str, io.ProcessManager)
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
        bamfile = star.align(
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
        return bamfile, upload_manager

    def create_read_array(bamfile, index, aws_upload_key, min_poly_t,
                          max_transcript_length):
        """Create or download a ReadArray object.

        :param max_transcript_length:
        :param str bamfile: filename of .bam file
        :param str index: directory containing index files
        :param str aws_upload_key: key where aws files should be uploaded
        :param int min_poly_t: minimum number of poly_t nucleotides for a read to be valid
        :returns ReadArray, UploadManager: ReadArray object, bamfile ProcessManager
        """
        log.info('Filtering aligned records and constructing record database.')
        # Construct translator
        translator = GeneIntervals(
            index + 'annotations.gtf', max_transcript_length=max_transcript_length)
        read_array = ReadArray.from_alignment_file(
            bamfile, translator, min_poly_t)

        # converting sam to bam and uploading to S3, else removing bamfile
        if aws_upload_key:
            log.info('Uploading bam file to S3.')
            upload_bam = 'aws s3 mv {fname} {s3link}{prefix}_Aligned.out.bam'.format(
                fname=bamfile, s3link=aws_upload_key, prefix=args.output_prefix)
            print(upload_bam)
            upload_manager = io.ProcessManager(upload_bam)
            upload_manager.run_all()
        else:
            log.info('Removing bamfile for memory management.')
            rm_bamfile = 'rm %s' % bamfile
            io.ProcessManager(rm_bamfile).run_all()
            upload_manager = None
        return read_array, upload_manager

    # ######################## MAIN FUNCTION BEGINS HERE ################################

    log.setup_logger(args.log_name)

    with ec2.instance_clean_up(
            email=args.email, upload=args.upload_prefix, log_name=args.log_name,
            debug=args.debug):
        pigz, mutt = verify.executables('pigz', 'mutt')
        if mutt:
            log.notify('mutt executable identified, email will be sent when run '
                       'terminates. ')
        else:
            log.notify('mutt was not found on this machine; an email will not be sent to '
                       'the user upon termination of SEQC run.')

        max_insert_size = args.max_insert_size
        if (args.platform == "ten_x") or (args.platform == "ten_x_v2"):
            max_insert_size = 10000
            log.notify("Full length transcripts are used for read mapping in 10x data.")
            args.filter_low_coverage = False

        log.args(args)

        output_dir, output_prefix = os.path.split(args.output_prefix)
        if not output_dir:
            output_dir = '.'

        # check if the platform name provided is supported by seqc
        # todo move into verify for run
        platform_name = verify.platform_name(args.platform)
        platform = platforms.AbstractPlatform.factory(platform_name)  # returns platform

        n_processes = multiprocessing.cpu_count() - 1  # get number of processors

        merge, align, process_bamfile = determine_start_point(args)

        args = download_input(output_dir, args)

        if args.platform == "in_drop_v5":
            platform = platform.build_cb2_barcodes(args.barcode_files)
            log.notify("Built cb2 barcode hash for v5 barcodes.")

        if merge:
            if args.min_poly_t is None:  # estimate min_poly_t if it was not provided
                args.min_poly_t = filter.estimate_min_poly_t(
                    args.barcode_fastq, platform)
                log.notify('Estimated min_poly_t={!s}'.format(args.min_poly_t))

            args.merged_fastq = merge_fastq_files(
                platform, args.barcode_fastq, args.output_prefix, args.genomic_fastq)

        # SEQC was started from input other than fastq files
        if args.min_poly_t is None:
            args.min_poly_t = 0
            log.notify('Warning: SEQC started from step other than unmerged fastq with '
                       'empty --min-poly-t parameter. Continuing with --min-poly-t 0.')

        if align:
            upload_merged = args.upload_prefix if merge else None
            args.alignment_file, manage_merged = align_fastq_records(
                args.merged_fastq, output_dir, args.star_args,
                args.index, n_processes, upload_merged)
        else:
            manage_merged = None

        if process_bamfile:
            upload_bamfile = args.upload_prefix if align else None

            ra, manage_bamfile, = create_read_array(
                args.alignment_file, args.index, upload_bamfile, args.min_poly_t,
                max_insert_size)

        else:
            manage_bamfile = None
            ra = ReadArray.load(args.read_array)

        # create the first summary section here
        status_filters_section = Section.from_status_filters(ra, 'initial_filtering.html')
        sections = [status_filters_section]

        # Skip over the corrections if read array is specified by the user
        if not args.read_array:

            # Correct barcodes
            log.info('Correcting barcodes and estimating error rates.')
            error_rate = platform.apply_barcode_correction(ra, args.barcode_files)

            # Resolve multimapping
            log.info('Resolving ambiguous alignments.')
            mm_results = ra.resolve_ambiguous_alignments()

            # correct errors
            log.info('Identifying RMT errors.')
            platform.apply_rmt_correction(ra, error_rate)

            # Apply low coverage filter
            if platform.filter_lonely_triplets:
                log.info('Filtering lonely triplet reads')
                ra.filter_low_coverage(alpha=args.low_coverage_alpha)

            log.info('Saving read array.')
            ra.save(args.output_prefix + '.h5')

            # Summary sections
            # create the sections for the summary object
            sections += [
                Section.from_cell_barcode_correction(ra, 'cell_barcode_correction.html'),
                Section.from_rmt_correction(ra, 'rmt_correction.html'),
                Section.from_resolve_multiple_alignments(mm_results, 'multialignment.html')]

        # create a dictionary to store output parameters
        mini_summary_d = dict()

        # filter non-cells
        log.info('Creating counts matrix.')
        sp_reads, sp_mols = ra.to_count_matrix(
            sparse_frame=True, genes_to_symbols=args.index + 'annotations.gtf')

        # Save sparse matrices
        log.info('Saving sparse matrices')
        scipy.io.mmwrite(args.output_prefix + '_sparse_read_counts.mtx', sp_reads.data)
        scipy.io.mmwrite(args.output_prefix + '_sparse_molecule_counts.mtx', sp_mols.data)
        # Indices
        df = np.array([np.arange(sp_reads.shape[0]), sp_reads.index]).T
        np.savetxt(
            args.output_prefix + '_sparse_counts_barcodes.csv', df,
            fmt='%d', delimiter=',')
        # Columns
        df = np.array([np.arange(sp_reads.shape[1]), sp_reads.columns]).T
        np.savetxt(
            args.output_prefix + '_sparse_counts_genes.csv', df,
            fmt='%s', delimiter=',')

        log.info('Creating filtered counts matrix.')
        cell_filter_figure = args.output_prefix +  '_cell_filters.png'

        # By pass low count filter for mars seq
        sp_csv, total_molecules, molecules_lost, cells_lost, cell_description = (
            filter.create_filtered_dense_count_matrix(
                sp_mols, sp_reads, mini_summary_d, plot=True, figname=cell_filter_figure,
                filter_low_count=platform.filter_low_count,
                filter_mitochondrial_rna=args.filter_mitochondrial_rna,
                filter_low_coverage=args.filter_low_coverage,
                filter_low_gene_abundance=args.filter_low_gene_abundance))

        # Output files
        files = [cell_filter_figure,
                 args.output_prefix + '.h5',
                 args.output_prefix + '_sparse_read_counts.mtx',
                 args.output_prefix + '_sparse_molecule_counts.mtx',
                 args.output_prefix + '_sparse_counts_barcodes.csv',
                 args.output_prefix + '_sparse_counts_genes.csv']

        # Summary sections
        # create the sections for the summary object
        sections += [
            Section.from_cell_filtering(cell_filter_figure, 'cell_filtering.html'),
            Section.from_run_time(args.log_name, 'seqc_log.html')]

        # get alignment summary
        if os.path.isfile(output_dir + '/alignments/Log.final.out'):
            os.rename(output_dir + '/alignments/Log.final.out',
                      output_dir + '/' + args.output_prefix + '_alignment_summary.txt')

            # Upload files and summary sections
            files += [output_dir + '/' + args.output_prefix + '_alignment_summary.txt']
            sections.insert(
                0, Section.from_alignment_summary(
                    output_dir + '/' + args.output_prefix + '_alignment_summary.txt',
                    'alignment_summary.html'))

        cell_size_figure = 'cell_size_distribution.png'
        index_section = Section.from_final_matrix(
            sp_csv, cell_size_figure, 'cell_distribution.html')
        seqc_summary = Summary(
            output_dir + '/' + args.output_prefix + '_summary', sections, index_section)
        seqc_summary.prepare_archive()
        seqc_summary.import_image(cell_filter_figure)
        seqc_summary.import_image(cell_size_figure)
        seqc_summary.render()
        summary_archive = seqc_summary.compress_archive()
        files += [summary_archive]

        # Create a mini summary section
        alignment_summary_file = output_dir + '/' + args.output_prefix + '_alignment_summary.txt'
        seqc_mini_summary = MiniSummary(
            args.output_prefix, mini_summary_d, alignment_summary_file, cell_filter_figure,
            cell_size_figure)
        seqc_mini_summary.compute_summary_fields(ra, sp_csv)
        seqc_mini_summary_json, seqc_mini_summary_pdf = seqc_mini_summary.render()
        files += [seqc_mini_summary_json, seqc_mini_summary_pdf]

        # Running MAST for differential analysis
        # file storing the list of differentially expressed genes for each cluster
        de_gene_list_file = run_mast(
            seqc_mini_summary.get_counts_filtered(), seqc_mini_summary.get_clustering_result(),
            args.output_prefix)
        files += [de_gene_list_file]

        # adding the cluster column and write down gene-cell count matrix
        dense_csv = args.output_prefix + '_dense.csv'
        sp_csv.insert(loc=0, column='CLUSTER', value=seqc_mini_summary.get_clustering_result())
        sp_csv.to_csv(dense_csv)
        files += [dense_csv]

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

        if manage_merged:
            manage_merged.wait_until_complete()
            log.info('Successfully uploaded %s to the specified S3 location "%s"' %
                     (args.merged_fastq, args.upload_prefix))
        if manage_bamfile:
            manage_bamfile.wait_until_complete()
            log.info('Successfully uploaded %s to the specified S3 location "%s"'
                     % (args.alignment_file, args.upload_prefix))

        log.info('SEQC run complete. Cluster will be terminated')

        # upload logs
        if args.upload_prefix:
            # Upload count matrices files, logs, and return
            bucket, key = io.S3.split_link(args.upload_prefix)
            for item in [args.log_name, './nohup.log']:
                try:
                    # Make a copy of the file with the output prefix
                    copyfile(item, args.output_prefix + '_' + item)
                    print(args.output_prefix + '_' + item)
                    ec2.Retry(retries=5)(io.S3.upload_file)(
                        args.output_prefix + '_' + item, bucket, key)
                    log.info('Successfully uploaded %s to the specified S3 location '
                             '"%s".' % (item, args.upload_prefix))
                except FileNotFoundError:
                    log.notify('Item %s was not found! Continuing with upload...' % item)

        # todo local test does not send this email
        if mutt:
            email_body = (
                '<font face="Courier New, Courier, monospace">'
                'SEQC RUN COMPLETE.\n\n'
                'The run log has been attached to this email and '
                'results are now available in the S3 location you specified: '
                '"%s"\n\n' % args.upload_prefix)
            email_body = email_body.replace('\n', '<br>').replace('\t', '&emsp;')
            email_user(summary_archive, email_body, args.email)
