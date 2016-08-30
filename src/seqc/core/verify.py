from seqc import io
import shutil


def arguments(args, basespace_token: str) -> float:
    """
    verify data input through the command line arguments and throws
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
        io.S3.check_links(seqc_input)

        # checking size of input file
        input_fastq = barcode_fastq + genomic_fastq
        for item in input_fastq:
            total += io.S3.obtain_size(item)
        total += (total * 14) + cushion
    if samfile:
        if any((merged, barcode_fastq, genomic_fastq, basespace, read_array)):
            raise ValueError(multi_input_error_message)
        seqc_input += [samfile]
        io.S3.check_links(seqc_input)

        # checking size of input file
        total += io.S3.obtain_size(samfile)
        total += (total * 2) + 2e10
    if merged:
        if any((samfile, barcode_fastq, genomic_fastq, basespace, read_array)):
            raise ValueError(multi_input_error_message)
        seqc_input += [merged]
        io.S3.check_links(seqc_input)

        # checking size of input file
        total += io.S3.obtain_size(merged)
        total += (total * 13) + cushion
    if basespace:
        if any((samfile, merged, barcode_fastq, genomic_fastq, read_array)):
            raise ValueError(multi_input_error_message)
        if not basespace_token or basespace_token == 'None':
            raise ValueError(
                'If the --basespace argument is used, the BaseSpace token must be '
                'specified in the seqc config file.')
        seqc_input += [basespace]
        io.S3.check_links(seqc_input)
    if read_array:
        if any((samfile, merged, barcode_fastq, genomic_fastq, basespace)):
            raise ValueError(multi_input_error_message)
        seqc_input += [read_array]
        seqc_input.remove(index_dir)
        io.S3.check_links(seqc_input)

        # checking size of input
        for item in seqc_input:
            total += io.S3.obtain_size(item)
        total += 1e10
    if basespace:
        io.BaseSpace.check_sample(basespace, basespace_token)
        # checking size of input file
        total = io.BaseSpace.check_size(basespace, basespace_token)
        total += (total * 14) + cushion

    # return total size needed for EBS volume
    return total


def executables(*execs):
    """
    checks whether executables are installed on the machine of the
    current seqc run.

    :param execs: Tuple of executables to check
    :returns : Tuple of boolean (True if a specific executable is installed).
    """
    return tuple(map(lambda exe: shutil.which(exe) is not None, execs))
