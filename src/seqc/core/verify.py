from seqc import io, platforms
import shutil
import inspect
from math import ceil
from seqc.core import config


def estimate_required_volume_size(args):
    """

    :param args:
    :return:
    """

    # keep track of how much space is needed given input
    seqc_input = args.barcode_files + [args.index]

    # using worst-case estimates to make sure we don't run out of space
    cushion = 9e10
    total = 0

    if args.barcode_fastq and args.genomic_fastq:
        seqc_input = seqc_input + args.barcode_fastq + args.genomic_fastq
        io.S3.check_links(seqc_input)

        # checking size of input file
        input_fastq = args.barcode_fastq + args.genomic_fastq
        for item in input_fastq:
            total += io.S3.obtain_size(item)
        total += (total * 14) + cushion
    if args.samfile:
        seqc_input += [args.samfile]
        io.S3.check_links(seqc_input)

        # checking size of input file
        total += io.S3.obtain_size(args.samfile)
        total += (total * 2) + 2e10
    if args.merged_fastq:
        seqc_input += [args.merged_fastq]
        io.S3.check_links(seqc_input)

        # checking size of input file
        total += io.S3.obtain_size(args.merged_fastq)
        total += (total * 13) + cushion
    if args.basespace:
        if not args.basespace_token or args.basespace_token == 'None':
            raise ValueError(
                'If the --basespace argument is used, the basespace token must be '
                'specified in the seqc config file or passed as --basespace-token')
        seqc_input += [args.basespace]
        io.S3.check_links(seqc_input)
    if args.read_array:
        seqc_input += [args.read_array]
        seqc_input.remove(args.index)
        io.S3.check_links(seqc_input)

        # checking size of input
        for item in seqc_input:
            total += io.S3.obtain_size(item)
        total += 1e10
    if args.basespace:
        io.BaseSpace.check_sample(args.basespace, args.basespace_token)
        # checking size of input file
        total = io.BaseSpace.check_size(args.basespace, args.basespace_token)
        total += (total * 14) + cushion

    # return total size needed for EBS volume
    return ceil(total * 1e-9)


def run(args) -> float:
    """
    verify data input through the command line arguments, fixes minor issues, and
    throws exceptions if invalid parameters are encountered

    additionally, this function obtains a rough estimate of how much
    volume storage is needed for a remote run.

    :param Namespace args: Namespace object, output from ArgumentParser.parse_args()
    :returns total: float, estimated Kb of Volume space needed to run SEQC remotely.
    """

    if args.output_prefix.endswith('/'):
        raise ValueError('output_stem should not be a directory.')
    if not args.index.endswith('/'):
        raise ValueError('index must be a directory, and must end with "/"')

    # check platform name; raises ValueError if invalid
    platform_name(args.platform)

    # check to make sure that --email-status is passed with remote run
    if args.remote and not args.email:
        raise ValueError('Please supply the --email-status flag for a remote SEQC run.')
    # if args.instance_type not in ['c3', 'c4', 'r3']:  # todo fix this instance check
    #     raise ValueError('All AWS instance types must be either c3, c4, or r3.')
    # if args.terminate not in ['True', 'true', 'False', 'false', 'on-success']:
    #     raise ValueError('the --no-terminate flag must be either True, False, '
    #                      'or on-success.')

    # make sure at least one input has been passed
    valid_inputs = (
        args.barcode_fastq, args.genomic_fastq, args.merged_fastq, args.samfile,
        args.basespace, args.read_array)
    if not any(valid_inputs):
        raise ValueError(
            'At least one input argument (-b/-g, -m, -s, -r, --basespace) must be passed '
            'to SEQC.')
    if not args.barcode_files:  # todo clean this up and fold into platform somehow
        if args.platform != 'drop_seq':
            raise ValueError('--barcode-files is required for this platform.')

    # make sure at most one input type has been passed
    num_inputs = 0
    if args.barcode_fastq or args.genomic_fastq:
        if not all((args.barcode_fastq, args.genomic_fastq)):
            raise ValueError(
                'if either genomic or barcode fastq are provided, both must be provided')
        num_inputs += 1
    num_inputs += sum(1 for i in (args.merged_fastq, args.samfile,
                                  args.basespace, args.read_array) if i)
    if num_inputs > 1:
        raise ValueError(
            'user should provide at most one input argument (-b/-g, -m, -s, -r, '
            '--basespace')

    # if basespace is being used, make sure there is a valid basespace token
    if args.basespace and not hasattr(args, 'basespace_token'):
        configuration = config.read_config('~/.seqc/config')
        setattr(
            args, 'basespace_token', configuration['BaseSpaceToken']['base_space_token'])
    else:
        setattr(args, 'basespace_token', None)

    # check that spot-bid is correct
    if args.spot_bid is not None:
        if args.spot_bid < 0:
            raise ValueError('bid %f must be a non-negative float.' % args.spot_bid)

    if args.upload_prefix and not args.upload_prefix.startswith('s3://'):
        raise ValueError('upload_prefix should be an s3 address beginning with s3://')

    if args.volume_size is None:
        setattr(args, 'volume_size', estimate_required_volume_size(args))

    return args


def index(args):
    """add a default volume_size if it was not otherwise passed to seqc.

    :param args:
    :return:
    """
    if args.volume_size is None:
        setattr(args, 'volume_size', 100)
    return args


def install():
    raise NotImplementedError  # todo implement; verify aws key, everything else.


def executables(*execs):
    """
    checks whether executables are installed on the machine of the
    current seqc run.

    :param execs: Tuple of executables to check
    :returns : Tuple of boolean (True if a specific executable is installed).
    """
    return tuple(map(lambda exe: shutil.which(exe) is not None, execs))


def platform_name(name: str):
    """
    checks whether the platform name supplied by the user is supported by the current
    iteration of seqc.
    :param name: string of platform name to check
    :return: name (if supported by seqc).
    """
    choices = [x[0] for x in inspect.getmembers(platforms, inspect.isclass) if
               issubclass(x[1], platforms.AbstractPlatform)][1:]
    if name not in choices:
        raise ValueError('Please specify a valid platform name for SEQC. The available '
                         'options are: {}'.format(choices))
    # throw error for mars1_seq since we don't have the appropriate primer length yet
    if name == 'mars1_seq':
        raise ValueError('Mars1-seq is currently not stable in this version of SEQC.')
    return name
