import os
import shutil
import inspect
from math import ceil
from seqc import io, platforms


def filesize(filename):
    """return filesize of filename in bytes

    :param str filename: full path to file
    :return int: number of bytes in filename
    """
    return os.stat(filename).st_size


def validate_and_return_size(filename):
    """return true if a link or filepath points to a valid file or directory

    :param str filename: filepath or s3 link
    :return None: raises errors if path or link is invalid.
    """
    if filename.startswith('s3://'):
        io.S3.check_links([filename])
        return io.S3.obtain_size(filename)
    else:
        if os.path.isfile(filename):
            return filesize(filename)
        elif os.path.isdir(filename.rstrip('/')):
            return sum(filesize(filename + f) for f in os.listdir(filename))
        else:
            print(filename)
            raise ValueError('%s does not point to a valid file')


def estimate_required_volume_size(args):
    """estimate the size of volume that should be attached to an aws instance to run SEQC

    :param args: namespace object containing filepaths or download links to input data
    :return int: size of volume in gb
    """

    total = sum(validate_and_return_size(f) for f in args.barcode_files)

    # using worst-case estimates to make sure we don't run out of space

    # todo stopped here; remove aws dependency
    if args.barcode_fastq and args.genomic_fastq:
        total += sum(validate_and_return_size(f) for f in args.barcode_fastq) * 14 + 9e10
        total += sum(validate_and_return_size(f) for f in args.genomic_fastq) * 14 + 9e10
        total += validate_and_return_size(args.index)

    elif args.samfile:
        total += (validate_and_return_size(args.samfile) * 2) + 2e10
        total += validate_and_return_size(args.index)

    elif args.merged_fastq:
        total += (validate_and_return_size(args.merged_fastq) * 13) + 9e10
        total += validate_and_return_size(args.index)

    elif args.read_array:
        total += validate_and_return_size(args.read_array)
    if args.basespace:
        if not args.basespace_token or args.basespace_token == 'None':
            raise ValueError(
                'If the --basespace argument is used, the basespace token must be '
                'specified in the seqc config file or passed as --basespace-token')

        io.BaseSpace.check_sample(args.basespace, args.basespace_token)
        total += io.BaseSpace.check_size(args.basespace, args.basespace_token)

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

    if args.rsa_key is None:
        raise ValueError('-k/--rsa-key does not point to a valid file object. ')
    if not os.path.isfile(args.rsa_key):
        raise ValueError('-k/--rsa-key does not point to a valid file object. ')

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
        raise RuntimeError('if --basespace input is selected, user must provide an OAuth '
                           'token using the --basespace-token parameter.')

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

    :param args: namespace object from argparse
    :return: updated namespace object with volume_size set.
    """
    if args.volume_size is None:
        setattr(args, 'volume_size', 100)
    return args


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
