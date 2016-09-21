import argparse
import sys
from seqc import remote
from seqc import version
from seqc.exceptions import ArgumentParserError


class NewArgumentParser(argparse.ArgumentParser):
    """
    Simple wrapper for ArgumentParser that allows flags to be caught before the presence
    of 'required' arguments are tested. Allows us to add flags, for example, for
    checking the progress of existing SEQC runs.
    """
    def error(self, message):
        # checks to see whether user wants to check remote experiment status
        if '--check-progress' in sys.argv[1:]:
            remote.check_progress()
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
    meta = argparse.ArgumentParser(
        description='Processing Tools for scRNA-seq Experiments')
    meta.add_argument('-v', '--version', action='version',
                      version='{} {}'.format(meta.prog, version.__version__))
    subparsers = meta.add_subparsers(help='sub-command help', dest='subparser_name')

    # subparser for running experiments
    p = subparsers.add_parser('run', help='run help')
    # p.add_argument('platform',
    #                choices=['in_drop', 'drop_seq', 'mars1_seq',
    #                         'mars2_seq', 'in_drop_v2', 'in_drop_v3', 'ten_x'],
    #                help='which platform are you merging annotations from?')

    # a = p.add_argument_group('required arguments')
    # a.add_argument('-o', '--output-stem', metavar='O', required=True,
    #                help='If remote=True (default), an s3 link to the directory where all '
    #                     'output data should be uploaded when the run completes (e.g. '
    #                     's3://my-bucket/seqc_output_run_x/). If remote=False, the '
    #                     'complete path and prefix for the output data '
    #                     '(e.g. /Users/me/data/seqc/run_x_output). Cannot have a terminal '
    #                     '/.')
    # a.add_argument('-i', '--index', metavar='I', required=True,
    #                help='Local folder or s3 link to a directory containing the STAR '
    #                     'index used for alignment.')
    # a.add_argument('--barcode-files', nargs='*', metavar='BF', default=[],
    #                help='Either (a) an s3 link to a folder containing only barcode '
    #                     'files, or (b) the full file path of each file on the local '
    #                     'machine.')
    # a.add_argument('--email-status', metavar='E', default=None,
    #                help='Email address to receive run summary or errors when running '
    #                     'remotely. Optional if running locally.')
    #
    # i = p.add_argument_group('input arguments')
    # i.add_argument('-g', '--genomic-fastq', nargs='*', metavar='G', default=[],
    #                help='List of fastq file(s) containing genomic information, or an s3 '
    #                     'link to a directory containing only genomic fastq file(s).')
    # i.add_argument('-b', '--barcode-fastq', nargs='*', metavar='B', default=[],
    #                help='List of fastq file(s) containing barcode information, or an s3 '
    #                     'link to a directory containing only barcode fastq file(s).')
    # i.add_argument('-m', '--merged-fastq', nargs='?', metavar='M', default='',
    #                help='Filename or s3 link to a fastq file containing genomic '
    #                     'information annotated with barcode data.')
    # i.add_argument('-s', '--samfile', nargs='?', metavar='S', default='',
    #                help='Filename or s3 link to a .sam file containing aligned, merged '
    #                     'sequence records.')
    # i.add_argument('-r', '--read-array', nargs='?', metavar='RA', default='',
    #                help='Filename or s3 link to a ReadArray (.h5) archive containing '
    #                     'processed sam records.')
    # i.add_argument('--basespace', metavar='BS',
    #                help='BaseSpace sample ID. The string of numbers indicating the id '
    #                     'of the BaseSpace sample. (e.g. if the link to the sample is '
    #                     'https://basespace.illumina.com/sample/34000253/0309, '
    #                     'then --basespace would be 34000253.')
    #
    # f = p.add_argument_group('filter arguments')
    # f.add_argument('--max-insert-size', metavar='F',
    #                help='maximum paired-end insert size that is considered a valid '
    #                     'record. For multialignment correction. Not currently used.',
    #                default=1000)
    # f.add_argument('--min-poly-t', metavar='T',
    #                help='minimum size of poly-T tail that is required for a barcode to '
    #                     'be considered a valid record (default=None, automatically '
    #                     'estimates the parameter from the sequence length)',
    #                default=None, type=int)
    # f.add_argument('--max-dust-score', metavar='D',
    #                help='maximum complexity score for a read to be considered valid. '
    #                     '(default=10, higher scores indicate lower complexity.)')
    # f.add_argument('--max-ed', metavar='ED',
    #                help='Maximum hamming distance for correcting barcode errors. Should '
    #                     'be set to 1 when running in_drop legacy barcodes. (Default=2).',
    #                default=2, type=int)
    # f.add_argument('--singleton-weight', metavar='SW',
    #                help='Weight to apply to singletons in the count matrix. Float '
    #                     'between 0 and 1, default=1 (all molecules get full weight)',
    #                default=1.0, type=float)
    # f.set_defaults(filter_mitochondrial_rna=True)
    # f.add_argument('--no-filter-mitochondrial-rna', action='store_false',
    #                dest='filter_mitochondrial_rna',
    #                help='Do not filter cells with greater than 20 percent mitochondrial '
    #                     'RNA ')
    #
    # s = p.add_argument_group('alignment arguments')
    # s.add_argument('--star-args', default=None, nargs='*',
    #                help='additional arguments that should be passed to the STAR '
    #                     'aligner. For example, to set the maximum allowable times for a '
    #                     'read to align to 20, one would set '
    #                     '--star-args outFilterMultimapNmax=20. Additional arguments can '
    #                     'be provided as a white-space separated list.')
    #
    # r = p.add_argument_group('remote run arguments')
    # r.set_defaults(remote=True)
    # r.set_defaults(check=False)
    # r.add_argument('--local', dest="remote", action="store_false",
    #                help='Run SEQC locally instead of initiating on AWS EC2 servers.')
    # r.add_argument('--aws', default=False, action='store_true',
    #                help='Automatic flag; no need for user specification.')
    # r.add_argument('--no-terminate', default='False',
    #                help='Do not terminate the EC2 instance after program completes. If '
    #                     '"on-success" is specified, the EC2 instance does not terminate '
    #                     'in case the SEQC run throws an error.')
    # r.add_argument('--check-progress', dest="check", action="store_true",
    #                help='Check progress of all currently running remote SEQC runs.')
    # r.add_argument('--instance-type', default='c4',
    #                help='AWS instance (c3, c4, r3) to run SEQC remotely. Default=c4.')
    # r.add_argument('--spot-bid', type=float, default=None,
    #                help='float, Amount to bid for a spot instance. Default=None (will '
    #                     'reserve a non-spot instance). WARNING: using spot instances '
    #                     'will cause your instance to terminate if instance prices exceed '
    #                     'your spot bid during runtime.')
    # r.add_argument('--log-name', type=str, default='seqc.log',
    #                help='Output log name (default=seqc.log)')

    # add subparser to check progress
    progress = subparsers.add_parser('check_progress', help='check SEQC run progress')

    # add subparser to create index
    index = subparsers.add_parser('create_index', help='create a SEQC index')

    index.add_argument(
        '-o', '--organism', required=True,
        help='organism to create index for. Must be genus_species in all lower-case. '
             'e.g. human is homo_sapiens.')
    index.add_argument('-f', '--folder', default='')
    index.add_argument('-u', '--upload-location', help='s3 link')
    index.add_argument('--ids', '--additional-id-types',
                       help='additional ids to check against')
    index.add_argument('-b', '--valid-biotypes', help='valid id biotypes')

    try:
        return p.parse_args(args)
    except ArgumentParserError:
        raise


def generate_remote_cmdline_args(argv: list) -> str:
    """recreate the command line arguments for a remote run, appending --local and --aws

    :param argv: the output of sys.argv[1:], in other words, all parameters passed to
      the script, omitting the script name (process_experiment.py)
    :return str: command line arguments for a remote run
    """
    return 'process_experiment.py ' + ' '.join(argv) + ' --local --aws'


def recreate_cmdline_args(argv: list) -> str:
    """recreate the command line arguments for a remote run, appending --local and --aws

    :param argv: the output of sys.argv[1:], in other words, all parameters passed to
      the script, omitting the script name (process_experiment.py)
    :return str: command line arguments for a remote run
    """
    return 'process_experiment.py ' + ' '.join(argv)
