import argparse
import sys
from subprocess import Popen, PIPE
from seqc import version
from seqc.exceptions import ArgumentParserError


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
    subparsers = meta.add_subparsers(dest='subparser_name')

    # subparser for running experiments
    # can use to make prettier: formatter_class=partial(argparse.HelpFormatter, width=200)
    p = subparsers.add_parser('run', help='initiate SEQC runs')

    p.add_argument('platform',
                   choices=['in_drop', 'drop_seq', 'mars1_seq',
                            'mars2_seq', 'in_drop_v2', 'in_drop_v3', 'ten_x'],
                   help='which platform are you merging annotations from?')

    a = p.add_argument_group('required arguments')
    a.add_argument('-o', '--output-prefix', metavar='O', required=True,
                   help='filename prefix for all seqc output. Should not be a directory.')
    a.add_argument('-i', '--index', metavar='I', required=True,
                   help='Local folder or s3 link to a directory containing the STAR '
                        'index used for alignment.')
    a.add_argument('--barcode-files', nargs='*', metavar='BF', default=[],
                   help='Either (a) an s3 link to a folder containing only barcode '
                        'files, or (b) the full file path of each file on the local '
                        'machine.')

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
                        'of the BaseSpace sample. (e.g. if the link to the sample is '
                        'https://basespace.illumina.com/sample/34000253/0309, '
                        'then --basespace would be 34000253.')
    i.add_argument('--basespace-token', metavar='BST', default=None,
                   help='OAuth token for basespace access. Required if BaseSpace input '
                        'is used.')

    f = p.add_argument_group('filter arguments')
    f.add_argument('--max-insert-size', metavar='F',
                   help='maximum distance from the TTS for a read to be considered valid'
                        '. For multialignment correction. Not currently used.',
                   default=1000)
    f.add_argument('--min-poly-t', metavar='T',
                   help='minimum size of poly-T tail that is required for a barcode to '
                        'be considered a valid record (default=None, automatically '
                        'estimates the parameter from the sequence length)',
                   default=None, type=int)
    f.add_argument('--max-dust-score', metavar='D', default=10, type=int,
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
    f.set_defaults(filter_mitochondrial_rna=True)
    f.add_argument('--no-filter-mitochondrial-rna', action='store_false',
                   dest='filter_mitochondrial_rna',
                   help='Do not filter cells with greater than 20 percent mitochondrial '
                        'RNA ')

    s = p.add_argument_group('alignment arguments')
    s.add_argument('--star-args', default=None, nargs='*',
                   help='additional arguments that should be passed to the STAR '
                        'aligner. For example, to set the maximum allowable times for a '
                        'read to align to 20, one would set '
                        '--star-args outFilterMultimapNmax=20. Additional arguments can '
                        'be provided as a white-space separated list.')

    progress = subparsers.add_parser('progress', help='check SEQC run progress')
    progress.set_defaults(remote=False)
    progress.add_argument(
        '-i', '--instance-ids', help='check the progress of run(s)', nargs='+')
    progress.add_argument(
        '-k', '--rsa-key', help='RSA key registered to your aws account', required=True)

    terminate = subparsers.add_parser('terminate', help='terminate SEQC runs')
    terminate.set_defaults(remote=False)
    terminate.add_argument(
        '-i', '--instance-ids', help='terminate these instance(s)', nargs='+')

    instances = subparsers.add_parser('instances', help='list all running instances')
    instances.set_defaults(remote=False)
    instances.add_argument(
        '-k', '--rsa-key', help='RSA key registered to your aws account', required=True)

    pindex = subparsers.add_parser('index', help='create a SEQC index')
    pindex.add_argument(
        '-o', '--organism', required=True,
        help='organism to create index for. Must be formatted as genus_species in all '
             'lower-case. e.g. human is homo_sapiens.')
    pindex.add_argument(
        '-f', '--folder', default=None,
        help='folder in which to create the index. Defaults to the name of the organism, '
             'which is created in the current directory.')
    pindex.add_argument(
        '--ids', '--additional-id-types', nargs='*',
        help='names of additional ids from other consortia to check against. If '
             'provided, each ENSEMBL gene id must also be annotated by at least one of '
             'these consortia to be considered valid and appear in the final SEQC count '
             'matrix.')
    pindex.add_argument(
        '-b', '--valid-biotypes', default=('protein_coding', 'lincRNA'),
        help='list of gene biotypes that are considered valid. Defaults are '
             'protein_coding and lincRNA. In most cases, other biotypes are not expected '
             'to be captured by SEQC, and should be excluded')

    for parser in [pindex, p]:
        r = parser.add_argument_group('Amazon Web Services arguments')
        r.set_defaults(remote=True)
        r.set_defaults(terminate=True)
        r.set_defaults(log_name='seqc_log.txt')  # changed to .txt for email
        r.add_argument(
            '--local', dest="remote", action="store_false",
            help='Run locally instead of on an aws instance')
        r.add_argument(
            '-u', '--upload-prefix', metavar='U', default=None,
            help='s3 location for data to be uploaded.')
        r.add_argument(
            '--no-terminate', action='store_false', dest='terminate',
            help='Do not terminate the EC2 instance after program completes. Normally, '
                 'the instance is terminated unless an error is encountered.')
        r.add_argument(
            '--instance-type', default='c4.8xlarge',
            help='AWS instance type to initialize for this job. '
                 'See https://aws.amazon.com/ec2/instance-types/ for valid types')
        r.add_argument(
            '--spot-bid', type=float, default=None,
            help='float, Amount to bid for a spot instance. Default=None (will reserve a '
                 'non-spot instance). WARNING: using spot instances will cause your '
                 'instance to terminate if instance prices exceed your spot bid during '
                 'runtime.')
        r.add_argument(
            '--volume-size', type=int, default=None,
            help='size in Gb required to execute the requested process. If not provided, '
                 'it will be estimated from passed parameters.')
        r.add_argument(
            '-e', '--email', metavar='E', default=None,
            help='Email address to receive run summary or errors when running remotely. '
                 'Optional only if running locally.')
        r.add_argument(
            '-k', '--rsa-key', metavar='K', default=None,
            help='RSA key registered to your aws account that allowed access to ec2 '
                 'resources. Required if running instance remotely.')

    # custom help handling
    if len(args) == 0:  # print help if no args are passed
        meta.print_help()
        sys.exit(1)
    if args == ['run', '-h']:  # send help for run to less, is too long
        pipe = Popen(['less'], stdin=PIPE)
        pipe.communicate(p.format_help().encode())
        sys.exit(1)

    try:
        return meta.parse_args(args)
    except ArgumentParserError:
        raise


def generate_remote_cmdline_args(argv: list) -> str:
    """recreate the command line arguments for a remote run, appending --local and --aws

    :param argv: the output of sys.argv[1:], in other words, all parameters passed to
      the script, omitting the script name (process_experiment.py)
    :return str: command line arguments for a remote run
    """
    return 'process_experiment.py ' + ' '.join(argv) + ' --local --aws'
