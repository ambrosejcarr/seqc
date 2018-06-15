import os
import argparse
import sys
import inspect
from subprocess import Popen, PIPE
from seqc import version, platforms


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

    # Platform choices
    choices = [x[0] for x in inspect.getmembers(platforms, inspect.isclass) if
           issubclass(x[1], platforms.AbstractPlatform)][1:]
    p.add_argument('platform',
                   choices=choices,
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
    i.add_argument('-a', '--alignment-file', nargs='?', metavar='A', default='',
                   help='Filename or s3 link to a .sam or .bam file containing aligned, '
                        'merged sequence records.')
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
    f.add_argument('--max-insert-size', metavar='F', type=int,
                   help='the maximum fragment size in bp. Aligments that are further '
                        'than this distance from a TTS are discarded. Default=1000',
                   default=1000)
    f.add_argument('--min-poly-t', metavar='T',
                   help='minimum size of poly-T tail that is required for a barcode to '
                        'be considered a valid record (default=None, automatically '
                        'estimates the parameter from the sequence length)',
                   default=None, type=int)
    # f.add_argument('--max-dust-score', metavar='D', default=10, type=int,
    #                help='maximum complexity score for a read to be considered valid. '
    #                     '(default=10, higher scores indicate lower complexity.)')
    f.add_argument('--singleton-weight', metavar='SW',
                   help='Weight to apply to singletons in the count matrix. Float '
                        'between 0 and 1, default=1 (all molecules get full weight)',
                   default=1.0, type=float)
    f.set_defaults(filter_mitochondrial_rna=True)
    f.add_argument('--no-filter-mitochondrial-rna', action='store_false',
                   dest='filter_mitochondrial_rna',
                   help='Do not filter cells with greater than 20 percent mitochondrial '
                        'RNA ')
    f.set_defaults(filter_low_coverage=True)
    f.add_argument('--no-filter-low-coverage', action='store_false',
                   dest='filter_low_coverage',
                   help='Do not filter cells with low coverage')
    f.set_defaults(filter_low_gene_abundance=True)
    f.add_argument('--no-filter-low-gene-abundance', action='store_false',
                   dest='filter_low_gene_abundance',
                   help='Do not filter cells with low coverage')
    f.add_argument('--low-coverage-alpha', metavar='LA',
                   help='FDR rate for low coverage reads filter in mars-seq datasets. '
                        'Float between 0 and 1, default=0.25',
                   default=0.25, type=float)

    s = p.add_argument_group('alignment arguments')
    s.add_argument('--star-args', default=None, nargs='*',
                   help='additional arguments that should be passed to the STAR '
                        'aligner. For example, to set the maximum allowable times for a '
                        'read to align to 20, one would set '
                        '--star-args outFilterMultimapNmax=20. Additional arguments can '
                        'be provided as a white-space separated list.')

    # PROGRESS PARSER
    progress = subparsers.add_parser('progress', help='check SEQC run progress')
    progress.set_defaults(remote=False)
    progress.add_argument(
        '-i', '--instance-ids', help='check the progress of run(s)', nargs='+')
    progress.add_argument(
        '-k', '--rsa-key', help='RSA key registered to your aws account',
        default=None)

    # TERMINATE PARSER
    terminate = subparsers.add_parser('terminate', help='terminate SEQC runs')
    terminate.set_defaults(remote=False)
    terminate.add_argument(
        '-i', '--instance-ids', help='terminate these instance(s)', nargs='+')

    # INSTANCES PARSER
    instances = subparsers.add_parser('instances', help='list all running instances')
    instances.set_defaults(remote=False)
    instances.add_argument(
        '-k', '--rsa-key', help='RSA key registered to your aws account',
        default=None)

    # START PARSER
    start = subparsers.add_parser(
        'start', help='initialize a seqc-ready instance')
    start.set_defaults(remote=False)
    start.add_argument(
        '-s', '--volume-size', help='size of volume (Gb) to attach to instance',
        default=5, type=int)
    start.add_argument(
        '-b', '--spot-bid', help='amount to bid for instance in fractions of dollars',
        type=float, default=None)
    start.add_argument(
        '-t', '--instance-type', default='c4.8xlarge',
        help='AWS instance type to initialize. '
             'See https://aws.amazon.com/ec2/instance-types/ for valid types')
    start.add_argument(
        '-k', '--rsa-key', help='RSA key registered to your aws account',
        default=None)

    # NOTEBOOK PARSERS
    notebook_sp = subparsers.add_parser('notebook', help='notebook tools')
    _nb_parser = notebook_sp.add_subparsers(dest='subsubparser_name')

    # NOTEBOOK MERGE PARSER
    merge = _nb_parser.add_parser(
        'merge', help='merge multiple datasets prior to running an analysis notebook')
    merge.add_argument(
        '-o', '--output-filename', help='name for merged fastq file', required=True)
    merge.add_argument(
        '-i', '--input-data', nargs='+', help='count matrices to merge', required=True)

    # NOTEBOOK GENERATE PARSER
    generate = _nb_parser.add_parser('generate', help='generate a notebook from a dataset')
    generate.add_argument(
        '-i', '--input-count-matrix', help='count matrix file', required=True)
    generate.add_argument(
        '-o', '--output-stem', help='directory and filestem for output', required=True)

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
        r.add_argument('--debug', default=False, action='store_true',
                       help='If debug is set, runs that throw errors do not '
                            'terminate the instance they were run on.')
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

    parsed = meta.parse_args(args)

    if hasattr(parsed, 'rsa_key'):
        if parsed.rsa_key is None:
            try:
                parsed.rsa_key = os.environ['AWS_RSA_KEY']
            except KeyError:
                pass

    return parsed
