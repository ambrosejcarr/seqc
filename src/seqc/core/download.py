import os
from subprocess import Popen
from seqc import io, log


def basespace(args, output_dir: str, basespace_token: str) -> (str, str):
    """
    If --basespace argument is passed, download BaseSpace data

    :param args: Namespace object, output of ArgumentParser.parse_args()
    :param output_dir: str, output directory
    :param basespace_token: str, OAuth token for BaseSpace authentication
    :return barcode_fastq, genomic_fastq:
    """
    log.info('BaseSpace link provided for fastq argument. Downloading input data.')

    # making extra directories for BaseSpace download, changing permissions
    bspace_dir = output_dir + '/Data/Intensities/BaseCalls/'
    bf = Popen(['sudo', 'mkdir', '-p', bspace_dir])
    bf.communicate()
    if args.aws:  # changing permissions is unnecessary if local run
        bf2 = Popen(['sudo', 'chown', '-c', 'ubuntu', bspace_dir])
        bf2.communicate()
    barcode_fastq, genomic_fastq = io.BaseSpace.download(
        args.platform, args.basespace, output_dir, basespace_token)
    return barcode_fastq, genomic_fastq


def s3_fastq(fastq_file: list, output_dir: str, ftype: str) -> list:
    """
    Checks if -g/--genomic-fastq was passed. If it was, downloads the necessary files if
    any passed arguments were s3 links.

    :param fastq_file: list, a list of fastq files or s3 links to fastq files
    :param output_dir: directory in which to download fastq files
    :param ftype: denotes whether fastq files are genomic or barcode
    :returns fastq_file: list, filename(s) of local genomic_fastq files.
    """
    if fastq_file:
        if not fastq_file[0].startswith('s3://'):
            for gf in fastq_file:
                if not os.path.isfile(gf):
                    raise ValueError('Provided %s fastq files: "[%s]" is '
                                     'neither an s3 link or a valid filepath' %
                                     (ftype, ', '.join(map(str, fastq_file))))
        else:
            log.info('Downloading {} fastq files from Amazon s3 link.'.format(ftype))
            if fastq_file[0].endswith('/'):
                # s3 directory specified, download all files recursively
                bucket, prefix = io.S3.split_link(fastq_file[0])
                cut_dirs = prefix.count('/')
                fastq_file = io.S3.download_files(
                    bucket=bucket, key_prefix=prefix, output_prefix=output_dir,
                    cut_dirs=cut_dirs)
            else:
                # individual s3 links provided, download each fastq file
                downloaded_files = []
                for s3link in fastq_file:
                    bucket, prefix = io.S3.split_link(s3link)
                    _, fname = os.path.split(prefix)
                    fname = output_dir + '/' + fname
                    io.S3.download_file(bucket, prefix, fname)
                    downloaded_files.append(fname)
                fastq_file = sorted(downloaded_files)
            # note that this is printed regardless of whether a file is downloaded
            log.info('%s fastq files [%s] successfully installed.' %
                     (ftype, ', '.join(map(str, fastq_file))))
    return fastq_file


def barcodes(platform: str, barcode_files: list, output_dir: str) -> list:
    """
    if an s3 link was passed and the platform was not drop-seq, download barcode files

    :param platform: str, type of experiment
    :param barcode_files: list, files containing barcode information or s3 links
    :param output_dir: str, directory for output files
    :returns barcode_files: list, files containing barcodes
    """
    if platform != 'drop_seq':
        if not barcode_files[0].startswith('s3://'):
            for cb in barcode_files:
                if not os.path.isfile(cb):
                    raise ValueError('provided barcode files: "[%s]" is neither '
                                     'an s3 link or a valid filepath' %
                                     ', '.join(map(str, barcode_files)))
        else:
            log.info('AWS s3 link provided for barcodes. Downloading files.')
            if not barcode_files[0].endswith('/'):
                barcode_files[0] += '/'
            bucket, prefix = io.S3.split_link(barcode_files[0])
            cut_dirs = prefix.count('/')
            barcode_files = io.S3.download_files(
                bucket=bucket, key_prefix=prefix, output_prefix=output_dir,
                cut_dirs=cut_dirs)
            log.info('Barcode files [%s] successfully downloaded.' %
                     ', '.join(map(str, barcode_files)))
    return barcode_files


def index(output_dir: str, index_dir: str, read_array: str, samfile: str) -> str:
    """
    Checks the index parameter, downloading files from s3 if an s3 link was passed

    :param output_dir: str, directory that the index will be downloaded to if necessary
    :param index_dir: value of -i/--index parameter
    :param read_array: str, value of -r/--read-array parameter
    :param samfile: str, value of -s/--samfile parameter
    :return index: str, path to downloaded or local index
    """
    if not index_dir.startswith('s3://'):
        if not os.path.isdir(index_dir):
            raise ValueError('Provided index: "%s" is neither an s3 link or a valid '
                             'filepath' % index_dir)
    else:
        log.info('AWS s3 link provided for index. Downloading index.')
        bucket, prefix = io.S3.split_link(index_dir)
        index_dir = output_dir + '/index/'  # set index  based on s3 download
        cut_dirs = prefix.count('/')
        # install whole index
        if not any([samfile, read_array]):
            io.S3.download_files(bucket=bucket, key_prefix=prefix,
                                 output_prefix=index_dir, cut_dirs=cut_dirs)
        else:  # either samfile or read array provided, only download annotations file
            annotations_file = index_dir + 'annotations.gtf'
            io.S3.download_file(bucket, prefix + 'annotations.gtf', annotations_file)
    return index_dir




