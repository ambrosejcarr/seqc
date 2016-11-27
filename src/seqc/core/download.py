import os
from subprocess import Popen
from seqc import io, log


# todo fold this module into io.S3 and io.Basespace
def basespace(args, output_dir: str, basespace_token: str) -> (str, str):
    """
    If --basespace argument is passed, download BaseSpace data

    :param args: Namespace object, output of ArgumentParser.parse_args()
    :param output_dir: str, output directory
    :param basespace_token: str, OAuth token for BaseSpace authentication
    :return barcode_fastq, genomic_fastq:
    """
    log.info('BaseSpace link provided for fastq argument. Downloading input data.')

    # making extra directories for BaseSpace download
    bspace_dir = output_dir + '/Data/Intensities/BaseCalls/'
    bf = Popen(['mkdir', '-p', bspace_dir])
    bf.communicate()
    barcode_fastq, genomic_fastq = io.BaseSpace.download(
        args.platform, args.basespace, output_dir, basespace_token)
    return barcode_fastq, genomic_fastq


def s3_data(files_or_links, output_prefix):
    """downloads any data provided by s3 links, otherwise gets list of files.

    :param list files_or_links: str files or str s3 links to files
    :param str output_prefix: prefix to prepend files
    :returns list files: filename(s) of downloaded files
    """
    files = []
    for f in files_or_links:
        if not f.startswith('s3://'):
            if f.endswith('/'):
                files.extend(f + subfile for subfile in os.listdir(f))
            else:
                files.append(f)
        else:
            recursive = True if f.endswith('/') else False
            files.extend(io.S3.download(f, output_prefix, overwrite=True,
                                        recursive=recursive))
    return files
