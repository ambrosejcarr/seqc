import os
from seqc import io


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
