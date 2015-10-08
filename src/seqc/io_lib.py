__author__ = 'ambrose'


# interfaces with ftp and s3 go here.
import boto3
import os
import ftplib
from threading import Thread
from queue import Queue, Empty

# these may be broken by os.makedirs, if that function throws any errors when directories
#exist
def download_file(bucket, key, fout=None, overwrite=False):
    """download the file key located in bucket, sending output to filename fout"""
    if not overwrite:
        if os.path.isfile(fout):
            raise FileExistsError('file "%s" already exists. Set overwrite=True to '
                                  're-download' % fout)

    s3_client = boto3.client('s3')

    if key.startswith('/'):
        key = key[1:]

    if fout is None:
        fout = './' + key

    *dirs, filename = fout.split('/')
    dirs = '/'.join(dirs)
    os.makedirs(dirs)
    s3_client.download_file(bucket, key, fout)


def recursive_download(bucket_name, key_prefix, output_prefix='./', cut_dirs=0,
                       overwrite=False):
    """recursively download objects from amazon s3

    recursively downloads objects from bucket_name starting with key_prefix. If desired,
    removes a number of leading directories equal to cut_dirs. Finally, can overwrite
    files lying in the download path if overwrite is True.
    """
    # get bucket and filenames
    s3 = boto3.resource('s3')
    s3_client = boto3.client('s3')
    bucket = s3.Bucket(bucket_name)
    keys = [k.key for k in bucket.objects.all() if k.key.startswith(key_prefix)]
    if not output_prefix.endswith('/'):
        output_prefix += '/'

    # download data
    for k in keys:
        # drop first directory from output name, place in data folder
        fout = output_prefix + '/'.join(k.split('/')[cut_dirs:])
        dirs = '/'.join(fout.split('/')[:-1])
        os.makedirs(dirs)
        if os.path.isfile(fout):
            if overwrite is False:
                continue  # don't replace existing files
        s3_client.download_file(bucket_name, k, fout)


def _download_sra_file(link_queue, prefix, clobber=False, verbose=True):
    """downloads ftp_file available at 'link' into the 'prefix' directory"""

    while True:
        try:
            link = link_queue.get()
        except Empty:
            break

        # check link validity
        if not link.startswith('ftp://'):
            raise ValueError(
                'link must start with "ftp://". Provided link is not valid: %s'
                % link)

        ip, *path, file_name = link.split('/')[2:]  # [2:] -- eliminate 'ftp://'
        path = '/'.join(path)

        # check if file already exists
        if os.path.isfile(prefix + file_name):
            if not clobber:
                continue  # try to download next file

        ftp = ftplib.FTP(ip)
        try:
            ftp.login()
            ftp.cwd(path)
            with open(prefix + file_name, 'wb') as fout:
                if verbose:
                    print('beginning download of file: "%s"' % link.split('/')[-1])
                ftp.retrbinary('RETR %s' % file_name, fout.write)
                if verbose:
                    print('download of file complete: "%s"' % link.split('/')[-1])
        finally:
            ftp.close()


def download_srp(srp, prefix, max_concurrent_dl, clobber=False):
    """in-parallel download of an srp experiment"""

    if not srp.startswith('ftp://'):
        raise ValueError(
            'link must start with "ftp://". Provided link is not valid: %s' % srp)

    # create prefix directory if it does not exist
    if not os.path.isdir(prefix):
        os.makedirs(prefix)

    # make sure prefix has a trailing '/':
    if not prefix.endswith('/'):
        prefix += '/'

    # parse the download link
    ip, *path = srp.split('/')[2:]  # [2:] -- eliminate leading 'ftp://'
    path = '/' + '/'.join(path) + '/'

    # get all SRA files from the SRP experiment
    ftp = ftplib.FTP(ip)
    files = []
    try:
        ftp.login()
        ftp.cwd(path)
        dirs = ftp.nlst()  # all SRA experiments are nested in directories of the SRP
        for d in dirs:
            files.extend([path + d + '/' + f for f in ftp.nlst()])
    finally:
        ftp.close()

    if not files:
        raise ValueError('no files found in ftp directory: "%s"' % path)

    # create set of links
    if not srp.endswith('/'):
        srp += '/'

    for_download = Queue()
    for f in files:
        for_download.put(srp + f)

    threads = []
    for i in range(max_concurrent_dl):
        threads.append(Thread(target=_download_sra_file,
                              args=([for_download, prefix, clobber])))

        threads[i].start()

    for t in threads:
        t.join()
