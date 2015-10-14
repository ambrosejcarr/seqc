__author__ = 'ambrose'


# interfaces with ftp and s3 go here.
from glob import glob
import boto3
import os
import ftplib
from threading import Thread
from queue import Queue, Empty
from subprocess import Popen, check_output
from itertools import zip_longest


class S3:
    """A series of methods to upload and download files from amazon s3"""

    @staticmethod
    def download_file(bucket, key, fout=None, overwrite=False):
        """download the file key located in bucket, sending output to filename fout"""
        if not overwrite:
            if os.path.isfile(fout):
                raise FileExistsError('file "%s" already exists. Set overwrite=True to '
                                      're-download' % fout)

        # key should not start with a forward slash
        if key.startswith('/'):
            key = key[1:]

        # if fout is not provided, download the key, cutting all directories
        if fout is None:
            fout = key.split('/')[-1]

        # check if all directories exist. if not, create them.
        *dirs, filename = fout.split('/')
        dirs = '/'.join(dirs)
        if not os.path.isdir(dirs):
            os.makedirs(dirs)

        # download the file
        client = boto3.client('s3')
        client.download_file(bucket, key, fout)

        return fout

    # todo return all downloaded files
    @classmethod
    def download_files(cls, bucket, key_prefix, output_prefix='./', cut_dirs=0,
                       overwrite=False):
        """
        recursively download objects from amazon s3

        recursively downloads objects from bucket starting with key_prefix.
        If desired, removes a number of leading directories equal to cut_dirs. Finally,
        can overwrite files lying in the download path if overwrite is True.
        """
        # get bucket and filenames
        client = boto3.client('s3')
        keys = cls.list(bucket, key_prefix)

        # output prefix needs to end in '/'
        if not output_prefix.endswith('/'):
            output_prefix += '/'

        # download data
        for k in keys:

            # drop first directory from output name, place in data folder
            fout = output_prefix + '/'.join(k.split('/')[cut_dirs:])
            dirs = '/'.join(fout.split('/')[:-1])

            # make directories if they don't exist
            if not os.path.isdir(dirs):
                os.makedirs(dirs)

            # check for overwriting
            if os.path.isfile(fout):
                if overwrite is False:
                    continue

            client.download_file(bucket, k, fout)

    @staticmethod
    def upload_file(filename, bucket, key):
        """upload filename to aws at s3://bucket/key/filename"""

        if key.startswith('/'):
            key = key[1:]

        if key.endswith('/'):
            file_id = filename.split('/')[-1]  # get file id, stripping directories
            key += file_id

        if not os.path.isfile(filename):
            raise FileNotFoundError('file "%s" is not a valid file identifier' % filename)

        client = boto3.client('s3')
        client.upload_file(filename, bucket, key)

    @staticmethod
    def upload_files(file_prefix, bucket, key_prefix, cut_dirs=True):
        """
        upload all files f found at file_prefix to s3://bucket/key_prefix/f

        This function eliminates any uninformative directories. For example, if uploading
        a file_prefix such as '/data/*' to bucket 'MyBucket', at key_prefix 'tmp/', which
        produces the following set of files:
        /data/useless/file1
        /data/useless/file2
        /data/useless/dir1/file3
        /data/useless/dir2/dir3/file4

        the upload with proceed as follows unless cut_dirs is False:
        s3://MyBucket/tmp/file1
        s3://MyBucket/tmp/file2
        s3://MyBucket/tmp/dir1/file3
        s3://MyBucket/tmp/dir2/dir3/file4
        """

        # if a wildcard was present, we will need to filter files in the last directory
        all_files = []
        if "*" in file_prefix:
            allowed_prefixes = set(glob(file_prefix))
            for file_or_dir in allowed_prefixes:
                if os.path.isdir(file_or_dir):
                    for path, subdirs, files in os.walk(file_or_dir):
                        for name in files:
                            all_files.append(os.path.join(path, name))
                else:
                    all_files.append(file_or_dir)
        else:  # if no wildcard, walk the directory to get all of the file ids
            for path, subdirs, files in os.walk(file_prefix):
                for name in files:
                    all_files.append(os.path.join(path, name))

        if not key_prefix.endswith('/'):
            key_prefix += '/'

        if cut_dirs:  # get uninformative directories.
            # zip together each directory level, starting from the root. If a file
            # length is too short, zip_longest inserts "None" and that directory is
            # therefore informative.
            directories = zip_longest(*(f.lstrip('/').split('/') for f in all_files))
            n_cut = 0
            for directory in directories:
                if len(set(directory)) == 1:
                    n_cut += 1
                else:
                    break
            upload_keys = [key_prefix + '/'.join(f.lstrip('/').split('/')[n_cut:])
                           for f in all_files]
        else:
            upload_keys = [key_prefix + f.lstrip('/') for f in all_files]

        client = boto3.client('s3')
        for file_, key in zip(all_files, upload_keys):
            client.upload_file(file_, bucket, key)


    @staticmethod
    def list(bucket, key_prefix):
        """
        list all objects beginning with key_prefix

        since amazon stores objects as keys, not as a filesystem, this is synonymous to
        listing the directory defined by key_prefix
        """
        s3 = boto3.resource('s3')
        bucket = s3.Bucket(bucket)
        keys = [k.key for k in bucket.objects.all() if k.key.startswith(key_prefix)]
        return keys

    @staticmethod
    def remove_file(bucket, key):
        """delete AWS s3 file at s3://bucket/key"""
        client = boto3.client('s3')
        _ = client.delete_object(Bucket=bucket, Key=key)

    @classmethod
    def remove_files(cls, bucket, key_prefix):
        keys = cls.list(bucket, key_prefix)
        client = boto3.client('s3')
        for k in keys:
            _ = client.delete_object(Bucket=bucket, Key=k)


def _download_sra_file(link_queue, prefix, clobber=False, verbose=True):
    """downloads ftp_file available at 'link' into the 'prefix' directory"""

    while True:
        try:
            link = link_queue.get_nowait()
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


def parallel_download_srp(srp, prefix, max_concurrent_dl, verbose=True, clobber=False):
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

    if not srp.endswith('/'):
        srp += '/'

    # parse the download link
    ip, *path = srp.split('/')[2:]  # [2:] -- eliminate leading 'ftp://'
    path = '/' + '/'.join(path)

    # get all SRA files from the SRP experiment
    ftp = ftplib.FTP(ip)
    files = []
    try:
        ftp.login()
        ftp.cwd(path)
        dirs = ftp.nlst()  # all SRA experiments are nested in directories of the SRP
        for d in dirs:
            files.append('%s/%s.sra' % (d, d))
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
                              args=([for_download, prefix, clobber, verbose])))

        threads[i].start()

    for t in threads:
        t.join()

    # get output files
    output_files = []
    for f in files:
        output_files.append(prefix + f.split('/')[-1])

    return files


def _extract_fastq(sra_queue, verbose=True):

    while True:
        try:
            file_ = sra_queue.get_nowait()
        except Empty:
            break

        *dir, file = file_.split('/')
        dir = '/'.join(dir)
        os.chdir(dir)  # set working directory as files are output in cwd

        # extract file
        if verbose:
            print('beginning extraction of file: "%s"' % file_)
        Popen(['fastq-dump', '--split-3', file_])
        if verbose:
            print('extraction of file complete: "%s"' % file_)


def parallel_extract_fastq(sra_files, max_concurrent, verbose):
    """requires fastq-dump from sra-tools"""

    # check that fastq-dump exists
    if not check_output(['which', 'fastq-dump']):
        raise EnvironmentError('fastq-dump not found. Please verify that fastq-dump is '
                               'installed and retry.')

    to_extract = Queue()
    for f in sra_files:
        to_extract.put(f)

    threads = []
    for i in range(max_concurrent):
        threads.append(Thread(target=_extract_fastq,
                              args=([to_extract, verbose])))
        threads[i].start()

    for t in threads:
        t.join()

    # get output files
    forward = []
    reverse = []
    for f in sra_files:
        forward.append(f.replace('.sra', '_1.fastq'))
        reverse.append(f.reaplce('.sra', '_2.fastq'))

    return forward, reverse
