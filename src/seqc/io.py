__author__ = 'ambrose'

from glob import glob
import gzip
import bz2
import os
import ftplib
import threading
from multiprocessing import Process, Queue
import socket
import threading
from queue import Empty
from subprocess import Popen, check_output, PIPE
from itertools import zip_longest
import seqc
import boto3
import logging
from pyftpdlib.handlers import FTPHandler
from pyftpdlib.servers import FTPServer
from pyftpdlib.authorizers import DummyAuthorizer
import requests
from multiprocessing import Pool
from functools import partial


# turn off boto3 non-error logging, otherwise it logs tons of spurious information
logging.getLogger('botocore').setLevel(logging.CRITICAL)
logging.getLogger('boto3').setLevel(logging.CRITICAL)


class S3:
    """A series of methods to upload and download files from amazon s3"""

    @staticmethod
    def download_file(bucket: str, key: str, fout: str=None, overwrite: bool=False):
        """download the file key located in bucket, sending output to filename fout"""

        # check argument types
        seqc.util.check_type(bucket, str, 'bucket must be type str, not %s'
                             % type(bucket))
        seqc.util.check_type(key, str, 'key must be type str, not %s' % key)
        seqc.util.check_type(fout, str, 'fout must be type str, not %s' % key)

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
        try:
            client = boto3.client('s3')
            client.download_file(bucket, key, fout)
        except FileNotFoundError:
            raise FileNotFoundError('No file was found at the specified s3 location: '
                                    '"%s".' % bucket + '/' + key)

        return fout

    # todo return all downloaded filenames
    # todo implement glob-based filtering
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
        keys = cls.listdir(bucket, key_prefix)

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
    def listdir(bucket, key_prefix):
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
        keys = cls.listdir(bucket, key_prefix)
        client = boto3.client('s3')
        for k in keys:
            _ = client.delete_object(Bucket=bucket, Key=k)

    @staticmethod
    def split_link(link_or_prefix):
        """
        take an amazon s3 link or link prefix and return the bucket and key for use with
        S3.download_file() or S3.download_files()

        args:
        -----
        link_or_prefix: an amazon s3 link, e.g. s3://dplab-data/genomes/mm38/chrStart.txt
          link prefix e.g. s3://dplab-data/genomes/mm38/

        returns:
        --------
        bucket: the aws bucket used in the link. Above, this would be dplab-data
        key_or_prefix: the aws key or prefix provided in link_or_prefix. for the above
          examples, either genomes/mm38/chrStart.txt (link) or genomes/mm38/ (prefix)
        """
        if not link_or_prefix.startswith('s3://'):
            raise ValueError('aws s3 links must start with s3://')
        link_or_prefix = link_or_prefix[5:]  # strip leading s3://
        bucket, *key_or_prefix = link_or_prefix.split('/')
        return bucket, '/'.join(key_or_prefix)


class GEO:

    @staticmethod
    def _ftp_login(ip, port=0, username='anonymous', password=''):
        ftp = ftplib.FTP()
        ftp.connect(ip, port)
        ftp.login(user=username, passwd=password)
        return ftp

    @classmethod
    def _download_sra_file(cls, link_queue, prefix, clobber=False, verbose=True, port=0):
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

            ftp = cls._ftp_login(ip, port)
            ftp.cwd(path)
            with open(prefix + file_name, 'wb') as fout:
                if verbose:
                    print('beginning download of file: "%s"' % link.split('/')[-1])
                ftp.retrbinary('RETR %s' % file_name, fout.write)
                if verbose:
                    print('download of file complete: "%s"' % link.split('/')[-1])
            ftp.close()

    @classmethod
    def download_sra_file(cls, link, prefix, clobber=False, verbose=True, port=0):

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
                return

        ftp = cls._ftp_login(ip, port)
        ftp.cwd(path)
        with open(prefix + file_name, 'wb') as fout:
            if verbose:
                print('beginning download of file: "%s"' % link.split('/')[-1])
            ftp.retrbinary('RETR %s' % file_name, fout.write)
            if verbose:
                print('download of file complete: "%s"' % link.split('/')[-1])
        ftp.close()

        return prefix + file_name.split('/')[-1]


    @classmethod
    def download_srp(cls, srp, prefix, max_concurrent_dl, verbose=True, clobber=False,
                     port=0):
        """download all files in an SRP experiment into directory 'prefix'."""

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

        ftp = cls._ftp_login(ip, port)
        files = []
        ftp.cwd(path)
        dirs = ftp.nlst()  # all SRA experiments are nested in directories of the SRP
        for d in dirs:
            files.append('%s/%s.sra' % (d, d))
        ftp.close()

        if not files:
            raise ValueError('no files found in ftp directory: "%s"' % path)

        # create set of links
        if not srp.endswith('/'):
            srp += '/'

        for_download = Queue()
        for f in files:
            for_download.put(srp + f)

        processes = []
        for i in range(max_concurrent_dl):
            processes.append(Process(target=cls._download_sra_file,
                                  args=([for_download, prefix, clobber, verbose])))

            processes[i].start()

        for t in processes:
            t.join()

        # get output files
        output_files = []
        for f in files:
            output_files.append(prefix + f.split('/')[-1])

        return output_files

    @staticmethod
    def _extract_fastq(sra_queue, working_directory, verbose=True, paired_end=False,
                       clobber=False):

        while True:
            try:
                file_ = sra_queue.get_nowait()
            except Empty:
                break

            if not clobber:
                if paired_end:
                    if all([os.path.isfile(file_.replace('.sra', '_1.fastq')),
                            os.path.isfile(file_.replace('.sra', '_2.fastq'))]):
                        continue
                else:
                    if os.path.isfile(file_.replace('.sra', '.fastq')):
                        continue

            # extract file
            if verbose:
                print('beginning extraction of file: "%s"' % file_)
            p = Popen(['fastq-dump', '--split-3', '--outdir', working_directory, file_],
                      stderr=PIPE, stdout=PIPE)
            _, err = p.communicate()
            if err:
                raise ChildProcessError(err)
            if verbose:
                print('extraction of file complete: "%s"' % file_)

    @classmethod
    def extract_fastq(cls, sra_files, max_concurrent, working_directory='.',
                      verbose=True, paired_end=False, clobber=False):
        """requires fastq-dump from sra-tools"""

        # check that fastq-dump exists
        if not check_output(['which', 'fastq-dump']):
            raise EnvironmentError(
                'fastq-dump not found. Please verify that fastq-dump is installed.')

        to_extract = Queue()
        for f in sra_files:
            to_extract.put(f)

        processes = []
        for i in range(max_concurrent):
            processes.append(Process(
                target=cls._extract_fastq,
                args=([to_extract, working_directory, verbose, paired_end, clobber])
            ))
            processes[i].start()

        for t in processes:
            t.join()

        # get output files
        if paired_end:
            forward = [f.replace('.sra', '_1.fastq') for f in sra_files]
            reverse = [f.replace('.sra', '_2.fastq') for f in sra_files]
            return forward, reverse
        else:
            forward = [f.replace('.sra', '.fastq') for f in sra_files]
            return forward


class BaseSpace:

    @classmethod
    def download_fastq(cls, sample_id: str, access_token: str, dest_path: str)\
            -> (list, list):
        """
        Downloads all files related to a sample from the basespace API

        args:
        -----
        sample_id: The sample id, taken directory from the basespace link for a
         sample (experiment). e.g. if the link is:
         "https://basespace.illumina.com/sample/30826030/Day0-ligation-11-17", then the
         sample_id is "30826030"
        access_token: a string access token that allows permission to access the ILLUMINA
         BaseSpace server and download the requested data. Access tokens can be obtained
         by (1) logging into https://developer.basespace.illumina.com, (2), creating a
         "new application", and (3) going to the credentials tab of that application to
         obtain the access token.
        dest_path: the location that the downloaded files should be placed.

        returns:
        forward, reverse: lists of fastq files
        """

        # check types
        seqc.util.check_type(sample_id, str, 'sample_id')
        seqc.util.check_type(access_token, str, 'access_token')
        seqc.util.check_type(dest_path, str, 'dest_path')

        response = requests.get('https://api.basespace.illumina.com/v1pre3/samples/' +
                                sample_id +
                                '/files?Extensions=gz&access_token=' +
                                access_token)

        # check that a valid request was sent
        if response.status_code != 200:
            raise ValueError('Invalid access_token or sample_id. BaseSpace could not '
                             'find the requested data. Response status code: %d'
                             % response.status_code)

        data = response.json()

        func = partial(cls._download_content, data['Response']['Items'], access_token,
                       dest_path)
        seqc.log.info('BaseSpace API link provided, downloading files from BaseSpace.')
        with Pool(len(data['Response']['Items'])) as pool:
            pool.map(func, range(len(data['Response']['Items'])))

        # get downloaded forward and reverse fastq files
        filenames = [f['Name'] for f in data['Response']['Items']]
        forward_fastq = [dest_path + '/' + f for f in filenames if '_R1_' in f]
        reverse_fastq = [dest_path + '/' + f for f in filenames if '_R2_' in f]

        return forward_fastq, reverse_fastq


    @staticmethod
    def _download_content(item_data, access_token, dest_path, index):
        """gets the content of a file requested from the BaseSpace REST API."""
        item = item_data[index]
        response = requests.get('https://api.basespace.illumina.com/v1pre3/files/' +
                                item['Id'] + '/content?access_token=' +
                                access_token, stream=True)
        path = dest_path + '/' + item['Path']
        with open(path, "wb") as fd:
            for chunk in response.iter_content(104857600):  # chunksize = 100MB
                fd.write(chunk)
            fd.close()


def open_file(filename):
    if filename.endswith('.gz'):
        return gzip.open(filename, 'rt')
    elif filename.endswith('.bz2'):
        return bz2.open(filename, 'rt')
    else:
        return open(filename)


class DummyFTPClient(threading.Thread):
    """A threaded FTP server used for running tests.

    This is basically a modified version of the FTPServer class which
    wraps the polling loop into a thread.

    The instance returned can be used to start(), stop() and
    eventually re-start() the server.

    The instance can also launch a client using ftplib to navigate and download files.
    it will serve files from home.
    """
    handler = FTPHandler
    server_class = FTPServer

    def __init__(self, addr=None, home=None):

        try:
            host = socket.gethostbyname('localhost')
        except socket.error:
            host = 'localhost'

        threading.Thread.__init__(self)
        self.__serving = False
        self.__stopped = False
        self.__lock = threading.Lock()
        self.__flag = threading.Event()
        if addr is None:
            addr = (host, 0)

        if not home:
            home = os.getcwd()

        authorizer = DummyAuthorizer()
        authorizer.add_anonymous(home, perm='erl')
        # authorizer.add_anonymous(home, perm='elr')
        self.handler.authorizer = authorizer
        # lower buffer sizes = more "loops" while transfering data
        # = less false positives
        self.handler.dtp_handler.ac_in_buffer_size = 4096
        self.handler.dtp_handler.ac_out_buffer_size = 4096
        self.server = self.server_class(addr, self.handler)
        self.host, self.port = self.server.socket.getsockname()[:2]
        self.client = None

    def __repr__(self):
        status = [self.__class__.__module__ + "." + self.__class__.__name__]
        if self.__serving:
            status.append('active')
        else:
            status.append('inactive')
        status.append('%s:%s' % self.server.socket.getsockname()[:2])
        return '<%s at %#x>' % (' '.join(status), id(self))

    def generate_local_client(self):
        self.client = ftplib.FTP()
        self.client.connect(self.host, self.port)
        self.client.login()
        return self.client

    @property
    def running(self):
        return self.__serving

    def start(self, timeout=0.001):
        """Start serving until an explicit stop() request.
        Polls for shutdown every 'timeout' seconds.
        """
        if self.__serving:
            raise RuntimeError("Server already started")
        if self.__stopped:
            # ensure the server can be started again
            DummyFTPClient.__init__(self, self.server.socket.getsockname(), self.handler)
        self.__timeout = timeout
        threading.Thread.start(self)
        self.__flag.wait()

    def run(self):
        logging.basicConfig(filename='testing.log', level=logging.DEBUG)
        self.__serving = True
        self.__flag.set()
        while self.__serving:
            self.__lock.acquire()
            self.server.serve_forever(timeout=self.__timeout, blocking=False)
            self.__lock.release()
        self.server.close_all()

    def stop(self):
        """Stop serving (also disconnecting all currently connected
        clients) by telling the serve_forever() loop to stop and
        waits until it does.
        """
        if not self.__serving:
            raise RuntimeError("Server not started yet")
        self.__serving = False
        self.__stopped = True
        self.join()
        self.client.close()
