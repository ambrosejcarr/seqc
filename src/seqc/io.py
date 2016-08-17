import sys
import os
import ftplib
import shlex
from glob import glob
from functools import partial
from multiprocessing import Process, Pool
from queue import Queue, Empty
from subprocess import Popen, check_output, PIPE, CalledProcessError
from itertools import zip_longest
import fcntl
import boto3
import logging
import requests
import time
from seqc import log

# turn off boto3 non-error logging, otherwise it logs tons of spurious information
logging.getLogger('botocore').setLevel(logging.CRITICAL)
logging.getLogger('boto3').setLevel(logging.CRITICAL)
logging.getLogger('requests.packages.urllib3.connectionpool').setLevel(logging.CRITICAL)
logging.getLogger('requests').setLevel(logging.CRITICAL)


class S3:
    """A series of methods to upload and download files from amazon s3"""

    @staticmethod
    def download_file(bucket: str, key: str, fout: str=None, overwrite: bool=False,
                      boto: bool=False):
        """
        download the file key located in bucket, sending output to filename fout

        :param overwrite: True if overwrite existing file
        :param fout: name of output file
        :param key: key of S3 bucket
        :param bucket: name of S3 bucket
        :param boto: True if download using boto3 (default=False, uses awscli)
        """

        if not overwrite:
            if os.path.isfile(fout):
                log.info('Skipped download of file "{fout}" because it already '
                         'exists.'.format(fout=fout))
                return fout

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
        if not boto:
            s3link = 's3://' + bucket + '/' + key
            cmd = 'aws s3 cp {s3link} {fout}'.format(s3link=s3link, fout=fout)
            download_cmd = shlex.split(cmd)
            Popen(download_cmd).wait()
        else:
            client = boto3.client('s3')
            client.download_file(bucket, key, fout)

        return fout

    @classmethod
    def download_files(cls, bucket, key_prefix, output_prefix='./', cut_dirs=0,
                       overwrite=False, boto=False, filters=None):
        """
        recursively download objects from amazon s3
        recursively downloads objects from bucket starting with key_prefix.
        If desired, removes a number of leading directories equal to cut_dirs. Finally,
        can overwrite files lying in the download path if overwrite is True.

        :param overwrite: True if overwrite existing file
        :param cut_dirs: number of leading directories to remove
        :param output_prefix: location to download file
        :param key_prefix: key of S3 bucket to download
        :param bucket: name of S3 bucket
        :param boto: use boto3 to download files from S3 (takes longer than awscli)
        :param filters: a list of file extensions to download without the period
        (ex) ['h5', 'log'] for .h5 and .log files
        :return: sorted list of file names that were downloaded
        """

        # get bucket and file names
        client = boto3.client('s3')
        keys = cls.listdir(bucket, key_prefix)

        if filters:
            to_download = []
            for f_ext in filters:
                to_download += [item for item in keys if f_ext in item]
            keys = to_download

        # output prefix needs to end in '/'
        if not output_prefix.endswith('/'):
            output_prefix += '/'

        # download data
        output_files = []
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
                    log.info('Skipped download of file "{fout}" because it already '
                             'exists.'.format(fout=fout))
                    output_files.append(fout)
                    continue

            # download file
            if not boto:
                s3link = 's3://' + bucket + '/' + k
                cmd = 'aws s3 cp {s3link} {fout}'.format(s3link=s3link, fout=fout)
                download_cmd = shlex.split(cmd)
                Popen(download_cmd).wait()
            else:
                client.download_file(bucket, k, fout)

            output_files.append(fout)
        return sorted(output_files)

    @staticmethod
    def upload_file(filename, bucket, key, boto=False):
        """upload filename to aws at s3://bucket/key/filename
        :param key: key of S3 bucket to download
        :param bucket: name of S3 bucket
        :param filename: name of file to download
        :param boto: True if download using boto3 (default=False, uses awscli)
        """

        if key.startswith('/'):
            key = key[1:]

        if key.endswith('/'):
            file_id = filename.split('/')[-1]  # get file id, stripping directories
            key += file_id

        if not os.path.isfile(filename):
            raise FileNotFoundError('file "%s" is not a valid file identifier' % filename)

        if not boto:
            s3link = 's3://' + bucket + '/' + key
            cmd = 'aws s3 cp {fname} {s3link}'.format(fname=filename, s3link=s3link)
            download_cmd = shlex.split(cmd)
            Popen(download_cmd).wait()
        else:
            client = boto3.client('s3')
            client.upload_file(filename, bucket, key)

    @staticmethod
    def upload_files(file_prefix, bucket, key_prefix, cut_dirs=True, boto=False):
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

        :param cut_dirs: number of leading directories to remove
        :param key_prefix: key of S3 bucket to upload
        :param bucket: name of S3 bucket
        :param file_prefix: name of files' prefix to upload
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
            if not boto:
                s3link = 's3://' + bucket + '/' + key
                cmd = 'aws s3 cp {fname} {s3link}'.format(fname=file_, s3link=s3link)
                download_cmd = shlex.split(cmd)
                Popen(download_cmd).wait()
            else:
                client.upload_file(file_, bucket, key)

    @staticmethod
    def listdir(bucket, key_prefix):
        """
        list all objects beginning with key_prefix
        since amazon stores objects as keys, not as a filesystem, this is synonymous to
        listing the directory defined by key_prefix
        :param key_prefix: prefix of S3 bucket to be listed
        :param bucket: name of S3 bucket
        :return: keys of items to be downloaded from S3
        """
        client = boto3.client('s3')
        objs = client.list_objects(Bucket=bucket, Prefix=key_prefix)['Contents']
        keys = [x['Key'] for x in objs]
        return keys

    @staticmethod
    def remove_file(bucket: str, key: str):
        """delete AWS s3 file at s3://bucket/key
        :param key: key of S3 bucket to remove
        :param bucket: name of S3 bucket
        """

        client = boto3.client('s3')
        _ = client.delete_object(Bucket=bucket, Key=key)

    @classmethod
    def remove_files(cls, bucket: str, key_prefix: str):
        """
        :param bucket: name of S3 bucket
        :param key_prefix: key of S3 bucket containing files to remove
        """

        keys = cls.listdir(bucket, key_prefix)
        client = boto3.client('s3')
        for k in keys:
            _ = client.delete_object(Bucket=bucket, Key=k)

    @staticmethod
    def split_link(link_or_prefix: str):
        """
        take an amazon s3 link or link prefix and return the bucket and key for use with
        S3.download_file() or S3.download_files()

        :param link_or_prefix: str, an amazon s3 link
        :returns: tuple, (bucket, key_or_prefix)
        """
        if not link_or_prefix.startswith('s3://'):
            raise ValueError('aws s3 links must start with s3://')
        link_or_prefix = link_or_prefix[5:]  # strip leading s3://
        bucket, *key_or_prefix = link_or_prefix.split('/')
        return bucket, '/'.join(key_or_prefix)

    @classmethod
    def check_links(cls, input_args: list) -> None:
        """determine if valid arguments were passed before initiating run,
        specifically whether s3 links exist

        :param input_args: list of files that should be checked
        """

        s3 = boto3.resource('s3')
        for infile in input_args:
            try:
                if infile.startswith('s3://'):
                    if not infile.endswith('/'):  # check that s3 link for file exists
                        bucket, key = cls.split_link(infile)
                        s3.meta.client.head_object(Bucket=bucket, Key=key)
                    else:
                        cmd = 'aws s3 ls ' + infile  # directory specified in s3 link
                        res = check_output(cmd.split())
                        if b'PRE ' in res:  # subdirectories present
                            raise ValueError
            except CalledProcessError:
                log.notify(
                    'Failed to access %s with "aws s3 ls", check your link' % infile)
                sys.exit(2)
            except ValueError:
                log.notify(
                    'Error: Provided s3 link "%s" does not contain the proper '
                    'input files to SEQC.' % infile)
                sys.exit(2)

    @staticmethod
    def obtain_size(item: str) -> int:
        """
        obtains the size of desired item, used to determine how much volume
        should be allocated for the remote AWS instance

        :param item: str, name of file input
        """

        cmd = 'aws s3 ls --summarize --recursive ' + item + ' | grep "Total Size"'
        obj_size = int(check_output(cmd, shell=True).decode().split()[-1])
        return obj_size


class GEO:
    """
    Group of methods for downloading files from NCBI GEO
    """

    @staticmethod
    def _ftp_login(ip, port=0, username='anonymous', password=''):
        """
        method to log in using ftf
        :param: port: name of port to use
        :param: username: default = anonymous
        :param: password: default = ''
        :return: ftp connection
        """
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
    def download_sra_file(cls, link: str, prefix: str, clobber=False, verbose=True,
                          port=0) -> str:
        """
        Downloads file from ftp server found at link into directory prefix

        :param link: ftp link to file
        :param prefix: directory into which file should be downloaded
        :param clobber: If False, will not download if a file is already present in
          prefix with the same name
        :param verbose: If True, status updates will be printed throughout file download
        :param port: Port for login. for NCBI, this should be zero (default).
        :return: downloaded filename
        """

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
    def download_srp(cls, srp: str, prefix: str, max_concurrent_dl: int, verbose=True,
                     clobber=False, port=0) -> [str]:
        """
        Download all files in an SRP experiment into directory prefix

        :param srp: the complete ftp link to the folder for the SRP experiment
        :param prefix: the name of the folder in which files should be saved
        :param max_concurrent_dl: number of processes to spawn for parallel downloading
        :param verbose: If True, status updates will be printed throughout file download
        :param clobber: If True, overwrite existing files
        :param port: Port for login. for NCBI, this should be zero (default).
        :return: list of downloaded files
        """

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
        """
        Private method. Serially extracts .SRA archive into paired fastq.

        :param sra_queue: SRA files to be extracted
        :param working_directory: str, directory for output files
        :param verbose: if True, print status updates
        :param paired_end: if True, extracts paired-ended data
        :param clobber: if True, will overwrite existing files that were previously
          extracted
        :return: None
        """

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
        """
        :param sra_files: list of str, .SRA files to be extracted
        :param max_concurrent: int, maximum number of processes to initiate
        :param working_directory: str, directory for fastq output
        :param verbose: bool, if True, prints status updates
        :param paired_end: bool, if True, generates paired-end fastq
        :param clobber: bool, if True, will overwrite existing fastq files with the
          same filestem as the archive that is being extracted.
        :return: None
        """

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

    @staticmethod
    def _download_basespace_content(item_data, access_token, dest_path, index):
        """
        obtains the content of a file requested from the BaseSpace REST API
        :param item_data:
        :param access_token: basespace access token
        :param dest_path: path to where basespace content to be downloaded
        :param index:
        """

        item = item_data[index]
        response = requests.get('https://api.basespace.illumina.com/v1pre3/files/' +
                                item['Id'] + '/content?access_token=' +
                                access_token, stream=True)
        path = dest_path + '/' + item['Path']
        with open(path, "wb") as fd:
            for chunk in response.iter_content(104857600):  # chunksize = 100MB
                fd.write(chunk)
            fd.close()

    @classmethod
    def check_sample(cls, sample_id: str, access_token: str):
        """
        checks whether provided basespace sample id is valid

        :param sample_id: id of sample to check
        :param access_token: basespace access token
        """

        # send request to see whether basespace sample exists
        response = requests.get('https://api.basespace.illumina.com/v1pre3/samples/' +
                                sample_id +
                                '/files?Extensions=gz&access_token=' +
                                access_token)

        # check that a valid request was sent
        if response.status_code != 200:
            raise ValueError('Invalid access_token or sample_id. BaseSpace could not '
                             'find the requested data. Response status code: %d'
                             % response.status_code)

    @classmethod
    def check_size(cls, sample_id, access_token):
        """
        checks size of basespace files to be downloaded.
        this function will be executed after check_sample(), so
        the API request will already have returned status code 200

        :param sample_id: sample id of basespace files to download
        :param access_token: basespace token required for sample download
         """

        # send request to obtain json metadata
        response = requests.get('https://api.basespace.illumina.com/v1pre3/samples/' +
                                sample_id +
                                '/files?Extensions=gz&access_token=' +
                                access_token)

        # obtain sizes
        total_size = 0
        resp = response.json()['Response']['Items']
        for item in resp:
            total_size += item['Size']
        return total_size

    @classmethod
    def download(
            cls, platform, sample_id: str, dest_path: str, access_token: str=None
    ) -> (list, list):
        """
        Downloads all files related to a sample from the basespace API

        :param platform: the type of data that is being downloaded
        :param sample_id: The sample id, taken directory from the basespace link for a
          sample (experiment). e.g. if the link is:
          "https://basespace.illumina.com/sample/30826030/Day0-ligation-11-17", then the
          sample_id is "30826030"
        :param access_token: a string access token that allows permission to access the
          ILLUMINA BaseSpace server and download the requested data. Access tokens can be
          obtained by (1) logging into https://developer.basespace.illumina.com, (2),
          creating a "new application", and (3) going to the credentials tab of that
          application to obtain the access token.
        :param dest_path: the location that the downloaded files should be placed.
        :returns: (list, list), forward, reverse: lists of fastq files
        """

        # validity of response will already have been checked
        response = requests.get('https://api.basespace.illumina.com/v1pre3/samples/' +
                                sample_id +
                                '/files?Extensions=gz&access_token=' +
                                access_token)
        data = response.json()

        func = partial(cls._download_basespace_content, data['Response']['Items'],
                       access_token, dest_path)
        log.info('BaseSpace API link provided, downloading files from BaseSpace.')
        with Pool(len(data['Response']['Items'])) as pool:
            pool.map(func, range(len(data['Response']['Items'])))

        # get downloaded forward and reverse fastq files
        filenames = [f['Name'] for f in data['Response']['Items']]

        # fixed location for how BaseSpace installs files
        dest_path += '/Data/Intensities/BaseCalls/'

        if 'mars' not in platform:
            barcode_fastq = [dest_path + f for f in filenames if '_R1_' in f]
            genomic_fastq = [dest_path + f for f in filenames if '_R2_' in f]
        else:
            genomic_fastq = [dest_path + f for f in filenames if '_R1_' in f]
            barcode_fastq = [dest_path + f for f in filenames if '_R2_' in f]

        return barcode_fastq, genomic_fastq


class ProcessManager:
    """
    Manages processes in the background to prevent blocking main loop.
    Processes can either be left running in the background with run_all(),
    or blocked until completion with wait_until_complete().
    """

    def __init__(self, wait=False, *args):
        """
        For sequential processes, pass individual args
        For piped processes, pass as one

        (ex)
        seqc.io.ProcessManager('df -h')
        seqc.io.ProcessManager('df -h', 'sleep 10')
        seqc.io.ProcessManager('echo hello world | grep hello')

        sample use:
        test = seqc.io.ProcessManager('df -h', 'sleep 10')
        test.run_all()
        test.wait_until_complete()  # optional, this blocks the process

        """
        self.args = args
        self.wait = wait
        self.nproc = len(args)
        self.processes = []

    @staticmethod
    def format_processes(proc: str):
        """
        :param proc: string argument to be executed
        :return cmd: Properly formatted command (list) for Popen
        """

        if '|' in proc:
            cmd = [shlex.split(item) for item in proc.split('|')]
        else:
            cmd = [shlex.split(proc)]
        return cmd

    @staticmethod
    def check_launch(proc):
        """
        checks whether a newly launched background process has any errors;
        otherwise error goes undetected unless wait_until_complete() is called.
        -> stderr is modified such that it can be read without blocking process.
        :param proc: an instance of subprocess.Popen
        :return:
        """

        fcntl.fcntl(proc.stderr.fileno(), fcntl.F_SETFL, os.O_NONBLOCK)
        error_msg = proc.stderr.read()
        if error_msg:
            # catch "error" from samtools in order to prevent premature exiting
            if 'SAM header is present' in error_msg.decode():
                pass
            else:
                raise ChildProcessError(error_msg)

    def run_background_processes(self, proc_list: list, wait=False):
        """
        Executes processes in proc_list and saves them in self.processes.
        All processes executed in this function are non-blocking.
        :param proc_list: Command to be executed (obtained from format_process).
        :param wait: return after process has finished executing
        """

        for i, cmd in enumerate(proc_list):
            if i != 0:
                proc = Popen(cmd, stdin=self.processes[i-1].stdout, stdout=PIPE,
                             stderr=PIPE)
            else:
                proc = Popen(cmd, stdout=PIPE, stderr=PIPE)
            # wait a few seconds, otherwise stderr may still return None
            time.sleep(2)
            self.check_launch(proc)
            self.processes.append(proc)
            finished = False
            if self.wait:  # if processes are chained, wait until completion w/o blocking
                while not finished:
                    return_code = proc.poll()
                    if return_code is not None:
                        finished = True

    def run_all(self):
        """
        This function must be called in order for the processes to be executed.
        All processes are non-blocking and executed in the background.
        They are independent processes and are spawned one after another.
        :return:
        """

        for arg in self.args:
            proc_list = self.format_processes(arg)
            self.run_background_processes(proc_list)

    def wait_until_complete(self):
        """
        This function blocks until all processes in self.processes are complete.
        Any error calls are raised to notify the user.
        :return: list of outputs from each executed process
        """

        output = []
        for proc in self.processes:
            out, err = proc.communicate()
            if err:
                raise ChildProcessError(err)
            out = out.decode().strip()
            output.append(out)
        return output
