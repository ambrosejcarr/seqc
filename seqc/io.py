from glob import glob
import os
import ftplib
from functools import partial
from multiprocessing import Process, Pool
import configparser
from queue import Queue, Empty
from subprocess import Popen, check_output, PIPE
from itertools import zip_longest
import boto3
import logging
import requests
import seqc

# turn off boto3 non-error logging, otherwise it logs tons of spurious information
logging.getLogger('botocore').setLevel(logging.CRITICAL)
logging.getLogger('boto3').setLevel(logging.CRITICAL)


class S3:
    """A series of methods to upload and download files from amazon s3"""

    @staticmethod
    def download_file(bucket: str, key: str, fout: str=None, overwrite: bool=False):
        """download the file key located in bucket, sending output to filename fout

        :param overwrite:
        :param fout:
        :param key:
        :param bucket:
        """

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

    # todo
    # implement glob-based filtering, return all downloaded filenames
    @classmethod
    def download_files(cls, bucket, key_prefix, output_prefix='./', cut_dirs=0,
                       overwrite=False):
        """
        recursively download objects from amazon s3
        recursively downloads objects from bucket starting with key_prefix.
        If desired, removes a number of leading directories equal to cut_dirs. Finally,
        can overwrite files lying in the download path if overwrite is True.

        :param overwrite:
        :param cut_dirs:
        :param output_prefix:
        :param key_prefix:
        :param bucket:
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
        """upload filename to aws at s3://bucket/key/filename
        :param key:
        :param bucket:
        :param filename:
        """

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

        :param cut_dirs:
        :param key_prefix:
        :param bucket:
        :param file_prefix:
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
        :param key_prefix:
        :param bucket:
        """
        s3 = boto3.resource('s3')
        bucket = s3.Bucket(bucket)
        keys = [k.key for k in bucket.objects.all() if k.key.startswith(key_prefix)]
        return keys

    @staticmethod
    def remove_file(bucket: str, key: str):
        """delete AWS s3 file at s3://bucket/key
        :type bucket: str
        :param key: str
        :param bucket:
        """
        client = boto3.client('s3')
        _ = client.delete_object(Bucket=bucket, Key=key)

    @classmethod
    def remove_files(cls, bucket: str, key_prefix: str):
        """
        :param bucket: str
        :param key_prefix: str
        :return:
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


class GEO:
    """
    Group of methods for downloading files from NCBI GEO
    """

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
        :returns: list of downloaded files
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
        :param sra_files:
        :param max_concurrent:
        :param working_directory:
        :param verbose:
        :param paired_end:
        :param clobber:
        :return:
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

        # if the access token isn't provided, look for it in config
        if not access_token:
            config = configparser.ConfigParser()
            config.read(os.path.expanduser('~/.seqc/config'))
            access_token = config['BaseSpaceToken']['base_space_token']
            print(access_token)

        response = requests.get('https://api.basespace.illumina.com/v1pre3/samples/' +
                                sample_id +
                                '/files?Extensions=gz&access_token=' +
                                access_token)
        print(response)

        # check that a valid request was sent
        if response.status_code != 200:
            raise ValueError('Invalid access_token or sample_id. BaseSpace could not '
                             'find the requested data. Response status code: %d'
                             % response.status_code)

        data = response.json()

        func = partial(cls._download_basespace_content, data['Response']['Items'],
                       access_token, dest_path)
        seqc.log.info('BaseSpace API link provided, downloading files from '
                          'BaseSpace.')
        with Pool(len(data['Response']['Items'])) as pool:
            pool.map(func, range(len(data['Response']['Items'])))

        # get downloaded forward and reverse fastq files
        filenames = [f['Name'] for f in data['Response']['Items']]
        if 'mars' not in platform:
            barcode_fastq = [dest_path + '/' + f for f in filenames if '_R1_' in f]
            genomic_fastq = [dest_path + '/' + f for f in filenames if '_R2_' in f]
        else:
            genomic_fastq = [dest_path + '/' + f for f in filenames if '_R1_' in f]
            barcode_fastq = [dest_path + '/' + f for f in filenames if '_R2_' in f]

        return barcode_fastq, genomic_fastq
