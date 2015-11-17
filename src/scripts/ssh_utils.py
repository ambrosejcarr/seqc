import paramiko
import sys
import time
import boto3 #trying this here
import seqc.io_lib as iolib

class SSHServer(object):

    def __init__(self, inst_id, keypath):
        #TODO think about class args
        #private_key = ~/.ssh/eyc2120.rsa, path to RSA key
        #do you want to make this part of a separate function?
        # self.ec2 = boto3.resource('ec2')
        # self.instance = self.ec2.Instance(inst_id)
        ec2 = boto3.resource('ec2')
        self.instance = ec2.Instance(inst_id)
        self.key = keypath
        self.ssh = paramiko.SSHClient()
        # if not private_key:
        #     print('private rsa key required!')
        #     sys.exit(2)
        #if not private key, raise some kind of descriptive error

    def connect(self):
        max_attempts = 5
        self.ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
        dns = self.instance.public_dns_name
        for attempt in range(max_attempts):
            try:
                self.ssh.connect(dns, username='ubuntu', key_filename=self.key)
            except Exception as e:
                print(e)
                print('instance not ready for connection, sleeping...')
                self.instance.reload()
                time.sleep(30)
            # except TimeoutError:
            #     print('the connection timed out')
            #     time.sleep(30)
            # except TypeError:
            #     print('weird type error?')
            #     time.sleep(30)
        # raise RuntimeError("maximum number of unsuccessful attempts reached")

    def is_connected(self):
        if self.ssh.get_transport() is None:
            return False
        else:
            return True

    def disconnect(self):
        if self.is_connected():
            self.ssh.close()

    # will be using iolib's S3 class to do file downloading/uploading
    def download_files(self):
        raise NotImplementedError

        # it says you need to create and connect to a transport:
        # transport = paramiko.Transport((host, port))
        # transport.connect(username = username, pkey = mykey)
        # which is also what starcluster did, but what's the difference between that
        # and an SFTPClient? (should look up)

    def get_file(self, localfile, remotefile):
        ftp = self.ssh.open_sftp()
        ftp.get(remotefile, localfile)
        ftp.close()
        # OSError if socket is closed, can get as many files as you want

    # SFTP is more restrictive, expects full path of file location --> no wildcards

    def put_file(self, localfile, remotefile):
        ftp = self.ssh.open_sftp()
        ftp.put(localfile, remotefile)
        ftp.close()

    def exec_command(self, args):
        if not self.is_connected():
            print('you are not connected!')
            sys.exit(2)
        stdin, stdout, stderr = self.ssh.exec_command(args)
        stdin.flush() #--> doesn't seem to be absolutely necessary
        # if stderr.read():
        #     print('ERROR: %s' % stderr.read().decode())
        #     sys.exit(2)
        data = stdout.read().decode().splitlines()  # response in bytes
        errs = stderr.read().decode().splitlines()
        # for line in data:
        #     print(line)
        #option to return or print (user-def)
        return data, errs

