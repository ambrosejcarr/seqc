import paramiko
import sys
import time
import boto3
import seqc


class SSHServer(object):

    def __init__(self, inst_id, keypath):
        ec2 = boto3.resource('ec2')
        self.instance = ec2.Instance(inst_id)
        self.key = keypath
        self.ssh = paramiko.SSHClient()

    def connect(self):
        max_attempts = 5
        self.ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
        dns = self.instance.public_dns_name
        for attempt in range(max_attempts):
            try:
                self.ssh.connect(dns, username='ubuntu', key_filename=self.key)
            # except paramiko.AuthenticationException:
            #     print('autherror')
            #     print('instance not ready for connection, sleeping...')
            #     self.instance.reload()
            #     time.sleep(30)
            # except paramiko.SSHException:
            #     print('ssherror')
            #     print('instance not ready for connection, sleeping...')
            #     self.instance.reload()
            #     time.sleep(30)
            # except FileNotFoundError:
            #     print('the key %s was not found!' %self.key)
            #     sys.exit(2)
            # except paramiko.BadHostKeyException:
            #     print('the host key %s could not be verified!' %self.key)
            #     sys.exit(2)
            except Exception as e:
                seqc.log.notify('Not yet connected, sleeping...')
                time.sleep(30)
                # continue
        # gname = self.instance.security_groups[0]['GroupName']
        # gid = self.instance.security_groups[0]['GroupId']
        # self.instance.terminate()
        # boto3.client('ec2').delete_security_group(GroupName=gname,GroupId=gid)
        # raise RuntimeError("connection failed: maximum number of unsuccessful attempts reached")

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

    def get_file(self, localfile, remotefile):
        if not self.is_connected():
            seqc.log.notify('You are not connected!')
            sys.exit(2)
        ftp = self.ssh.open_sftp()
        ftp.get(remotefile, localfile)
        ftp.close()

    def put_file(self, localfile, remotefile):
        if not self.is_connected():
            seqc.log.notify('You are not connected!')
            sys.exit(2)
        ftp = self.ssh.open_sftp()
        ftp.put(localfile, remotefile)
        print('successfully placed %s in %s!' %(localfile,remotefile))
        ftp.close()

    def exec_command(self, args):
        if not self.is_connected():
            seqc.log.notify('You are not connected!')
            sys.exit(2)
        stdin, stdout, stderr = self.ssh.exec_command(args)
        stdin.flush()
        data = stdout.read().decode().splitlines()  # response in bytes
        errs = stderr.read().decode().splitlines()
        return data, errs