import time
import string
import sys
import os
import configparser
import random
from subprocess import Popen, PIPE
import shutil
import paramiko
import boto3
from botocore.exceptions import ClientError
import seqc


class EC2RuntimeError(Exception):
    pass


class VolumeCreationError(Exception):
    pass


class ClusterServer(object):
    """Connects to AWS instance using paramiko and a private RSA key,
    allows for the creation/manipulation of EC2 instances and executions
    of commands on the remote server"""

    def __init__(self):

        self.keyname = None
        self.keypath = None
        self.image_id = None
        self.inst_type = None
        self.subnet = None
        self.zone = None
        self.ec2 = boto3.resource('ec2')
        self.inst_id = None
        self.n_tb = None
        self.sg = None
        self.serv = None
        self.aws_id = None
        self.aws_key = None

    def create_security_group(self, name=None):
        """Creates a new security group for the cluster
        :param name: cluster name if provided by user
        """
        if name is None:
            name = 'SEQC-%07d' % random.randint(1, int(1e7))
            seqc.log.notify('No instance name provided, assigned %s.' % name)
        try:
            sg = self.ec2.create_security_group(GroupName=name, Description=name)
            sg.authorize_ingress(IpProtocol="tcp", CidrIp="0.0.0.0/0", FromPort=22,
                                 ToPort=22)
            sg.authorize_ingress(SourceSecurityGroupName=name)
            self.sg = sg.id

            seqc.log.notify('Created security group %s (%s).' % (name, sg.id))
        except ClientError:
            seqc.log.notify('Instance %s already exists! Exiting.' % name)
            sys.exit(2)

    def configure_cluster(self, config_file):
        """configures the newly created cluster according to config
        :param config_file: /path/to/seqc/config
        """
        config = configparser.ConfigParser()
        config.read(config_file)
        template = config['global']['default_template']
        self.keyname = config['key']['rsa_key_name']
        self.keypath = os.path.expanduser(config['key']['rsa_key_location'])
        self.image_id = config[template]['node_image_id']
        self.inst_type = config[template]['node_instance_type']
        self.subnet = config['c4']['subnet_id']
        self.zone = config[template]['availability_zone']
        self.n_tb = config['raid']['n_tb']
        self.aws_id = config['aws_info']['aws_access_key_id']
        self.aws_key = config['aws_info']['aws_secret_access_key']

    def create_cluster(self):
        """creates a new AWS cluster with specifications from config"""
        if 'c4' in self.inst_type:
            if not self.subnet:
                raise ValueError('A subnet-id must be specified for C4 instances!')
            else:
                clust = self.ec2.create_instances(ImageId=self.image_id, MinCount=1,
                                                  MaxCount=1,
                                                  KeyName=self.keyname,
                                                  InstanceType=self.inst_type,
                                                  Placement={
                                                      'AvailabilityZone': self.zone},
                                                  SecurityGroupIds=[self.sg],
                                                  SubnetId=self.subnet)
        elif 'c3' in self.inst_type:
            clust = self.ec2.create_instances(ImageId=self.image_id, MinCount=1,
                                              MaxCount=1,
                                              KeyName=self.keyname,
                                              InstanceType=self.inst_type,
                                              Placement={'AvailabilityZone': self.zone},
                                              SecurityGroupIds=[self.sg])
        else:
            raise ValueError('self.inst_type must be a c3 or c4 instance')
        instance = clust[0]
        seqc.log.notify('Created new instance %s. Waiting until instance is running' %
                        instance)
        instance.wait_until_exists()
        instance.wait_until_running()
        seqc.log.notify('Instance %s now running.' % instance)
        self.inst_id = instance

    def cluster_is_running(self):
        """checks whether a cluster is running"""
        if self.inst_id is None:
            raise EC2RuntimeError('No inst_id assigned. Instance was not successfully '
                                  'created!')
        self.inst_id.reload()
        if self.inst_id.state['Name'] == 'running':
            return True
        else:
            return False

    def restart_cluster(self):
        """restarts a stopped cluster"""
        if self.inst_id.state['Name'] == 'stopped':
            self.inst_id.start()
            self.inst_id.wait_until_running()
            seqc.log.notify('Stopped instance %s has restarted.' % self.inst_id.id)
        else:
            seqc.log.notify('Instance %s is not in a stopped state!' %
                            self.inst_id.id)

    def stop_cluster(self):
        """stops a running cluster"""
        if self.cluster_is_running():
            self.inst_id.stop()
            self.inst_id.wait_until_stopped()
            seqc.log.notify('Instance %s is now stopped.' % self.inst_id)
        else:
            seqc.log.notify('Instance %s is not running!' % self.inst_id)

    def create_volume(self):
        """creates a volume of size vol_size and returns the volume's id"""
        # todo: change this back to 1024 after testing
        vol_size = 50 #1024
        vol = self.ec2.create_volume(Size=vol_size, AvailabilityZone=self.zone,
                                     VolumeType='standard')
        vol_id = vol.id
        vol_state = vol.state
        max_tries = 40
        i = 0
        while vol_state != 'available':
            time.sleep(3)
            vol.reload()
            i += 1
            if i >= max_tries:
                raise VolumeCreationError('Volume could not be created.')
            vol_state = vol.state
        seqc.log.notify('Volume %s created successfully.' % vol_id)
        return vol_id

    def attach_volume(self, vol_id, dev_id):
        """attaches a vol_id to inst_id at dev_id
        :param dev_id: where volume will be mounted
        :param vol_id: ID of volume to be attached
        """
        vol = self.ec2.Volume(vol_id)
        self.inst_id.attach_volume(VolumeId=vol_id, Device=dev_id)
        max_tries = 40
        i = 0
        while vol.state != 'in-use':
            time.sleep(.5)
            vol.reload()
            i += 1
            if i >= max_tries:
                raise VolumeCreationError('Volume could not be attached.')
        resp = self.inst_id.modify_attribute(
            BlockDeviceMappings=[
                {'DeviceName': dev_id, 'Ebs': {'VolumeId': vol.id,
                                               'DeleteOnTermination': True}}])
        if resp['ResponseMetadata']['HTTPStatusCode'] != 200:
            EC2RuntimeError('Something went wrong modifying the attribute of the Volume!')

        # wait until all volumes are attached
        device_info = self.inst_id.block_device_mappings
        for i in range(1, len(device_info)):
            status = device_info[i]['Ebs']['Status']
            i = 0
            while status != 'attached':
                time.sleep(.5)
                self.inst_id.reload()
                device_info = self.inst_id.block_device_mappings
                status = device_info[i]['Ebs']['Status']
                i += 1
                if i >= max_tries:
                    raise VolumeCreationError('All Volumes could not be attached')
        seqc.log.notify('Volume %s attached to %s at %s.' %
                        (vol_id, self.inst_id.id, dev_id))

    def connect_server(self):
        """connects to the aws instance"""
        ssh_server = SSHServer(self.inst_id.id, self.keypath)
        seqc.log.notify('Connecting to instance %s...' % self.inst_id.id)
        ssh_server.connect()
        if ssh_server.is_connected():
            seqc.log.notify('Connection successful!')
        self.serv = ssh_server

    def create_raid(self):
        """creates a raid array of a specified number of volumes on /data"""
        seqc.log.notify('Creating and attaching storage volumes.')
        dev_base = "/dev/xvd"
        alphabet = string.ascii_lowercase[5:]  # starts at f
        dev_names = []
        for i in range(int(self.n_tb)):
            seqc.log.notify("Creating volume %s of %s..." % (i + 1, self.n_tb))
            vol_id = self.create_volume()
            dev_id = dev_base + alphabet[i]
            dev_names.append(dev_id)
            self.attach_volume(vol_id, dev_id)
        seqc.log.notify("Successfully attached %s TB in %s volumes." %
                        (self.n_tb, self.n_tb))
        seqc.log.notify("Creating logical RAID device...")
        all_dev = ' '.join(dev_names)

        num_retries = 30
        mdadm_failed = True
        for i in range(num_retries):
            _, res = self.serv.exec_command(
                "sudo mdadm --create --verbose /dev/md0 --level=0 --name=my_raid "
                "--raid-devices=%s %s" % (self.n_tb, all_dev))
            if 'started' in ''.join(res):
                mdadm_failed = False
                break
            else:
                time.sleep(2)
        if mdadm_failed:
            EC2RuntimeError('Error creating raid array md0 with mdadm function.')

        ls_failed = True
        for i in range(num_retries):
            out, err = self.serv.exec_command('ls /dev | grep md0')
            if out:
                ls_failed = False
                break
            else:
                time.sleep(1)
        if ls_failed:
            EC2RuntimeError('Error creating raid array md0 with mdadm function.')
        self.serv.exec_command("sudo mkfs.ext4 -L my_raid /dev/md0")
        self.serv.exec_command("sudo mkdir -p /data")
        self.serv.exec_command("sudo mount LABEL=my_raid /data")

        # checking for proper raid mounting
        md0_failed = True
        for i in range(num_retries):
            output, err = self.serv.exec_command('df -h | grep /dev/md0')
            if output and output[0].endswith('/data'):
                seqc.log.notify("Successfully created RAID array in /data.")
                md0_failed = False
                break
            else:
                seqc.log.notify('retrying df -h')
                time.sleep(1)
        if md0_failed:
            EC2RuntimeError("Error occurred in mounting RAID array to /data! Exiting.")

    def git_pull(self):
        """installs the SEQC directory in /data/software"""
        # todo: replace this with git clone once seqc repo is public

        folder = '/data/software/'
        seqc.log.notify('Installing SEQC on remote instance.')
        self.serv.exec_command("sudo mkdir %s" % folder)
        self.serv.exec_command("sudo chown -c ubuntu /data")
        self.serv.exec_command("sudo chown -c ubuntu %s" % folder)

        location = folder + 'seqc.tar.gz'
        # todo: change this back to v0.1.6 after committing
        self.serv.exec_command(
            'curl -H "Authorization: token a22b2dc21f902a9a97883bcd136d9e1047d6d076" -L '
            'https://api.github.com/repos/ambrosejcarr/seqc/tarball/remote | '
            'sudo tee %s > /dev/null' % location)
        self.serv.exec_command('cd %s; mkdir seqc && tar -xvf seqc.tar.gz -C seqc '
                               '--strip-components 1' % folder)
        self.serv.exec_command('cd %s; sudo pip3 install -e ./' % folder + 'seqc')
        num_retries = 30
        install_fail = True
        for i in range(num_retries):
            out, err = self.serv.exec_command('process_experiment.py -h | grep RNA')
            if not out:
                time.sleep(2)
            else:
                install_fail = False
                break
        if not install_fail:
            seqc.log.notify('SEQC successfully installed in %s.' % folder)
        else:
            raise EC2RuntimeError('Error installing SEQC on the cluster.')

    def set_credentials(self):
        self.serv.exec_command('aws configure set aws_access_key_id %s' % self.aws_id)
        self.serv.exec_command(
            'aws configure set aws_secret_access_key %s' % self.aws_key)
        self.serv.exec_command('aws configure set region %s' % self.zone[:-1])

    def cluster_setup(self, name=None):
        config_file = os.path.expanduser('~/.seqc/config')
        self.configure_cluster(config_file)
        self.create_security_group(name)
        self.create_cluster()
        self.connect_server()
        self.create_raid()
        self.git_pull()
        self.set_credentials()
        seqc.log.notify('Remote instance successfully configured.')


def terminate_cluster(instance_id):
    """terminates a running cluster
    :param instance_id:
    """
    ec2 = boto3.resource('ec2')
    instance = ec2.Instance(instance_id)

    try:
        if instance.state['Name'] == 'running':
            instance.terminate()
            instance.wait_until_terminated()
            seqc.log.notify('termination complete!')
        else:
            seqc.log.notify('instance %s is not running!' % instance_id)
    except ClientError:
        seqc.log.notify('instance %s does not exist!' % instance_id)


def remove_sg(sg_id):
    ec2 = boto3.resource('ec2')
    sg = ec2.SecurityGroup(sg_id)
    sg_name = sg.group_name
    try:
        sg.delete()
        seqc.log.notify('security group %s (%s) successfully removed' % (
            sg_name, sg_id))
    except ClientError:
        seqc.log.notify('security group %s (%s) is still in use!' % (sg_name, sg_id))


def email_user(attachment: str, email_body: str, email_address: str) -> None:
    """
    sends an email to email address with text contents of email_body and attachment
    attached. Email will come from "Ubuntu@<ec2-instance-ip-of-aws-instance>

    :param attachment: the file location of the attachment to append to the email
    :param email_body: text to send in the body of the email
    :param email_address: the address to which the email should be sent.
    """
    if isinstance(email_body, str):
        email_body = email_body.encode()
    # Note: exceptions used to be logged here, but this is not the right place for it.
    email_args = ['mutt', '-a', attachment, '-s', 'Remote Process', '--', email_address]
    email_process = Popen(email_args, stdin=PIPE)
    email_process.communicate(email_body)


def upload_results(output_stem: str, email_address: str, aws_upload_key: str) -> None:
    """
    :param output_stem: specified output directory in cluster
    :param email_address: e-mail where run summary will be sent
    :param aws_upload_key: tar gzipped files will be uploaded to this S3 bucket
    """
    prefix, directory = os.path.split(output_stem)

    samfile = prefix + '/alignments/Aligned.out.sam'
    h5_archive = output_stem + '.h5'
    merged_fastq = output_stem + '_merged.fastq'
    counts = output_stem + '_read_and_count_matrices.p'
    shutil.copyfile(prefix + '/alignments/Log.final.out', output_stem +
                    '_alignment_summary.txt')
    alignment_summary = output_stem + '_alignment_summary.txt'
    log = prefix + '/seqc.log'
    files = [samfile, h5_archive, merged_fastq, counts, alignment_summary, log]

    # gzip everything and upload to aws_upload_key
    archive_name = output_stem + '.tar.gz'
    gzip_args = ['tar', '-czf', archive_name] + files
    gzip = Popen(gzip_args)
    gzip.communicate()
    bucket, key = seqc.io.S3.split_link(aws_upload_key)
    seqc.io.S3.upload_file(archive_name, bucket, key)

    # todo @AJC put this back in
    # generate a run summary and append to the email
    # exp = seqc.Experiment.from_npz(counts)
    # run_summary = exp.summary(alignment_summary)
    run_summary = ''

    # get the name of the output file
    archive_suffix = archive_name.split('/')[-1]

    # email results to user
    body = ('SEQC RUN COMPLETE.\n\n'
            'The run log has been attached to this email and '
            'results are now available in the S3 location you specified: '
            '"%s%s"\n\n'
            'RUN SUMMARY:\n\n%s' % (aws_upload_key, archive_suffix, repr(run_summary)))
    email_user(log, body, email_address)


class SSHServer(object):
    def __init__(self, inst_id, keypath):
        ec2 = boto3.resource('ec2')
        self.instance = ec2.Instance(inst_id)

        if not os.path.isfile(keypath):
            raise ValueError('ssh key not found at provided keypath: %s' % keypath)
        self.key = keypath
        self.ssh = paramiko.SSHClient()

    def connect(self):
        max_attempts = 25
        attempt = 1
        self.ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
        dns = self.instance.public_dns_name
        while True:
            try:
                self.ssh.connect(dns, username='ubuntu', key_filename=self.key)
                break
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
            except FileNotFoundError:
                seqc.log.notify('The key %s was not found!' % self.key)
                sys.exit(2)
            # except paramiko.BadHostKeyException:
            #     print('the host key %s could not be verified!' %self.key)
            #     sys.exit(2)
            except:
                seqc.log.notify('Not yet connected, sleeping (try %d of %d)' % (
                    attempt, max_attempts))
                time.sleep(4)
                attempt += 1
                if attempt > max_attempts:
                    raise

    def is_connected(self):
        if self.ssh.get_transport() is None:
            return False
        else:
            return True

    def disconnect(self):
        if self.is_connected():
            self.ssh.close()

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
        seqc.log.info('Successfully placed {local_file} in {remote_file}.'.format(
            local_file=localfile, remote_file=remotefile))
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
