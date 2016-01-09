import boto3
import time
import string
import sys
import os
import configparser
import random
import seqc
from subprocess import Popen, PIPE
from botocore.exceptions import ClientError

class EC2RuntimeError(Exception):
    pass


# TODO check if any errors in the arguments here
# maybe keep track of all the volumes associated with it too
class ClusterServer(object):
    """Connects to AWS instance using paramiko and a private RSA key,
    allows for the creation/manipulation of EC2 instances and executions
    of commands on the remote server"""

    # todo many of these parameters aren't used. Either integrate or remove.
    def __init__(self, private_key=None, security_group=None,
                 image_id=None, zone=None, server=None,
                 instance_type=None, subnet_id=None):
        self.keyname = None
        self.keypath = None
        self.image_id = None
        self.inst_type = None
        self.subnet = None
        self.zone = None
        self.ec2 = boto3.resource('ec2')
        self.inst_id = None
        self.dir_name = None
        self.n_tb = None
        self.sg = None
        self.serv = None
        self.aws_id = None
        self.aws_key = None

    def create_security_group(self, name=None):
        """Creates a new security group for the cluster"""
        if name is None:
            name = 'SEQC-%07d' % random.randint(1, int(1e7))
            seqc.log.notify('No instance name provided, assigned %s.' % name)
        try:
            sg = self.ec2.create_security_group(GroupName=name, Description=name)
            sg.authorize_ingress(IpProtocol="tcp", CiIp="0.0.0.0/0", FromPort=22,
                                 ToPort=22)
            sg.authorize_ingress(SourceSecurityGroupName=name)
            self.sg = sg.id

            seqc.log.notify('Created security group %s (%s).' % (name, sg.id))
        except ClientError:
            seqc.log.notify('Instance %s already exists! Exiting.' % name)
            sys.exit(2)

    # todo catch errors in cluster configuration
    def configure_cluster(self, config_file):
        """configures the newly created cluster according to config"""
        config = configparser.ConfigParser()
        config.read(config_file)
        template = config['global']['default_template']
        self.keyname = config['key']['rsa_key_name']
        self.keypath = os.path.expanduser(config['key']['rsa_key_location'])
        self.image_id = config[template]['node_image_id']
        self.inst_type = config[template]['node_instance_type']
        if template == 'c4':
            self.subnet = config[template]['subnet_id']
        self.zone = config[template]['availability_zone']
        self.n_tb = config['raid']['n_tb']
        self.dir_name = config['gitpull']['dir_name']
        self.aws_id = config['aws_info']['aws_access_key_id']
        self.aws_key = config['aws_info']['aws_secret_access_key']

    def create_cluster(self):
        """creates a new AWS cluster with specifications from config"""
        if 'c4' in self.inst_type:
            if not self.subnet:
                ValueError('A subnet-id must be specified for C4 instances!')
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
            seqc.log.notify('Instance %s is not in a stopped state!' % self.inst_id.id)

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
        vol_size = 1024
        vol = self.ec2.create_volume(Size=vol_size, AvailabilityZone=self.zone,
                                     VolumeType='standard')
        vol_id = vol.id
        vol_state = vol.state
        while vol_state != 'available':
            time.sleep(3)
            vol.reload()
            vol_state = vol.state
        seqc.log.notify('Volume %s created successfully.' % vol_id)
        return vol_id

    # todo deal with volume creation errors
    def attach_volume(self, vol_id, dev_id):
        """attaches a vol_id to inst_id at dev_id"""
        vol = self.ec2.Volume(vol_id)
        self.inst_id.attach_volume(VolumeId=vol_id, Device=dev_id)
        # todo potential infinite loop
        while vol.state != 'in-use':
            time.sleep(3)
            vol.reload()
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
            # todo possibly infinite loop
            while status != 'attached':
                time.sleep(1)
                self.inst_id.reload()
                device_info = self.inst_id.block_device_mappings
                status = device_info[i]['Ebs']['Status']
        seqc.log.notify('Volume %s attached to %s at %s.' %
                        (vol_id, self.inst_id.id, dev_id))

    def connect_server(self):
        """connects to the aws instance"""
        ssh_server = seqc.ssh_utils.SSHServer(self.inst_id.id, self.keypath)
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

        # check if this sleep is necessary for successful execution of mdadm function
        num_retries = 30
        mdadm_failed = True
        for i in range(num_retries):
            _, res = self.serv.exec_command(
                "sudo mdadm --create --verbose /dev/md0 --level=0 --name=my_raid "
                "--raid-devices=%s %s" % (
                    self.n_tb, all_dev))
            if 'started' in ''.join(res):
                mdadm_failed = False
                break
            else:
                seqc.log.notify('Retrying sudo mdadm.')
                time.sleep(1)
        if mdadm_failed:
            EC2RuntimeError('Error creating raid array md0 with mdadm function.')

        # todo | may not need to do repeated checks
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
        # todo replace this with git clone once repo is public
        # works with normal public git repository
        # install seqc on AMI to simplify

        if not self.dir_name.endswith('/'):
            self.dir_name += '/'
        folder = self.dir_name
        seqc.log.notify('Installing SEQC on remote instance.')
        self.serv.exec_command("sudo mkdir %s" % folder)
        self.serv.exec_command("sudo chown -c ubuntu /data")
        self.serv.exec_command("sudo chown -c ubuntu %s" % folder)

        location = folder + "seqc.tar.gz"
        # todo | get rid of "nuke_sc" branch here, just for testing
        self.serv.exec_command(
            'curl -H "Authorization: token a22b2dc21f902a9a97883bcd136d9e1047d6d076" -L '
            'https://api.github.com/repos/ambrosejcarr/seqc/tarball/develop | '
            'sudo tee %s > /dev/null' % location)
        # todo implement some sort of ls grep check system here
        self.serv.exec_command('sudo pip3 install %s' % location)
        seqc.log.notify('SEQC successfully installed in %s.' % folder)

    def set_credentials(self):
        self.serv.exec_command('aws configure set aws_access_key_id %s' % self.aws_id)
        self.serv.exec_command(
            'aws configure set aws_secret_access_key %s' % self.aws_key)
        self.serv.exec_command('aws configure set region %s' % self.zone[:-1])

    def cluster_setup(self, name):
        # config_file = '/'.join(seqc.__file__.split('/')[:-3]) + '/src/plugins/config'
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
    """terminates a running cluster"""
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
        seqc.log.notify('security group %s (%s) successfully removed' % (sg_name, sg_id))
    except ClientError:
        seqc.log.notify('security group %s (%s) is still in use!' % (sg_name, sg_id))


def email_user(attachment: str, email_body: str, email_address: str) -> None:
    """
    sends an email to email address with text contents of email_body and attachment
    attached. Email will come from "Ubuntu@<ec2-instance-ip-of-aws-instance>

    args:
    -----
    attachment: the file location of the attachment to append to the email
    email_body: text to send in the body of the email
    email_address: the address to which the email should be sent.

    returns:
    --------
    None
    """
    if isinstance(email_body, str):
        email_body = email_body.encode()
    # Note: exceptions used to be logged here, but this is not the right place for it.
    email_args = ['mutt', '-a', attachment, '-s', 'Remote Process', '--', email_address]
    email_process = Popen(email_args, stdin=PIPE)
    email_process.communicate(email_body)


def upload_results(output_prefix: str, email_address: str, aws_upload_key: str) -> None:
    """
    todo document me!

    args:
    -----

    returns:
    --------
    None
    """
    prefix, directory = seqc.core.fix_output_paths(output_prefix)

    samfile = directory + 'Aligned.out.sam'
    h5_archive = prefix + '.h5'
    merged_fastq = directory + 'merged.fastq'
    counts = prefix + '_sp_counts.npz'
    id_map = prefix + '_gene_id_map.p'
    alignment_summary = prefix + '_alignment_summary.txt'
    log = 'seqc.log'
    files = [samfile, h5_archive, merged_fastq, counts, id_map, alignment_summary, log]

    # gzip everything and upload to aws_upload_key
    archive_name = prefix + '.tar.gz'
    gzip_args = ['tar', '-czf', archive_name] + files
    gzip = Popen(gzip_args)
    gzip.communicate()
    bucket, key = seqc.io.S3.split_link(aws_upload_key)
    seqc.io.S3.upload_file(archive_name, bucket, key)

    # generate a run summary and append to the email
    exp = seqc.Experiment.from_npz(counts)
    run_summary = exp.summary(alignment_summary)

    # get the name of the output file
    archive_suffix = archive_name.split('/')[-1]

    # email results to user
    body = ('SEQC RUN COMPLETE.\n\n'
            'The run log has been attached to this email and '
            'results are now available in the S3 location you specified: '
            '"%s%s"\n\n'
            'RUN SUMMARY:\n\n%s' % (aws_upload_key, archive_suffix, repr(run_summary)))
    email_user(log, body, email_address)
