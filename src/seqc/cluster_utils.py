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

#TODO check if any errors in the arguments here
# maybe keep track of all the volumes associated with it too
class ClusterServer(object):
    """Connects to AWS instance using paramiko and a private RSA key,
    allows for the creation/manipulation of EC2 instances and executions
    of commands on the remote server"""

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
        if name == None:
            name = 'seqc_' + str(random.randint(1, int(1e12)))
            print('no name assigned, chose %s' %name)
        try:
            sg = self.ec2.create_security_group(GroupName=name, Description=name)
            sg.authorize_ingress(IpProtocol="tcp", CidrIp="0.0.0.0/0", FromPort=22, ToPort=22)
            sg.authorize_ingress(SourceSecurityGroupName=name)
            self.sg = sg.id
            print('created security group %s (%s)' % (name,sg.id))
        except ClientError:
            print('the cluster %s already exists!' %name)
            sys.exit(2)

    # todo catch errors in cluster configuration
    def configure_cluster(self, config_file):
        """configures the newly created cluster according to aws.config"""
        config = configparser.ConfigParser()
        config.read(config_file)
        template = config['global']['DEFAULT_TEMPLATE']
        self.keyname = config['key']['KEY_NAME']
        self.keypath = os.path.expanduser(config['key']['KEY_LOCATION'])
        self.image_id = config[template]['NODE_IMAGE_ID']
        self.inst_type = config[template]['NODE_INSTANCE_TYPE']
        if template == 'c4':
            self.subnet = config[template]['SUBNET_ID']
        self.zone = config[template]['AVAILABILITY_ZONE']
        self.n_tb = config['raid']['n_tb']
        self.dir_name = config['gitpull']['dir_name']
        self.aws_id = config['aws_info']['AWS_ACCESS_KEY_ID']
        self.aws_key = config['aws_info']['AWS_SECRET_ACCESS_KEY']

    def create_cluster(self):
        """creates a new AWS cluster with specifications from aws.config"""
        if 'c4' in self.inst_type:
            if not self.subnet:
                print('You must specify a subnet-id for C4 instances!')
                sys.exit(2)
            else:
                clust = self.ec2.create_instances(ImageId=self.image_id, MinCount=1, MaxCount=1,
                                                  KeyName=self.keyname, InstanceType=self.inst_type,
                                                  Placement={'AvailabilityZone': self.zone},
                                                  SecurityGroupIds=[self.sg], SubnetId=self.subnet)
        elif 'c3' in self.inst_type:
            clust = self.ec2.create_instances(ImageId=self.image_id, MinCount=1, MaxCount=1,
                                              KeyName=self.keyname, InstanceType=self.inst_type,
                                              Placement={'AvailabilityZone': self.zone},
                                              SecurityGroupIds=[self.sg])
        instance = clust[0]
        print('created new instance %s' % instance)
        instance.wait_until_exists()
        instance.wait_until_running()
        self.inst_id = instance

    def cluster_is_running(self):
        """checks whether a cluster is running"""
        if self.inst_id is None:
            print('No instance created!')
            sys.exit(2)
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
            print('stopped instance %s has restarted' % self.inst_id.id)
        else:
            print('instance %s is not in a stopped state!' % self.inst_id.id)

    def stop_cluster(self):
        """stops a running cluster"""
        if self.cluster_is_running():
            self.inst_id.stop()
            self.inst_id.wait_until_stopped()
            print('instance %s is now stopped' % self.inst_id)
        else:
            print('instance %s is not running!' % self.inst_id)

    def create_volume(self):
        """creates a volume of size vol_size and returns the volume's id"""
        vol_size = 50 #1024 --> just testing
        vol = self.ec2.create_volume(Size=vol_size, AvailabilityZone=self.zone,
                                     VolumeType='standard')
        vol_id = vol.id
        vol_state = vol.state
        while vol_state != 'available':
            time.sleep(3)
            vol.reload()
            vol_state = vol.state
        print('volume %s created successfully' % vol_id)
        return vol_id

    # todo deal with volume creation errors
    def attach_volume(self, vol_id, dev_id):
        """attaches a vol_id to inst_id at dev_id"""
        vol = self.ec2.Volume(vol_id)
        self.inst_id.attach_volume(VolumeId=vol_id, Device=dev_id)
        while vol.state != 'in-use':
            time.sleep(3)
            vol.reload()
        resp = self.inst_id.modify_attribute(
            BlockDeviceMappings=[{'DeviceName': dev_id, 'Ebs': {'VolumeId': vol.id, 'DeleteOnTermination': True}}])
        if resp['ResponseMetadata']['HTTPStatusCode'] != 200:
            print('Something went wrong modifying the attribute of the Volume!')
            sys.exit(2)
        print('volume %s attached to %s at %s' % (vol_id, self.inst_id.id, dev_id))

    def connect_server(self):
        """connects to the aws instance"""
        ssh_server = seqc.ssh_utils.SSHServer(self.inst_id.id, self.keypath)
        print('connecting to instance %s...' % self.inst_id.id)
        ssh_server.connect()
        if ssh_server.is_connected():
            print('connection successful!')
        self.serv = ssh_server

    def create_raid(self):
        """creates a raid array of a specified number of volumes on /data"""
        dev_base = "/dev/xvd"
        alphabet = string.ascii_lowercase[5:]  # starts at f
        dev_names = []
        for i in range(int(self.n_tb)):
            print("creating volume %s of %s..." % (i + 1, self.n_tb))
            vol_id = self.create_volume()
            dev_id = dev_base + alphabet[i]
            dev_names.append(dev_id)
            self.attach_volume(vol_id, dev_id)
        print("successfully attached %s TB in %s volumes!" % (self.n_tb, self.n_tb))
        print("creating logical RAID device...")
        all_dev = ' '.join(dev_names)

        # check if this sleep is necessary for successful execution of mdadm function
        time.sleep(5)
        self.serv.exec_command(
            "sudo mdadm --create --verbose /dev/md0 --level=0 --name=my_raid --raid-devices=%s %s" % (
                self.n_tb, all_dev))
        # grep for md0 as a check here in dev --> function
        out, err = self.serv.exec_command('ls /dev | grep "md0"')
        if not out:
            print('error with mdadm function')
            print(err)
            sys.exit(2)

        self.serv.exec_command("sudo mkfs.ext4 -L my_raid /dev/md0")
        self.serv.exec_command("sudo mkdir -p /data")
        self.serv.exec_command("sudo mount LABEL=my_raid /data")

        # checking for proper raid mounting
        output, err = self.serv.exec_command('df -h | grep "/dev/md0"')
        if output:
            if output[0].endswith('/data'):
                print("successfully created RAID array in /data!")
        else:
            print("error occurred in mounting RAID array to /data")
            print(err)
            sys.exit(2)

    def git_pull(self):
        """installs the SEQC directory in /data/software"""
        # todo replace this with git clone once repo is public
        # works with normal public git repository
        # install seqc on AMI to simplify

        if not self.dir_name.endswith('/'):
            self.dir_name += '/'
        folder = self.dir_name
        print('installing seqc.tar.gz...')
        self.serv.exec_command("sudo mkdir %s" % folder)
        self.serv.exec_command("sudo chown -c ubuntu /data")
        self.serv.exec_command("sudo chown -c ubuntu %s" % folder)

        location = folder + "seqc.tar.gz"
        # todo | get rid of "nuke_sc" branch here, just for testing
        self.serv.exec_command(
            'curl -H "Authorization: token a22b2dc21f902a9a97883bcd136d9e1047d6d076" -L '
            'https://api.github.com/repos/ambrosejcarr/seqc/tarball/nuke_sc | sudo tee %s > /dev/null' % location)
        # todo implement some sort of ls grep check system here
        self.serv.exec_command('sudo pip3 install %s' % location)
        print('successfully installed seqc.tar.gz in %s on the cluster!' %folder)

    def set_credentials(self):
        self.serv.exec_command('aws configure set aws_access_key_id %s' %self.aws_id)
        self.serv.exec_command('aws configure set aws_secret_access_key %s' %self.aws_key)
        self.serv.exec_command('aws configure set region %s' %self.zone[:-1])

    def cluster_setup(self, name):
        print('setting up cluster %s...' % name)
        config_file = '/'.join(seqc.__file__.split('/')[:-3]) + '/src/plugins/aws.config'
        self.configure_cluster(config_file)
        self.create_security_group(name)
        self.create_cluster()
        self.connect_server()
        self.create_raid()
        self.git_pull()
        self.set_credentials()
        print('sucessfully set up the remote cluster environment!')

# todo test this function
def terminate_cluster(instance_id):
    """terminates a running cluster"""
    ec2 = boto3.resource('ec2')
    instance = ec2.Instance(instance_id)
    if instance.state['Name'] == 'running':
        gname = instance.security_groups[0]['GroupName']
        sg = instance.security_groups[0]['GroupId']
        instance.terminate()
        instance.wait_until_terminated()
        print('instance %s has successfully terminated' % instance_id)

        print('removing security group %s...' %gname)
        boto3.client('ec2').delete_security_group(GroupName=gname,GroupId=sg)
        print('termination complete!')
    else:
        print('instance %s is not running!' % instance_id)

def email_user(attachment, email_body, email_address: str) -> None:
    """
    todo document me!

    args:
    -----

    returns:
    --------
    None
    """
    seqc.log.exception()
    email_args = ['mutt', '-a', attachment, '-s', 'Remote Process', '--', email_address]
    message = Popen(['echo',email_body],stdout=PIPE,stderr=PIPE)
    email_process = Popen(email_args, stdin=message.stdout,stdout=PIPE,stderr=PIPE)
    email_process.communicate()

def upload_results(output_prefix: str, email_address: str, aws_upload_key) -> None:
    """
    todo document me!

    args:
    -----

    returns:
    --------
    None
    """
    prefix, directory = seqc.core.fix_output_paths(output_prefix)

    # todo
    # at some point, script should write a metadata summary and upload that in place of
    # the alignment summary

    samfile = directory + 'Aligned.out.sam'
    h5_archive = prefix + '.h5'
    merged_fastq = directory + 'merged.fastq'
    counts = prefix + '_sp_counts.npz'
    id_map = prefix + '_gene_id_map.p'
    summary = prefix + '_alignment_summary.txt'
    files = [samfile, h5_archive, merged_fastq, counts, id_map, summary]

    # gzip everything and upload to aws_upload_key
    archive_name = prefix + '.tar.gz'
    gzip_args = ['tar', '-czf', archive_name] + files
    gzip = Popen(gzip_args)
    gzip.communicate()
    bucket, key = seqc.io.S3.split_link(aws_upload_key)
    seqc.io.S3.upload_file(archive_name, bucket, key)

    # gzip small files for email
    attachment = prefix + '_counts_and_metadata.tar.gz'
    gzip_args = ['tar', '-czf', attachment, counts, id_map, summary]
    gzip = Popen(gzip_args)
    gzip.communicate()

    # email results to user
    body = ('SEQC run complete -- see attached .npz file. The rest of the output is '
            'available as %s in your specified S3 bucket.' % archive_name)
    email_user(attachment, body, email_address)
