import boto3
import time
import string
import ssh_utils as sshutils
import sys
import os
import configparser


# instance.update() to refresh
# workflow: 1) run code w/ paramters 2) this thing does everything and uploads file onto S3
# 3) sends back to user
# TODO: currently designed to work only with the instance created by the class method!!!
class ClusterServer(object):
    """Connects to AWS instance using paramiko and a private RSA key,
    allows for the creation/manipulation of EC2 instances and executions
    of commands on the remote server"""

    # check if you need to alter these parameters
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
        # maybe keep track of all the volumes associated with it too

    def create_security_group(self, name):
        try:
            sg = self.ec2.create_security_group(GroupName=name, Description=name)
            sg.authorize_ingress(IpProtocol="tcp", CidrIp="0.0.0.0/0", FromPort=22, ToPort=22)
            sg.authorize_ingress(SourceSecurityGroupName=name)
            self.sg = sg.id
            print('created security group seqc (%s)' % sg.id)
        except self.ec2.ClientError:
            print('the cluster %s already exists!')

    # TODO catch errors in cluster configuration
    def configure_cluster(self, config_file):
        config = configparser.ConfigParser()
        config.read(config_file)
        self.keyname = config['key']['KEY_NAME']
        self.keypath = os.path.expanduser(config['key']['KEY_LOCATION'])
        self.image_id = config['defaultcluster']['NODE_IMAGE_ID']
        self.inst_type = config['defaultcluster']['NODE_INSTANCE_TYPE']
        # self.subnet = config['defaultcluster']['SUBNET_ID']
        self.zone = config['defaultcluster']['AVAILABILITY_ZONE']
        self.n_tb = config['raid']['n_tb']
        self.dir_name = config['gitpull']['dir_name']

    def create_cluster(self, clust_type):
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
        if self.inst_id is None:
            print('No instance created!')
            sys.exit(2)
        self.inst_id.reload()
        if self.inst_id.state['Name'] == 'running':
            return True
        else:
            return False

    def restart_cluster(self):
        if self.inst_id.state['Name'] == 'stopped':
            self.inst_id.start()
            self.inst_id.wait_until_running()
            print('stopped instance %s has restarted' % self.inst_id.id)
        else:
            print('instance %s is not in a stopped state!' % self.inst_id.id)

            # documenting test code
            # resp = instance.start()
            # state = resp['Name']
            # while state != 'running':
            #     time.sleep(5)
            #     instance = ec2.Instance(inst_id)
            #     resp = instance.start()
            #     state = resp['Name']
            # nice function here
            # instance.wait_until_running()
            # print('stopped instance %s has restarted' % self.inst_id.id)

            # def stop_cluster(self, inst_id):
            # instance = self.ec2.Instance(inst_id)

    def stop_cluster(self):
        if self.cluster_is_running():
            self.inst_id.stop()
            self.inst_id.wait_until_stopped()
            print('instance %s is now stopped' % self.inst_id)
        else:
            print('instance %s is not running!' % self.inst_id)
            # resp = instance.stop()
            # stopping seems to be instantaneous, see if you want something like this here
            # if resp['StoppingInstances'][0]['CurrentState']['Name'] == 'stopped':
            #     print('Instance %s is now stopped' %inst_id)
            # instance.wait_until_stopped()
            # print('Instance %s is now stopped' % inst_id)

            # def terminate_cluster(self, inst_id):
            #     instance = self.ec2.Instance(inst_id)
            #     resp = instance.terminate()
            #     instance.wait_until_terminated()
            # terminating seems to be instantaneous, see if you want something like this here
            # if resp['Terminating Instances'][0]['CurrentState']['Name'] == 'terminated':
            #     print('Instance %s is now stopped' %inst_id)
            # maybe also do some error handling
            # print('Instance %s has successfully terminated' % inst_id)

    def terminate_cluster(self):
        if self.cluster_is_running():
            self.inst_id.terminate()
            self.inst_id.wait_until_terminated()
            print('instance %s has successfully terminated' % self.inst_id)
        else:
            print('instance %s is not running!' % self.inst_id)

    # should test out this code
    def create_volume(self, vol_size):
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

    # resp['Attachments']['DeleteOnTermination'] = True
    # def attach_volume(self, inst_id, vol_id, dev_id):
    def attach_volume(self, vol_id, dev_id):
        """attaches a vol_id to inst_id at dev_id"""
        # instance = self.ec2.Instance(inst_id)
        vol = self.ec2.Volume(vol_id)
        # instance.attach_volume(VolumeId=vol_id, Device=dev_id)
        self.inst_id.attach_volume(VolumeId=vol_id, Device=dev_id)
        # vol_state = vol.attachments[0]['State']
        # while vol_state != 'attached':
        while vol.state != 'in-use':
            time.sleep(3)
            vol.reload()
            # vol_state = vol.attachments[0]['State']
        resp = self.inst_id.modify_attribute(
            BlockDeviceMappings=[{'DeviceName': dev_id, 'Ebs': {'VolumeId': vol.id, 'DeleteOnTermination': True}}])
        if resp['ResponseMetadata']['HTTPStatusCode'] != 200:
            print('Something went wrong modifying the attribute of the Volume!')
        # deal with error somehow
        print('volume %s attached to %s at %s' % (vol_id, self.inst_id.id, dev_id))

    # should test this code out
    # def create_raid(self, inst_id, vol_size):
    def connect_server(self):
        ssh_server = sshutils.SSHServer(self.inst_id.id, self.keypath)
        print('connecting to instance %s...' % self.inst_id.id)
        ssh_server.connect()
        if ssh_server.is_connected():
            print('connection successful!')
        self.serv = ssh_server

    def create_raid(self, vol_size):
        dev_base = "/dev/xvd"
        alphabet = string.ascii_lowercase[5:]  # starts at f
        dev_names = []
        for i in range(int(self.n_tb)):
            print("creating volume %s of %s..." % (i + 1, self.n_tb))
            vol_id = self.create_volume(vol_size)
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

    def git_pull(self):
        # TODO replace this with public stuff
        # works with normal public git repository
        # install seqc on AMI to simplify
        if not self.dir_name.endswith('/'):
            self.dir_name += '/'
        folder = self.dir_name
        self.serv.exec_command("sudo mkdir %s" % folder)
        location = folder + "seqc.tar.gz"
        self.serv.exec_command(
            'sudo curl -H "Authorization: token a22b2dc21f902a9a97883bcd136d9e1047d6d076" -L '
            'https://api.github.com/repos/ambrosejcarr/seqc/tarball | sudo tee %s > /dev/null' % location)
        # implement some sort of ls grep check system here
        self.serv.exec_command('sudo pip3 install %s' % location)
