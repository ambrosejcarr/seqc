import time
import sys
import os
import configparser
import random
from subprocess import Popen, PIPE, check_output
import shutil
import paramiko
import boto3
from botocore.exceptions import ClientError
from seqc.exceptions import (retry_boto_call, SpotBidError, BotoCallError,
                             VolumeCreationError, EC2RuntimeError)
from seqc.stats import ExperimentalYield
from seqc import io
from seqc import log
import socket
import logging

# turn off paramiko non-error logging
logging.getLogger('paramiko').setLevel(logging.CRITICAL)


class ClusterServer(object):
    """
    Connects to AWS instance using paramiko and a private RSA key,
    allows for the creation/manipulation of EC2 instances and executions
    of commands on the remote server
    """

    def __init__(self):

        self.keyname = None
        self.keypath = None
        self.image_id = None
        self.inst_type = None
        self.subnet = None
        self.zone = None
        self.ec2 = boto3.resource('ec2')
        self.inst_id = None
        self.sg = None
        self.serv = None
        self.aws_id = None
        self.aws_key = None
        self.spot_bid = None

    def create_security_group(self):
        """creates a new security group for the cluster"""

        def create_security_group_attempt():
            name = 'SEQC-%07d' % random.randint(1, int(1e7))
            sg = self.ec2.create_security_group(GroupName=name, Description=name)
            log.notify('Assigned instance name %s.' % name)
            sg.authorize_ingress(IpProtocol="tcp", CidrIp="0.0.0.0/0", FromPort=22,
                                 ToPort=22)
            sg.authorize_ingress(SourceSecurityGroupName=name)
            # check to make sure that security group exists
            time.sleep(2)
            client = boto3.client('ec2')
            client.describe_security_groups(GroupIds=[sg.id])
            self.sg = sg.id
            log.notify('Created security group %s (%s).' % (name, sg.id))

        retry_boto_call(
            func=create_security_group_attempt, retries=20,
            exceptions_to_catch=ClientError, delay_retry=2,
            exception_to_raise=TimeoutError("Failed to create a security group"))()

    def configure_cluster(self, config_file, aws_instance, spot_bid=None):
        """
        configures the newly created cluster according to config
        :param config_file: /path/to/seqc/config
        :param aws_instance: [c3, c4, r3] for config template
        :param spot_bid: Desired amount to bid for spot instance
        """

        config = configparser.ConfigParser()
        config.read(config_file)
        template = aws_instance
        self.keypath = str(os.path.expanduser(config['key']['path_to_rsa_key']))
        self.keyname = self.keypath.split('/')[-1].split('.')[0]
        self.image_id = config[template]['node_image_id']
        self.inst_type = config[template]['node_instance_type']
        self.subnet = config['c4']['subnet_id']
        self.zone = config[template]['availability_zone']
        self.aws_id = config['aws_info']['aws_access_key_id']
        self.aws_key = config['aws_info']['aws_secret_access_key']
        self.spot_bid = spot_bid

    def generate_launch_specs(self) -> dict:
        """
        Generates the launch specifications requested by boto3.Client
        Shared by create_cluster and create_spot_cluster.
        :return: a dict describing an instance to be launched
        """
        launch_specs = {
            'ImageId': self.image_id,
            'KeyName': self.keyname,
            'InstanceType': self.inst_type,
            'Placement': {'AvailabilityZone': self.zone},
            'SecurityGroupIds': [self.sg]}
        if 'c3' in self.inst_type:
            pass
        elif 'c4' in self.inst_type or 'r3' in self.inst_type:
            if not self.subnet:
                raise ValueError('A subnet-id must be specified for R3/C4 instances!')
            else:
                launch_specs['SubnetId'] = self.subnet
        else:
            raise ValueError('Seq-c is unable to configure the following instance type: '
                             '%s' % str(self.inst_type))
        return launch_specs

    def create_spot_cluster(self, volume_size):
        """
        launches an instance using the specified spot bid
        and cancels bid in case of error or timeout
        :param volume_size: size of volume (GB) to be attached to instance
        """

        if self.spot_bid is None:
            raise RuntimeError('Cannot create a spot instance without a spot_bid value')
        client = boto3.client('ec2')
        log.notify('Launching cluster with volume size {volume_size} GB at spot bid '
                   '{spot_bid}.'.format(volume_size=volume_size, spot_bid=self.spot_bid))
        launch_specs = self.generate_launch_specs()
        launch_specs['BlockDeviceMappings'] = [{'DeviceName': '/dev/xvdf',
                                                'Ebs': {'VolumeSize': volume_size,
                                                        'DeleteOnTermination': True}
                                                }]
        resp = client.request_spot_instances(DryRun=False, SpotPrice=str(self.spot_bid),
                                             LaunchSpecification=launch_specs)
        # check status of spot bid request
        all_resp = client.describe_spot_instance_requests()['SpotInstanceRequests']
        sec_groups = []
        for i in range(len(all_resp)):
            item = all_resp[i]
            try:
                sgid = item['LaunchSpecification']['SecurityGroups'][0]['GroupId']
                sec_groups.append(sgid)
            except KeyError:
                sec_groups.append('NA')
                continue
        idx = sec_groups.index(self.sg)
        request_id = resp['SpotInstanceRequests'][0]['SpotInstanceRequestId']
        log.notify("Starting bidding process...")
        instance_id = self._make_bid(client, idx, request_id)
        self.inst_id = self.ec2.Instance(instance_id)
        # sleep for 5s just in case boto call needs a bit more time
        time.sleep(5)
        retry_boto_call(self.wait_for_cluster)(instance_id)

    @staticmethod
    def _make_bid(client, idx, request_id):
        """
        Manages the bidding process by repeatedly asking for a spot, and handling bid
        failures

        :param client: Boto3 EC2 client2 client
        :param idx: index of the specific spot request being made
        :param request_id: Spot request id, used to shut down the bid process in case of
         problems
        :return: ID of instance that is being launched
        """
        bad_status = ['price-too-low', 'capacity-oversubscribed',
                      'capacity-not-available', 'launch-group-constraint',
                      'az-group-constraint', 'placement-group-constraint',
                      'constraint-not-fulfillable', 'schedule-expired',
                      'bad-parameters', 'system-error', 'canceled-before-fulfillment']

        def attempt_bidding():
            spot_resp = (
                client.describe_spot_instance_requests()['SpotInstanceRequests'][idx])
            status_code = spot_resp['Status']['Code']
            if status_code in bad_status:
                raise SpotBidError('Please adjust your spot bid request.')
            log.notify('The current status of your request is: {status}'.format(
                status=status_code))
            if spot_resp['State'] == 'active':
                return spot_resp
            raise BotoCallError()

        try:
            bid_response = retry_boto_call(
                func=attempt_bidding, retries=40, delay_retry=15,
                exceptions_to_catch=BotoCallError,
                exception_to_raise=SpotBidError(
                    'Timeout: spot bid could not be fulfilled'))()

            log.notify('Spot bid request was successfully fulfilled!')
            return bid_response['InstanceId']
        except SpotBidError:
            client.cancel_spot_instance_requests(
                DryRun=False, SpotInstanceRequestIds=[request_id])
            raise

    def create_cluster(self):
        """creates a new AWS cluster with specifications from config"""
        launch_specs = self.generate_launch_specs()
        launch_specs['MaxCount'] = 1
        launch_specs['MinCount'] = 1
        clust = self.ec2.create_instances(**launch_specs)
        instance = clust[0]
        self.inst_id = instance
        log.notify('Created new instance %s. Waiting until instance is running' %
                   instance)
        # sleep for 5s just in case boto call needs a bit more time
        time.sleep(5)
        retry_boto_call(self.wait_for_cluster)(instance.id)

    @staticmethod
    def wait_for_cluster(inst_id: str):
        """
        waits until newly created cluster exists and is running
        changing default waiter settings to avoid waiting forever
        :param inst_id: instance id of AWS cluster
        """

        client = boto3.client('ec2')
        exist_waiter = client.get_waiter('instance_exists')
        run_waiter = client.get_waiter('instance_running')
        run_waiter.config.delay = 10
        run_waiter.config.max_attempts = 20
        exist_waiter.config.max_attempts = 30
        exist_waiter.wait(InstanceIds=[inst_id])
        run_waiter.wait(InstanceIds=[inst_id])
        log.notify('Instance %s now running.' % inst_id)

    def is_cluster_running(self):
        """checks whether a cluster is running"""

        if self.inst_id is None:
            raise EC2RuntimeError('No inst_id assigned. Instance was not successfully '
                                  'created!')

        self.inst_id.reload()
        if (self.inst_id.state['Name'] == 'running' or
                self.inst_id.state['Name'] == 'pending'):
            return True
        else:
            return False

    def restart_cluster(self):
        """restarts a stopped cluster"""

        if self.inst_id.state['Name'] == 'stopped':
            self.inst_id.start()
            self.inst_id.wait_until_running()
            log.notify('Stopped instance %s has restarted.' % self.inst_id.id)
        else:
            log.notify('Instance %s is not in a stopped state!' % self.inst_id.id)

    def stop_cluster(self):
        """stops a running cluster"""

        if self.is_cluster_running():
            self.inst_id.stop()
            self.inst_id.wait_until_stopped()
            log.notify('Instance %s is now stopped.' % self.inst_id)
        else:
            log.notify('Instance %s is not running!' % self.inst_id)

    def create_and_attach_volume(self, vol_size, dev_id):
        """
        creates a volume of size vol_size and returns the volume's id
        :param vol_size: size(GB) of volume to be created
        :param dev_id: device where to attach the volume
        :return: volume id of newly created volume
        """

        if vol_size <= 0:
            raise ValueError("Cannot create a volume of size {vol}".format(vol=vol_size))

        def reload_vol_until(desired_state):
            vol.reload()
            if vol.state != desired_state:
                raise BotoCallError()

        # Create
        vol = self.ec2.create_volume(Size=vol_size, AvailabilityZone=self.zone,
                                     VolumeType='gp2')
        retry_boto_call(
            func=reload_vol_until, retries=40, delay_retry=3,
            exception_to_raise=VolumeCreationError('Volume could not be created.'),
            exceptions_to_catch=BotoCallError)(desired_state='available')
        log.notify('Volume %s created successfully.' % vol.id)
        # Attach
        self.inst_id.attach_volume(VolumeId=vol.id, Device=dev_id)
        retry_boto_call(
            func=reload_vol_until, retries=40, delay_retry=.5,
            exception_to_raise=VolumeCreationError('Volume could not be attached.'),
            exceptions_to_catch=BotoCallError)(desired_state='in-use')
        # Apply changes
        resp = self.inst_id.modify_attribute(
            BlockDeviceMappings=[
                {'DeviceName': dev_id,
                 'Ebs': {'VolumeId': vol.id, 'DeleteOnTermination': True}}])
        if resp['ResponseMetadata']['HTTPStatusCode'] != 200:
            EC2RuntimeError('Something went wrong modifying the attribute of the Volume!')
        # Wait for volume to be attached

        def query_device_info():
            self.inst_id.reload()
            device_info = self.inst_id.block_device_mappings
            status = device_info[1]['Ebs']['Status']
            if status != 'attached':
                raise BotoCallError()

        retry_boto_call(
            func=query_device_info, exceptions_to_catch=(IndexError, BotoCallError),
            exception_to_raise=VolumeCreationError('New volume could not be attached'),
            retries=40, delay_retry=.5,)()
        log.notify('Volume %s attached to %s at %s.' % (vol.id, self.inst_id.id, dev_id))

    def connect_server(self):
        """
        connects to the aws instance
        :return: id of newly connected AWS instance
        """

        ssh_server = SSHServer(self.inst_id.id, self.keypath)
        log.notify('Connecting to instance %s...' % self.inst_id.id)
        ssh_server.connect()
        if ssh_server.is_connected():
            log.notify('Connection successful!')
        self.serv = ssh_server
        return self.inst_id.id

    def allocate_space(self, vol_size: int, spot: bool=False, dev_id: str="/dev/xvdf"):
        """
        dynamically allocates the specified amount of space on /data
        :param spot: True if spot bid, False otherwise
        :param vol_size: size (GB) of volume to be attached to instance
        :param dev_id: location where to attach the volume
        """

        if not spot:
            log.notify("Creating volume of size %d GB..." % vol_size)
            self.create_and_attach_volume(vol_size, dev_id)  # called for side-effect
            log.notify("Successfully attached %d GB in 1 volume." % vol_size)

        if not self.serv or not self.serv.is_connected:
            self.connect_server()

        log.notify("Formatting and mounting %s" % dev_id)
        self.serv.exec_command("sudo mkfs -t ext4 %s" % dev_id)
        self.serv.exec_command("sudo mkdir -p /data")
        self.serv.exec_command("sudo mount %s /data" % dev_id)
        log.notify("Successfully mounted new volume onto /data.")

    def git_pull(self):
        """installs the SEQC directory in /data/software"""
        # todo: replace this with git clone once seqc repo is public

        folder = '/data/software/'
        log.notify('Installing SEQC on remote instance.')
        self.serv.exec_command("sudo mkdir %s" % folder)
        self.serv.exec_command("sudo chown -c ubuntu /data")
        self.serv.exec_command("sudo chown -c ubuntu %s" % folder)

        location = folder + 'seqc.tar.gz'
        # todo: this currently pulls from the develop branch of SEQC.
        # to install a specific version, change version='develop' to the tag (ex: 0.1.6)
        # or seqc.__version__
        self.serv.exec_command(
            'curl -H "Authorization: token a22b2dc21f902a9a97883bcd136d9e1047d6d076" -L '
            'https://api.github.com/repos/ambrosejcarr/seqc/tarball/{version} | '
            'sudo tee {location} > /dev/null'.format(
                location=location, version='master'))
        self.serv.exec_command('cd %s; mkdir seqc && tar -xvf seqc.tar.gz -C seqc '
                               '--strip-components 1' % folder)
        log.notify("Sources are installed and decompressed, setting up seqc")
        self.serv.exec_command('cd %s; sudo pip3 install -e ./' % folder + 'seqc')

        def test_install():
            out_test_install, _ = self.serv.exec_command(
                'process_experiment.py -h | grep RNA')
            if not out_test_install:
                raise BotoCallError()

        retry_boto_call(
            func=test_install, exceptions_to_catch=BotoCallError,
            exception_to_raise=EC2RuntimeError('Error installing SEQC on the cluster.'),
            retries=30, delay_retry=2)()
        log.notify('SEQC successfully installed in %s.' % folder)

    def upload_and_install_seqc(self):
        folder = '/data/software/'
        location = folder + 'seqc.tar.gz'
        log.notify('Installing SEQC on remote instance.')
        self.serv.exec_command("sudo mkdir %s" % folder)
        self.serv.exec_command("sudo chown -c ubuntu /data")
        self.serv.exec_command("sudo chown -c ubuntu %s" % folder)

        seqc_distribution = os.path.expanduser('~/.seqc/seqc.tar.gz')
        self.serv.put_file(seqc_distribution, location)

        self.serv.exec_command('cd %s; mkdir seqc && tar -xvf seqc.tar.gz -C seqc '
                               '--strip-components 1' % folder)

        log.notify("Sources are installed and decompressed, setting up seqc")
        self.serv.exec_command('cd %s; sudo pip3 install -e ./' % folder + 'seqc')

        def test_install():
            out_test_install, _ = self.serv.exec_command(
                'process_experiment.py -h | grep RNA')
            if not out_test_install:
                raise BotoCallError()

        retry_boto_call(
            func=test_install, exceptions_to_catch=BotoCallError,
            exception_to_raise=EC2RuntimeError('Error installing SEQC on the cluster.'),
            retries=30, delay_retry=2)()
        log.notify('SEQC successfully installed in %s.' % folder)

    def set_credentials(self):
        """sets aws credentials on remote instance from user's local config file"""

        self.serv.exec_command('aws configure set aws_access_key_id %s' % self.aws_id)
        self.serv.exec_command(
            'aws configure set aws_secret_access_key %s' % self.aws_key)
        self.serv.exec_command('aws configure set region %s' % self.zone[:-1])

    def setup_cluster(self, volsize, aws_instance, spot_bid=None):
        """
        creates a new cluster, attaches the appropriate volume, configures
        :param spot_bid: amount to use for spot bid, default is None
        :param volsize: size (GB) of volume to be attached
        :param aws_instance: instance type (c3, c4, r3)
        """

        config_file = os.path.expanduser('~/.seqc/config')
        self.configure_cluster(config_file, aws_instance, spot_bid)
        self.create_security_group()
        # modified cluster creation for spot bid
        if self.spot_bid is not None:
            self.create_spot_cluster(volsize)
            self.connect_server()
            self.allocate_space(volsize, True)
        else:
            self.create_cluster()
            self.connect_server()
            self.allocate_space(volsize, False)
        log.notify('Use the following command to log in:\t\t ssh -i {rsa_path} '
                   'ubuntu@{dns_name}'.format(
                       rsa_path=self.keypath, dns_name=self.inst_id.public_dns_name))
        self.git_pull()
        self.set_credentials()
        log.notify('Remote instance successfully configured.')


def terminate_cluster(instance_id):
    """
    terminates a running cluster
    :param instance_id: id of instance to terminate
    """

    ec2 = boto3.resource('ec2')
    instance = ec2.Instance(instance_id)
    try:
        if instance.state['Name'] != 'terminated' and instance.state['Name'] != \
                'shutting-down':
            log.notify('Requesting termination of instance {id}'.format(
                id=instance_id))
            instance.terminate()
            instance.wait_until_terminated()
            log.notify('Termination of instance {id} complete!'.format(
                id=instance_id))
        else:
            log.notify('Instance {id} is not running!'.format(id=instance_id))
    except ClientError:
        log.notify('Instance {id} does not exist!'.format(id=instance_id))


def remove_sg(sg_id):
    """
    removes security group to avoid accumulation
    :param sg_id: security group id number to be deleted
    """

    ec2 = boto3.resource('ec2')
    sg = ec2.SecurityGroup(sg_id)
    sg_name = sg.group_name
    try:
        sg.delete()
        log.notify('security group %s (%s) successfully removed' % (
            sg_name, sg_id))
    except ClientError:
        log.notify('security group %s (%s) is still in use!' % (sg_name, sg_id))


def email_user(attachment: str, email_body: str, email_address: str) -> None:
    """
    sends an email to email address with text contents of email_body and attachment
    attached. Email will come from "Ubuntu@<ec2-instance-ip-of-aws-instance>

    :param attachment: the file location of the attachment to append to the email
    :param email_body: text to send in the body of the email
    :param email_address: the address to which the email should be sent"""

    if isinstance(email_body, str):
        email_body = email_body.encode()

    email_args = ['mutt', '-e', 'set content_type="text/html"', '-a', attachment, '-s',
                  'Remote Process', '--', email_address]
    email_process = Popen(email_args, stdin=PIPE)
    email_process.communicate(email_body)


def upload_results(output_stem: str, email_address: str, aws_upload_key: str,
                   start_pos: str, summary: dict, log_name: str, email: bool) -> None:
    """
    uploads remaining files from the SEQC run
    :param output_stem: specified output directory in cluster
    :param email_address: e-mail where run summary will be sent
    :param aws_upload_key: tar gzipped files will be uploaded to this S3 bucket
    :param start_pos: determines where in the script SEQC started
    :param summary: dictionary of summary statistics from SEQC run
    :param log_name: log name of SEQC run provided by user
    :param email: True if mutt is installed, False otherwise.
    """

    prefix, directory = os.path.split(output_stem)
    counts = output_stem + '_read_and_count_matrices.p'
    log_text = prefix + '/' + log_name

    # generate a run summary and append to log + email
    run_summary = ExperimentalYield.construct_run_summary(summary)
    log.info('A copy of the SEQC run summary can be found below.\nRUN SUMMARY:\n'
             '{run_summary}'.format(run_summary=run_summary))
    files = [counts, log_text]  # counts and log will always be uploaded

    if start_pos == 'start' or start_pos == 'merged':
        alignment_summary = output_stem + '_alignment_summary.txt'
        # copying over alignment summary for upload
        shutil.copyfile(prefix + '/alignments/Log.final.out', output_stem +
                        '_alignment_summary.txt')
        files.append(alignment_summary)

    if aws_upload_key:
        bucket, key = io.S3.split_link(aws_upload_key)
        for item in files:
            try:
                retry_boto_call(io.S3.upload_file)(item, bucket, key)
                item_name = item.split('/')[-1]
                log.info('Successfully uploaded %s to the specified S3 location '
                         '"%s%s".' % (item, aws_upload_key, item_name))
            except FileNotFoundError:
                log.notify('Item %s was not found! Continuing with upload...' % item)

        # get the name of the output file
        log.info('Upload complete. An e-mail will be sent to %s.' % email_address)

    # email results to user
    body = ('<font face="Courier New, Courier, monospace">'
            'SEQC RUN COMPLETE.\n\n'
            'The run log has been attached to this email and '
            'results are now available in the S3 location you specified: '
            '"%s"\n\n'
            'RUN SUMMARY:\n\n%s'
            '</font>' % (aws_upload_key, run_summary))
    log.info('SEQC run complete. Cluster will be terminated unless --no-terminate '
             'flag was specified.')
    body = body.replace('\n', '<br>')
    body = body.replace('\t', '&emsp;')
    if email:
        email_user(log_text, body, email_address)


def check_progress():
    """
    flag that can be used with process_experiment.py --check-progress
    to retrieve remote log and track the progress of the remote run
    """

    # reading in configuration file
    config_file = os.path.expanduser('~/.seqc/config')
    config = configparser.ConfigParser()
    if not config.read(config_file):
        raise ValueError('Please run ./configure (found in the seqc directory) before '
                         'attempting to run process_experiment.py.')

    # obtaining rsa key from configuration file
    rsa_key = os.path.expanduser(config['key']['path_to_rsa_key'])

    # checking for instance status
    inst_file = os.path.expanduser('~/.seqc/instance.txt')
    try:
        with open(inst_file, 'r') as f:
            for line in f:
                entry = line.strip('\n')
                inst_id, run_name = entry.split(':')

                # connecting to the remote instance
                s = SSHServer(inst_id, rsa_key)
                try:
                    inst_state = s.instance.state['Name']
                    if inst_state != 'running':
                        print('Cluster (%s) for run "%s" is currently %s.' %
                              (inst_id, run_name, inst_state))
                        continue
                except:  # todo what should be raised here?
                    print('Cluster (%s) for run "%s" has been terminated.'
                          % (inst_id, run_name))
                    continue

                s.connect()
                out, err = s.exec_command('cd /data; ls *.log')
                if not out:
                    print('ERROR: SEQC log file not found in cluster (%s) for run "%s." '
                          'Something went wrong during remote run.' % (inst_id, run_name))
                    continue
                logfile = out[0]
                out, err = s.exec_command('cat {fname}'.format(fname='/data/'+logfile))
                print('-'*80)
                print('Printing contents of the remote SEQC log file for run "%s":'
                      % run_name)
                print('-'*80)
                for x in out:
                    print(x)
                print('-'*80 + '\n')
    except FileNotFoundError:
        print('You have not started a remote instance -- exiting.')
        sys.exit(0)


def cluster_cleanup():
    """
    checks for all security groups that are unused and deletes
    them prior to each SEQC run. Instance ids that are created by seqc
    will be documented in instances.txt and can used to clean up unused
    security groups/instances after SEQC runs have terminated.
    """

    cmd = 'aws ec2 describe-instances --output text'
    cmd2 = 'aws ec2 describe-security-groups --output text'
    in_use = set([x for x in check_output(cmd.split()).split() if x.startswith(b'sg-')])
    all_sgs = set([x for x in check_output(cmd2.split()).split() if x.startswith(b'sg-')])
    to_delete = list(all_sgs - in_use)

    # set up ec2 resource from boto3
    ec2 = boto3.resource('ec2')

    # iteratively remove unused SEQC security groups
    for sg in to_delete:
        sg_id = sg.decode()
        sg_name = ec2.SecurityGroup(sg_id).group_name
        if 'SEQC-' in sg_name:
            remove_sg(sg_id)

    # check which instances in instances.txt are still running
    inst_file = os.path.expanduser('~/.seqc/instance.txt')

    if os.path.isfile(inst_file):
        with open(inst_file, 'r') as f:
            seqc_list = [line.strip('\n') for line in f]
        if seqc_list:
            with open(inst_file, 'w') as f:
                for i in range(len(seqc_list)):
                    try:
                        entry = seqc_list[i]
                        inst_id = entry.split(':')[0]
                        instance = ec2.Instance(inst_id)
                        if (instance and instance.state and
                                instance.state['Name'] == 'running'):
                            f.write('%s\n' % entry)
                    except (ClientError, KeyError):
                        continue
    else:
        pass  # instances.txt file has not yet been created


class SSHServer(object):
    """Class that wraps the newly launched AWS instance to allow for connecting
    and execution of commands remotely"""

    _error_msg = ('You need to specify a valid RSA key to connect to Amazon EC2 '
                  'instances, see https://github.com/ambrosejcarr/seqc#create-an-rsa-key'
                  '-to-allow-you-to-launch-a-cluster')

    def __init__(self, inst_id, keypath):
        ec2 = boto3.resource('ec2')
        self.instance = ec2.Instance(inst_id)
        self.key = keypath
        self.check_key_file()
        self.ssh = paramiko.SSHClient()

    def check_key_file(self):
        """Checks the rsa file is present"""
        if not os.path.isfile(self.key):
            log.notify('The key %s was not found!' % self.key)
            raise FileNotFoundError(self._error_msg, 'The key file %s does not exist' %
                                    self.key)

    def connect(self):
        """connects to a remote instance"""
        self.ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
        try:
            retry_boto_call(
                func=self.ssh.connect, retries=25, exceptions_to_catch=socket.error,
                exception_to_raise=RuntimeError('Authentication to EC2 failed with key '
                                                '%s.' % self.key),
                delay_retry=5)(self.instance.public_dns_name, username='ubuntu',
                               key_filename=self.key, timeout=3.0)
        except (paramiko.AuthenticationException, paramiko.SSHException):
            log.notify('Authentication failure using key %s, try to log in manually with '
                       'ssh' % self.key)
            raise

    def is_connected(self):
        """checks whether the user is connected to the desired instance"""

        if self.ssh.get_transport() is None:
            return False
        else:
            return True

    def disconnect(self):
        """closes the ssh connection to a remote instance"""

        if self.is_connected():
            self.ssh.close()

    def get_file(self, local_file, remote_file):
        """
        obtains a file located on a remote AWS instance and downloads it locally
        :param local_file: name of file to be downloaded on the local machine
        :param remote_file: name of file on AWS instance to be fetched
        """

        if not self.is_connected():
            log.notify('You are not connected!')
            raise ConnectionError("")
        ftp = self.ssh.open_sftp()
        ftp.get(remote_file, local_file)
        ftp.close()

    def put_file(self, local_file, remote_file):
        """
        places a file from the local machine onto a remote instance
        :param local_file: name of file to be copied to remote instance
        :param remote_file: name of file placed remotely
        """

        if not self.is_connected():
            log.notify('You are not connected!')
            sys.exit(2)
        ftp = self.ssh.open_sftp()
        ftp.put(local_file, remote_file)
        log.info('Successfully placed {local_file} in {remote_file}.'.format(
            local_file=local_file, remote_file=remote_file))
        ftp.close()

    def exec_command(self, args):
        """
        executes the specified arguments remotely on an AWS instance
        :param args: args to be executed remotely
        :return: decoded stdout and stderr parsed into lines
        """

        if not self.is_connected():
            log.notify('You are not connected!')
            sys.exit(2)
        stdin, stdout, stderr = self.ssh.exec_command(args)
        stdin.flush()
        data = stdout.read().decode().splitlines()  # response in bytes
        errs = stderr.read().decode().splitlines()
        return data, errs
