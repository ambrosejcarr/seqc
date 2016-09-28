import os
import sys
import logging
import time
import random
import configparser
import traceback
from functools import wraps
from contextlib import closing
import paramiko
import boto3
from seqc import log
from seqc.exceptions import (
    RetryLimitExceeded, ConfigurationError, InstanceNotRunningError)
from botocore.exceptions import ClientError

logging.getLogger('paramiko').setLevel(logging.CRITICAL)

config_file = os.path.expanduser('~/.seqc/config')


def _set_default_ec2_configuration_values():
    """parses configuration file for default aws configuation"""
    config = configparser.ConfigParser()
    config.read(config_file)
    return config

_defaults = _set_default_ec2_configuration_values()


def _get_default_ec2_parameter(*args, config_dict=_defaults):
    if isinstance(args, tuple) and len(args) == 1:
        try:
            return config_dict[args[0]]
        except KeyError:
            raise ConfigurationError(
                'no value for {k} found in config file {c}, and {k} was '
                'not provided to AWSInstance constructor.'.format(
                    c=config_file, k=args[0]))
    else:
        return _get_default_ec2_parameter(
            *args[1:], config_dict=_defaults[args[0]])


class retry:

    def __init__(
            self,
            retries: int=10,
            catch=(ClientError,),
            delay: int=1,
            verbose=False):
        self.retries = retries
        self.exceptions_to_catch = catch
        self.delay_retry = delay
        self.verbose = verbose

    def __call__(self, function):

        @wraps(function)
        def wrapper(*args, **kwargs):
            retries = self.retries
            while True:
                try:
                    return function(*args, **kwargs)
                except self.exceptions_to_catch:
                    if retries > 0:
                        retries -= 1
                        if self.verbose:
                            log.notify(
                                'Non fatal error in function {} (retrying in '
                                '{!s}s):\n{}'.format(
                                    function.__qualname__, self.delay_retry,
                                    traceback.format_exc()))
                        time.sleep(self.delay_retry)
                    else:
                        raise RetryLimitExceeded(
                            'fatal error in function {} occurred {} times at {!s}s call '
                            'interval:\n{}'.format(
                                function.__qualname__, self.retries, self.delay_retry,
                                traceback.format_exc()))
        return wrapper


class AWSInstance(object):
    """
    Connects to AWS instance using paramiko and a private RSA key,
    allows for the creation/manipulation of EC2 instances and executions
    of commands on the remote server
    """

    ec2 = boto3.resource('ec2')
    client = boto3.client('ec2')

    # todo make security_group and instance properties so they support string arguments
    # (do the typechecking at the property setter level)
    def __init__(
            self, instance_alias, **kwargs):
        """
        required:
        ---------
        :param str instance_alias: type of instance to create. options=('c4, 'c3', 'r3')

        passed keyword arguments or present in config:
        ---------------------------
        :param str rsa_key_path:
        :param str image_id:
        :param str subnet_id:
        :param str availability_zone:
        :param str aws_public_access_key:
        :param str aws_secret_access_key:

        optional arguments:
        -------------------
        :param boto3.resources.factory.ec2.Instance instance: instance object
        :param boto3.resources.factory.ec2.SecurityGroup security_group: security group
          object
        :param float spot_bid:
        :param int volume_size: size in Gb of volume to attach to /data. Must be smaller
          than 16 and greater than 0.

        """

        # dummy arguments set by kwargs or configuration; for inspection only.
        self.aws_public_access_key = None
        self.aws_secret_access_key = None
        self.availability_zone = None
        self.rsa_key_path = None
        self.subnet_id = None

        self.image_id = None
        self.instance_type = None

        self.instance = None
        self.security_group = None
        self.spot_bid = None
        self.volume_size = 5  # by default, create a 5gb volume

        # additional properties
        self._ssh_connection = None

        # list of arguments to search config for
        security = ('aws_public_access_key', 'aws_secret_access_key', 'availability_zone',
                    'rsa_key_path', 'subnet_id')
        instance = ('image_id', 'instance_type')
        optional = ('instance', 'security_group', 'spot_bid', 'volume_size')

        # set defaults from config file if they weren't passed
        for k in security:
            if k not in kwargs.keys():
                setattr(self, k, _get_default_ec2_parameter('security', k))

        for k in instance:
            if k not in kwargs.keys():
                setattr(self, k, _get_default_ec2_parameter(instance_alias, k))

        for k in optional:
            if k in kwargs.keys():
                setattr(self, k, kwargs[k])

    @classmethod
    def create_security_group(cls, name=None):
        """creates a new security group

        :param str name: optional, name of the security group to create. Note that this
          name must be unique or an error will be thrown. Default: SEQC + random 7-integer
          number.
        """

        # todo get list of existing groups; check against
        if name is None:
            name = 'SEQC-%07d' % random.randint(1, int(1e7))
        sg = cls.ec2.create_security_group(GroupName=name, Description=name)
        log.notify('Created new security group: %s (name=%s).' % (sg.id, name))
        return sg

    @classmethod
    @retry()
    def enable_ssh(cls, security_group):
        security_group.authorize_ingress(
            IpProtocol="tcp", CidrIp="0.0.0.0/0", FromPort=22, ToPort=22)
        security_group.authorize_ingress(
            SourceSecurityGroupName=security_group.description)
        log.notify('Enabled ssh access via port 22 for security group %s' %
                   security_group.id)

    @classmethod
    @retry(retries=20, delay=0.5)
    def verify_security_group(cls, security_group) -> None:
        """If the security group has not been recognized server-side, this will thrown an
        exception.
        """
        _ = security_group.description

    @classmethod
    @retry(retries=10, delay=0.5)
    def remove_security_group(cls, security_group) -> None:
        security_group.delete()
        log.notify('security group %s successfully removed.' % (
            security_group.id))

    @property
    def launch_specification(self) -> dict:
        """return the specification for launching an instance with parameters defined
        by the class constructor
        """
        if self.security_group is None:
            sg = self.create_security_group()
            self.verify_security_group(sg)
            self.enable_ssh(sg)
        else:
            sg = self.security_group

        spec = {
            'ImageID': self.image_id,
            'KeyName': self.rsa_key_path.split('/')[-1],
            'InstanceType': self.instance_type,
            'Placement': {'AvailabilityZone': self.availability_zone},
            'SecurityGroupIds': [sg],
            'BlockDeviceMappings': [{'DeviceName': '/dev/xvdf',
                                     'Ebs': {'VolumeSize': self.volume_size,
                                             'DeleteOnTermination': True}}],
        }
        if self.instance_type in ['c4', 'r3']:
            spec['SubnetId'] = self.subnet_id
        return spec

    @staticmethod
    @retry(retries=40, catch=(ClientError, InstanceNotRunningError), delay=5)
    def verify_instance_running(instance):
        """wait for instance to reach 'running' state, then return"""
        if not instance.state['Name'] == 'running':
            raise InstanceNotRunningError
        log.notify('instance %s in running state' % instance.id)

    def create_instance(self) -> None:
        if self.instance is not None:
            specification = self.launch_specification
            specification['MinCount'] = specification['MaxCount'] = 1
            instance = self.ec2.create_instances(**specification)
            self.verify_instance_running(instance)
            self.instance = instance

    def setup_seqc(self):
        if self.instance is None:
            self.create_instance()  # todo split here if spot
        with SSHConnection(self.instance) as ssh:
            log.notify('Installing SEQC on remote instance.')
            seqc_distribution = os.path.expanduser('~/.seqc/seqc.tar.gz')
            ssh.exec_command('sudo mkdir -p /data/software')
            ssh.exec_command('sudo chown -c ubuntu /data')
            ssh.exec_command('sudo chown -c ubuntu /data/software')
            ssh.put_file(seqc_distribution, '/data/software/seqc.tar.gz')
            ssh.exec_command(
                'cd /data/software/; '
                'mkdir seqc && tar -xvf seqc.tar.gz -C seqc --strip-components 1')

            log.notify("Sources are installed and decompressed, setting up seqc")
            ssh.exec_command('sudo pip3 install -e data/software/seqc/')

            # test the installation
            try:
                ssh.exec_command('SEQC.py -h')
            except:
                log.notify('SEQC installation failed.')
                log.exception()
                raise

    # todo implement these.
    def start(self):
        raise NotImplementedError

    def stop(self):
        raise NotImplementedError

    def terminate(self):
        raise NotImplementedError

    def create_spot_bid(self):
        response = self.client.request_spot_instances(
            DryRun=False, SpotPrice=str(self.spot_bid),
            LaunchSpecification=self.launch_specification)
        raise NotImplementedError

    def __enter__(self):
        raise NotImplementedError

    def __exit__(self, exc_type, exc_val, exc_tb):
        raise NotImplementedError

    def __call__(self, *args, **kwargs):
        raise NotImplementedError


class SSHConnection:
    _error_msg = ('You need to specify a valid RSA key to connect to Amazon EC2 '
                  'instances, see https://github.com/ambrosejcarr/seqc#create-an-rsa-key'
                  '-to-allow-you-to-launch-a-cluster')

    ec2 = boto3.resource('ec2')

    def __init__(self, instance, rsa_key_path=None):
        self._instance = self.ec2.Instance(instance)
        if rsa_key_path is None:
            self.rsa_key_path = _get_default_ec2_parameter('rsa_key_path')
        self.check_key_file()
        self.ssh = paramiko.SSHClient()

    @property
    def instance(self):
        return self._instance

    @instance.setter
    def instance(self, item):
        if isinstance(item, str):
            self._instance = self.ec2.Instance(item)
        elif isinstance(item, boto3.resources.factory.ec2.Instance):
            self._instance = item
        else:
            raise TypeError('instance must either be a Boto3 instance object or a string '
                            'instance id')

    def check_key_file(self):
        """Checks the rsa file is present"""
        if not os.path.isfile(self.rsa_key_path):
            log.notify('The key %s was not found!' % self.rsa_key_path)
            raise FileNotFoundError(self._error_msg, 'The key file %s does not exist' %
                                    self.rsa_key_path)

    @retry(retries=50, delay=2.5)
    def connect(self):
        """connects to a remote instance"""
        self.ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
        try:
            self.ssh.connect(self.instance.public_dns_name, username='ubuntu',
                             key_filename=self.rsa_key_path, timeout=3.0)
        except (paramiko.AuthenticationException, paramiko.SSHException):
            log.notify('Authentication failure using key %s, try to log in manually with '
                       'ssh' % self.rsa_key_path)
            log.exception()
            raise

    def is_connected(self):
        """checks whether the user is connected to the desired instance"""
        return False if self.ssh.get_transport() is None else True

    def disconnect(self):
        """closes the ssh connection to a remote instance"""
        if self.is_connected():
            self.ssh.close()

    def get_file(self, remote_file, local_file) -> None:
        """obtains a file located on a remote AWS instance and downloads it locally

        :param remote_file: name of file on AWS instance to be fetched
        :param local_file: name of file to be downloaded on the local machine
        """
        if not self.is_connected():
            self.connect()
        with closing(self.ssh.open_sftp()) as ftp:
            ftp.get(remote_file, local_file)  # todo check order

    def put_file(self, local_file, remote_file):
        """places a file from the local machine onto a remote instance

        :param local_file: name of file to be copied to remote instance
        :param remote_file: name of file placed remotely
        """

        if not self.is_connected():
            self.connect()
        with closing(self.ssh.open_sftp()) as ftp:
            ftp.put(local_file, remote_file)  # todo check order
            log.info('placed {lfile} at {rfile}.'.format(
                lfile=local_file, rfile=remote_file))

    def exec_command(self, args):
        """executes the specified arguments remotely on an AWS instance

        :param args: args to be executed remotely
        :return: decoded stdout and stderr parsed into lines
        """

        if not self.is_connected():
            self.connect()
        stdin, stdout, stderr = self.ssh.exec_command(args)
        stdin.flush()
        data = stdout.read().decode().splitlines()
        errs = stderr.read().decode().splitlines()
        return data, errs

    def login_command(self):
        return ('ssh -i {rsa_path} ubuntu@{dns_name}'.format(
            rsa_path=self.rsa_key_path, dns_name=self.instance.public_dns_name))

    def __enter__(self):
        self.connect()
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.disconnect()
        if not exc_type:
            return True
