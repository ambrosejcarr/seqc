import os
import logging
import time
import random
import configparser
import traceback
import types
import dill
from functools import wraps
from contextlib import closing
import paramiko
import boto3
from seqc import log
from seqc.exceptions import (
    RetryLimitExceeded, ConfigurationError, InstanceNotRunningError, EC2RuntimeError)
from botocore.exceptions import ClientError

# change some logging defaults
logging.getLogger('paramiko').setLevel(logging.CRITICAL)
logging.getLogger('boto3').setLevel(logging.CRITICAL)


# set some configuration values
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


class Retry:

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
    client = boto3.client('ec2'
                          )

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
        :param bool terminate_on_complete: flag for use in contextmanager or function
          decorator context. If true, the instance will be automatically shut down when
          either (1) the contextmanager exits or (2) the decorated function completes

        """

        # todo change these dummy arguments into properties
        # dummy arguments set by kwargs or configuration; for inspection only.
        self.aws_public_access_key = None
        self.aws_secret_access_key = None
        self.availability_zone = None
        self._rsa_key_path = None
        self.subnet_id = None

        self.image_id = None
        self.instance_type = None

        self._instance = None
        self.security_group = None
        self.spot_bid = None
        self.terminate_on_complete = False
        self.volume_size = 5

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

    # todo define def __repr__(self):

    @classmethod
    def _check_instance(cls, instance):
        if isinstance(instance, str):
            instance = cls.ec2.Instance(instance)
            try:
                instance.description
            except ClientError:
                raise ValueError(
                    'provided instance id does not map to an existing instance')
        return instance

    @property
    def instance(self):
        return self._instance

    @instance.setter
    def instance(self, instance_or_id):
        self._instance = self._check_instance(instance_or_id)

    @property
    def rsa_key_path(self):
        return self._rsa_key_path

    @rsa_key_path.setter
    def rsa_key_path(self, value):
        if not isinstance(value, str):
            raise ValueError('rsa_key_path must be type str')
        self._rsa_key_path = os.path.expanduser(value)

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
    @Retry(retries=5)
    def enable_ssh(cls, security_group):
        try:
            security_group.authorize_ingress(
                IpProtocol="tcp", CidrIp="0.0.0.0/0", FromPort=22, ToPort=22)
            security_group.authorize_ingress(
                SourceSecurityGroupName=security_group.description)
        except ClientError as e:
            if 'InvalidPermission.Duplicate' not in e.args[0]:  # todo figure out why this is happening
                raise
        log.notify('Enabled ssh access via port 22 for security group %s' %
                   security_group.id)

    @classmethod
    @Retry(retries=20, delay=0.5)
    def verify_security_group(cls, security_group) -> None:
        """If the security group has not been recognized server-side, this will thrown an
        exception.
        """
        _ = security_group.description

    @classmethod
    @Retry(retries=10, delay=0.5)
    def remove_security_group(cls, security_group) -> None:
        security_group.delete()
        log.notify('security group %s successfully removed.' % (
            security_group.id))

    # @property
    def launch_specification(self) -> dict:
        """return the specification for launching an instance with parameters defined
        by the class constructor
        """
        if self.security_group is None:
            sg = self.create_security_group()
            self.verify_security_group(sg)
            self.enable_ssh(sg)
            self.security_group = sg
        else:
            sg = self.security_group

        spec = {
            'ImageId': self.image_id,
            'KeyName': self.rsa_key_path.split('/')[-1].split('.')[0],
            'InstanceType': self.instance_type,
            'Placement': {'AvailabilityZone': self.availability_zone},
            'SecurityGroupIds': [sg.id],
            'BlockDeviceMappings': [{'DeviceName': '/dev/xvdf',
                                     'Ebs': {'VolumeSize': self.volume_size,
                                             'VolumeType': 'gp2',
                                             'DeleteOnTermination': True}}],
        }
        if self.instance_type in ['c4', 'r3']:
            spec['SubnetId'] = self.subnet_id
        return spec

    @staticmethod
    @Retry(retries=40, catch=(ClientError, InstanceNotRunningError), delay=5)
    def verify_instance_running(instance):
        """wait for instance to reach 'running' state, then return"""
        if not instance.state['Name'] == 'running':
            raise InstanceNotRunningError
        log.notify('instance %s in running state' % instance.id)

    def create_instance(self) -> None:
        if self.instance is None:
            if self.spot_bid:
                self.create_spot_instance()
            else:
                specification = self.launch_specification()
                specification['MinCount'] = specification['MaxCount'] = 1
                self.instance = self.ec2.create_instances(**specification)[0]
                log.notify('instance %s created' % self.instance.id)
                self.instance.wait_until_running()

    def setup_seqc(self):
        if self.instance is None:
            self.create_instance()  # todo split here if spot
        with SSHConnection(self.instance) as ssh:
            log.notify('Installing SEQC on remote instance.')
            seqc_distribution = os.path.expanduser('~/.seqc/seqc.tar.gz')
            ssh.execute('sudo mkdir -p /data/software')
            ssh.execute('sudo chown -c ubuntu /data')
            ssh.execute('sudo chown -c ubuntu /data/software')
            ssh.put_file(seqc_distribution, '/data/software/seqc.tar.gz')
            ssh.execute(
                'cd /data/software/; '
                'mkdir seqc && tar -xvf seqc.tar.gz -C seqc --strip-components 1')

            log.notify("Sources are installed and decompressed, setting up seqc")
            ssh.execute('sudo pip3 install -e data/software/seqc/')

            # test the installation
            try:
                ssh.execute('SEQC.py -h')
            except:
                log.notify('SEQC installation failed.')
                log.exception()
                raise
            log.notify('SEQC setup complete.')

    def start(self):
        self.setup_seqc()
        log.notify('Instance set-up complete.')

    def stop(self):
        """stops a running instance"""
        if self.instance is None:
            raise RuntimeError('Instance not yet created, nothing to be stopped.')
        elif self.instance.state['Name'] not in ('stopped', 'terminated', 'shutting-down'):
            iid = self.instance.id
            log.notify('requesting termination of instance {id}'.format(id=iid))
            self.instance.stop()
            self.instance.wait_until_stopped()
            log.notify('instance {id} stopped.'.format(id=iid))
        else:
            log.notify('instance is not running')

    def restart(self):
        """restarts a stopped instance"""
        if self.instance is None:
            raise RuntimeError('Instance not yet created, nothing to be restarted.')
        if self.instance.state['Name'] == 'stopped':
            self.instance.start()
            self.instance.wait_until_running()
            log.notify('Stopped instance %s has restarted.' % self.instance.id)
        else:
            log.notify('Instance %s in state "%s" must be in a stopped state to be '
                       'restarted.' % (self.instance.id, self.instance.state['Name']))

    def terminate(self):
        """terminates an instance in any state (including stopped)"""
        if self.instance is None:
            raise RuntimeError('Instance not yet created, nothing to be restarted.')
        elif self.instance.state['Name'] not in ('terminated', 'shutting-down'):
            iid = self.instance.id
            log.notify('requesting termination of instance {id}'.format(id=iid))
            self.instance.terminate()
            self.instance.wait_until_terminated()
            log.notify('instance {id} terminated.'.format(id=iid))
        else:
            log.notify('Instance %s in state "%s" must be running to be stopped.' %
                       (self.instance.id, self.instance.state['Name']))

    @classmethod
    @Retry(retries=40, delay=5, catch=InstanceNotRunningError)
    def verify_spot_bid_fulfilled(cls, sir_id):
        result = cls.client.describe_spot_instance_requests(
            SpotInstanceRequestIds=[sir_id])
        status = result['SpotInstanceRequests'][0]['Status']['Code']
        if status not in ['pending-evaluation', 'pending-fulfillment', 'fulfilled']:
            raise EC2RuntimeError('spot request bad-status: %s' % status)
        elif status != 'fulfilled':
            raise InstanceNotRunningError
        return cls.ec2.Instance(result['SpotInstanceRequests'][0]['InstanceId'])

    def create_spot_instance(self):
        if not self.spot_bid:
            raise ValueError('must pass constructor spot_bid price (float) to create a '
                             'spot bid request.')
        response = self.client.request_spot_instances(
            DryRun=False, SpotPrice=str(self.spot_bid),
            LaunchSpecification=self.launch_specification())
        sir_id = response['SpotInstanceRequests'][0]['SpotInstanceRequestId']
        log.notify(
            'spot instance requested (%s), waiting for bid to be accepted.' % sir_id)
        instance = self.verify_spot_bid_fulfilled(sir_id)
        if instance is None:
            raise InstanceNotRunningError(
                'spot bid of %f was not fulfilled, please try a higher bid or ')
        log.notify('spot bid accepted, waiting for instance to attain running state.')
        self.instance = instance
        self.instance.wait_until_running()
        log.notify('spot instance (id=%s) in running state' % self.instance.id)
        # self.verify_instance_running(instance)



    def __enter__(self):
        self.setup_seqc()
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        if self.terminate_on_complete:
            self.instance.terminate()
        if not exc_type:
            return True

    @staticmethod
    def pickle_function(function: object, args, kwargs) -> str:
        """ pickle and function and its arguments

        :param object function: function to be pickled
        :param tuple args: positional arguments for the function
        :param dict kwargs: keyword arguments for the function
        :return str: filename of the pickled function
        """
        filename = '{}{}.p'.format(os.environ['TMPDIR'], function.__name__)
        with open(filename, 'wb') as f:
            dill.dump(dict(function=function, args=args, kwargs=kwargs), f)
        return filename

    @staticmethod
    def get_imports():
        for alias, val in globals().items():
            if isinstance(val, types.ModuleType):
                yield (val.__name__, alias)

    @classmethod
    def format_importlist(cls):
        importlist = ''
        for name, alias in cls.get_imports():
            if name != alias:
                importlist += 'import {name} as {alias}\n'.format(name=name, alias=alias)
            else:
                importlist += 'import {name}\n'.format(name=name)
        return importlist

    @classmethod
    def write_script(cls, function) -> str:
        """generate a python script that calls function after importing required modules

        :param object function: function to be called
        :return str: filename of the python script
        """
        script_name = '{}{}.py'.format(os.environ['TMPDIR'], function.__name__)
        script_body = (
            '{imports}'
            'with open("/data/func.p", "rb") as fin:\n'
            '    data = dill.load(fin)\n'
            'results = data["function"](*data["args"], **data["kwargs"])\n'
            'with open("/data/results.p", "wb") as f:\n'
            '    dill.dump(results, f)\n'
        )
        script_body = script_body.format(imports=cls.format_importlist())

        with open(script_name, 'w') as f:
            print('writing script to file:\n%s' % script_body)
            f.write(script_body)
        return script_name

    def __call__(self, function):

        def function_executed_on_aws(*args, **kwargs):

            # dump original function to file
            script = self.write_script(function)
            func = self.pickle_function(function, args, kwargs)

            # create an instance
            self.setup_seqc()

            with SSHConnection(self.instance, self.rsa_key_path) as ssh:
                ssh.put_file(script, '/data/script.py')
                ssh.put_file(func, '/data/func.p')
                if self.terminate_on_complete:
                    ssh.execute('python3 /data/script.py')
                    results_name = os.environ['TMPDIR'] + function.__name__ + '_results.p'
                    ssh.get_file('/data/results.p', results_name)
                    with open(results_name, 'rb') as f:
                        results = dill.load(f)
                else:
                    ssh.execute('nohup python3 /data/script.py > /dev/null/ 2>&1 &')
                    results = None

            if self.terminate_on_complete:
                self.terminate(self.instance)

            return results

        return function_executed_on_aws


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

    # todo define def __repr__(self):

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

    @Retry(retries=50, delay=2.5)
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

    def execute(self, args):
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
