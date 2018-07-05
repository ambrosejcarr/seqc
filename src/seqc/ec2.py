import os
import time
import random
import configparser
import traceback
import types
import dill
from functools import wraps
from contextlib import closing
from paramiko.ssh_exception import NoValidConnectionsError
import paramiko
import boto3
import socket
from subprocess import Popen, PIPE
from seqc import log, io
from seqc.core import verify
from seqc.exceptions import (
    RetryLimitExceeded, InstanceNotRunningError, EC2RuntimeError)
from botocore.exceptions import ClientError

# change some logging defaults
log.logging.getLogger('paramiko').setLevel(log.logging.CRITICAL)
log.logging.getLogger('boto3').setLevel(log.logging.CRITICAL)

# set default values for a few parameters
IMAGE_ID = 'ami-8927f1f3'


def _get_ec2_configuration():
    """assumes you have awscli and that you have configured it. If so, the default values
    for credentials will be searchable!"""
    defaults = {}
    config = configparser.ConfigParser()
    config.read(os.path.expanduser('~/.aws/config'))
    defaults['region'] = config['default']['region']
    config.read(os.path.expanduser('~/.aws/credentials'))
    defaults['aws_access_key_id'] = config['default']['aws_access_key_id']
    defaults['aws_secret_access_key'] = config['default']['aws_secret_access_key']
    return defaults


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
    client = boto3.client('ec2')

    def __init__(self, rsa_key, instance_type, instance_id=None, security_group_id=None,
                 spot_bid=None, synchronous=False, volume_size=5, **kwargs):
        """

        passed keyword arguments or present in config:
        ---------------------------
        :param str rsa_key:
        :param str image_id:
        :param str instance_type:

        optional arguments:
        -------------------
        :param str instance_id: id for ec2 instance
        :param str security_group_id: id for ec2 security group object
        :param float spot_bid: amount of money to pay per hour in dollars
        :param int volume_size: size in Gb of volume to attach to home directory. Must be
          between 1 and 16384gb.
        :param bool synchronous: flag for use in contextmanager or function
          decorator context. If true, the instance will be automatically shut down when
          either (1) the contextmanager exits or (2) the decorated function completes

        """
        # todo allow overwriting of these arguments with **kwargs

        defaults = _get_ec2_configuration()
        self.aws_public_access_key = defaults['aws_access_key_id']
        self.aws_secret_access_key = defaults['aws_secret_access_key']
        self.region = defaults['region']
        self._rsa_key = rsa_key
        self.image_id = IMAGE_ID
        self.instance_type = instance_type

        self._instance_id = instance_id
        self._security_group_id = security_group_id
        self.spot_bid = spot_bid
        self.synchronous = synchronous

        if not isinstance(volume_size, int) or not 1 <= volume_size < 16384:
            raise ValueError('volume size must be an integer ')
        self.volume_size = volume_size

        # additional properties
        self._ssh_connection = None

    # todo define def __repr__(self):

    @property
    def instance_id(self):
        return self._instance_id

    @instance_id.setter
    def instance_id(self, value):
        if not isinstance(value, str):
            raise ValueError('instance must be a string instance id')
        if not value.startswith('i-'):
            raise ValueError('valid instance identifiers must start with "i-"')
        self._instance_id = value

    @property
    def security_group_id(self):
        return self._security_group_id

    @security_group_id.setter
    def security_group_id(self, value):
        if not isinstance(value, str):
            raise ValueError('instance must be a string instance id')
        if not value.startswith('sg-'):
            raise ValueError('valid instance identifiers must start with "i-"')
        self._security_group_id = value

    @property
    def rsa_key(self):
        return self._rsa_key

    @rsa_key.setter
    def rsa_key(self, value):
        if not isinstance(value, str):
            raise ValueError('rsa_key_path must be type str')
        self._rsa_key = os.path.expanduser(value)

    @property
    def Instance(self):
        return self.ec2.Instance(self.instance_id)

    # todo set this security group name according to the seqc run
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
        return sg.id

    @classmethod
    @Retry(retries=5)
    def enable_ssh(cls, security_group_id):
        security_group = cls.ec2.SecurityGroup(security_group_id)
        try:
            security_group.authorize_ingress(
                IpProtocol="tcp", CidrIp="0.0.0.0/0", FromPort=22, ToPort=22)
            security_group.authorize_ingress(
                SourceSecurityGroupName=security_group.description)
        except ClientError as e:  # todo figure out why this is happening
            if 'InvalidPermission.Duplicate' not in e.args[0]:
                raise
        log.notify('Enabled ssh access via port 22 for security group %s' %
                   security_group_id)

    @classmethod
    @Retry(retries=20, delay=0.5)
    def verify_security_group(cls, security_group_id) -> None:
        """If the security group has not been recognized server-side, this will thrown an
        exception.
        """
        _ = cls.ec2.SecurityGroup(security_group_id).description

    @classmethod
    @Retry(retries=10, delay=0.5)
    def remove_security_group(cls, security_group_id) -> None:
        cls.ec2.SecurityGroup(security_group_id).delete()
        log.notify('security group %s successfully removed.' % (
            security_group_id))

    def launch_specification(self) -> dict:
        """return the specification for launching an instance with parameters defined
        by the class constructor
        """
        if self.security_group_id is None:
            sg_id = self.create_security_group()
            self.verify_security_group(sg_id)
            self.enable_ssh(sg_id)
            self.security_group_id = sg_id
        else:
            sg_id = self.security_group_id

        spec = {
            'ImageId': self.image_id,
            'KeyName': self.rsa_key.split('/')[-1].split('.')[0],
            'InstanceType': self.instance_type,
            'SecurityGroupIds': [sg_id],
            'BlockDeviceMappings': [{'DeviceName': '/dev/xvdf',
                                     'Ebs': {'VolumeSize': self.volume_size,
                                             'VolumeType': 'gp2',
                                             'DeleteOnTermination': True}}],
        }
        return spec

    @Retry(retries=40, catch=(ClientError, InstanceNotRunningError), delay=5)
    def verify_instance_running(self, instance_id):
        """wait for instance to reach 'running' state, then return"""
        instance = self.ec2.Instance(id=instance_id)
        if not instance.state['Name'] == 'running':
            raise InstanceNotRunningError
        log.notify('instance %s in running state' % instance_id)

    def create_instance(self) -> None:
        if self.instance_id is not None:
            raise RuntimeError('instance %s already exists.' % self.instance_id)
        if self.spot_bid:
            self.create_spot_instance()
        else:
            specification = self.launch_specification()
            specification['MinCount'] = specification['MaxCount'] = 1
            instance = self.ec2.create_instances(**specification)[0]
            self.instance_id = instance.id
            log.notify('instance %s created, waiting until running' %
                       self.instance_id)
            instance.wait_until_running()
            log.notify('instance %s in running state' % self.instance_id)

    @staticmethod
    def mount_volume(ssh, directory='/home/ec2-user'):
        """mount /dev/xvdf to /data given an ssh client with access to an instance

        :param str directory: directory to mount the drive to. Note that odd behavior may
          be encountered if the home directory is not the mount target, since paramiko
          automatically uses the home directory as the point of execution.
        """
        try:  # pass errors related to the drive already being mounted
            log.notify("Formatting and mounting /dev/xvdf to %s" % directory)
            ssh.execute("sudo mkfs -t ext4 /dev/xvdf 2>&1")  # redir; errors invisible
            ssh.execute("sudo cp -a %s/. /tmp/directory/" % directory)  # copy original
            ssh.execute("sudo mkdir -p %s" % directory)
            ssh.execute("sudo mount /dev/xvdf %s && sudo cp -a /tmp/directory/. %s/"
                        % (directory, directory))
            ssh.execute("sudo chown ec2-user:ec2-user %s/lost+found && "
                        "chmod 755 %s/lost+found" % (directory, directory))
            log.notify("Successfully mounted new volume onto %s." % directory)
        except ChildProcessError as e:
            if not ('mount: according to mtab, /dev/xvdf is already mounted on %s'
                        % directory in ' '.join(e.args[0])):
                raise

    def set_credentials(self, ssh):
        """sets aws credentials on remote instance from user's local config file"""

        ssh.execute('aws configure set aws_access_key_id %s' % self.aws_public_access_key)
        ssh.execute(
            'aws configure set aws_secret_access_key %s' % self.aws_secret_access_key)
        ssh.execute('aws configure set region %s' % self.region)

    def setup_seqc(self):
        if self.instance_id is None:
            self.create_instance()
        with SSHConnection(instance_id=self.instance_id,
                           rsa_key=self.rsa_key) as ssh:
            self.mount_volume(ssh)
            log.notify('setting aws credentials.')
            self.set_credentials(ssh)
            log.notify('uploading local SEQC installation to remote instance.')
            seqc_distribution = os.path.expanduser('~/.seqc/seqc.tar.gz')
            ssh.execute('mkdir -p software/seqc')
            ssh.put_file(seqc_distribution, 'software/seqc.tar.gz')
            ssh.execute(
                'tar -m -xvf software/seqc.tar.gz -C software/seqc --strip-components 1')
            log.notify("Sources are uploaded and decompressed, installing seqc.")
            try:
                ssh.execute('sudo -H pip3 install -e software/seqc/')
            except ChildProcessError as e:
                if 'pip install --upgrade pip' in str(e):
                    pass
                else:
                    raise

            try:  # test the installation
                ssh.execute('SEQC -h')
            except:
                log.notify('SEQC installation failed.')
                log.exception()
                raise
            log.notify('SEQC setup complete.')
            log.notify('instance login: %s' % ssh.login_command())

    def start(self):
        self.setup_seqc()
        log.notify('Instance set-up complete.')

    def stop(self):
        """stops a running instance"""
        if self.instance_id is None:
            raise RuntimeError('Instance not yet created, nothing to be stopped.')
        instance = self.ec2.Instance(self.instance_id)
        if instance.state['Name'] not in (
                'stopped', 'terminated', 'shutting-down'):
            log.notify('requesting termination of instance {id}'.format(
                id=self.instance_id))
            instance.stop()
            instance.wait_until_stopped()
            log.notify('instance {id} stopped.'.format(id=self.instance_id))
        else:
            log.notify('instance is not running')

    def restart(self):
        """restarts a stopped instance"""
        if self.instance_id is None:
            raise RuntimeError('Instance not yet created, nothing to be restarted.')
        instance = self.ec2.Instance(self.instance_id)
        if instance.state['Name'] == 'stopped':
            instance.start()
            instance.wait_until_running()
            log.notify('Stopped instance %s has restarted.' % self.instance_id)
        else:
            log.notify('Instance %s in state "%s" must be in a stopped state to be '
                       'restarted.' % (self.instance_id, instance.state['Name']))

    def terminate(self):
        """terminates an instance in any state (including stopped)"""
        if self.instance_id is None:
            raise RuntimeError('Instance not yet created, nothing to be restarted.')
        instance = self.ec2.Instance(self.instance_id)
        if instance.state['Name'] not in ('terminated', 'shutting-down'):
            log.notify('requesting termination of instance {id}'.format(
                id=self.instance_id))
            instance.terminate()
            instance.wait_until_terminated()
            log.notify('instance {id} terminated.'.format(id=self.instance_id))
        else:
            log.notify('Instance %s in state "%s" must be running to be stopped.' %
                       (self.instance_id, instance.state['Name']))

    @classmethod
    @Retry(retries=40, delay=5, catch=(InstanceNotRunningError, ClientError))
    def verify_spot_bid_fulfilled(cls, sir_id):
        result = cls.client.describe_spot_instance_requests(
            SpotInstanceRequestIds=[sir_id])
        status = result['SpotInstanceRequests'][0]['Status']['Code']
        if status not in ['pending-evaluation', 'pending-fulfillment', 'fulfilled']:
            raise EC2RuntimeError('spot request bad-status: %s' % status)
        elif status != 'fulfilled':
            raise InstanceNotRunningError
        return result['SpotInstanceRequests'][0]['InstanceId']

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
        self.instance_id = self.verify_spot_bid_fulfilled(sir_id)
        if self.instance_id is None:
            raise InstanceNotRunningError(
                'spot bid of %f was not fulfilled, please try a higher bid or ')
        log.notify('spot bid accepted, waiting for instance (id=%s) to attain running '
                   'state.' % self.instance_id)
        self.ec2.Instance(self.instance_id).wait_until_running()
        log.notify('spot instance (id=%s) in running state' % self.instance_id)

    def __enter__(self):
        try:
            self.setup_seqc()
        except:
            if self.synchronous and self.instance_id:
                log.notify('error occurred during setup, attemption instance termination')
                log.exception()
                try:
                    self.terminate()
                except ClientError:
                    pass
            raise
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        if self.synchronous:
            self.terminate()
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
        filename = '{}{!s}_{}.p'.format(
            os.environ['TMPDIR'], random.randint(0, 1e9), function.__name__)

        with open(filename, 'wb') as f:
            dill.dump(dict(function=function, args=args, kwargs=kwargs), f)
        return filename

    # todo this doesn't work; it gets THIS module's imports, but not the calling module!
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
        script_name = '{}{!s}_{}.py'.format(
            os.environ['TMPDIR'], random.randint(0, 1e9), function.__name__)
        script_body = (
            '{imports}'
            'with open("func.p", "rb") as fin:\n'
            '    data = dill.load(fin)\n'
            'results = data["function"](*data["args"], **data["kwargs"])\n'
            'with open("results.p", "wb") as f:\n'
            '    dill.dump(results, f)\n'
        )
        script_body = script_body.format(imports=cls.format_importlist())

        with open(script_name, 'w') as f:
            # log.notify('writing script to file:\n%s' % script_body)
            f.write(script_body)
        return script_name

    def __call__(self, function):

        def function_executed_on_aws(*args, **kwargs):

            # dump original function to file
            script = self.write_script(function)
            func = self.pickle_function(function, args, kwargs)

            # create an instance, or ensure the passed instance has the necessary
            # packages installed
            self.setup_seqc()

            with SSHConnection(self.instance_id, self.rsa_key) as ssh:
                ssh.put_file(script, 'script.py')
                ssh.put_file(func, 'func.p')
                if self.synchronous:
                    ssh.execute('python3 script.py')
                    results_name = os.environ['TMPDIR'] + function.__name__ + '_results.p'
                    ssh.get_file('results.p', results_name)
                    with open(results_name, 'rb') as f:
                        results = dill.load(f)
                else:
                    ssh.execute('nohup python3 script.py > nohup.log 2>&1 &')
                    results = None

            if self.synchronous:
                self.terminate()

            return results

        return function_executed_on_aws


class SSHConnection:
    _error_msg = ('You need to specify a valid RSA key to connect to Amazon EC2 '
                  'instances, see https://github.com/ambrosejcarr/seqc#create-an-rsa-key'
                  '-to-allow-you-to-launch-a-cluster')

    ec2 = boto3.resource('ec2')

    def __init__(self, instance_id, rsa_key):
        if not isinstance(instance_id, str):
            raise ValueError('instance must be a string instance id')
        if not instance_id.startswith('i-'):
            raise ValueError('valid instance identifiers must start with "i-"')
        self._instance_id = instance_id
        self.rsa_key = os.path.expanduser(rsa_key)

        self.check_key_file()
        self.ssh = paramiko.SSHClient()
        self.ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())

    # todo define def __repr__(self):

    @property
    def instance_id(self):
        return self._instance_id

    @instance_id.setter
    def instance_id(self, value):
        if isinstance(value, str):
            raise ValueError('instance must be a string instance id')
        if not value.startswith('i-'):
            raise ValueError('valid instance identifiers must start with "i-"')
        self._instance_id = value

    def check_key_file(self):
        """Checks the rsa file is present"""
        if not self.rsa_key:
            log.notify('The key %s was not found!' % self.rsa_key)
            raise FileNotFoundError(self._error_msg, 'The key file %s does not exist' %
                                    self.rsa_key)

    @Retry(retries=40, delay=2.5, catch=(NoValidConnectionsError, socket.error))
    def connect(self):
        """connects to a remote instance"""
        instance = self.ec2.Instance(self.instance_id)
        try:
            self.ssh.connect(instance.public_dns_name, username='ec2-user',
                             key_filename=self.rsa_key, timeout=3.0)
        except NoValidConnectionsError:
            state = instance.state['Name']
            if state not in ['running', 'pending']:
                raise InstanceNotRunningError(
                    'instance %s in state %s. Only running instances can be connected to.'
                    % (self.instance_id, state))
            else:
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
            ftp.get(remote_file, local_file)

    def put_file(self, local_file, remote_file):
        """places a file from the local machine onto a remote instance

        :param local_file: name of file to be copied to remote instance
        :param remote_file: name of file placed remotely
        """

        if not self.is_connected():
            self.connect()
        with closing(self.ssh.open_sftp()) as ftp:
            ftp.put(local_file, remote_file)
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
        if errs:
            raise ChildProcessError('\n'.join(errs))
        return data, errs

    def login_command(self):
        instance = self.ec2.Instance(self.instance_id)
        return ('ssh -i {rsa_path} ec2-user@{dns_name}'.format(
            rsa_path=self.rsa_key, dns_name=instance.public_ip_address))

    def __enter__(self):
        log.notify('connecting to instance %s via ssh' % self.instance_id)
        self.connect()
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.disconnect()
        if not exc_type:
            return True


class instance_clean_up:

    # todo in the pipeline, self.terminate is always True
    def __init__(self, email=None, upload=None, log_name='seqc.log', terminate=True,
                 debug=False):
        """Execution context for on-server code execution with defined clean-up practices.

        This is the cognate context manager to aws_setup, and when used with aws_setup(),
        ensures that all errors are captured and, if desirable, instances can be properly
        terminated.

        usage:
        ------
        with aws_execute(email=addr, upload=s3://my_bucket/my_key/, terminate=True):
            <execute code here>

        :param email: email address to send the log containing execution summary and any
          errors
        :param upload: s3 location to upload the log to
        :param log_name: name of the log object
        :param terminate: if True, terminate the instance upon completion, provided that
          no errors occurred.
        :param debug: if True, instance is not terminated when an error is raised
        """
        self.email = email
        self.log_name = log_name  # changed mount point to ~, cwd is good.
        self.terminate = terminate  # only terminate if no errors occur
        self.aws_upload_key = upload
        self.err_status = False
        self.mutt = verify.executables('mutt')[0]  # unpacking necessary for singleton
        self.debug = debug

    @staticmethod
    def email_user(attachment: str, email_body: str, email_address: str) -> None:
        """
        sends an email to email address with text contents of email_body and attachment
        attached. Email will come from "ec2-user@<ec2-instance-ip-of-aws-instance>

        :param attachment: the file location of the attachment to append to the email
        :param email_body: text to send in the body of the email
        :param email_address: the address to which the email should be sent"""

        if isinstance(email_body, str):
            email_body = email_body

        email_args = (
            'echo "{b}" | mutt -e "set content_type="text/html"" -s '
            '"Remote Process" {e} -a "{a}"'.format(
                b=email_body, a=attachment, e=email_address))
        email_process = Popen(email_args, shell=True, stderr=PIPE, stdout=PIPE)
        out, err = email_process.communicate(email_body)
        if err:
            raise ChildProcessError(err.decode())

    def __enter__(self):
        pass
        # log.setup_logger(self.log_name)
        # log.notify('Beginning protected execution')  # todo only run if verbose

    @staticmethod
    def _get_instance_id():
        """get an aws instances id from it's private ip address"""
        p = Popen('/sbin/ifconfig eth0 | grep "inet addr"', shell=True, stdout=PIPE,
                  stderr=PIPE)
        ip, err = p.communicate()
        if err:  # not an ec2 linux instance, nothing to terminate
            return
        else:
            ip = ip.decode().strip().split()[1].replace('addr:', '')
        ec2 = boto3.resource('ec2')
        try:  # if a non-ec2 instance, no instance will pass filter
            instances = ec2.instances.filter(
                Filters=[{'Name': 'private-ip-address', 'Values': [ip]}])
        except StopIteration:
            return
        return next(iter(instances)).id

    def __exit__(self, exc_type, exc_val, exc_tb):
        """If an exception occurs, log the exception, email if possible, then terminate
        the aws instance if requested by the user

        :param exc_type: type of exception encountered
        :param exc_val: value of exception
        :param exc_tb: exception traceback
        """

        # log any exceptions, set email body based on error / terminate status

        if exc_type is not None:
            log.exception()
            email_body = 'Process interrupted -- see attached error message'
        elif self.terminate:
            email_body = 'Process completed successfully -- see attached log'
            log.info('Execution completed successfully, instance will be terminated.')
        else:
            email_body = 'Process completed successfully -- see attached log'
            log.info('Execution completed successfully, but user requested no '
                     'termination. Instance will continue to run.')

        # todo this is the source of the second email for successful runs
        # email user if possible; catch exceptions if email fails.
        if self.email and self.mutt:
            log.notify('Emailing user.')
            try:
                self.email_user(
                    attachment=self.log_name, email_body=email_body,
                    email_address=self.email)
            except ChildProcessError:
                log.exception()

        # upload data if requested
        if self.aws_upload_key:
            log.notify('Uploading log to {}'.format(self.aws_upload_key))
            bucket, key = io.S3.split_link(self.aws_upload_key)

            @Retry(catch=Exception)
            def upload_file():
                io.S3.upload_file(self.log_name, bucket, key)
            upload_file()

        # terminate if no errors and debug is False
        if self.terminate:
            if exc_type and self.debug:
                return  # don't terminate if an error was raised and debug was set
            instance_id = self._get_instance_id()
            if instance_id is None:
                return  # todo notify if verbose
            ec2 = boto3.resource('ec2')
            instance = ec2.Instance(instance_id)
            log.notify('instance %s termination requested. If successful, this is the '
                       'final log entry.' % instance_id)
            instance.terminate()
            instance.wait_until_terminated()


def remove_inactive_security_groups():
    """security groups cannot be automatically cleaned up when instances self-terminate.

    This function finds all inactive security groups. Note that it is NOT limited to your
    user account
    """
    ec2 = boto3.resource('ec2')
    for s in ec2.security_groups.all():
        try:
            s.delete()
        except ClientError:  # is in use
            pass


def check_bucket(s3_uri):
    """Check whether a bucket exists and you have access.

    :param str s3_uri: name of uri in a bucket to check
    """
    if not s3_uri.startswith('s3://'):
        raise ValueError('%s is not a valid s3 URI' % s3_uri)
    bucket = s3_uri[5:].split('/')[0]
    s3 = boto3.resource('s3')
    try:
        s3.meta.client.head_bucket(Bucket=bucket)
    except ClientError:
        raise ValueError('Bucket %s for s3 URI %s does not exist' % (bucket, s3_uri))