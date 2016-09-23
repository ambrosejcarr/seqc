import os
import dill  # can pickle lambdas
import types
import boto3
from subprocess import check_output
from seqc import log, io, exceptions, remote
from functools import wraps
from seqc.core import verify


class aws_execute:

    def __init__(self, email=None, upload=None, log_name='seqc.log', terminate=True):
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
        """
        self.email = email
        self.log_name = '/data/' + log_name  # todo think about this!
        self.terminate = terminate  # only terminate if error occurs
        self.aws_upload_key = upload
        self.err_status = False
        self.mutt = verify.executables('mutt')[0]  # unpacking necessary for singleton

    def __enter__(self):
        log.setup_logger(self.log_name)

    @staticmethod
    def _get_instance_id():
        """get an aws instances id from it's private ip address"""
        ip = check_output(
            '/sbin/ifconfig eth0 | grep "inet addr"'
            ).decode().strip().split()[0].replace('inet addr:', '')
        ec2 = boto3.resource('ec2')
        instances = ec2.instances.filter(
            Filters=[{'Name': 'private-ip-address', 'Values': [ip]}])
        return next(iter(instances)).id  # todo add assertion that len(instances == 1)

    def __exit__(self, exc_type, exc_val, exc_tb):
        """If an exception occurs, log the exception, email if possible, then terminate
        the aws instance if requested by the user

        :param exc_type: type of exception encountered
        :param exc_val: value of exception
        :param exc_tb: exception traceback
        """

        # log any exceptions, set email body based on error / terminate status
        if issubclass(exc_type, BaseException):
            log.exception()
            email_body = 'Process interrupted -- see attached error message'
        elif self.terminate:
            email_body = 'Process completed successfully -- see attached log'
            log.info('Execution completed successfully, instance terminated.')
        else:
            email_body = 'Process completed successfully -- see attached log'
            log.info('Execution completed successfully, but user requested no '
                     'termination. Instance will continue to run.')

        # email user if possible
        if self.email and self.mutt:
            remote.email_user(attachment=self.log_name, email_body=email_body,
                              email_address=self.email)

        # upload data if requested
        if self.aws_upload_key:
            bucket, key = io.S3.split_link(self.aws_upload_key)
            exceptions.retry_boto_call(io.S3.upload_file)(
                self.log_name, bucket, key)

        # terminate if requested and no errors
        if self.terminate and exc_type:
            remote.terminate_cluster(self._get_instance_id())


class aws_setup:

    def __init__(self, inst_type=None, spot_bid=None, volsize=None, terminate=True):
        """Create an aws instance for the remote execution of passed command strings.
        Unless requested, the instance is terminated when this context is exited. Note
        that this is NOT desirable in instances where asynchronous commands are being
        passed

        usage:
        ------
        with aws_setup(inst_type='c4', spot_bid=1.0, volsize=5, terminate=True) as inst:
            inst.put('important_local_file', '/data/important_file_now_on_remote')
            inst.execute('script_name.sh')

        :param inst_type: type of aws instance (options: 'c4', 'c3', 'r3')
        :param spot_bid: amount of money to bid per hour in dollars for instance. If None,
          the instance is reserved
        :param volsize: size of volume to be mounted to the instance at '/data'
        """

        self.instance_type = inst_type
        self.spot_bid = spot_bid
        self.volsize = int(volsize)
        self.terminate = terminate
        self.cluster = None

    def __enter__(self):
        """create an instance and make it accessible within the context"""

        # create an aws instance
        cluster = remote.ClusterServer()
        try:
            cluster.setup_cluster(  # todo split up this call for more expressive errors during setup
                self.volsize, self.instance_type, spot_bid=self.spot_bid)
            self.cluster = cluster
        except Exception as e:
            log.notify('Exception occurred during instance setup!'.format(e=e))
            log.exception()
            raise

        # connect to the instance
        try:
            cluster.serv.connect()
        except:
            log.notify('Could not connect to instance!')
            log.exception()
            raise

        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """terminate the instance if requested, log all exceptions"""
        if exc_type:  # error occurred during setup
            log.notify('Exception occurred during remote execution, within aws_setup '
                       'environment.')
            log.exception()

        if self.terminate:
            log.notify('Terminating the instance as requested.')
            self.cluster.serv.disconnect()
            remote.terminate_cluster(self.cluster.inst_id.instance_id)

    def execute(self, command_string):
        """run the remote function on the server, capturing output and errors"""
        data, errs = self.cluster.serv.exec_command(command_string)
        if errs:
            raise ChildProcessError('Error captured from remote execution: %s' % errs)
        return data

    def put_file(self, local_file, remote_file):
        self.cluster.serv.put_file(local_file, remote_file)

    def get_file(self, remote_file, local_file):
        self.cluster.serv.get_file(local_file, remote_file)


class Remote:

    # todo add verbosity (don't print all the log outputs)
    def __init__(self, inst_type='c4', volsize=100, spot_bid=None, retrieve=False,
                 log_name='seqc.log'):
        """function decorator to synchronously execute the decorated function on an aws
        instance. Instance is terminated upon exit from the context environment

        :param inst_type: type of instance to start (default c4.8xlarge, options: 'c4',
          'c3', 'r3')
        :param volsize: size of volume in GB to mount to the instance's /data directory
        :param spot_bid: amount of money to bid per hour in dollars for instance. If None,
          the instance is reserved
        :param retrieve: whether or not the output of the called function should be
          retrieved and returned. Default False (many functions upload the results to
          aws s3.)
        :param log_name: name of the log file to generate
        """
        self.inst_type = inst_type
        self.volsize = volsize
        self.spot_bid = spot_bid
        self.retrieve = retrieve
        self.log_name = log_name

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
        log.setup_logger(self.log_name)

        @wraps(function)
        def wrapped(*args, **kwargs):
            script = self.write_script(function)
            func = self.pickle_function(function, args, kwargs)
            with aws_setup(self.inst_type, self.spot_bid, self.volsize,
                           terminate=True) as s:
                s.put_file(script, '/data/script.py')
                s.put_file(func, '/data/func.p')
                s.execute('python3 /data/script.py')
                if self.retrieve:
                    results_name = os.environ['TMPDIR'] + function.__name__ + '_results.p'
                    s.get_file('/data/results.p', results_name)
                    with open(results_name, 'rb') as f:
                        return dill.load(f)
        return wrapped


class AsyncRemote(Remote):

    def __init__(self, email, terminate, upload, *args,
                 **kwargs):
        """function decorator to asynchronously execute the decorated function on an aws
        instance.

        :param inst_type: type of instance to start (default c4.8xlarge, options: 'c4',
          'c3', 'r3')
        :param volsize: size of volume in GB to mount to the instance's /data directory
        :param spot_bid: amount of money to bid per hour in dollars for instance. If None,
          the instance is reserved
        :param retrieve: whether or not the output of the called function should be
          retrieved and returned. Default False (many functions upload the results to
          aws s3.)
        :param log_name: name of the log file to generate
        :param email: email address to email results
        :param terminate: if True, instance will be terminated when passed code completes
          unless errors occur
        :param upload: an aws s3 link where the remote log will be uploaded
        """
        super().__init__(*args, **kwargs)
        self.email = email
        self.terminate = terminate
        self.upload = upload

    def write_script(self, function) -> str:
        """generate a python script that calls function after importing
        required modules

        :param object function: function to be called
        :return str: filename of the python script
        """
        script_name = '{}{}.py'.format(os.environ['TMPDIR'], function.__name__)
        script_body = (
            '{imports}'
            'from seqc.core import execution_context\n'
            'with open("/data/func.p", "rb") as fin:\n'
            '    data = dill.load(fin)\n'
            'with execution_context.aws_execute('
            'email={email}, upload={upload}, '
            'log_name={log_name}, terminate={terminate}):\n'
            '    results = data["function"](*data["args"], **data["kwargs"])\n'
            '    with open("/data/results.p", "wb") as f:\n'
            '        dill.dump(results, f)\n'
        )
        script_body = script_body.format(
            imports=self.format_importlist(),
            email=self.email,
            upload=self.upload,
            log_name=self.log_name,
            terminate=self.terminate
        )

        with open(script_name, 'w') as f:
            print('writing script to file:\n%s' % script_body)
            f.write(script_body)
        return script_name

    def __call__(self, function):
        log.setup_logger(self.log_name)

        @wraps(function)
        def wrapped(*args, **kwargs):
            script = self.write_script(function)
            func = self.pickle_function(function, args, kwargs)
            with aws_setup(self.inst_type, self.spot_bid, self.volsize,
                           terminate=False) as s:
                s.put_file(script, '/data/script.py')
                s.put_file(func, '/data/func.p')
                s.execute('nohup python3 /data/script.py > /dev/null 2>&1 &')
            log.notify('Code executing on remote instance. You will be emailed at {} '
                       'when function completes. Instance will be automatically '
                       'terminated. '.format(self.email))

        return wrapped
