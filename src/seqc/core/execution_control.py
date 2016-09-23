import os
import dill  # can pickle lambdas
import types
from seqc import log, io, exceptions, remote
from functools import wraps
from seqc.core import verify
from math import ceil
from warnings import warn


class local_instance_cleanup:

    def __init__(self, args):
        """Execute the seqc code on an instance with defined cleanup practices"""
        self.args = args
        self.err_status = False
        self.email = verify.executables('mutt')[0]  # unpacking necessary for singleton
        self.aws_upload_key = args.output_stem

    def __enter__(self):
        """No entrance behavior is necessary to wrap the main function"""
        return

    def __exit__(self, exc_type, exc_val, exc_tb):
        """If an exception occurs, log the exception, email if possible, then terminate
        the aws instance if requested by the user

        :param exc_type: type of exception encountered
        :param exc_val: value of exception
        :param exc_tb: exception traceback
        """

        # log any non-SystemExit exceptions if they were encountered, and email user if
        # specified in args.
        if issubclass(exc_type, Exception):
            log.exception()
            self.err_status = True  # we have encountered an error
            if self.args.email_status and self.email:
                email_body = 'Process interrupted -- see attached error message'
                if self.args.aws:
                    attachment = '/data/' + self.args.log_name
                else:
                    attachment = self.args.log_name
                if self.aws_upload_key:
                    bucket, key = io.S3.split_link(self.aws_upload_key)
                    exceptions.retry_boto_call(io.S3.upload_file)(
                        attachment, bucket, key)
                remote.email_user(attachment=attachment, email_body=email_body,
                                  email_address=self.args.email_status)

        # determine if instance needs clean up
        if self.args.remote:
            return True  # this is a remote run, there is no instance to terminate.

        # determine how to deal with termination in the event of errors
        if self.args.no_terminate == 'on-success':
            if self.err_status:
                no_terminate = 'True'
            else:
                no_terminate = 'False'
        else:
            no_terminate = self.args.no_terminate

        if no_terminate in ['False', 'false']:
            fpath = '/data/instance.txt'
            if os.path.isfile(fpath):
                with open(fpath, 'r') as f:
                    inst_id = f.readline().strip('\n')
                remote.terminate_cluster(inst_id)
            else:
                log.info('File containing instance id is unavailable!')
        else:
            log.info('no-terminate={}; cluster not terminated. User is responsible for '
                     'clean-up'.format(self.args.no_terminate))

        return True  # signals successful cleanup for contextmanager


class remote_execute:

    def __init__(self, instance_type=None, spot_bid=None, volsize=None):
        """Create a temporary cluster for the remote execution of passed command strings

        :param instance_type:
        :param spot_bid:
        :param volsize:
        """

        self.instance_type = instance_type
        self.spot_bid = spot_bid
        self.volsize = int(ceil(volsize / 1e9))
        self.cluster = None
        self.async_process = False

    def __enter__(self):
        """create a cluster and make it accessible within the context"""
        try:
            cluster = remote.ClusterServer()
            cluster.setup_cluster(
                self.volsize, self.instance_type, spot_bid=self.spot_bid)
            cluster.serv.connect()
            self.cluster = cluster
        except Exception as e:
            log.notify('Exception {e} occurred during cluster setup!'.format(e=e))
            raise
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """terminate the instance unless a clean asynchronous call was made"""

        if exc_type or not self.async_process:
            log.notify('%s: %s\n%s' % (exc_type, exc_val, exc_tb))
            self.cluster.serv.disconnect()
            remote.terminate_cluster(self.cluster.inst_id.instance_id)
        return True  # signal clean exit

    def execute(self, command_string):
        """run the remote function on the server, capturing output and errors"""
        data, errs = self.cluster.serv.exec_command(command_string)
        if errs:
            raise ChildProcessError(errs)
        return data

    def async_execute(self, command_string):
        """run the remote function on the server, signaling to remote_execute that there
        may be a process that remains running in the background after this function
        returns.

        If an async_execute function returns without error, instance termination is
        halted, and the user must manually terminate the instance with
        "SEQC.py instance terminate -i <instance id>"

        :param command_string: command to be remote executed. Assumed to be called with
          nohup, to prevent hangup when the ssh shell is disconnected, and to be placed
          in the background with "&", so that async_execute() returns immediately. If
          these requirements are not met, there may be undefined results, and the
          process will warn the user.
        """
        if not all((s in command_string for s in ['nohup', '&'])):
            warn('Excecuting command that may not be asynchronous. User is '
                 'responsible for including commands that place the called function '
                 'into the background. Missing "nohup" or "&". If a synchronous command '
                 'is desired, please use the execute() method.')
        self.cluster.serv.exec_command(command_string)
        self.async_process = True

    def put_file(self, local_file, remote_file):
        self.cluster.serv.put_file(local_file, remote_file)

    def get_file(self, remote_file, local_file):
        self.cluster.serv.get_file(remote_file, local_file)

""" list of requirements for a remote execution method/decorator/execution context:

must:
1. initiate a logger locally or remotely
2. call the function locally if remote is false
3. call the function remotely if remote is true
4. catch errors and email user if executing remotely
5. shut down instance if errors occur remotely
6. shut down instance if errors occur locally
7. be able to execute synchronously or asynchronously
8. track (with locked access!!!!) the instance ids in a local file such that the
   security group can later be cleaned up.
9. others?

a function decorator @remote could handle some of the above. It could:

1. initialize a logger when @remote is called
2. check if remote is True. If true, re-call the function with remote=False. if false,
   execute the function
3. The function itself could be called within an execution_context wherein errors are
   caught, and dealt with depending whether the function is executing locally or remotely

   if remote, start an instance, initiate the remote instance, call the function remotely
   using MPI/pickle somehow?

   problem: how to pickle a local function (and environment) and execute it in a remote
   context?
   - if possible to check imports, could pickle the function and all the arguments, then
     write a file which imports all the imports and calls the pickled function with the
     pickled arguments.

"""

# needs to not be used as a decorator, but rather as a function wrapping another function
# defined at the module level. THAT should work.

class Remote:

    def __init__(self, instance_type, volsize, spot_bid, retrieve=False,
                 logname='seqc.log'):
        """
        1. sets up a logger
        2. pickles the passed function

        :param instance_type:
        :param volsize:
        :param spot_bid:
        :param retrieve:
        :param logname:
        """
        self.instance_type = instance_type
        self.volsize = volsize
        self.spot_bid = spot_bid
        self.retrieve = retrieve
        self.logname = logname

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
        :param list args: positional arguments for the function
        :param dict kwargs: keyword arguments for the function
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
        script_body.format(imports=cls.format_importlist())

        with open(script_name, 'w') as f:
            f.write(script_body)
        return script_name

    def __call__(self, function):

        log.setup_logger(self.logname)

        @wraps(function)
        def wrapped(*args, **kwargs):
            script = self.write_script(function)
            func = self.pickle_function(function, args, kwargs)
            with remote_execute(self.instance_type, self.volsize, self.spot_bid) as s:
                s.put_file(script, '/data/script.py')
                s.put_file(func, '/data/func.p')
                s.execute('python3 /data/script.py')
                if self.retrieve:
                    results_name = os.environ['TMPDIR'] + function.__name__ + '_results.p'
                    s.get_file('/data/results.p', results_name)
                    with open(results_name, 'rb') as f:
                        return dill.load(f)
        return wrapped


# class asyncRemote(Remote):
#
#     @classmethod
#     def write_script(cls, function) -> str:
#         """generate a python script that calls function after importing required modules
#
#         :param object function: function to be called
#         :param list args: positional arguments for the function
#         :param dict kwargs: keyword arguments for the function
#         :return str: filename of the python script
#         """
#         script_name = '{}{}.py'.format(os.environ['TMPDIR'], function.__name__)
#         script_body = (
#             '{imports}'
#             'with open("/data/func.p", "rb") as fin:\n'
#             '    data = pickle.load(fin)\n'
#             'results = data["function"](*data["args"], **data["kwargs"])\n'
#             'with open("results.p", "wb") as f:\n'
#             '    pickle.dump(results, f)\n'
 #         )
#         script_body.format(imports=cls.format_importlist())
#
#         with open(script_name, 'w') as f:
#             f.write(script_body)
#         return script_name
#
#     def __call__(self, f):
#         def wrapped(*args, **kwargs):
#             script = self.write_script(f)
#             func = self.pickle_function(f, args, kwargs)
#             with remote_execute(self.instance_type, self.volsize, self.spot_bid) as s:
#                 s.put_file(script, '/data/script.py')
#                 s.put_file(func, '/data/func.p')
#                 if self.asynchronous:
#                     s.async_execute('python3 /data/script.py')
#
#
#             pass
#         return wrapped_f



# write a test before proceeding that tests:

# def wrapper(function):
#     pickle.dump(function)
#
# def function(x):
#     return x
#
# p = wrapper(function)
# p()




class Wrapper:

    def __init__(self, arg1, arg2):
        self.arg1 = arg1
        self.arg2 = arg2

    def __call__(self, f):
        with open('test_dill.p', 'wb') as fout:
            dill.dump(f, fout)

        def wrapper(f):
            print(f)
            return 10

        return wrapper
