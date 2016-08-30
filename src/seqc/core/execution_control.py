import os
from seqc import log, io, exceptions, remote
from seqc.core import verify


class cleanup:

    def __init__(self, args):
        """Execute the seqc code with defined cleanup practices"""
        self.args = args
        self.err_status = False
        self.email = verify.executables('mutt')[0]  # unpacking necessary for singleton
        self.aws_upload_key = args.output_stem
        print('initialized cleanup execution context')

    def __enter__(self):
        """No entrance behavior is necessary to wrap the main function"""
        print('entering cleanup execution context')
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
