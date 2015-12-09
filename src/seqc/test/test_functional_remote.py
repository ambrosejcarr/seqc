
import unittest
import seqc
import nose2
import os
from copy import copy
import sys
import json
import logging
from seqc.test.test_unit import config

logging.getLogger('nose2').setLevel(logging.CRITICAL)

# set the input args as they would be retrieved from argv


class TestFunctionalRemote(unittest.TestCase):

    args = [
        'in-drop',
        '--reverse', 's3://dplab-data/seqc/test_seqc/in_drop/fastq/reverse/'
                     'test_seqc_r2.fastq',
        '--forward', 's3://dplab-data/seqc/test_seqc/in_drop/fastq/forward/'
                     'test_seqc_r1.fastq',
        '--aws-upload-key', 's3://dplab-data/seqc/test_seqc/',
        '--n-threads', '30',
        '--index', 's3://dplab-data/genomes/mm38_chr19/',
        '--output-prefix', '/data/software/output',
        '--frag-len', '1000',
        '--email-status', 'ambrosejcarr@gmail.com',
        '--barcodes', 's3://dplab-data/barcodes/in_drop/serial/barcodes.p',
        '--remote'
    ]

    if args[16] is None:
        raise ValueError('Please input user email for testing on line 30 of '
                         'test_functional_remote.py, then rerun the test.')

    def test_seqc_raw_fastq_input_remote_test(self):

        def initialize_logging():
            arg_copy = copy(kwargs)
            del arg_copy['func']  # function is not serializable

            seqc.log.info('SEQC version: %s' % seqc.__version__)
            seqc.log.info('SEQC working directory: %s' % os.getcwd())
            seqc.log.info('Passed command line arguments: %s' %
                          json.dumps(arg_copy, separators=(',', ': '), indent=4))

        parser = seqc.core.create_parser()
        kwargs = seqc.core.parse_args(parser, self.args)
        seqc.log.setup_logger()

        if kwargs['remote']:
            try:
                # log command line arguments for debugging
                initialize_logging()
                seqc.log.info('Starting REMOTE run')

                # run remotely
                del kwargs['remote']
                seqc.core.run_remote(kwargs)
            except:
                seqc.log.exception()
                raise

        elif kwargs['email_status']:
            try:
                # log command line arguments for debugging
                arg_copy = copy(kwargs)
                initialize_logging()

                # run locally
                func = kwargs['func']
                func(**kwargs)
                seqc.cluster_utils.upload_results(
                    kwargs['output_prefix'], kwargs['email_status'], kwargs['aws_upload_key'])
            except:
                seqc.log.exception()
                email_body = 'Process interrupted -- see attached error message'
                seqc.cluster_utils.email_user(attachment='seqc.log', email_body=email_body,
                                              email_address=kwargs['email_status'])
                sys.exit(1)
            finally:
                # we are currently within a cluster:
                if not kwargs['no_terminate']:
                    if os.path.isfile('/data/software/instance.txt'):
                        with open('/data/software/instance.txt','r') as f:
                            inst_id = f.readline().strip('\n')
                        seqc.cluster_utils.terminate_cluster(inst_id)
                    else:
                        seqc.log.info('file containing instance id is unavailable!')
                else:
                    seqc.log.info('not terminating cluster -- user responsible for cleanup')
        else:  # not email_status and not remote
            try:
                # log command line arguments for debugging
                arg_copy = copy(kwargs)
                initialize_logging()

                # run locally
                func = kwargs['func']
                func(**kwargs)
            except:
                seqc.log.exception()
                raise


if __name__ == "__main__":
    seqc.log.setup_logger('seqc_functest_remote.log')
    try:
        nose2.main()#exit=False, module=__name__, verbosity=5)
    except:
        seqc.log.exception()
        raise
