from subprocess import Popen, PIPE
from seqc import ec2
from paramiko.ssh_exception import AuthenticationException
from botocore.exceptions import ClientError


def progress(args):
    """print progress of requested seqc run(s) to less

    :param args: namespace object from argparse, must include rsa-key and instance-id
    :return None:
    """
    if args.rsa_key is None:
        raise ValueError('User must supply -k/--rsa-key or set the environment variable '
                         'AWS_RSA_KEY')

    if args.instance_ids is None:
        raise ValueError('No instances specified. Please supply an instance using the -i '
                         'parameter.')

    for id_ in args.instance_ids:
        connection = ec2.SSHConnection(id_, args.rsa_key)
        try:
            out, err = connection.execute('cat ./seqc_log.txt')
        except AuthenticationException:
            raise ValueError('instance %s cannot be found.' % repr(id_))
        except ClientError:
            raise ValueError('instance %s cannot be found.' % repr(id_))
        p = Popen(['less'], stdin=PIPE)
        p.communicate(input='\n'.join(out).encode())
