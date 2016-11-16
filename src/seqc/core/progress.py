from subprocess import Popen, PIPE
from seqc import ec2


def progress(args):
    """print progress of requested seqc run(s) to less

    :param args: namespace object from argparse, must include rsa-key and instance-id
    :return None:
    """
    for id_ in args.instance_ids:
        connection = ec2.SSHConnection(id_, args.rsa_key)
        out, err = connection.execute('cat ./seqc_log.txt')
        p = Popen(['less'], stdin=PIPE)
        p.communicate(input='\n'.join(out).encode())
