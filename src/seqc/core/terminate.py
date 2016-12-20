import boto3
from botocore.exceptions import ClientError


def terminate(args):
    """print progress of requested seqc run to top

    :param args: namespace object from argparse, must include rsa-key and instance-id
    :return None:
    """
    ec2 = boto3.resource('ec2')
    for id_ in args.instance_ids:
        instance = ec2.Instance(id=id_)
        try:
            response = instance.terminate()
            print('termination signal sent:\n%s' % response)
        except ClientError:
            print('instance %s does not exist')
