import boto3


def instances(args):
    """list instances and return

    :param args: namespace object from argparse, must contain args.rsa_key, the path to
      the rsa-key used to start the instances you want to list
    :return None:
    """

    keyname = args.rsa_key.rpartition('.')[0].rpartition('/')[-1]

    ec2 = boto3.resource('ec2')
    all_instances = ec2.instances.filter(
        Filters=[
            {'Name': 'key-name',
             'Values': [keyname]}])
    for i in all_instances.all():
        print('id: %s, type: %s, launch-time: %s, state: %s, ip %s' % (
            i.id, i.instance_type, str(i.launch_time), i.state, i.public_ip_address))
