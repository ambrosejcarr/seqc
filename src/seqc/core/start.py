from seqc import ec2
import os


def start(args):
    """start an aws instance"""

    if args.rsa_key is None:
        raise ValueError('-k/--rsa-key does not point to a valid file object. ')
    if not os.path.isfile(args.rsa_key):
        raise ValueError('-k/--rsa-key does not point to a valid file object. ')

    instance = ec2.AWSInstance(
        rsa_key=args.rsa_key, instance_type=args.instance_type, spot_bid=args.spot_bid,
        volume_size=args.volume_size)
    instance.start()
