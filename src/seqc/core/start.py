from seqc import ec2


def start(args):
    """start an aws instance"""
    instance = ec2.AWSInstance(
        rsa_key=args.rsa_key, instance_type=args.instance_type, spot_bid=args.spot_bid,
        volume_size=args.volume_size)
    instance.start()
