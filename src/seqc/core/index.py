
def index(args):
    """create an index for SEQC.

    :param args: parsed arguments. This function is only called if subprocess_name is
      'index'
    """

    # functions to be pickled and run remotely must import all their own modules
    from seqc import ec2, log
    from seqc.sequence.index import Index

    log.setup_logger(args.log_name)
    with ec2.instance_clean_up(args.email, args.upload, log_name=args.log_name):
        idx = Index(args.organism, args.additional_id_types)
        idx.create_index(args.upload_location)
