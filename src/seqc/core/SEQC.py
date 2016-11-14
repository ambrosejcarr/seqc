#!/usr/local/bin/python3

import sys
from seqc import core
from seqc.core import parser, verify
from seqc import ec2


def main(argv):
    """Check arguments, then call the appropriate sub-module

    Created to allow the main pipeline to be tested from the earliest entry point
    (command-line arguments).

    :param argv: output of sys.argv[1:]
    """
    arguments = parser.parse_args(argv)
    func = getattr(core, arguments.subparser_name)
    assert func is not None
    if arguments.remote:
        # todo improve how verification works; it's not really necessary, what is needed
        # is a method to determine volume size for remote.
        verification_func = getattr(verify, arguments.subparser_name)
        verified_args = verification_func(arguments)
        remote_args = {
            k: getattr(verified_args, k) for k in
            ('rsa_key', 'instance_type', 'spot_bid', 'volume_size') if
            getattr(verified_args, k)}
        ec2.AWSInstance(synchronous=False, **remote_args)(func)(verified_args)
    else:
        func(arguments)


if __name__ == '__main__':
    main(sys.argv[1:])
