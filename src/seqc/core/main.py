#!/usr/local/bin/python3

import sys
from seqc import core
from seqc.core import parser, verify
from seqc import ec2
import boto3

def clean_up_security_groups():
    '''
    Cleanning all the unused security groups that were created/started using SEQC
    when the number of unused ones is greater than 300 
    '''
    ec2 = boto3.resource('ec2') 
    sgs = list(ec2.security_groups.all())
    insts = list(ec2.instances.all())
    
    all_sgs = set([sg.group_name for sg in sgs])                    # get all security groups                        
    all_inst_sgs = set([sg['GroupName'] 
        for inst in insts for sg in inst.security_groups])          # get security groups associated with instances
    unused_sgs = all_sgs - all_inst_sgs                             # get ones without instance association
    
    if len(unused_sgs) >= 300:
        print("Cleaning up the unused security groups:")
        client = boto3.client('ec2')
        for g in unused_sgs:
            all_inst_sgs = set([sg['GroupName'] for inst in insts for sg in inst.security_groups])      # since deleting ones takes a while, doublecheck whether 
            if g.startswith("SEQC") and (g not in all_inst_sgs):    # only cleaning ones associated with SEQC                                        # the security group is still unused
                client.delete_security_group(GroupName=g)
                print(g+" deleted")


def main(argv):
    """Check arguments, then call the appropriate sub-module

    Created to allow the main pipeline to be tested from the earliest entry point
    (command-line arguments).

    :param argv: output of sys.argv[1:]
    """
    arguments = parser.parse_args(argv)

    func = getattr(core, arguments.subparser_name)
    assert func is not None

    # notebooks execute local
    if arguments.subparser_name == 'notebook':
        return func(arguments)

    if arguments.remote:
        # todo improve how verification works; it's not really necessary, what is needed
        # is a method to determine volume size for remote.
        verification_func = getattr(verify, arguments.subparser_name)
        verified_args = verification_func(arguments)
        remote_args = {
            k: getattr(verified_args, k) for k in
            ('rsa_key', 'instance_type', 'spot_bid', 'volume_size') if
            getattr(verified_args, k)}
        clean_up_security_groups()
        ec2.AWSInstance(synchronous=False, **remote_args)(func)(verified_args)
    else:
        func(arguments)


if __name__ == '__main__':
    main(sys.argv[1:])
