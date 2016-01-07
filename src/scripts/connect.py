#!/usr/local/bin/python3

import boto3
from subprocess import call
import seqc
import configparser
import os
import sys

inst_file = '/'.join(seqc.__file__.split('/')[:-3]) + '/src/scripts/instance.txt'
try:
    with open(inst_file, 'r') as f:
        inst_id = f.readline().strip('\n')
except FileNotFoundError:
    print('Error: Please check that the cluster was set up appropriately!')
    sys.exit(1)

# extract hostname and key information
ec2 = boto3.resource('ec2')
instance = ec2.Instance(inst_id)
hostname = 'ubuntu@' + instance.public_dns_name
config_file = '/'.join(seqc.__file__.split('/')[:-3]) + '/config'
config = configparser.ConfigParser()
config.read(config_file)
keypath = os.path.expanduser(config['key']['rsa_key_location'])

# connect to instance
call(['ssh', '-i', keypath, hostname])
