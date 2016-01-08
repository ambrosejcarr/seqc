#!/usr/local/bin/python3
__author__ = 'kristychoi'

import boto3
from subprocess import call, check_output
import re
import seqc
import configparser
import os

ec2 = boto3.resource('ec2')
line = ' '.join(check_output(['grep', '-i', 'i-[a-z0-9]', 'seqc.log']).decode().split())
result = re.search('i-[a-z0-9]+', line)
if not result:
    print('')
else:
    inst_id = result.group(0)
    instance = ec2.Instance(inst_id)
    hostname = 'ubuntu@' + instance.public_dns_name
    config_file = '/'.join(seqc.__file__.split('/')[:-3]) + '/config'
    config = configparser.ConfigParser()
    config.read(config_file)
    keypath = os.path.expanduser(config['key']['rsa_key_location'])
    call(['ssh', '-i', keypath, hostname])
