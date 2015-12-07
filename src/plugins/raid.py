from starcluster.clustersetup import ClusterSetup
from starcluster.logger import log
import string
import re
import time
import os
from subprocess import check_output, call


class RaidAutomator(ClusterSetup):
    def __init__(self, n_tb=4, availability_zone='us-east-1e'):
        self.n_tb = n_tb
        self.availability_zone = availability_zone
        log.debug('total number of terabytes = %s' % n_tb)

    def run(self, nodes, master, user, user_shell, volumes):
        log.info("setting up RAID array on master node...")
        vol_num = self.n_tb
        vol_size = 1024
        inst_id = master.id
        dev_base = "/dev/xvd"
        alphabet = string.ascii_lowercase[5:]  # starts at f
        dev_names = []
        vol_names = []
        availability_zone = self.availability_zone

        for i in range(int(vol_num)):
            log.info("creating volume %s of %s..." % (i + 1, vol_num))
            vol = check_output(
                ["aws", "ec2", "create-volume", "--size", str(vol_size), "--region",
                 "us-east-1", "--availability-zone", availability_zone, "--volume-type",
                 "gp2"])
            vol_list = vol.decode().split()
            vol_id = [s for s in vol_list if 'vol-' in s][0]
            vol_names.append(vol_id)

            # generating new device id to mount
            dev_id = dev_base + alphabet[i]
            dev_names.append(dev_id)

            log.info("waiting for volume to become available...")
            vol_stat = check_output(
                ["aws", "ec2", "describe-volumes", "--volume-ids", vol_id])
            vol_stat = vol_stat.decode().split()

            while vol_stat[6] != 'available':
                time.sleep(5)
                vol_stat = check_output(["aws", "ec2", "describe-volumes", "--volume-ids", vol_id])
                vol_stat = vol_stat.decode().split()
                log.info('...')

            log.info("attaching new volume...")
            call(["aws", "ec2", "attach-volume", "--volume-id", vol_id, "--instance-id", inst_id, "--device", dev_id])

            log.info("waiting for volume to finish attaching...")
            vol_stat = check_output(["aws", "ec2", "describe-volumes", "--volume-ids", vol_id])
            vol_stat = vol_stat.decode().split()

            while vol_stat[14] != 'attached':
                time.sleep(5)
                vol_stat = check_output(["aws", "ec2", "describe-volumes", "--volume-ids", vol_id])
                vol_stat = vol_stat.decode().split()
                log.info("...")
            # configure volumes here to delete on termination :(
            call(["aws", "ec2", "modify-instance-attribute", "--instance-id", inst_id, "--block-device-mappings",
                  "[{\"DeviceName\": \"" + dev_id + "\",\"Ebs\":{\"DeleteOnTermination\":true}}]"])
            i += 1

        log.info("successfully attached %s TB in %s volumes!" % (self.n_tb, vol_num))
        log.info("creating logical RAID device...")
        all_dev = ' '.join(dev_names)

        #writing names to file for future volume cleanup
        home = os.path.expanduser("~")
        fpath = home + "/.starcluster/plugins/vol_names.txt"
            with open(fpath,'w') as f:
            for x in vol_names:
                f.write('%s\n' %x)

        master.ssh.execute(
            "sudo mdadm --create --verbose /dev/md0 --level=0 --name=my_raid "
            "--raid-devices=%s %s" % (vol_num, all_dev))
        master.ssh.execute("sudo mkfs.ext4 -L my_raid /dev/md0")
        master.ssh.execute("sudo mkdir -p /data")
        master.ssh.execute("sudo mount LABEL=my_raid /data")
        log.info("successfully created RAID array in /data!")
