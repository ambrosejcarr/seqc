from starcluster.clustersetup import ClusterSetup
from starcluster.logger import log
from subprocess import check_output, call


class DownloadRepo(ClusterSetup):
    def __init__(self, dir_name='/data/software/'):
        if not dir_name.endswith('/'):
            dir_name += '/'
        self.dir_name = dir_name

    def run(self, nodes, master, user, user_shell, volumes):
        folder = self.dir_name
        log.info('installing seqc repo onto %s' % folder)

        master.ssh.execute("mkdir %s" % folder)
        location = folder + "seqc.tar.gz"
        master.ssh.execute(
            'curl -H "Authorization: token a22b2dc21f902a9a97883bcd136d9e1047d6d076" -L '
            'https://api.github.com/repos/ambrosejcarr/seqc/tarball > %s' % location)
        log.info("seqc.tar.gz has been downloaded in /data/software directory")
        master.ssh.execute('pip3 install %s' % location)
