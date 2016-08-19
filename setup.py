import os
import sys
import shutil
from subprocess import call
from setuptools import setup
from warnings import warn

if sys.version_info.major != 3:
    raise RuntimeError('SEQC requires Python 3')
if sys.version_info.minor < 5:
    warn('Multiprocessing analysis methods may not function on Python versions < 3.5')

# install phenograph
call(['pip3', 'install', 'git+https://github.com/jacoblevine/phenograph.git'])

# get version
with open('src/seqc/version.py') as f:
    exec(f.read())

setup(name='seqc',
      version=__version__,  # read in from the exec of version.py; ignore error
      description='Single Cell Sequencing Processing and QC Suite',
      author='Ambrose J. Carr',
      author_email='mail@ambrosejcarr.com',
      package_dir={'': 'src'},
      packages=['seqc', 'seqc.sequence', 'seqc.alignment'],
      install_requires=[
          'numpy>=1.10.0',
          'awscli',
          'cython>0.14',
          'numexpr>=2.4',
          'pandas>=0.18.1',
          'paramiko>=2.0.2',
          'regex',
          'requests',
          'nose2',
          'scipy>=0.14.0',
          'boto3',
          'intervaltree',
          'matplotlib',
          'tinydb',
          'tables',
          'seaborn',
          'fastcluster',
          'statsmodels',
          'ecdsa',
          'pycrypto',
          'scikit_learn>=0.17'],
      scripts=['src/seqc/process_experiment.py'],
      extras_require={
          'GSEA_XML': ['html5lib', 'lxml', 'BeautifulSoup4'],
      }
      )

# look for star
if not shutil.which('STAR'):
    warn('SEQC: STAR is not installed. SEQC will not be able to align files.')

# get location of setup.py
setup_dir = os.path.dirname(os.path.realpath(__file__))
tools_dir = os.path.expanduser('~/.seqc/tools')
seqc_dir = os.path.expanduser('~/.seqc/seqc')

print('setup_dir: {}'.format(setup_dir))
print('tools_dir: {}'.format(tools_dir))
print('seqc_dir: {}'.format(seqc_dir))

if os.path.isdir(tools_dir):
    shutil.rmtree(tools_dir)
if os.path.isdir(seqc_dir):
    shutil.rmtree(seqc_dir)


def ignore_test_and_tools(dir_, files):
    """do not copy the test directory"""
    return [f for f in files if (f == 'test' or f == 'tools')]

# install tools and a local copy of seqc.
shutil.copytree(setup_dir + '/tools/', tools_dir)
shutil.copytree(setup_dir, seqc_dir, ignore=ignore_test_and_tools)  # copy seqc repository
shutil.make_archive(seqc_dir, 'gztar', seqc_dir)
shutil.unpack_archive(tools_dir + '/mouse_gene_sets.tar.gz', tools_dir)
shutil.unpack_archive(tools_dir + '/human_gene_sets.tar.gz', tools_dir)
