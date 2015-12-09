__author__ = 'Ambrose J. Carr'

from setuptools import setup
from warnings import warn
import os
import shutil
import configparser

__version__ = '0.1.5'

# pip3 cannot install external dependencies for python; warn user if external dependencies
# are missing; do this at the end so that the users are more likely to see it.

# look in /usr/local/ and /usr/local/hdf5/ for hdf5 libraries;
# if found in /usr/local/hdf5/, set an environment variable to help pip3 install it.
h5fail = True
if os.path.isfile('/usr/lib/libhdf5.so'):
    h5file = False
elif os.path.isfile('/usr/local/lib/libhdf5.so'):
    h5fail = False
elif os.path.isfile('/usr/hdf5/lib/libhdf5.so'):
    os.environ['HDF5_DIR'] = '/usr/hdf5/'
elif os.path.isfile('/usr/local/hdf5/lib/libhdf5.so'):
    os.environ['HDF5_DIR'] = '/usr/local/hdf5/'
    h5fail = False

setup(name='seqc',
      version=__version__,
      description='Single Cell Sequencing Processing and QC Suite',
      author='Ambrose J. Carr',
      author_email='mail@ambrosejcarr.com',
      package_dir={'': 'src'},
      # test_suite=nose2.collector.collector,
      # note: requires numpy > 1.10.0
      packages=['seqc', 'seqc.sa_postprocess', 'seqc.sa_preprocess', 'seqc.sa_process'],
      install_requires=[
          'numpy>=1.10.0',
          'cython>0.14',  # tables requirement
          'numexpr>=2.4',  # tables requirement
          'pandas>=0.16.0',
          'matplotlib>=1.5.0',
          'seaborn',
          'scipy>=0.14.0',
          'boto3',
          'pyftpdlib',
          'intervaltree',
          'sklearn>=0.17',
          # 'tsne',
          'tables',
          'nose2',
          'memory_profiler',
          'more_itertools',
          'paramiko', 'requests'],
      scripts=['src/scripts/SEQC',
               'src/scripts/PROCESS_BARCODES',
               'src/scripts/TEST_BARCODES',
               'src/scripts/process_multi_file_scseq_experiment.py',
               'src/scripts/process_single_file_scseq_experiment.py'],
      )

# print any warnings
if h5fail:
    warn("""
SEQC: libhdf5 shared library "libhdf5.so" not found in /usr/local/lib/,
/usr/lib/, /usr/hdf5/lib/, or /usr/local/lib/hdf5/.

"tables" will not find h5lib and installation will likely fail unless the HDF5_DIR
environment variable has been set to the location that HDF5 was installed into. If HDF5
is not installed, please install it prior to installing SEQC.
""")

# look for star
if not shutil.which('STAR'):
    warn('SEQC: STAR is not installed. SEQC will not be able to align files.')


















