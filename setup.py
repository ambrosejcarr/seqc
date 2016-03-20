from setuptools import setup
from warnings import warn
import os
import shutil

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
      version='0.1.6',
      description='Single Cell Sequencing Processing and QC Suite',
      author='Ambrose J. Carr',
      author_email='mail@ambrosejcarr.com',
      package_dir={'': 'src'},
      packages=['seqc', 'seqc.sequence', 'seqc.alignment'],
      install_requires=[
          'numpy>=1.10.0',
          'cython>0.14',  # tables requirement
          'numexpr>=2.4',  # tables requirement
          'pandas>=0.16.0',
          'paramiko',
          'regex',
          'requests',
          'scipy>=0.14.0',
          'boto3',
          'intervaltree',
          'tables'],
      scripts=['src/seqc/process_experiment.py']
      )

# print any warnings
if h5fail:
    warn('SEQC: libhdf5 shared library "libhdf5.so" not found in /usr/local/lib/, '
         '/usr/lib/, /usr/hdf5/lib/, or /usr/local/lib/hdf5/. '
         'tables will not find h5lib and installation will likely fail unless the '
         'HDF5_DIR environment variable has been set to the location that HDF5 was '
         'installed into. If HDF5 is not installed, please install it prior to '
         'installing SEQC. ')

# look for star
if not shutil.which('STAR'):
    warn('SEQC: STAR is not installed. SEQC will not be able to align files.')
