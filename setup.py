__author__ = 'Ambrose J. Carr'

from setuptools import setup

setup(name='seqc',
      version='0.1',
      description='Single Cell Sequencing Processing and QC Suite',
      author='Ambrose J. Carr',
      author_email='mail@ambrosejcarr.com',
      package_dir={'': 'src'},
      # note: requires numpy > 1.10.0
      packages=['seqc'], requires=['numpy', 'pandas', 'matplotlib', 'seaborn'],
      scripts=['src/scripts/SEQC', 'src/scripts/PROCESS_BARCODES'],
      )

