__author__ = 'ambrose'

__version__ = '0.1.5'


# supress warnings after importing seqc
import warnings
warnings.filterwarnings("ignore")

import pyximport; pyximport.install()
from . import gtf_new

from . import align
from . import arrays
from . import analyze
from . import barcodes
from . import convert_features
from . import fasta
from . import fastq
from . import gtf
from . import h5
from . import io
from . import log
from . import parallel
from . import plot
from . import sam
from . import sa_process
from . import sa_preprocess
from . import sa_postprocess
from . import encodings
from . import core
from . import util
from . import sra
from . import cluster_utils
from . import ssh_utils

# import some key classes into the general namespace
from .analyze import Experiment, DenseCounts, SparseCounts
from .arrays import ReadArray, UniqueReadArray