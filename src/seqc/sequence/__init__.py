from . import fastq
from . import fasta
from . import merge_functions


# _revcomp = {b'A': b'T', b'T': b'A', b'C': b'G', b'G': b'C', b'N': b'N',
#             b'a': b't', b't': b'a', b'c': b'g', b'g': b'c', b'n': b'n'}

_revcomp = {84: b'A', 71: b'C', 67: b'G', 65: b'T', 78: b'N'}

def revcomp(b: bytes) -> bytes:
    """
    returns the reverse-complement of s

    args:
    -----
    b: nucleotide sequence to be reverse complemented

    returns:
    --------
    reverse complement of b
    """
    try:
        return b''.join(_revcomp[n] for n in b[::-1])
    except KeyError:
        raise ValueError('{0} contains invalid nucleotide character. Supported '
                         'characters are A, C, G, T and N.'.format(b.decode()))
