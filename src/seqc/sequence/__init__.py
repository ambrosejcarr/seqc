from . import fastq
from . import merge_functions
from . import gtf
from . import encodings

_revcomp = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N',
            'a': 't', 't': 'a', 'c': 'g', 'g': 'c', 'n': 'n'}

_revcomp_bytes = {b'A': b'T', b'T': b'A', b'C': b'G', b'G': b'C', b'N': b'N',
                  b'a': b't', b't': b'a', b'c': b'g', b'g': b'c', b'n': b'n'}


def revcomp(s: str) -> str:
    """
    returns the reverse-complement of nucleotide sequence s

    :param s: str; nucleotide sequence to be reverse complemented
    :return: str; reverse complemented input string s
    """
    if not isinstance(s, str):
        raise TypeError('Nucleotide sequence "s" must be string, not %s' % type(s))
    try:
        return ''.join(_revcomp[n] for n in s[::-1])
    except KeyError:
        raise ValueError('%s contains invalid nucleotide character. Supported characters '
                         'are A, C, G, T and N.' % s)


def revcomp_bytes(b: bytes) -> bytes:
    """
    returns the reverse-complement of s

    :param b: bytes; nucleotide sequence to be reverse complemented
    :return: bytes; reverse complemented input string s
    """
    if not isinstance(b, bytes):
        raise TypeError('Nucleotide sequence "s" must be bytes, not %s' % type(b))
    try:
        return ''.join(_revcomp[n] for n in b[::-1])
    except KeyError:
        raise ValueError('%s contains invalid nucleotide character. Supported characters '
                         'are A, C, G, T and N.' % b)
