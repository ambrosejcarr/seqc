from seqc.sequence.encodings import DNA3Bit
from sys import maxsize

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
        return ''.join(_revcomp_bytes[n] for n in b[::-1])
    except KeyError:
        raise ValueError('%s contains invalid nucleotide character. Supported characters '
                         'are A, C, G, T and N.' % b)


def find_correct_barcode(code, barcodes_list):
    """
    For a given barcode find the closest correct barcode to it from the list (limited to
    one ED), a string representing the error and the edit distance
    NOTE: for now this function looks for a barcode with ED==1 and does not bother
    looking for the minimum

    :param barcodes_list:
    :param code:
    :returns:
    """
    if code in barcodes_list:
        return code, 0

    min_ed = maxsize
    cor_code = 0
    for bc in barcodes_list:
        hamm_d = hamming_dist_bin(code, bc)
        if hamm_d == 1:
            min_ed = 1
            cor_code = bc
            break
        if hamm_d < min_ed:
            min_ed = hamm_d
            cor_code = bc

    return cor_code, min_ed
        
        
def hamming_dist_bin(c1, c2):
    """Return the hamming distance between two numbers representing a sequence (3 bits
    per base)

    :param c1:
    :param c2:
    :return:
    """
    if DNA3Bit.seq_len(c1) != DNA3Bit.seq_len(c2):
        return maxsize
    d = 0
    while c1 > 0:
        if c1 & 0b111 != c2 & 0b111:
            d += 1
        c1 >>= 3
        c2 >>= 3
    return d
    
def list_errors(s1, s2):
    """
    Return the list of nucleotide transformations that turn s1 to s2.
    An error is a six bit int representing a two chr string of type "AG","CT", etc.

    :param s2:
    :param s1:

    :returns:
    """

    # return the actual error
    err_list = []
    while s1 > 0:
        if s1 & 0b111 != s2 & 0b111:
            err_list.append(DNA3Bit.ints2int([s1 & 0b111, s2 & 0b111]))
        s1 >>= 3
        s2 >>= 3
    return err_list

