from seqc.sequence.encodings import DNA3Bit
from sys import maxsize


def find_correct_barcode(code, barcodes_list, exact_match=False):
    """
    For a given barcode find the closest correct barcode to it from the list (limited to
    one ED), a string representing the error and the edit distance
    NOTE: for now this function looks for a barcode with ED==1 and does not bother
    looking for the minimum

    :param exact_match:
    :param barcodes_list:
    :param code:
    :returns:
    """

    # Return the barcode if it exists
    if code in barcodes_list:
        return code, 0

    # If perfect match is required, return an error since the barcode does not appear
    # in the correct barcode list
    if exact_match:
        return 0, maxsize

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
            err_list.append((s1 & 0b111, s2 & 0b111))
        s1 >>= 3
        s2 >>= 3
    return err_list
