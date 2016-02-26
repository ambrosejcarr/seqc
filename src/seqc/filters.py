import numpy as np
import seqc


def quality(read_array: seqc.arrays.ReadArray, minimum_average_forward_quality=25,
            minimum_average_reverse_quality=30):
    """
    returns a boolean vector indicating which reads have adequate forward and reverse
    qualities to be treated as valid sequencing outputs
    """
    passes_fwd = read_array['fwd_quality'] >= minimum_average_forward_quality
    passes_rev = read_array['rev_quality'] >= minimum_average_reverse_quality
    return passes_fwd & passes_rev


def complexity(read_array: seqc.arrays.ReadArray, maximum_complexity_score=5):
    """
    returns a boolean vector indicating which reads have adequate complexity
    """
    return read_array['dust_score'] <= maximum_complexity_score


def unique(read_array: seqc.arrays.ReadArray):  # todo implement
    """
    returns a boolean vector indicating which reads have unique features

    args:
    -----
    read_array: seqc.arrays.ReadArray object
    """
    pass


def valid_capture_primer(read_array: seqc.arrays.ReadArray, min_poly_t=3):
    """
    returns a boolean vector indicating which reads have valid capture primers

    args:
    -----
    read_array: seqc.arrays.ReadArray object
    min_poly_t: minimum number of T nucleotides for the capture primer to be valid
    """
    return read_array['n_poly_t'] >= min_poly_t


def genomic_contains_barcode():  # todo implement
    """
    """
    pass


def valid_cell(read_array: seqc.arrays.ReadArray, barcodes: seqc.barcodes.CellBarcodes):
    """
    return a boolean vector indicating whether or not the cell barcode is within 1 ED of
    an expected barcode.
    """
    codes = np.array(barcodes.error_codes)
    return np.in1d(read_array['cell'], codes)


def transcriptomic_false_positive(multi_alignment: tuple):
    """
    returns True if read is aligned to the transcriptome but has a better alignment to the
    genome

    args:
    -----
    multi_alignment
    """
    pass


def distance():  # todo implement
    """
    filter based on distance from nearest TTS
    """
    pass

