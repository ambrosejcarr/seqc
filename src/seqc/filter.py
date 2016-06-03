from seqc.sequence.fastq import Reader
from math import floor

_primer_lengths = dict(
    in_drop=47,
    in_drop_v2=49,
    drop_seq=20,
    mars1_seq=None,  # mars-seq has not been provided to us as pre-demultiplexed data
    mars2_seq=None
)


def estimate_min_poly_t(fastq_files: list, platform: str) -> int:
    """
    estimate the minimum size of poly-t tail that should be present on a properly captured molecule's
     forward read. If multiple fastq files are passed, the minimum value across all files will be
     returned

    :param fastq_files: list of fastq filenames
    :param platform: the platform used to generate this library
    :return: int minimum number of poly-t expected from a valid capture primer
    """
    min_vals = []
    try:
        primer_length = _primer_lengths[platform]
    except KeyError:
        raise ValueError('provided value {} is not a valid argument for platform.'.format(platform))
    if primer_length is None:
        raise RuntimeError('provided platform does not have a defined primer length, and thus the '
                           'min_poly_t parameter cannot be estimated. Please provide --min-poly-t '
                           'explicitly in process_experiment.py.')
    for f in fastq_files:
        mean = Reader(f).estimate_sequence_length()[0]
        available_nucleotides = max(0, mean - primer_length)
        min_vals.append(floor(min(available_nucleotides * .8, 20)))
    return min(min_vals)

