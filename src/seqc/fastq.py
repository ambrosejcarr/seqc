import gzip
import bz2
import numpy as np
import re
from seqc.encodings import DNA3Bit
from seqc.barcodes import CellBarcodes
import io
import pickle


__author__ = 'Ambrose J. Carr'


_revcomp = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}


def revcomp(s):
    return ''.join(_revcomp[n] for n in s[::-1])


def truncate_sequence_length(reverse_fastq, n, fname):
    """
    truncate the read length of a fastq file, often to equalize comparisons between
    different sequencing experiments
    """
    if not fname.endswith('.fastq'):
        fname += '.fastq'
    fin = open_file(reverse_fastq)
    try:
        with open(fname, 'w') as fout:
            for record in iter_records(fin):
                seq = seq.strip()[:n] + '\n'
                qual = qual.strip()[:n] + '\n'
                new_record = ''.join((record[0], seq, record[2], qual))
                fout.write(new_record)
    finally:
        fin.close()


def sequence_length_description(fastq):
    """get the sequence length of a fastq file"""
    with open(fastq, 'r') as fin:
        i = 0
        records = iter_records(fin)
        data = np.empty(2500, dtype=int)
        while i < 2500:
            try:
                seq = next(records)[1]
            except StopIteration:
                data = data[:i]
                break
            data[i] = len(seq) - 1
            i += 1
    return np.mean(data), np.std(data), np.unique(data, return_counts=True)


def paired_fastq_records(f, r):
    fit = iter(f)
    rit = iter(r)
    frecord = zip(*[fit] * 4)
    rrecord = zip(*[rit] * 4)
    return zip(iter(frecord), iter(rrecord))


def group_paired(f, r, n):
    fit = iter(f)
    rit = iter(r)
    frecord = zip(*[fit] * 4)
    rrecord = zip(*[rit] * 4)
    paired = zip(iter(frecord), iter(rrecord))
    return zip(*[iter(paired)] * n)


def iter_records(open_fastq_file):
    """return fastq records 1-by-1"""
    args = [iter(open_fastq_file)] * 4
    return zip(*args)


def open_file(filename):
    if filename.endswith('.gz'):
        return gzip.open(filename, 'rt')
    elif filename.endswith('.bz2'):
        return bz2.open(filename, 'rt')
    else:
        return open(filename, 'r')


def remove_homopolymer(r):
    """remove homopolymer sequences

    check for homopolymer sequences, trimming them from the start and end of each read as
    long as the percentage of homopolymers are greater than tolerance.
    """

    seq = r[1].strip()
    qual = r[3].strip()
    original_length = len(seq)

    # get first nucleotide, check for forward homopolymers
    first = seq[0]
    for i, n in enumerate(seq[1:]):
        if n == first:
            continue
        else:
            break
    if i >= 5:
        seq = seq[i + 1:]
        qual = qual[i + 1:]

    # get last nucleotide, check for reverse homopolymers
    last = seq[-1]
    for i, n in enumerate(seq[-2::-1]):
        if n == last:
            continue
        else:
            break
    if i >= 5:
        seq = seq[:-i - 1]
        qual = qual[:-i - 1]

    trimmed_bases = original_length - len(seq)
    return (r[0], seq + '\n', r[2], qual + '\n'), trimmed_bases


def dust_low_complexity_score(record):
    # Sequence
    seq = record[1].strip()

    # Counts of 3-mers in the sequence
    counts = {}
    for i in range(len(seq) - 2):
        kmer = seq[i:i + 3]
        counts[kmer] = counts.get(kmer, 0) + 1

    # Calculate dust score
    score = np.sum([i * (i - 1) / 2 for i in counts.values()]) / (len(seq) - 3)

    # Scale score (Max score possible is no. of 3mers/2)
    score = np.int8(score / ((len(seq) - 2) / 2) * 100)

    return score


def is_primer_dimer(cell_barcode, r):
    """
    determine if the reverse sequence r is a primer_dimer if r contains a cell barcode

    primer_kmer_map is a map containing all kmer pieces of cell barcodes. What size for k?

    For filtering primers:

    generate a suffix array that contains all k-length paths of nucleotides that are
    present in our primers. Paths of length k or longer should yield True. Paths that
    terminate before length k imply that there is no primer with that sequence. Therefore,
    each primer should be read from primer[-k:] and FURTHER.

    Then, a sliding window can be use to check each sequencing read for primer dimers.
    In order to positively call a contaminating primer, the entire sequence should be
    present in the sample

    UPDATE:
    - examining primer dimers in in-drop revealed that there was not a tremendous
    amount in each experiment. Specifically, approximately 0.5% of the data were
    primer dimers. There are primer-specific signals with the cell barcode being present
    in both the forward and reverse reads. I'll try to remove these, but otherwise I'll
    leave things be.

    What follows is the simplest possible check for primer dimers. it will not work for
    in-drop because cell_barcode returns a concatenated form. Can look only for first
    or second barcodes though, if desired.
    """
    return 1 if cell_barcode in r or revcomp(cell_barcode) in r else 0


def annotate_fastq_record(r, cell, rmt, n_poly_t, valid_cell, trimmed_bases, fwd_quality):
    name = (
        '@' + ':'.join(str(v) for v in [cell, rmt, n_poly_t, valid_cell, trimmed_bases,
                                        fwd_quality]) +
        ';' + r[0][1:])
    return ''.join([name] + list(r[1:]))


def average_quality(quality_string):
    """calculate the average quality of a sequencing read from and ASCII quality string"""
    quality_sum = sum(ord(q) for q in quality_string) - len(quality_string) * 33
    n_bases = len(quality_string)
    return quality_sum // n_bases


def process_record(forward, reverse, tbp, cb):
    """process a forward and reverse read pair; eventually, will want to add checks for
    primer dimers"""
    cell, rmt, n_poly_t = tbp.process_forward_sequence(forward[1][:-1])  # exclude \n
    valid_cell = cb.close_match(cell)
    # r, trimmed_bases = remove_homopolymer(reverse)
    # if len(r[1]) < 20:  # don't return short reads
    #     return
    dust_score = dust_low_complexity_score(reverse)
    fwd_quality = average_quality(forward[3][:-1])
    r = annotate_fastq_record(
        reverse, cell, rmt, n_poly_t, valid_cell, dust_score, fwd_quality)
    return r


# todo this n_low_complexity argument is very dubious
def merge_fastq(forward: list, reverse: list, exp_type, temp_dir, cb, n_low_complexity=0):
    """direct threadless fastq processing"""

    # open the files
    fout = '%s/merged_temp.fastq' % temp_dir.rstrip('/')
    merged_file = open(fout, 'w')
    try:
        for ffile, rfile in zip(forward, reverse):

            if not isinstance(ffile, io.TextIOBase):
                ffile = open_file(ffile)
            if not isinstance(rfile, io.TextIOBase):
                rfile = open_file(rfile)
            try:
                # initialize the processor
                tbp = DNA3Bit.default_processors(exp_type)
                if not isinstance(cb, CellBarcodes):
                    with open(cb, 'rb') as f:
                        cb = pickle.load(f)

                for f, r in paired_fastq_records(ffile, rfile):
                    merged_record = process_record(f, r, tbp, cb)
                    # if not merged_record:
                    #     n_low_complexity += 1
                    # else:
                    merged_file.write(merged_record)
            finally:
                ffile.close()
                rfile.close()
    finally:
        merged_file.close()

    return fout  #, n_low_complexity


def auto_detect_processor(experiment_name):
    processors = {
        r'[i|I]n.?[d|D]rop': 'in-drop',
        r'[a|A][v|V][o|O].?[s|S]eq': 'avo-seq',
        r'[d|D]rop.?[s|S]eq': 'drop-seq',
        r'[m|M][a|A][r|R][s|S].?[s|S]eq': 'mars-seq',
        r'[c|C][e|E][l|L].?[s|S]eq': 'cel-seq',
    }
    for p in processors:
        if re.search(p, experiment_name):
            return processors[p]
    raise NameError('pre-processor could not be auto-detected from experiment name, '
                    'please pass the pre-processor name using [-p, --processor]. '
                    'available processors can be found in seqdb.fastq.py')
