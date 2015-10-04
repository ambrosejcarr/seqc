__author__ = 'Ambrose J. Carr'

import gzip
import bz2
from threading import Thread
from queue import Queue, Empty, Full
from time import sleep, time
import re
from seqc.three_bit import ThreeBit, CellBarcodes
import io
import pickle


_revcomp = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}


def revcomp(s):
    return ''.join(_revcomp[n] for n in s[::-1])


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
    if i < 5:
        pass
    else:
        seq = seq[i + 1:]
        qual = qual[i + 1:]

    # get last nucleotide, check for reverse homopolymers
    last = seq[-1]
    for i, n in enumerate(seq[-2::-1]):
        if n == last:
            continue
        else:
            break
    if i < 5:
        pass
    else:
        seq = seq[:-i - 1]
        qual = qual[:-i - 1]
    trimmed_bases = original_length - len(seq)
    return (r[0], seq + '\n', r[2], qual + '\n'), trimmed_bases


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
    name = ('@' + ':'.join(str(v) for v in [cell, rmt, n_poly_t, valid_cell, trimmed_bases,
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
    r, trimmed_bases = remove_homopolymer(reverse)
    if len(r[1]) < 20:  # don't return short reads
        return ''
    fwd_quality = average_quality(forward[3][:-1])
    r = annotate_fastq_record(
        r, cell, rmt, n_poly_t, valid_cell, trimmed_bases, fwd_quality)
    return r


def merge_fastq(forward, reverse, exp_type, temp_dir, cb, n_low_complexity=0):
    """direct threadless fastq processing"""

    # open the files
    fout = '%s/merged_temp.fastq' % temp_dir
    merged_file = open(fout, 'w')
    try:
        for ffile, rfile in zip(forward, reverse):

            if not isinstance(ffile, io.TextIOBase):
                ffile = open_file(ffile)
            if not isinstance(rfile, io.TextIOBase):
                rfile = open_file(rfile)
            try:
                # initialize the processor
                tbp = ThreeBit.default_processors(exp_type)
                if not isinstance(cb, CellBarcodes):
                    with open(cb, 'rb') as f:
                        cb = pickle.load(f)

                for f, r in paired_fastq_records(ffile, rfile):
                    merged_record = process_record(f, r, tbp, cb)
                    if merged_record == '':
                        n_low_complexity += 1
                        continue
                    merged_file.write(merged_record)
            finally:
                ffile.close()
                rfile.close()
    finally:
        merged_file.close()

    return fout, n_low_complexity


# todo not working right now
def merge_fastq_threaded(forward, reverse, n_threads, exp_type, temp_dir,
                         cell_barcode_pickle):
    """multi-threaded merge of forward and reverse records into a single alignable file"""

    def read(forward_, reverse_, in_queue):
        """read a forward and reverse record and place on the queue"""
        for forward_file, reverse_file in zip(forward_, reverse_):
            # todo for production, boost this number up to 1e6 reads at a time
            paired = group_paired(forward_file, reverse_file, 10000)
            while True:  # loop over all input
                try:
                    data = next(paired)
                except StopIteration:
                    break
                while True:
                    try:
                        in_queue.put(data)
                        break
                    except Full:
                        sleep(0.1)
                        continue

    def process(in_queue, out_queue):
        while True:
            try:
                # forward_data, reverse_data = in_queue.get_nowait()
                data = in_queue.get_nowait()
            except Empty:
                if not read_thread.is_alive():
                    break
                else:
                    sleep(0.1)  # put on the brakes a bit
                    continue

            # PROCESS READS HERE
            merged = '\n'.join(process_record(f, r, tbp, cb) for (f, r) in data)

            while True:
                try:
                    out_queue.put(merged)
                    break
                except Full:
                    sleep(0.1)
                    continue

    def write(in_queue):
        merged = open('%s/merged_temp.fastq' % temp_dir, 'w')
        try:
            while True:
                try:
                    record = in_queue.get_nowait()
                    merged.write(record)
                except Empty:
                    if not any(t.is_alive() for t in process_threads):
                        break
                    else:
                        sleep(0.1)
                        continue
        finally:
            merged.close()

    # initialize the processor
    tbp = ThreeBit.default_processors(exp_type)
    # cb = CellBarcodes(processor, *cell_barcode_files)
    with open(cell_barcode_pickle, 'rb') as f:
        cb = pickle.load(f)

    # open the files
    if not isinstance(forward, io.TextIOBase):
        forward = open_file(forward)
    if not isinstance(reverse, io.TextIOBase):
        reverse = open_file(reverse)

    try:  # wrap to make sure files get closed
        # read the files
        paired_records = Queue(maxsize=2048)
        read_thread = Thread(target=read, args=([forward, reverse, paired_records]))
        read_thread.start()

        # process the data
        merged_records = Queue(maxsize=2048)
        # max --> make sure at least one thread starts
        process_threads = [Thread(target=process, args=([paired_records, merged_records]))
                           for _ in range(max(n_threads - 3, 1))]
        assert(len(process_threads) > 0)
        for t in process_threads:
            t.start()

        # write the results
        write_thread = Thread(target=write, args=[merged_records])
        write_thread.start()

        # wait for each process to finish
        read_thread.join()
        print('reading done: %f' % time())
        for t in process_threads:
            t.join()
        print('processing done: %f' % time())
        write_thread.join()
        print('writing done: %f' % time())
    finally:
        forward.close()
        reverse.close()

    return '%s/merged_temp.fastq' % temp_dir


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
