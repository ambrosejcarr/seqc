__author__ = 'ambrose'

from threading import Thread
from queue import Queue, Full, Empty
from time import sleep
import numpy as np
import collections
from itertools import islice, tee


def multi_delete(sorted_deque, *lists):
    while sorted_deque:
        i = sorted_deque.pop()
        for l in lists:
            del l[i]
    return lists


def mask_failing_cells(df_or_array, n_poly_t_required):
    return ((df_or_array['cell'] != 0) &
            (df_or_array['rmt'] != 0) &
            (df_or_array['n_poly_t'] >= n_poly_t_required) &
            (df_or_array['is_aligned'])
            )


# """
# the sole purpose of the two classes below is to defeat the confounding way that numpy
# treats tuples with len == 1 when constructing arrays by obscuring the existence of tuple
# data by 1 level. Pass a tuple as a feature or it will break
# """
#
#
# class ObfuscatedTuple:
#
#     def __init__(self, tuple_):
#         if not isinstance(tuple_, tuple):
#             tuple_ = tuple(tuple_)
#         self._tuple = tuple_
#
#     def __repr__(self):
#         return '<ObfuscatedTuple: %s>' % self._tuple.__repr__()
#
#     def __lt__(self, other):
#         return self._tuple.__lt__(other)
#
#     def __gt__(self, other):
#         return self._tuple.__gt__(other)
#
#     def __eq__(self, other):
#         return self._tuple.__eq__(other)
#
#     def __ne__(self, other):
#         return self._tuple.__ne__(other)
#
#     def __le__(self, other):
#         return self._tuple.__le__(other)
#
#     def __ge__(self, other):
#         return self._tuple.__ge__(other)
#
#     def __len__(self):
#         return len(self._tuple)
#
#     def __contains__(self, item):
#         return self._tuple.__contains__(item)
#
#     def __iter__(self):
#         return iter(self._tuple)
#
#     def __hash__(self):
#         return hash(self._tuple)
#
#     def to_tuple(self):
#         return self._tuple


class Peekable(collections.Iterator):

    def __init__(self, it):
        self.it, self.nextit = tee(iter(it))
        self._advance()

    def _advance(self):
        self.peek = next(self.nextit, None)

    def __next__(self):
        self._advance()
        return next(self.it)


def average_quality(quality_string):
    """calculate the average quality of a sequencing read from and ASCII quality string"""
    n_bases = len(quality_string)
    return (sum(ord(q) for q in quality_string) - n_bases * 33) // n_bases


def translate_feature(reference_name, strand, true_position, feature_table,
                      feature_positions):
    rounded_position = true_position // 100 * 100
    try:
        potential_scids = feature_table[(strand, reference_name, rounded_position)]
    except KeyError:
        return 0

    # only one feature in stranded libraries
    for scid in potential_scids:
        if any(s < true_position < e for (s, e) in feature_positions[scid]):
            return scid  # todo | test if ever > 1 feature
        else:
            return 0


def group_multialignments(alignments):
    iterator = Peekable(alignments)

    # get first non-header alignment
    multialignment = next(iterator)
    while multialignment.startswith('@'):
        multialignment = next(iterator)
    multialignment = [multialignment]

    while True:
        next_alignment = iterator.peek
        if not next_alignment:
            yield multialignment
            break
        if multialignment[0].split('\t')[0] == next_alignment.split('\t')[0]:
            multialignment.append(next(iterator))
        else:
            yield multialignment
            multialignment = [next(iterator)]


def process_multialignment(alignments, feature_positions, feature_table):
    """translate a sam record into a recarray row"""

    # all fields are identical except feature; get from first alignment
    first = alignments[0].strip().split('\t')

    rev_quality = average_quality(first[10])
    alignment_score = int(first[13].split(':')[-1])

    # parse data from name field, previously extracted from forward read
    forward_metadata = (int(f) for f in first[0].strip().split(';')[0].split(':'))
    cell, rmt, n_poly_t, valid_cell, trimmed_bases, fwd_quality = forward_metadata

    # get all features and positions
    if len(alignments) == 1:
        true_position = int(first[3])
        flag = int(first[1])
        strand = '-' if (flag & 16) else '+'
        reference_name = first[2]
        features = [translate_feature(reference_name, strand, true_position,
                                      feature_table, feature_positions)]
        positions = [true_position]
    else:
        features = []
        positions = []
        for alignment in alignments:
            alignment = alignment.strip().split('\t')
            true_position = int(alignment[3])
            flag = int(alignment[1])
            strand = '-' if (flag & 16) else '+'
            reference_name = alignment[2]
            features.append(translate_feature(reference_name, strand, true_position,
                                              feature_table, feature_positions))
            positions.append(true_position)

    delete = collections.deque()
    for i in range(len(features)):
        if features[i] == 0:
            delete.append(i)

    features, positions = multi_delete(delete, features, positions)

    positions = ObfuscatedTuple(positions)
    features = ObfuscatedTuple(features)
    is_aligned = True if features else False

    rec = (cell, rmt, n_poly_t, valid_cell, trimmed_bases, rev_quality, fwd_quality,
           features, positions, is_aligned, alignment_score)

    return rec


def create_structured_array(n):
    """pre-allocate a recarray of size n for sam processing

    note that the complete array will often not be utilized, since multimapping fastq
    records will be compressed to a single recarray row"""

    dtype = [
        ('cell', np.int64),
        ('rmt', np.int32),
        ('n_poly_t', np.uint8),
        ('valid_cell', np.bool),
        ('trimmed_bases', np.uint8),
        ('rev_quality', np.uint8),
        ('fwd_quality', np.uint8),
        ('is_aligned', np.bool),
        ('alignment_score', np.uint8)
    ]

    # trade-off: 50-95% of features will be unique. Therefore, it makes more sense to
    # store them as ints. lists have an overhead of ~ 64b; dictionaries have overhead of
    # 288 bytes (?) so it looks like we should have a VL array.
    # each record in the sam array corresponds to an entry in the VL array; we could know
    # the fixed size. To use a VL array would require two pointer lists, each with int8
    # size, and then the data itself.

    return np.zeros((n,), dtype=dtype)


def process_alignments(samfile, n_threads, gtf, fragment_len):

    """create a record array containing alignment, barcode, and filtering information"""

    def read(in_queue):
        s = open(samfile)
        iterator = Peekable(s)
        try:
            while True:
                try:
                    # get data
                    data = list(islice(iterator, int(1e6)))
                    # make sure we haven't bisected a multialignment
                    final_record = data[-1]
                    while True:
                        next_record = iterator.peek
                        if not next_record:
                            raise StopIteration
                        if final_record.split('\t')[0] == next_record.split('\t')[0]:
                            data.append(next(iterator))
                        else:
                            break
                except StopIteration:
                    if data:
                        while True:  # put the final data chunk on the queue
                            try:
                                in_queue.put(data)
                                break
                            except Full:
                                sleep(0.1)
                                continue
                    break

                # put the read data on the queue
                while True:
                    try:
                        in_queue.put(data)
                        break
                    except Full:
                        sleep(0.1)
                        continue
        finally:
            s.close()

    def process(in_queue, out_queue, feature_positions, feature_table):
        while True:
            # get data and process it into records
            try:
                data = in_queue.get_nowait()

                # process the alignment group
                arr = create_structured_array(len(data))
                for i, ma in enumerate(group_multialignments(data)):
                    row = process_multialignment(ma, feature_positions, feature_table)
                    arr[i] = row
                arr = arr[0:i + 1]  # cut any unfilled rows
            except Empty:  # wait a bit for the next set of data
                if not read_thread.is_alive():
                    break
                else:
                    sleep(0.1)
                    continue

            # put the data on the concatenation queue
            while True:
                try:
                    out_queue.put(arr)
                    break
                except Full:
                    sleep(0.1)
                    continue

    # create a feature table todo | load from index
    feature_table_, feature_positions_ = construct_feature_table(gtf, fragment_len)

    # read the files
    alignment_groups = Queue()
    read_thread = Thread(target=read, args=[alignment_groups])
    read_thread.start()

    # process the data
    processed_alignments = Queue()
    process_threads = [Thread(target=process,
                              args=[alignment_groups, processed_alignments,
                                    feature_positions_, feature_table_])
                       for _ in range(n_threads - 3)]

    for t in process_threads:
        t.start()

    # wait for each process to finish
    read_thread.join()
    for t in process_threads:
        t.join()

    # once all processing threads are dead, concatenate all of the arrays

    # get the first array
    while True:
        try:
            arr = processed_alignments.get_nowait()
            break
        except Empty:
            if not any(t.is_alive() for t in process_threads):
                raise ValueError('no processed input received')
            else:
                sleep(0.1)
                continue

    # concatenate all the other arrays
    while True:
        try:
            next_arr = processed_alignments.get_nowait()
            arr = np.hstack([arr, next_arr])
        except Empty:
            if not any(t.is_alive() for t in process_threads):
                break
            else:
                sleep(0.1)
                continue
    return arr


def process_fluidigm_data():
    """best we can do with this data is reads per gene

    process each gene with htseq-count
    """


def assess_outliers(samfile):
    """
    parse samfile, generate reads per cell, & molecules per cell lists and total
    alignments.

    Normalize two ways:
     1: total input reads
     2: aligned reads
    """
    raise NotImplementedError


def iterate(samfile):
    """iterate over the lines of a .sam or .bam file"""
    if samfile.endswith('.bam'):

        from subprocess import PIPE, Popen, check_output
        import sys

        samtools = check_output(['which', 'samtools'])
        if not samtools:
            print('Samtools binary not found on PATH. Exiting')
            sys.exit(1)
        else:
            samtools = samtools.decode('UTF-8').strip()
        args = [samtools, 'view', samfile]
        p = Popen(args, stdout=PIPE, bufsize=256)

        while True:
            line = p.stdout.readline()

            if not line:
                break

            yield line.decode('UTF-8')

    elif samfile.endswith('.sam'):
        with open(samfile, 'r') as f:
            for line in f:
                if line.startswith('@'):
                    continue
                else:
                    yield line
    else:
        raise ValueError('file extension not recognized')


def iterate_chunk(samfile, n):
    """read n gigabytes of sequencing information at a time."""

    if samfile.endswith('.bam'):

        from subprocess import PIPE, Popen, check_output
        import sys

        samtools = check_output(['which', 'samtools'])
        if not samtools:
            print('Samtools binary not found on PATH. Exiting')
            sys.exit(1)
        else:
            samtools = samtools.decode('UTF-8').strip()
        args = [samtools, 'view', samfile]
        p = Popen(args, stdout=PIPE, bufsize=-1)

        while True:
            lines = p.stdout.readlines(n * 1048576)

            if not lines:
                break
            lines = [line.decode('UTF-8') for line in lines]

            yield lines

    elif samfile.endswith('.sam'):
            with open(samfile, 'r') as f:
                first_lines = f.readlines(n * 1048576)
                first_lines = [line for line in first_lines if not
                               line.startswith('@')]

                while True:
                    lines = f.readlines(n * 1048576)

                    if first_lines:
                        lines = first_lines + lines
                        first_lines = None

                    if not lines:
                        break

                    yield lines

    else:
        raise ValueError('file extension not recognized')
