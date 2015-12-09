__author__ = 'Ambrose J. Carr'

import gzip
import bz2
import numpy as np
from multiprocessing import Process, Queue
from queue import Empty, Full
from time import sleep
import shutil
import os
import re
from itertools import islice
import seqc
import io
import random
from io import StringIO
from collections import namedtuple


_revcomp = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N',
            'a': 't', 't': 'a', 'c': 'g', 'g': 'c', 'n': 'n'}


def revcomp(s: str) -> str:
    """
    returns the reverse-complement of s

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


def truncate_sequence_length(reverse_fastq: str, n: int, fname: str) -> None:
    """
    truncate the read length of a fastq file, often to equalize comparisons between
    different sequencing experiments

    :param reverse_fastq: open file or str; fastq file to be truncated
    :param n: integer; number of bases to remain in truncated output file
    :param fname: str; name of the output file

    :return: None
    """

    # check n, fname types
    if not isinstance(n, int):
        raise TypeError('n must be an integer, not %s' % type(n))
    if not isinstance(fname, str):
        raise TypeError('fname must be a string, not %s' % type(fname))
    if not fname.endswith('.fastq'):
        fname += '.fastq'

    # check fastq type
    if isinstance(reverse_fastq, io.TextIOBase):
        fin = reverse_fastq
    elif isinstance(reverse_fastq, str):
        fin = open_file(reverse_fastq)
    else:
        raise TypeError('reverse_fastq must be a string, not %s' % type(reverse_fastq))

    # truncate the file
    try:
        with open(fname, 'w') as fout:
            for record in iter_records(fin):
                seq = record[1].strip()[:n] + '\n'
                qual = record[3].strip()[:n] + '\n'
                new_record = ''.join((record[0], seq, record[2], qual))
                fout.write(new_record)
    finally:
        fin.close()


def estimate_sequence_length(fastq: str) -> (int, int, (np.array, np.array)):
    """
    Rapidly estimate the mean, and standard deviation of the sequence length of a
    potentially large fastq file. Also returns the observed lengths and their frequencies.

    :param fastq: str or file object; fastq file whose sequence length is to be estimated

    :return: (int, int, (np.array, np.array); mean, standard deviation, (observed sizes,
      frequencies)
    """

    # check fastq type
    if isinstance(fastq, io.TextIOBase):
        fin = fastq
    elif isinstance(fastq, str):
        fin = open_file(fastq)
    else:
        raise TypeError('reverse_fastq must be a string, not %s' % type(fastq))

    # estimate sequence length from first 2500 records in the file
    try:
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
    finally:
        fin.close()

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
        return open(filename)


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


def merge_fastq(forward: list, reverse: list, exp_type: str, output_dir: str,
                cb: seqc.barcodes.CellBarcodes, n_processes: int) -> str:
    """
    Merges forward and reverse fastq files while extracting barcoding and quality
    data from forward read

    args:
    -----
    forward: list of forward fastq files
    reverse: list of reverse fastq read pairs
    exp_type: type of experiment (e.g. 'in-drop', 'drop-seq'...)
    output_dir: the directory in which to write the merged fastq file
    cb: cell barcode object, used to test whether each read has a valid cell barcode
    n_processes: the number of processes to apply to merging the files.
     Must be at least 3.

    returns:
    --------
    name of merged fastq file, equal to output_dir + merged.fastq

    """

    def read(forward_: list, reverse_: list, in_queue):
        """
        read chunks from fastq files and place them on the processing queue.
        It seems this should take < 1s per 1M read chunk

        # use bytes.index() to find '@' -- the beginning of the next record
        # could also likely find the final newline; this will prevent the method
        # from throwing out corner cases where reading a number of bytes breaks the reads
        """

        # set the number of reads in each chunk
        n = int(1e6)
        i = 0  # index for chunks
        # iterate over all input files
        for ffastq, rfastq in zip(forward_, reverse_):

            # open forward and reverse files for this iteration
            if not isinstance(ffastq, io.TextIOBase):
                ffastq = open_file(ffastq)
            if not isinstance(rfastq, io.TextIOBase):
                rfastq = open_file(rfastq)

            # get slices of reads and put them on the consume queue
            while True:
                seqc.log.info('%d Reading.' % i)
                data = (tuple(islice(ffastq, n * 4)), tuple(islice(rfastq, n * 4)))

                # check that files aren't exhausted
                if not any(d for d in data):  # don't check the index
                    break  # these files are exhausted

                # put chunk on the queue
                while True:
                    try:
                        in_queue.put_nowait((i, data))
                        seqc.log.info('%d Read. Putting on process Queue.' % i)
                        i += 1
                        break
                    except Full:
                        sleep(.1)

            # close fids
            ffastq.close()
            rfastq.close()

    def process(in_queue, out_queue):

        # method to group chunks into records
        def grouped_records(f_, r_):
            # note that this assert will exhaust iterators. pass ITERABLES not ITERATORS
            assert len(f_) == len(r_)
            fit, rit = iter(f_), iter(r_)
            while True:
                frecord, rrecord = tuple(islice(fit, 4)), tuple(islice(rit, 4))
                if not frecord:
                    return
                yield frecord, rrecord

        # get a chunk from queue until all chunks are processed
        while True:
            try:
                index, (forward_, reverse_) = in_queue.get_nowait()
            except Empty:
                try:
                    os.kill(read_pid, 0)  # does nothing if thread is alive
                    sleep(.1)
                    continue
                except ProcessLookupError:  # process is dead.
                    break

            # process records
            merged_filename = '%stemp_%d.fastq' % (output_dir, index)
            with open(merged_filename, 'w') as fout:
                for f, r in grouped_records(forward_, reverse_):
                    fout.write(process_record(f, r, tbp, cb))
            fout.close()

            # put filename out the output queue
            while True:
                try:
                    out_queue.put_nowait(merged_filename)
                    break
                except Full:
                    sleep(.1)

    def any_alive(pids):
        for id_ in pids:
            try:
                os.kill(id_, 0)
                return True  # at least one process alive; return True
            except ProcessLookupError:
                pass
        # no process was alive; return False
        return False

    def merge(out_queue):
        # set a destination file to write everythign into
        seed = open('%smerged.fastq' % output_dir, 'wb')

        # merge all remaining files into it.
        while True:
            try:
                next_file = out_queue.get_nowait()
                shutil.copyfileobj(open(next_file, 'rb'), seed)
                os.remove(next_file)
            except Empty:
                if any_alive(process_pids):
                    sleep(.1)
                    continue
                else:
                    break

        seed.close()

    # do some type checking of inputs
    seqc.util.check_type(forward, list, 'forward')
    seqc.util.check_type(reverse, list, 'reverse')
    seqc.util.check_type(exp_type, str, 'exp_type')
    seqc.util.check_type(output_dir, str, 'output_dir')
    if not isinstance(cb, (seqc.barcodes.CellBarcodes, str)):
        raise TypeError('cb must be one of seqc.barcodes.CellBarcodes or str, not %s' % (
                        type(cb)))

    tbp = seqc.three_bit.ThreeBit.default_processors(exp_type)
    if not isinstance(cb, seqc.barcodes.CellBarcodes):
        cb = seqc.barcodes.CellBarcodes.from_pickle(cb)

    # check that forward and reverse are equal length
    if not len(forward) == len(reverse):
        raise ValueError('Equal number of forward and reverse files must be passed. '
                         '%d != %d' % (len(forward), len(reverse)))

    # check that output dir has the terminal slash
    if not output_dir.endswith('/'):
        output_dir += '/'

    # set the number of processing threads
    n_proc = max(n_processes - 2, 1)

    # read the files
    paired_records = Queue(maxsize=n_proc)  # don't need more waiting items than threads
    read_proc = Process(target=read, args=([forward, reverse, paired_records]))
    read_proc.start()
    read_pid = read_proc.pid

    # process the data
    output_filenames = Queue()
    processors = [Process(target=process, args=([paired_records, output_filenames]))
                  for _ in range(n_proc)]
    for p in processors:
        p.start()

    process_pids = [p.pid for p in processors]

    # write the results
    merge_process = Process(target=merge, args=[output_filenames])
    merge_process.start()

    # wait for each process to finish
    read_proc.join()
    for p in processors:
        p.join()
    merge_process.join()

    return '%smerged.fastq' % output_dir


class GenerateFastq:

    # define some general constants
    _alphabet = ['A', 'C', 'G', 'T']

    def __init__(self):
        pass

    @classmethod
    def _forward_in_drop(cls, n, barcodes_):
        barcodes_ = seqc.barcodes.CellBarcodes.from_pickle(barcodes_)
        read_length = 50
        names = range(n)
        name2 = '+'
        quality = 'I' * read_length
        records = []
        umi_len = 6
        codes = list(barcodes_.perfect_codes)
        for name in names:
            # for now, translate barcode back into string code
            cb = random.choice(codes)
            c1, c2 = seqc.three_bit.ThreeBitInDrop.split_cell(cb)
            c1, c2 = [seqc.three_bit.ThreeBit.bin2str(c) for c in [c1, c2]]
            w1 = 'GAGTGATTGCTTGTGACGCCTT'
            cb = ''.join([c1, w1, c2])
            umi = ''.join(np.random.choice(cls._alphabet, umi_len))
            poly_a = (read_length - len(cb) - len(umi)) * 'T'
            records.append('\n'.join(['@%d' % name, cb + umi + poly_a, name2, quality]))
        forward_fastq = StringIO('\n'.join(records) + '\n')
        return forward_fastq

    @classmethod
    def _forward_drop_seq(cls, n, *args):  # args is for unused barcodes parameters
        names = range(n)
        name2 = '+'
        quality = 'I' * 20
        records = []
        for name in names:
            cb = ''.join(np.random.choice(cls._alphabet, 12))
            umi = ''.join(np.random.choice(cls._alphabet, 8))
            records.append('\n'.join(['@%d' % name, cb + umi, name2, quality]))
        forward_fastq = StringIO('\n'.join(records) + '\n')
        return forward_fastq

    @staticmethod
    def _reverse(n: int, read_length: int, fasta: str, gtf: str, tag_type='gene_id'):

        # read gtf
        reader = seqc.gtf.Reader(gtf)
        intervals = []
        for r in reader.iter_exons():
            end = int(r.end) - read_length
            start = int(r.start)
            if end > start:
                intervals.append((r.attribute[tag_type], start, end, r.strand))

        # pick intervals
        exon_selections = np.random.randint(0, len(intervals), (n,))

        # fasta information:
        with open(fasta) as f:
            fasta = f.readlines()[1:]
            fasta = ''.join(fasta)

        # generate sequences
        sequences = []
        tags = []
        for i in exon_selections:
            tag, start, end, strand = intervals[i]
            # get position within interval
            start = random.randint(start, end)
            end = start + read_length
            seq = fasta[start:end]
            if strand == '-':
                seq = revcomp(seq)
            sequences.append(seq)
            tags.append(tag)

        prefixes = range(n)
        name2 = '+'
        quality = 'I' * read_length
        records = []
        for name, tag, seq in zip(prefixes, tags, sequences):
            records.append('\n'.join(['@%d:%s' % (name, tag), seq, name2, quality]))
        reverse_fastq = StringIO('\n'.join(records) + '\n')
        return reverse_fastq

    @staticmethod
    def _reverse_three_prime(n: int, read_length: int, fasta: str, gtf: str,
                             tag_type='gene_id', fragment_length=1000):

        # read gtf
        reader = seqc.gtf.Reader(gtf)
        intervals = []
        for r in reader.iter_genes_final_nbases(fragment_length):
            for iv in r.intervals:
                start, end = int(iv[0]), int(iv[1])
            if (end - read_length) > start:
                intervals.append((r.attribute[tag_type], start, end, r.strand))

        # pick intervals
        exon_selections = np.random.randint(0, len(intervals), (n,))

        # fasta information:
        with open(fasta) as f:
            fasta = f.readlines()[1:]
            fasta = ''.join(fasta)

        # generate sequences
        sequences = []
        tags = []
        for i in exon_selections:
            tag, start, end, strand = intervals[i]
            # get position within interval
            start = random.randint(start, end)
            end = start + read_length
            seq = fasta[start:end]
            if strand == '-':
                seq = revcomp(seq)
            sequences.append(seq)
            tags.append(tag)

        prefixes = range(n)
        name2 = '+'
        quality = 'I' * read_length
        records = []
        for name, tag, seq in zip(prefixes, tags, sequences):
            records.append('\n'.join(['@%d:%s' % (name, tag), seq, name2, quality]))
        reverse_fastq = StringIO('\n'.join(records) + '\n')
        return reverse_fastq

    @classmethod
    def in_drop(cls, n, prefix, fasta, gtf, barcodes, tag_type='gene_id', replicates=3,
                *args, **kwargs):

        if not replicates >= 0:
            raise ValueError('Cannot generate negative replicates')

        # create directory if it doesn't exist
        directory = '/'.join(prefix.split('/')[:-1])
        if not os.path.isdir(directory):
            os.makedirs(directory)

        fwd_len = 50
        rev_len = 100
        forward = cls._forward_in_drop(n, barcodes)
        forward = forward.read()  # consume the StringIO object
        reverse = cls._reverse_three_prime(n, rev_len, fasta, gtf, tag_type=tag_type)
        reverse = reverse.read()  # consume the StringIO object
        with open(prefix + '_r1.fastq', 'w') as f:
            f.write(''.join([forward] * (replicates + 1)))
        with open(prefix + '_r2.fastq', 'w') as r:
            r.write(''.join([reverse] * (replicates + 1)))
        return prefix + '_r1.fastq', prefix + '_r2.fastq'

    @classmethod
    def drop_seq(cls, n, prefix, fasta, gtf, tag_type='gene_id', replicates=3, *args,
                 **kwargs):

        if not replicates >= 0:
            raise ValueError('Cannot generate negative replicates')

        # create directory if it doesn't exist
        directory = '/'.join(prefix.split('/')[:-1])
        if not os.path.isdir(directory):
            os.makedirs(directory)

        rev_len = 100
        forward = cls._forward_drop_seq(n)
        forward = forward.read()  # consume the StringIO object
        reverse = cls._reverse_three_prime(n, rev_len, fasta, gtf, tag_type=tag_type)
        reverse = reverse.read()  # consume the StringIO object
        with open(prefix + '_r1.fastq', 'w') as f:
            f.write(''.join([forward] * (replicates + 1)))
        with open(prefix + '_r2.fastq', 'w') as r:
            r.write(''.join([reverse] * (replicates + 1)))
        return prefix + '_r1.fastq', prefix + '_r2.fastq'

    @classmethod
    def simple_fastq(cls, n, length):
        sequences = np.random.choice(list('ACGT'), n * length)
        sequences = np.reshape(sequences, (n, length))
        qualities = np.random.choice(np.arange(30, 40), n * length)
        qualities = np.reshape(qualities, (n, length))

        fastq_data = ''
        for i in range(n):
            name = '@simple_fastq:i\n'
            seq = ''.join(sequences[i, :]) + '\n'
            name2 = '+\n'
            qual = ''.join(chr(s) for s in qualities[i, :]) + '\n'
            fastq_data += ''.join([name, seq, name2, qual])

        fastq_data = StringIO(fastq_data)
        fastq_data.seek(0)

        return fastq_data


FastqRecord = namedtuple('FastqRecord', ['name', 'seq', 'name2', 'qual'])


class Reader:
    """simple fastq reader, optimized for utility rather than speed"""

    def __init__(self, fastq_file):

        seqc.util.check_type(fastq_file, str, 'fastq_file')
        seqc.util.check_file(fastq_file, 'fastq_file')

        self._fastq = fastq_file
        try:
            fastq_iterator = iter(self)
            first = next(fastq_iterator)
        except:
            raise ValueError('%s is an invalid fastq_file. Please check file formatting.' %
                             fastq_file)
        if not first.name.startswith('@'):
            raise ValueError('Invalid formatting: all name lines must start with "@"')
        if not first.name2.startswith('+'):
            raise ValueError('Invalid formatting: all name2 lines must start with "+"')

    @property
    def fastq(self):
        return self._fastq

    def _open(self) -> io.TextIOBase:
        """
        seamlessly open self.fastq, whether gzipped or uncompressed

        returns:
        --------
        fobj: open file object
        """
        if self.fastq.endswith('.gz'):
            fobj = gzip.open(self.fastq, 'rt')
        else:
            fobj = open(self.fastq)
        return fobj

    def __len__(self):
        return sum(1 for _ in self)

    def __iter__(self):
        """return an iterator over all records in fastq file"""
        fobj = self._open()
        try:
            for record in iter_records(fobj):
                yield FastqRecord(*record)
        finally:
            fobj.close()
