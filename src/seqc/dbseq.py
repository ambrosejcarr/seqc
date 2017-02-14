import os
from functools import partial
from contextlib import closing
from collections import deque
from collections import namedtuple
from itertools import product
import multiprocessing as mp
import numpy as np
import tables as tb
import pysam
from seqc.sequence.encodings import DNA3Bit
from seqc.sequence import gtf
import time


# global values
GENOME_CHUNKSIZE = int(1e7)
ARRAY_CHUNKSIZE = int(1e6)
LOCK = mp.RLock()
DT = np.dtype([
        ('status', np.uint8),
        ('cell', np.int64),
        ('rmt', np.int32),
        ('n_poly_t', np.uint8),
        ('gene', np.uint32),
        ('position', np.uint64),
        ('n_aligned', np.uint8)])


def chunk_genome(chrom_name_length_file):
    """parse chromsome information embedded in a STAR index to iterate over genome chunks
    :param str chrom_name_length_file:
    :return Iterator(str, int, int):
    """
    with open(chrom_name_length_file, 'r') as f:
        chr_sizes = dict(line.strip().split() for line in f.readlines())

    for chromosome, size in chr_sizes.items():
        size = int(size)
        start = 0
        end = GENOME_CHUNKSIZE
        while True:
            if end > size:
                yield (chromosome, start, size)
                break
            else:
                yield (chromosome, start, end)
                start = end
                end += GENOME_CHUNKSIZE


def bam_iterator(bam_name, chromosome, start, end):
    """iterate over a chunk of a bamfile
    :param str bam_name:
    :param str chromosome:
    :param int start:
    :param int end:
    :return:
    """
    with closing(pysam.AlignmentFile(bam_name, 'rb')) as f:
        for aligned_segment in f.fetch(chromosome, start, end, multiple_iterators=True, until_eof=True):
            yield aligned_segment


def create_storage(archive_name):
    """
    :param str archive_name:
    :param int expectedrows: expected number of rows in the table
    :return:
    """
    

    blosc5 = tb.Filters(complib='blosc', complevel=5)
    with closing(tb.open_file(archive_name, mode='w', filters=blosc5)) as h5:
        records = h5.create_table(h5.root, 'records', DT, "Records")


def window(seq, n=2):
    it = iter(seq)
    win = deque((next(it, None) for _ in range(n)), maxlen=n)
    yield tuple(win)
    append = win.append
    for e in it:
        append(e)
        yield tuple(win)


def worker(chunk, archive_name, bam_name, translator):
    """
    :param (str, int, int) chunk:
    :param str archive_name:
    :param str bam_name:
    :return:
    """

    def store(h5, table):
        with LOCK:
            table.flush()
            h5.close()

    print('processing', *chunk)
    h5 = tb.open_file(archive_name, mode='a')
    table = h5.root.records
    i = 0
    for aligned_segment in bam_iterator(bam_name, *chunk):
        if i == ARRAY_CHUNKSIZE:
            store(h5, table)
            i = 0
            h5 = tb.open_file(archive_name, mode='a')
            table = h5.root.records
        else:
            rec = table.row
            #rec['status'] = 0
            rec['cell'], rec['rmt'], rec['n_poly_t'], rec['gene'], rec['position'], rec['n_aligned'] = sam_record(aligned_segment, translator)
            rec.append()
        i += 1

    store(h5, table)


def worker_multi(chunk, archive_name, bam_name, translator):
    """
    :param (str, int, int) chunk:
    :param str archive_file:
    :param str bam_name:
    :return:
    """
    storage = np.zeros((ARRAY_CHUNKSIZE,), dtype=DT)
    print('processing', *chunk)
    i = 0
    for aligned_segment in bam_iterator(bam_name, *chunk):
        if i == ARRAY_CHUNKSIZE:
            with LOCK:
                with closing(tb.open_file(archive_name, mode='a')) as h5:
                    table = h5.root.records
                    table.append(storage)
                    table.flush()
            i = 0
            storage = np.zeros((ARRAY_CHUNKSIZE, len(DT)), dtype=DT)
        else:
            storage[i]= from_sam_record(aligned_segment, translator)
        i += 1

    storage = storage[:i % ARRAY_CHUNKSIZE]
    with LOCK:
        with closing(tb.open_file(archive_name, mode='a')) as h5:
                table = h5.root.records
                table.append(storage)
                table.flush()


def from_sam_record(seg, translator):
    NameField = namedtuple(
        'NameField', ['pool', 'cell', 'rmt', 'poly_t', 'name'])
    qname = seg.query_name

    # Parse the name field
    fields, name = qname.split(';')
    processed_fields = fields.split(':')
    processed_fields.append(name)
    processed_fields = NameField(*processed_fields)

    pos = seg.reference_start
    minus_strand = seg.is_reverse
    if minus_strand:
        strand = '-'
    else:
        strand = '+'

    gene = translator.translate(seg.reference_id, strand, pos)
    if not gene:
        gene = 0

    if seg.is_unmapped:
        n_aligned = 0
    else:
        n_aligned = seg.get_tag('NH')

    cell = DNA3Bit.encode(processed_fields.cell)
    rmt = DNA3Bit.encode(processed_fields.rmt)
    n_poly_t = processed_fields.poly_t.count('T') + processed_fields.poly_t.count('N')

    return np.array([(0, cell, rmt, n_poly_t, gene, pos, n_aligned)], dtype=DT)



def main_single(archive_name, alignment_file, chrom_name_length_file, translator):
    create_storage(archive_name)
    print('storage created')
    iterator = chunk_genome(chrom_name_length_file)
    func = partial(worker, archive_name=archive_name, bam_name=alignment_file, translator=translator)
    for chunk in iterator:
        func(chunk)


def main_multiple(archive_name, alignment_file, chrom_name_length_file, translator):
    create_storage(archive_name)
    print('storage created')
    iterator = chunk_genome(chrom_name_length_file)
    func = partial(worker_multi, archive_name=archive_name, bam_name=alignment_file, translator=translator)
    with closing(mp.Pool()) as pool:
        # don't need lazy mapping, just passing chunk ids
        pool.map(func, iterator)


# MAKE SURE THERE IS BAM INDEX FILE IN SAME LOCATION
def test(proc_type):
    """
    :param str proc_type: [single | multiple]
    :return:
    """
    bamfile = os.path.expanduser('~/Dropbox/Research Peer/dbseq/sorted.bam')
    chrnamefile = os.path.expanduser('~/Dropbox/Research Peer/dbseq/chrNameLength.txt')
    annotation = os.path.expanduser('~/Dropbox/Research Peer/dbseq/annotations.gtf')
    translator = gtf.GeneIntervals(annotation, 10000)
    dir_ = os.environ['TMPDIR']
    archive_name = dir_ + 'test.h5'
    if os.path.isfile(archive_name):
        os.remove(archive_name)
    if proc_type == 'single':
        main_single(
            archive_name, bamfile, chrnamefile, translator)
    elif proc_type == 'multiple':
        main_multiple(
            archive_name, bamfile, chrnamefile, translator)



if __name__ == "__main__":
    start = time.time()
    #test('single')
    #end1 = time.time()
    test('multiple')
    end1 = time.time()
    print(end1-start)

