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


class DBSeq:

    GENOME_CHUNKSIZE = int(1e7)
    ARRAY_CHUNKSIZE = int(1e6)
    LOCK = mp.RLock()
    _DT = np.dtype([
        ('status', np.uint8),
        ('cell', np.int64),
        ('rmt', np.int32),
        ('n_poly_t', np.uint8),
        ('gene', np.uint32),
        ('position', np.uint64),
        ('n_aligned', np.uint8)])

    def __init__(self, archive_name):

        self.archive_name = archive_name
        self.archive = tb.open_file(self.archive_name, "r")
        self.table = self.archive.root.records

    def close(self):

        self.archive.close()

    @classmethod
    def from_alignment_file(cls, archive_name, bam_name, chrom_name_length_file, translator):
        cls.create_storage(archive_name)
        print('storage created')
        iterator = cls.chunk_genome(chrom_name_length_file)
        func = partial(cls.worker_multi, archive_name=archive_name,
                       bam_name=bam_name, translator=translator)
        with closing(mp.Pool()) as pool:
            # don't need lazy mapping, just passing chunk ids
            pool.map(func, iterator)

        return cls(archive_name)

    @classmethod
    def create_storage(cls, archive_name):
        """
        :param str archive_name:
        :return:
        """     
        blosc5 = tb.Filters(complib='blosc', complevel=5)
        with closing(tb.open_file(archive_name, mode='w', filters=blosc5)) as h5:
            records = h5.create_table(h5.root, 'records', cls._DT, "Records")

    @classmethod
    def chunk_genome(cls, chrom_name_length_file):
        """parse chromsome information embedded in a STAR index to iterate over genome chunks
        :param str chrom_name_length_file:
        :return Iterator(str, int, int):
        """
        with open(chrom_name_length_file, 'r') as f:
            chr_sizes = dict(line.strip().split() for line in f.readlines())

        for chromosome, size in chr_sizes.items():
            size = int(size)
            start = 0
            end = cls.GENOME_CHUNKSIZE
            while True:
                if end > size:
                    yield (chromosome, start, size)
                    break
                else:
                    yield (chromosome, start, end)
                    start = end
                    end += cls.GENOME_CHUNKSIZE
    @classmethod
    def bam_iterator(cls, bam_name, chromosome, start, end):
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

    @classmethod
    def worker_multi(cls, chunk, archive_name, bam_name, translator):
        """
        :param (str, int, int) chunk:
        :param str archive_file:
        :param str bam_name:
        :return:
        """
        storage = np.zeros((cls.ARRAY_CHUNKSIZE,), dtype=cls._DT)
        print('processing', *chunk)
        i = 0
        for aligned_segment in cls.bam_iterator(bam_name, *chunk):
            if i == cls.ARRAY_CHUNKSIZE:
                with cls.LOCK:
                    with closing(tb.open_file(archive_name, mode='a')) as h5:
                        table = h5.root.records
                        table.append(storage)
                        table.flush()
                i = 0
                storage = np.zeros((cls.ARRAY_CHUNKSIZE, len(cls._DT)), dtype=cls._DT)
            else:
                storage[i] = cls.from_sam_record(aligned_segment, translator)
            i += 1

        storage = storage[:i % cls.ARRAY_CHUNKSIZE]
        with cls.LOCK:
            with closing(tb.open_file(archive_name, mode='a')) as h5:
                table = h5.root.records
                table.append(storage)
                table.flush()

    @classmethod
    def from_sam_record(cls, seg, translator):
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
        n_poly_t = processed_fields.poly_t.count(
            'T') + processed_fields.poly_t.count('N')

        return np.array([(0, cell, rmt, n_poly_t, gene, pos, n_aligned)], dtype=cls._DT)



def test_class():

    bamfile = os.path.expanduser('~/Dropbox/Research Peer/dbseq/sorted.bam')
    chrnamefile = os.path.expanduser(
        '~/Dropbox/Research Peer/dbseq/chrNameLength.txt')
    annotation = os.path.expanduser(
        '~/Dropbox/Research Peer/dbseq/annotations.gtf')
    translator = gtf.GeneIntervals(annotation, 10000)
    dir_ = os.environ['TMPDIR']
    archive_name = dir_ + 'test.h5'

    d = DBSeq.from_alignment_file(
        archive_name, bamfile, chrnamefile, translator)
    entries = d.table.read()
    print(entries[0])

    d.close()


if __name__ == "__main__":
    start = time.time()
    # test('single')
    #end1 = time.time()
    test_class()
    end1 = time.time()
    print(end1 - start)
