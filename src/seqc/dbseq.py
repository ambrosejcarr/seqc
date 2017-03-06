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
import itertools


class DBSeq:

    GENOME_CHUNKSIZE = int(1e7)
    ARRAY_CHUNKSIZE = int(1e6)
    LOCK = mp.RLock()
    __DT = np.dtype([
        ('status', np.uint8),
        ('cell', np.int64),
        ('rmt', np.int32),
        ('n_poly_t', np.uint8),
        ('gene', np.uint32),
        ('position', np.uint64),
        ('n_aligned', np.uint8),
        ('read_id', np.uint32)])

    def __init__(self, archive_name):
        self.archive_name = archive_name
        

    def __enter__(self):
        self.archive = tb.open_file(self.archive_name, "r")
        self.table = self.archive.root.records
        return self

    def __exit__(self, *args):
        self.archive.close()

    #def close(self):
        #self.archive.close()

    def cell_selector(self, row):
        return row['cell']

    def cell_counts(self,mapping):
        cells = {}

        if mapping == "multimapped":
            for cell, rows_grouped_by_cell in itertools.groupby(self.table, self.cell_selector):
                cells[cell] = sum(1 for r in rows_grouped_by_cell if r['n_aligned'] > 1)
        if mapping == "mapped":
            for cell, rows_grouped_by_cell in itertools.groupby(self.table, self.cell_selector):
                cells[cell] = sum(1 for r in rows_grouped_by_cell if r['n_aligned'] == 1)
        if mapping == "unmapped":
            for cell, rows_grouped_by_cell in itertools.groupby(self.table, self.cell_selector):
                cells[cell] = sum(1 for r in rows_grouped_by_cell if r['n_aligned'] == 0)
        return cells


    @classmethod
    def from_alignment_file(cls, archive_name, bam_name, chrom_name_length_file, translator):
        cls.__create_storage(archive_name)
        print('storage created')
        iterator = cls.__chunk_genome(chrom_name_length_file)
        func = partial(cls.worker_multi, archive_name=archive_name,
                       bam_name=bam_name, translator=translator)
        with closing(mp.Pool()) as pool:
            # don't need lazy mapping, just passing chunk ids
            pool.map(func, iterator)

        return cls(archive_name)

    @classmethod
    def __create_storage(cls, archive_name):
        """
        :param str archive_name:
        :return:
        """     
        blosc5 = tb.Filters(complib='blosc', complevel=5)
        with closing(tb.open_file(archive_name, mode='w', filters=blosc5)) as h5:
            records = h5.create_table(h5.root, 'records', cls.__DT, "Records")

    @classmethod
    def __chunk_genome(cls, chrom_name_length_file):
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
    @staticmethod
    def __bam_iterator(bam_name, chromosome, start, end):
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
        storage = np.zeros((cls.ARRAY_CHUNKSIZE,), dtype=cls.__DT)
        print('processing', *chunk)
        i = 0
        for aligned_segment in cls.__bam_iterator(bam_name, *chunk):
            if i == cls.ARRAY_CHUNKSIZE:
                with cls.LOCK:
                    with closing(tb.open_file(archive_name, mode='a')) as h5:
                        table = h5.root.records
                        table.append(storage)
                        table.flush()
                i = 0
                storage = np.zeros((cls.ARRAY_CHUNKSIZE, len(cls.__DT)), dtype=cls.__DT)
            else:
                storage[i] = cls.__from_sam_record(aligned_segment, translator)
            i += 1

        storage = storage[:i % cls.ARRAY_CHUNKSIZE]
        with cls.LOCK:
            with closing(tb.open_file(archive_name, mode='a')) as h5:
                table = h5.root.records
                table.append(storage)
                table.flush()

    @classmethod
    def __from_sam_record(cls, seg, translator):
        NameField = namedtuple(
            'NameField', ['pool', 'cell', 'rmt', 'poly_t', 'name'])
        qname = seg.query_name

        # Parse the name field
        fields, name = qname.split(';')
        processed_fields = fields.split(':')
        processed_fields.append(name)
        processed_fields = NameField(*processed_fields)

        # Get read ID
        name = processed_fields.name
        name = name.split(':')
        r_ID = name[len(name)-3] + name[len(name)-2] + name[len(name)-1]
        r_ID = int(r_ID)

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

        return np.array([(0, cell, rmt, n_poly_t, gene, pos, n_aligned, r_ID)], dtype=cls.__DT)


#OLD DON'T USE 
def test_class():
    #check for index file
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

def test_class_2():
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
    
    with d as database:
        cells_unmapped = database.cell_counts("unmapped")
        cells_mapped = database.cell_counts("mapped")
        cells_multimapped = database.cell_counts("multimapped")


if __name__ == "__main__":
    start = time.time()
    # test('single')
    #end1 = time.time()
    test_class_2()
    end1 = time.time()
    print(end1 - start)
