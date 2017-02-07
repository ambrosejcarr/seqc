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
import gtf


# global values
K = 3
GENOME_CHUNKSIZE = int(1e7)
ARRAY_CHUNKSIZE = int(1e6)
KMER_MAP = {
    # change 3 to modify k
    kmer: i for i, kmer in enumerate(product('ACGTN', repeat=K))
}
LOCK = mp.RLock()


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
    dt = np.dtype([
        ('status', np.uint8),
        ('cell', np.int64),
        ('rmt', np.int32),
        ('n_poly_t', np.uint8),
        ('gene', np.uint32),
        ('position', np.uint64),
        ('n_aligned', np.uint8)])

    blosc5 = tb.Filters(complib='blosc', complevel=5)
    with closing(tb.open_file(archive_name, mode='w', filters=blosc5)) as h5:
        records = h5.create_table(h5.root, 'records', dt, "Records")


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

    def store(h5):
        with LOCK:
            h5.flush()
            h5.close()

    print('processing', *chunk)
    h5 = tb.open_file(archive_name, mode='a')
    table = h5.root.records
    i = 0
    for aligned_segment in bam_iterator(bam_name, *chunk):
        if i == ARRAY_CHUNKSIZE:
            store(h5)
            i = 0
            h5 = tb.open_file(archive_name, mode='a')
            table = h5.root.records
        else:
            rec = table.row
            #rec['status'] = 0
            rec['cell'], rec['rmt'], rec['n_poly_t'], rec['gene'], rec['position'], rec['n_aligned'] = sam_record(aligned_segment, translator)
            rec.append()
        i += 1

    store(h5)

def worker2(chunk, archive_name, bam_name, translator):
    """
    :param (str, int, int) chunk:
    :param str archive_name:
    :param str bam_name:
    :return:
    """

    
    with LOCK:

        print('processing', *chunk)
        h5 = tb.open_file(archive_name, mode='a')
        table = h5.root.records
        i = 0
        for aligned_segment in bam_iterator(bam_name, *chunk):
            if i == ARRAY_CHUNKSIZE:
                h5.flush()
                i = 0
            else:
                rec = table.row
                #rec['status'] = 0
                rec['cell'], rec['rmt'], rec['n_poly_t'], rec['gene'], rec['position'], rec['n_aligned'] = sam_record(aligned_segment, translator)
                rec.append()
            i += 1

        store(h5)


def sam_record(seg, translator):
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

    # FIX FOR GITHUB
    cell = DNA3Bit.encode(processed_fields.cell)
    rmt = DNA3Bit.encode(processed_fields.rmt)
    n_poly_t = processed_fields.poly_t.count('T') + processed_fields.poly_t.count('N')

    return cell, rmt, n_poly_t, gene, pos, n_aligned



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
    func = partial(worker2, archive_name=archive_name, bam_name=alignment_file, translator=translator)
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


class DNA3Bit(object):
    """
    Compact 3-bit encoding scheme for sequence data.
    """
    
    @staticmethod
    def bits_per_base():
        return 3

# TODO: The sam reader needs to be fixed so text files are read as text not binary
    str2bindict = {65: 0b100, 67: 0b110, 71: 0b101, 84: 0b011, 78: 0b111,
                   97: 0b100, 99: 0b110, 103: 0b101, 116: 0b011, 110: 0b111,
                   'A': 0b100, 'C': 0b110, 'G': 0b101, 'T': 0b011, 'N': 0b111,
                   'a': 0b100, 'c': 0b110, 'g': 0b101, 't': 0b011, 'n': 0b111}
    bin2strdict = {0b100: b'A', 0b110: b'C', 0b101: b'G', 0b011: b'T', 0b111: b'N'}
    
    @staticmethod
    def encode(b) -> int:
        """
        Convert string nucleotide sequence into binary, note: string is stored so
        that the first nucleotide is in the MSB position
        :param bytes|str b: sequence containing nucleotides to be encoded
        """
        res = 0
        for c in b:
            res <<= 3
            res += DNA3Bit.str2bindict[c]
        return res
        
    @staticmethod
    def decode(i: int) -> bytes:
        """
        Convert binary nucleotide sequence into string
        :param i: int, encoded sequence to be converted back to nucleotides
        """
        if i < 0:
            message = 'i must be an unsigned (positive) integer, not {0!s}'.format(i)
            raise ValueError(message)
        r = b''
        while i > 0:
            r = DNA3Bit.bin2strdict[i & 0b111] + r
            i >>= 3
        return r
        
    # TODO: another ooption is to use i.bit_length and take into account preceding 0's
    @staticmethod
    def seq_len(i: int) -> int:
        """
        Return the length of an encoded sequence based on its binary representation
        :param i: int, encoded sequence
        """
        l = 0
        while i > 0:
            l += 1
            i >>= 3
        return l
        
    @staticmethod
    def contains(s: int, char: int) -> bool:
        """
        return true if the char (bin representation) is contained in seq (binary
        representation)
        :param char: int, encoded character (one must be only one nucleotide)
        :param s: int, sequence of encoded nucleotides
        """
        while s > 0:
            if char == (s & 0b111):
                return True
            s >>= 3
        return False
    
    @staticmethod
    def ints2int(ints):
        """
        convert an iterable of sequences [i1, i2, i3] into a concatenated single integer
        0bi1i2i3. In cases where the sequence is longer than 64 bits, python will
        transition seamlessly to a long int representation, however the user must be
        aware that downsteam interaction with numpy or other fixed-size representations
        may not function
        :param ints: iterable of encoded sequences to concatenate
        """

        res = 0
        for num in ints:
            tmp = num
            # Get length of next number to concatenate (with enough room for leading 0's)
            while tmp > 0:
                res <<= 3
                tmp >>= 3
            res += num
        return res
    
    @staticmethod
    def count(seq, char_bin):
        """
        count how many times char is in seq.
        char needs to be an encoded value of one of the bases.
        """
        if char_bin not in DNA3Bit.bin2strdict.keys():
            raise ValueError("DNA3Bit.count was called with an invalid char code - "
                             "{}".format(char_bin))
        res = 0
        while seq > 0:
            if seq & 0b111 == char_bin:
                res += 1
            seq >>= 3
        return res

if __name__ == "__main__":
    test('multiple')
