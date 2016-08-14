import numpy as np
import tables as tb
from seqc.alignment import sam
from seqc.sequence.gtf import GeneIntervals
from seqc.sequence.encodings import DNA3Bit


class ReadArray:

    _dtype = [
        ('pool', np.int8),
        ('cell', np.int64),
        ('rmt', np.int32),
        ('n_poly_t', np.uint8),
        ('dust_score', np.uint8),
        ('gene', np.uint32),
        ('position', np.uint64)]

    def __init__(self, data):
        """
        Enhanced np.ndarray (structured array) with several companion functions to hold
        filter, sort, and access compressed fastq read information.

        :param data: np.ndarray with dtype of self._dtype

        :property data: stored np.ndarray
        :method save: saves the ReadArray in compressed .h5 format
        :method load: load a saved compressed .h5 representation of a ReadArray
        :method from_samfile: constructs a ReadArray object from a samfile (uniquely
          aligned records only)
        :method reads_passing_filters: Return a ReadArray containing only reads
          that pass all filters (n_poly_t, dust_complexity_score, presence of cell and
          rmt)
        """
        self._data = data

    @property
    def data(self):
        return self._data

    @classmethod
    def from_samfile(cls, samfile: str, gtf: str):
        """
        construct a ReadArray object from a samfile containing only uniquely aligned
        records

        :param gtf: str, filename of annotations.gtf file
        :param samfile: str, filename of alignment file.
        :return:
        """
        reader = sam.Reader(samfile)
        translator = GeneIntervals(gtf)
        num_records = 0

        data = np.recarray((len(reader),), cls._dtype)
        for i, alignment in enumerate(reader):
            num_records += 1
            gene = translator.translate(
                    alignment.rname, alignment.strand, alignment.pos)
            if gene is None:
                gene = 0
            pool = DNA3Bit.encode(alignment.pool)
            cell = DNA3Bit.encode(alignment.cell)
            rmt = DNA3Bit.encode(alignment.rmt)
            n_poly_t = alignment.poly_t.count(b'T') + alignment.poly_t.count(b'N')
            dust_score = alignment.dust_low_complexity_score
            data[i] = (pool, cell, rmt, n_poly_t, dust_score, gene, alignment.pos)

        return cls(data), num_records

    def reads_passing_filters(self, min_poly_t: int, max_dust_score: int):
        """
        Subset the ReadArray returning a new ReadArray containing eads that passed all
        filters

        :param min_poly_t: int, minimum number of T nucleotides that defined a valid
          capture primer
        :param max_dust_score: int, higher scores indicate increasingly degenerate
          sequences. Typically sequences with dust_score > 10 may be discarded.
        :return: ReadArray
        """
        phix_genes = np.array(range(1, 7)) * 111111111
        data = self.data[((self.data['n_poly_t'] >= min_poly_t) &
                          (self.data['dust_score'] <= max_dust_score) &
                          (self.data['gene'] != 0) &
                          (self.data['cell'] != 0) &
                          (self.data['rmt'] != 0))]

        # filter out phiX genes
        not_phix = ~np.in1d(data['gene'], phix_genes)
        data = data[not_phix]

        # filter out N's in rmt
        res = np.zeros(len(data), dtype=np.bool)
        rmt = data['rmt'].copy()
        while np.any(rmt):
            n_filter = rmt & 0b111 == 0b111
            res[n_filter] = True
            rmt >>= 3
        data = data[~res]
        return ReadArray(data)

    def save(self, archive_name: str) -> None:
        """save a ReadArray in .h5 format

        :param archive_name: filename of a new .h5 archive in which to save the ReadArray
        :return: None
        """

        # create table
        blosc5 = tb.Filters(complevel=5, complib='blosc')
        f = tb.open_file(archive_name, mode='w', title='Data for seqc.ReadArray',
                         filters=blosc5)

        # store data
        f.create_table(f.root, 'data', self._data)
        f.close()

    @classmethod
    def load(cls, archive_name: str):
        """load a ReadArray from a .h5 archive

        :param archive_name: name of a .h5 archive containing a saved ReadArray object
        :return: ReadArray
        """

        f = tb.open_file(archive_name, mode='r')
        data = f.root.data.read()
        f.close()
        return cls(data)

    def __len__(self):
        return len(self.data)
