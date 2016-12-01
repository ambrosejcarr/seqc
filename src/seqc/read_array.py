from multiprocessing import Pool
from contextlib import closing
import numpy as np
from seqc.alignment import sam
from seqc.sequence.encodings import DNA3Bit
from scipy.sparse import coo_matrix, csr_matrix
import seqc.sequence.barcodes
import tables as tb
from itertools import chain


class ReadArray:

    _dtype = [
        ('status', np.uint8),  # if > 8 tests, change to int16
        ('cell', np.int64),
        ('rmt', np.int32),
        ('n_poly_t', np.uint8)]

    def __init__(self, data, genes, positions):
        """
        Enhanced np.ndarray (structured array) with several companion functions to hold
        filter, sort, and access compressed fastq read information.

        :param np.recarray data: np.ndarray with dtype of self._dtype
        :param csr_matrix|np.ndarray genes:
        :param csr_matrix|np.ndarray positions:

        """

        if not isinstance(genes, (csr_matrix, np.ndarray)):
            raise TypeError(
                'genes must be a scipy csr_matrix or np.array, not %s'
                % repr(type(genes)))
        self._genes = genes
        if not isinstance(positions, (csr_matrix, np.ndarray)):
            raise TypeError(
                'positions must be a scipy csr_matrix or np.array, not %s'
                % repr(type(positions)))
        self._positions = positions

        if not isinstance(data, np.ndarray):
            raise TypeError(
                'data must be a structured np.array object, not %s' % repr(type(data)))
        self._data = data

        if isinstance(genes, csr_matrix):  # track if genes/positions are csr or np.array
            self._ambiguous_genes = True
        else:
            self._ambiguous_genes = False

    @property
    def data(self):
        return self._data

    @property
    def genes(self):
        return self._genes

    @property
    def positions(self):
        return self._positions

    def __len__(self):
        return len(self.data)

    def __iter__(self):
        """Returns an interator over the RecordArray

        If gene or position are zero, this indicates that those fields are empty. If the
        ReadArray has not yet been disambiguated, this function returns the first gene and
        position in each multialignment, if one is present.

        :return Iterator: iterator over (recarray_row, gene, position)
        """
        # todo test out speed loss with chain(), compare to numpy.concatenate, None
        if self._ambiguous_genes:  # genes, positions are scipy.sparse.csr_matrix
            for i in np.arange(self.data.shape[0]):
                yield chain(self.data[i], (self.genes[i, 0], self.positions[i, 0]))
        else:  # genes, positions are np.array
            for i in np.arange(self.data.shape[0]):
                yield chain(self.data[i], (self.genes[i], self.positions[i]))

    filter_codes = {
        'no_gene': 0b1,
        'no_rmt': 0b10,
        'no_cell': 0b100,  # todo no filter in this module for this method yet
        'low_polyt': 0b1000,
        'gene_not_unique': 0b10000,
        'no_spacer': 0b100000,
    }

    def calculate_filtering_results(self, required_poly_t=1):
        """Apply different filters to the read array. If a read fails a filter, it is not
        passed to the others and so the counts for each filter is the number of reads that
        failed that one but passed all previous ones Filters are not ordered in any
        particular way.

        :param required_poly_t: minimum number of T nucleotides in the primer tail of a
          valid read
        :return None: method sets the status vector according to the result of filtering.
        """

        failing = np.zeros(self.data.shape[0], dtype=np.int8)

        # genes are dealt with differently depending on the state of the array
        if self._ambiguous_genes:
            nnz = self.genes.getnnz(axis=1)
            failing[nnz == 0] |= self.filter_codes['no_gene']
            failing[nnz > 1] |= self.filter_codes['gene_not_unique']
        else:  # multiple gene filter is empty
            failing[self.genes == 0] |= self.filter_codes['no_gene']

        failing[self.data['rmt'] == 0] |= self.filter_codes['no_rmt']
        failing[self.data['cell'] == 0] |= self.filter_codes['no_cell']
        failing[self.data['n_poly_t'] < required_poly_t] |= self.filter_codes['low_polyt']

        self.data['status'] = np.bitwise_or(self.data['status'], failing)

    def filtering_mask(self, *ignore):
        """return a filtering mask that, when compared to status, will return False for
        reads failing any filters that are not specified in ignore

        :param *str ignore: fields to ignore when constructing filtering mask
        :return int: mask
        """
        mask = 2 ** 8 - 1  # 0b111... x16, permissive mask
        for filter_ in ignore:
            try:
                mask ^= self.filter_codes[filter_]  # mask filter
            except KeyError:
                raise KeyError('%s is not a valid filter. Please select from %s' %
                               (filter_, repr(self.filter_codes.keys())))
        return mask

    def iter_active(self, *ignore):
        """Iterator over active reads, ignoring any filter passed in ignore

        :param *str ignore: values to ignore when parsing active reads. choices:
          [no_gene, no_rmt, no_cell, low_polyt, gene_not_unique, no_spacer]
        :yields int, (np.array, int, int): iterator yields active records in form
          (index, (data_row, gene, pos))
        """

        if not ignore:  # save a bit of work by not &ing with mask
            for i, data in enumerate(self):
                if data['status'] == 0:
                    yield i, data

        else:  # create the appropriate mask for filters we want to ignore
            mask = self.filtering_mask(*ignore)
            for i, data in enumerate(self):
                if not data['status'] & mask:  # ignores fields in ignore
                    yield i, data

    @classmethod
    def from_alignment_file(cls, alignment_file, translator, required_poly_t):
        """
        construct a ReadArray object from a samfile containing only uniquely aligned
        records

        :param required_poly_t: number of poly_t required for a read to be considered
          a valid alignment
        :param GeneIntervals translator: translator created from the .gtf annotation
          file corresponding to the genome against which the reads in sam_file were
          aligned
        :param str alignment_file: filename of alignment file.
        :return:
        """

        # todo add a check for @GO query header (file matches sorting assumptions)

        reader = sam.Reader(alignment_file)  # todo swap to pysam reader, probably faster

        # todo allow reading of this from alignment summary
        num_reads = 0
        num_unique = 0
        prev_alignment_name = ''
        for alignment in reader:
            num_reads += 1
            if alignment.qname != prev_alignment_name:
                num_unique += 1
                prev_alignment_name = alignment.qname

        # pre-allocate arrays
        data = np.recarray((num_unique,), cls._dtype)
        row = np.zeros(num_reads, dtype=np.int32)
        col = np.zeros(num_reads, dtype=np.int8)
        position = np.zeros(num_reads, dtype=np.int32)
        gene = np.zeros(num_reads, dtype=np.int32)

        # loop over multialignments
        row_idx = 0  # identifies the read index
        arr_idx = 0  # identifies the alignment index across all reads
        max_ma = 0  # maximum number of alignments observed for a read
        for ma in reader.iter_multialignments():
            col_idx = 0  # identifies the alignment index within a read
            for a in ma:
                genes = translator.translate(a.rname, a.strand, a.pos)

                # if gene passes filter, store passing genes and increment appropriate
                # arr and col indices
                if genes is not None and len(genes) == 1:
                    row[arr_idx] = row_idx
                    col[arr_idx] = col_idx
                    gene[arr_idx] = genes[0]
                    position[arr_idx] = a.pos
                    arr_idx += 1
                    col_idx += 1
            max_ma = max(max_ma, col_idx)

            cell = seqc.sequence.encodings.DNA3Bit.encode(a.cell)
            rmt = seqc.sequence.encodings.DNA3Bit.encode(a.rmt)
            n_poly_t = a.poly_t.count('T') + a.poly_t.count('N')
            data[row_idx] = (0, cell, rmt, n_poly_t)
            row_idx += 1

        # some reads will not have aligned, throw away excess allocated space before
        # creating the ReadArray
        row, col, position, gene = (
            row[:arr_idx], col[:arr_idx], position[:arr_idx], gene[:arr_idx])

        gene = coo_matrix((gene, (row, col)), shape=(row_idx, max_ma), dtype=np.int32)
        position = coo_matrix((position, (row, col)), shape=(row_idx, max_ma),
                              dtype=np.int32)

        ra = cls(data, gene.tocsr(), position.tocsr())
        ra.calculate_filtering_results(required_poly_t=required_poly_t)
        return ra

    def group_filtered_indices(self, by, multimapping=False):
        """group the reads in ra.data by columns specified in "by", then return a tuple
        containing one array per group that consists of all the indicies corresponding
        to the positions of reads in self.data which pass filters.

        example usage:

        # group over molecules for multimapping
        self.group_filtered_indices(by=('cell, 'rmt'), multimapping=True)

        # group over gene/cell pairs for rmt correction; additional sort must be done
        # over genes within the error correction functions
        self.group_filtered_indices(by='cell')

        :param bool multimapping: if True, then reads are not filtered if they have
          alignments to more than one gene
        :param tuple|str by: group over these columns of self.data
        :return [np.array]: list of numpy arrays, each containing all of the indices for
          reads that correspond to a group, defined as a unique combination of the
          columns specified in parameter by.
        """

        idx = np.argsort(self.data, order=by)  # get index that executes sort

        # filter the index for reads that
        if multimapping:
            mask = self.filtering_mask('gene_not_unique')
            passing = (self.data['status'] & mask) == 0
        else:
            passing = self.data['status'] == 0
        idx = idx[passing]

        # determine which positions in idx are the start of new groups (boolean, True)
        breaks = np.zeros(idx.shape[0] - 1)
        for column in by:
            breaks |= np.diff(self.data[column][idx])

        # convert boolean positions to indices, add start and end points.
        molecule_breaks = np.where(breaks)[0] + 1

        # use these break points to split the filtered index according to "by"
        return np.split(idx, molecule_breaks)

    def _resolve_ambiguous_alignments(self, indices):
        """private function to resolve multi-aligning reads that is mapped over a Process
        pool.

        :param np.array indices: array of indices for all reads corresponding to a
          (cell, rmt) pair, defining one or more molecules.
        :return None: This function is side-effecting. It sets each ambiguous row
          identified by indices equal to the corrected gene and position.
        """
        raise NotImplementedError  # todo implement me

    def resolve_ambiguous_alignments(self):
        """Resolves ambiguous alignments according to the following logic:

        # todo more detail
        1. split reads that align to disjoint sets of molecules
        2. if now unique, return
        3. else determine if there is a single gene that explains all alignments. If so,
           set all alignments to be equal to that gene.

        :return None: this is a side-effecting function which produces the following
           modifications:
           1. self.genes, self.positions are transformed from scipy.sparse.csr_matrix to
              np.array vectors according to the above logic. If no unambiguous gene is
              found that explains the data, that read's gene and position are set to 0.
           2. self._ambiguous_genes is set to False
        """
        indices_grouped_by_molecule = self.group_filtered_indices(
            by=('cell', 'rmt'), multimapping=True)

        with closing(Pool()) as pool:
            pool.map(self._resolve_ambiguous_alignments, indices_grouped_by_molecule)

        self._genes = np.ravel(self.genes.tocsc()[:, 0])  # todo add setter to self.genes
        self._positions = np.ravel(self._positions.tocsc()[:, 0])

        # update filtering results
        self.data['status'] &= self.filtering_mask('gene_not_unique')
        self.data['status'][self.genes == 0] |= self.filter_codes['no_gene']
        self._ambiguous_genes = False

    def save(self, archive_name: str) -> None:
        """save a ReadArray object as an hdf5 archive

        :param archive_name: filestem for the new archive
        :return: None
        """

        def store_carray(archive, array, name):
            atom = tb.Atom.from_dtype(array.dtype)
            store = archive.create_carray(archive.root, name, atom, array.shape)
            store[:] = array
            store.flush()

        if not archive_name.endswith('.h5'):
            archive_name += '.h5'

        # construct container
        blosc5 = tb.Filters(complevel=5, complib='blosc')
        f = tb.open_file(archive_name, mode='w', title='Data for seqc.ReadArray',
                         filters=blosc5)

        f.create_table(f.root, 'data', self.data)

        if self._ambiguous_genes:
            # each array is data, indices, indptr
            store_carray(f, self.genes.indices, 'indices')
            store_carray(f, self.genes.indptr, 'indptr')
            store_carray(f, self.genes.data, 'gene_data')
            store_carray(f, self.positions.data, 'positions_data')
        else:
            store_carray(f, self.genes, 'genes')
            store_carray(f, self.positions, 'positions')

        f.close()

    @classmethod
    def load(cls, archive_name: str):
        """load a ReadArray from an hdf5 archive, note that ma_pos and ma_genes are
        discarded.

        :param archive_name: name of a .h5 archive containing a saved ReadArray object
        :return: ReadArray
        """

        f = tb.open_file(archive_name, mode='r')
        data = f.root.data.read()

        try:
            f.get_node('/genes')
            genes = f.root.genes.read()
            positions = f.root.positions.read()
        except tb.NoSuchNodeError:
            indptr = f.root.indptr.read()
            indices = f.root.indices.read()
            genes = f.root.gene_data.read()
            positions = f.root.positions_data.read()
            # shape = (indices.shape[0] - 1, np.max(indptr))
            genes = csr_matrix((genes, indices, indptr))
            positions = csr_matrix((positions, indices, indptr))

        f.close()

        return cls(data, genes, positions)
