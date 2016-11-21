import os
import numpy as np
from scipy.sparse import coo_matrix
from collections import OrderedDict
from seqc.sequence.gtf import create_gene_id_to_official_gene_symbol_map
from seqc.sequence.gtf import ensembl_gene_id_to_official_gene_symbol


class SparseFrame:

    def __init__(self, data, index, columns):
        """
        lightweight wrapper of scipy.stats.coo_matrix to provide pd.DataFrame-like access
        to index, column, and shape properties.

        :param data: scipy.stats.coo_matrix
        :param index: np.ndarray: row index
        :param columns: np.ndarray: column index

        :property data: scipy.stats.coo_matrix
        :property index: np.ndarray row index
        :property columns: np.ndarray column index
        :property shape: (int, int), number of rows and columns
        :method sum: wrapper of np.sum()
        """

        if not isinstance(data, coo_matrix):
            raise TypeError('data must be type coo_matrix')
        if not isinstance(index, np.ndarray):
            raise TypeError('index must be type np.ndarray')
        if not isinstance(columns, np.ndarray):
            raise TypeError('columns must be type np.ndarray')

        self._data = data
        self._index = index
        self._columns = columns

    @property
    def data(self):
        return self._data

    @data.setter
    def data(self, item):
        if not isinstance(item, coo_matrix):
            raise TypeError('data must be type coo_matrix')
        self._data = item

    @property
    def index(self):
        return self._index

    @index.setter
    def index(self, item):
        try:
            self._index = np.array(item)
        except:
            raise TypeError('self.index must be convertible into a np.array object')

    @property
    def columns(self):
        return self._columns

    @columns.setter
    def columns(self, item):
        try:
            self._columns = np.array(item)
        except:
            raise TypeError('self.columns must be convertible into a np.array object')

    @property
    def shape(self):
        return len(self.index), len(self.columns)

    def sum(self, axis=0):
        """
        sum over provided axis

        :param axis: options: 0 (rows) or 1 (columns)
        :return: np.ndarray vector of column or row sums
        """
        return self.data.sum(axis=axis)

    @classmethod
    def from_dict(cls, dictionary, genes_to_symbols=False):
        """create a SparseFrame from a dictionary

        :param dict dictionary: dictionary in form (cell, gene) -> count
        :param str|bool genes_to_symbols: convert genes into symbols. If not False, user
          must provide the location of a .gtf file to carry out conversion. Otherwise the
          column index will retain the original integer ids
        :return SparseFrame: SparseFrame containing dictionary data
        """

        i, j = (np.array(v, dtype=int) for v in zip(*dictionary.keys()))
        data = np.fromiter(dictionary.values(), dtype=int)

        # map cells to small values
        uniq_i = np.unique(i)
        imap = OrderedDict(zip(uniq_i, np.arange(uniq_i.shape[0])))

        uniq_j = np.unique(j)
        jmap = OrderedDict(zip(uniq_j, np.arange(uniq_j.shape[0])))

        i_inds = np.fromiter((imap[v] for v in i), dtype=int)
        j_inds = np.fromiter((jmap[v] for v in j), dtype=int)

        coo = coo_matrix((data, (i_inds, j_inds)), shape=(len(imap), len(jmap)),
                         dtype=np.int32)

        index = np.fromiter(imap.keys(), dtype=int)
        columns = np.fromiter(jmap.keys(), dtype=int)

        if genes_to_symbols:
            if not os.path.isfile(genes_to_symbols):
                raise ValueError('genes_to_symbols argument %s is not a valid annotation '
                                 'file' % repr(genes_to_symbols))
            gmap = create_gene_id_to_official_gene_symbol_map(genes_to_symbols)
            columns = np.array(ensembl_gene_id_to_official_gene_symbol(
                columns, gene_id_map=gmap))

        return cls(coo, index, columns)
