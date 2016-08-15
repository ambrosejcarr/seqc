import numpy as np
from scipy.sparse import coo_matrix


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
