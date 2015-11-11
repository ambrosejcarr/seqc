__author__ = 'ambrose'


import numpy as np
import pandas as pd
from collections import defaultdict
from scipy.sparse import coo_matrix
import seqc.plot

# check if display variable is set. If not, use a limited plotting suite.


class CompareExperiments:
    pass


class ReadsPerMolecule:

    def __init__(self, molecule_counts):
        if not isinstance(molecule_counts, pd.Series):
            raise ValueError('molecule_counts must be a pandas Series')
        self._counts = molecule_counts

    @property
    def counts(self):
        return self._counts

    @classmethod
    def from_read_array(cls, read_array):
        # filter failing molecules
        record_mask = read_array.mask_failing_records()
        support_mask = read_array.mask_low_support_molecules()
        view = read_array.data[record_mask & support_mask][['cell', 'rmt']]

        counts = pd.DataFrame(view).groupby(['cell', 'rmt']).size()
        return cls(counts)

    def plot_distribution(self, smooth=False, log=False):

        if log:
            vals = np.log(self.counts.values)
        else:
            vals = self.counts.values

        with sns.set_style('ticks'):
            f, ax = plt.subplots()
            if smooth:
                sns.kdeplot(vals, ax=ax)
            else:
                ax.hist(vals, bins=50)

        # label
        plt.xlabel('molecule counts')
        plt.ylabel('relative frequency')
        plt.title('Reads per Molecule Distribution')
        sns.despine()
        plt.tight_layout()
        return f, ax


class SparseCounts:

    def __init__(self, sparse_counts, index, columns):
        self._counts = sparse_counts
        self._index = index
        self._columns = columns

    @property
    def counts(self):
        return self._counts

    @property
    def index(self):
        return self._index

    @property
    def columns(self):
        return self._columns

    @classmethod
    def from_read_array(cls, read_array, collapse_molecules=True, n_poly_t_required=3):

        # mask failing cells and molecules with < 2 reads supporting them.
        read_mask = read_array.mask_failing_records(n_poly_t_required)
        low_coverage_mask = read_array.mask_low_support_molecules()
        unmasked_inds = np.arange(read_array.data.shape[0])[read_mask & low_coverage_mask]
        if unmasked_inds.shape[0] == 0:
            raise ValueError('Zero reads passed filters. SparseCounts would be empty')
        molecule_counts = defaultdict(dict)

        if collapse_molecules:
            # get molecule counts

            for i in unmasked_inds:
                feature = read_array.features[i]
                if len(feature) > 1:
                    continue
                cell = read_array.data['cell'][i]
                rmt = read_array.data['rmt'][i]
                try:
                    molecule_counts[int(feature)][cell].add(rmt)
                except KeyError:
                    molecule_counts[int(feature)][cell] = {rmt}

            # convert to molecule counts
            for f in molecule_counts.keys():
                for c, rmts in molecule_counts[f].items():
                    molecule_counts[f][c] = len(rmts)
        else:
            for i in unmasked_inds:
                feature = read_array.features[i]
                if len(feature) > 1:
                    continue
                cell = read_array.data['cell'][i]
                try:
                    molecule_counts[int(feature)][cell] += 1
                except KeyError:
                    molecule_counts[int(feature)][cell] = 1

        # convert to values, row, col form for scipy.coo
        # pre-allocate arrays
        size = sum(len(c) for c in molecule_counts.values())
        values = np.empty(size, dtype=int)
        row = np.empty(size, dtype=int)
        col = np.empty(size, dtype=int)
        i = 0
        for feature in molecule_counts:
            for cell, count in molecule_counts[feature].items():
                values[i] = count
                row[i] = cell
                col[i] = feature
                i += 1

        # get max count to shrink dtype if possible
        maxcount = np.max(values)

        # set dtype
        if 0 < maxcount < 2 ** 8:
            dtype = np.uint8
        elif maxcount < 2 ** 16:
            dtype = np.uint16
        elif maxcount < 2 ** 32:
            dtype = np.uint32
        elif maxcount < 2 ** 64:
            dtype = np.uint64
        elif maxcount < 0:
            raise ValueError('Negative count value encountered. These values are not'
                             'defined and indicate a probable upstream bug')
        else:
            raise ValueError('Count values too large to fit in int64. This is very '
                             'unlikely, and will often cause Memory errors. Please check '
                             'input data.')

        # map row and cell to integer values for indexing
        unq_row = np.unique(row)  # these are the ids for the new rows / cols of the array
        unq_col = np.unique(col)
        row_map = dict(zip(unq_row, np.arange(unq_row.shape[0])))
        col_map = dict(zip(unq_col, np.arange(unq_col.shape[0])))
        row_ind = np.array([row_map[i] for i in row])
        col_ind = np.array([col_map[i] for i in col])

        # change dtype, set shape
        values = values.astype(dtype)
        shape = (unq_row.shape[0], unq_col.shape[0])

        # return a sparse array
        coo = coo_matrix((values, (row_ind, col_ind)), shape=shape, dtype=dtype)
        return cls(coo, unq_row, unq_col)

    # def plot_observations_per_cell(self, ax=None):
    #     cell_sums = self.counts.sum(axis=1)
    #     f, ax = seqc.plot.get_axes(ax=ax)

    def to_dense(self, convert_ids=True):
        """
        Implement once we have a better sense of how to threshold cells

        return a pandas dataframe with rows = cells and columns = scids or genes if
        convert_ids is True
        """
        # mask cells that fail filter
        # convert to pandas dataframe
        # merge scids that associate with the same genes
        # return df or DenseCounts (not sure if necessary to have a separate class)
        raise NotImplementedError


