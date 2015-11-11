__author__ = 'ambrose'


import numpy as np
import pandas as pd
from collections import defaultdict, Mapping
from scipy.sparse import coo_matrix
from seqc import plot, gtf
import pickle
import os


class CompareExperiments:
    pass


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

    def cells_above_threshold(self, t):
        """return the number of cells with > t observations"""
        return np.sum(self.counts.sum(axis=1) > t)

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

    def convert_ids(self, fgtf=None, scid_map=None):
        """
        Convert scids to gene identifiers either by parsing the gtf file (slow) or by
        directly mapping scids to genes (fast). In the latter case, the gtf_map must be
        pre-processed and saved to the index folder.

        see seqc.gtf.Reader.scid_to_gene(gtf, save=<output_file>) to save this gtf_map for
        repeat-use.
        """
        if not any([fgtf, scid_map]):
            raise ValueError('one of gtf or scid_map must be passed for conversion to'
                             'occur')

        # load the gene map (gmap)
        if not scid_map:
            r = gtf.Reader(fgtf)
            gmap = r.scid_to_gene()
        elif isinstance(scid_map, str):
            # try to load pickle
            if not os.path.isfile(scid_map):
                raise FileNotFoundError('scid_map not found: %s' % scid_map)
            with open(scid_map, 'rb') as f:
                gmap = pickle.load(f)
                if not isinstance(gmap, Mapping):
                    raise TypeError('scid_map file object did not contain a '
                                    'dictionary')
        elif isinstance(scid_map, Mapping):
            gmap = scid_map
        else:
            raise TypeError('scid_map must be the location of the scid_map pickle file '
                            'or the loaded dictionary object, not type: %s' %
                            type(scid_map))

        # convert ids
        self._index = np.array([gmap[i] for i in self._index])

    def plot_observations_per_cell(self, smooth=False, xlabel='', ylabel='',
                                   title=''):
        cellsums = self.counts.sum(axis=1)
        f, ax = plot.get_axes()
        if not xlabel:
            xlabel = 'observations per cell'
        if not ylabel:
            ylabel = 'number of cells'
        if not title:
            title = 'distribution of observations per cell'
        if smooth:
            plot.kde(cellsums, fig=f, ax=ax, xlabel=xlabel, ylabel=ylabel,
                     title=title)
        else:
            plot.histogram(cellsums, fig=f, ax=ax, xlabel=xlabel, ylabel=ylabel,
                           title=title)
        return f, ax

    def to_npz(self, npz_file):
        """
        Save an .npz archive containing all information necessary to reconstruct the
        SparseCounts object
        """
        raise NotImplementedError

    @classmethod
    def from_npz(cls, npz_file):
        """
        Load a SparseCounts object from a .npz file.
        """
        raise NotImplementedError


# class ReadsPerMolecule:
#
#     def __init__(self, molecule_counts):
#         if not isinstance(molecule_counts, pd.Series):
#             raise ValueError('molecule_counts must be a pandas Series')
#         self._counts = molecule_counts
#
#     @property
#     def counts(self):
#         return self._counts
#
#     @classmethod
#     def from_read_array(cls, read_array):
#         # filter failing molecules
#         record_mask = read_array.mask_failing_records()
#         support_mask = read_array.mask_low_support_molecules()
#         view = read_array.data[record_mask & support_mask][['cell', 'rmt']]
#
#         counts = pd.DataFrame(view).groupby(['cell', 'rmt']).size()
#         return cls(counts)
#
#     def plot_distribution(self, smooth=False, log=False):
#
#         if log:
#             vals = np.log(self.counts.values)
#         else:
#             vals = self.counts.values
#
#         with sns.set_style('ticks'):
#             f, ax = plt.subplots()
#             if smooth:
#                 sns.kdeplot(vals, ax=ax)
#             else:
#                 ax.hist(vals, bins=50)
#
#         # label
#         plt.xlabel('molecule counts')
#         plt.ylabel('relative frequency')
#         plt.title('Reads per Molecule Distribution')
#         sns.despine()
#         plt.tight_layout()
#         return f, ax
