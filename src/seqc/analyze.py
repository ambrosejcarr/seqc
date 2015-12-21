__author__ = 'ambrose'


import seqc
import numpy as np
import pandas as pd
from collections import defaultdict, Mapping, Iterable
from scipy.sparse import coo_matrix
try:
    from tsne import bh_sne
except ImportError:
    bh_sne = None
from sklearn.preprocessing import PolynomialFeatures
from sklearn.linear_model import LinearRegression
from sklearn.pipeline import Pipeline
from sklearn.mixture import GMM
import pickle
import matplotlib
import os
try:
    os.environ['DISPLAY']
except KeyError:
    matplotlib.use('Agg')
import matplotlib.pyplot as plt


# Below are functions that are shared by sparse and dense counts. However, it is
# not appropriate to have one class inherit from the other, so they are defined once here
# and referenced in both classes
def _plot_cell_gc_content_bias(cell_counts, sequences, fig=None, ax=None,
                               molecules=False, reads=False):
    # figure out the y label
    if all([molecules, reads]):
        raise ValueError('cell_counts must be either molecules or reads, not both')
    elif molecules:
        ctype = '(molecules)'
    elif reads:
        ctype = '(reads)'
    else:
        ctype = ''

    xlabel = 'Barcode GC Content (%)'
    ylabel = 'Library Size %s' % ctype
    title = 'Cell Sequencing Coverage vs. GC Content'

    # get gc content of sequences
    gc_content = np.array([seqc.encodings.ThreeBit.gc_content(s) for s in sequences])
    fig, ax = seqc.plot.scatter_density(gc_content, cell_counts, fig=fig, ax=ax,
                                        xlabel=xlabel, ylabel=ylabel, title=title)
    return fig, ax


def _plot_fraction_mitochondrial_rna(cell_sums, mt_sums, fig=None, ax=None,
                                     molecules=False, reads=False):
    # determine y label:
    if all([molecules, reads]):
        raise ValueError('cell_counts must be either molecules or reads, not both')
    elif molecules:
        ctype = '(molecules)'
    elif reads:
        ctype = '(reads)'
    else:
        ctype = ''

    xlabel = 'Library Size %s' % ctype
    ylabel = 'Mitochondrial fraction size'
    title = 'Mitochondrial fraction separates dying cells'

    mt_fraction = mt_sums / cell_sums
    mt_fraction[np.isinf(mt_fraction)] = 1
    mt_fraction[np.isnan(mt_fraction)] = 1
    seqc.plot.scatter_density(cell_sums, mt_fraction, fig=fig, ax=ax,
                              xlabel=xlabel, ylabel=ylabel, title=title)
    return fig, ax


def _discard_high_value_outliers():
    raise NotImplementedError


class SparseCounts:
    """
    SparseCounts is a thin wrapper around a scipy.sparse.coo_matrix object that provides
    additional functionality for the analysis of sc-RNA-seq data
    """

    def __init__(self, sparse_counts, index, columns):
        if not isinstance(sparse_counts, coo_matrix):
            raise TypeError('invalid input %s: sparse_counts must be of type coo_matrix '
                            ' not %s' % (repr(sparse_counts), type(sparse_counts)))
        self._counts = sparse_counts
        # weird bug with numpy generating masked arrays, but this fixes it...
        if isinstance(index, np.ma.MaskedArray):
            index = np.ma.compressed(index)
        self._index = index
        self._columns = columns

    def __repr__(self):
        return ('<SparseCounts wapper of coo_matrix with shape (%d, %d)>' %
                self.counts.shape)

    @property
    def counts(self):
        return self._counts

    @property
    def index(self):
        return self._index

    @property
    def columns(self):
        return self._columns

    def _select(self, mask, axis=1):
        """return a new array containing only objects passing mask"""
        if axis == 1:
            csr = self.counts.tocsr()
            selected_rows = csr[mask, :]
            csc = selected_rows.tocsc()
            selected_cols = np.ravel(csc.sum(axis=0) > 0)
            selected = csc[:, selected_cols].tocoo()
            index = self.index[selected_rows]
            columns = self.columns[selected_cols]
            return SparseCounts(selected, index, columns)
        else:
            raise NotImplementedError

    @classmethod
    def _deprecated_from_read_array(cls, read_array, collapse_molecules=True,
                                    n_poly_t_required=0):
        """
        Construct a SparseCounts object from a ReadArray object.

        args:
        -----
        read_array: ReadArray input object containing alignment filter results.
        collapse_molecules (default: True): if True, counts molecules, else counts reads.
        n_poly_t_required (default: 0): The number of T nucleotides that must be observed
          on the 3' of the forward read. If a poly-T tail is expected, the lack of this
          can cause the primer to randomly prime, producing unexpected results. Setting
          this value to 0 will count all observations.

        returns:
        --------
        SparseCounts object
        """

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

    def cells_above_threshold_molecules(self, t):
        """return the number of cells with > t observations"""
        return np.sum(self.counts.sum(axis=1) > t)

    def cells_coexpressing(self, markers):
        """
        returns the number of cells coexpressing all markers

        can be passed either a scid or gene name, but will return 0 if you give the
        wrong encoding, so pay attention!
        """
        csc = self._counts.tocsc()
        if isinstance(markers, Iterable) and not isinstance(markers, str):
            mask = np.ones_like(self.index, dtype=np.bool)
            for m in markers:
                mask &= np.ravel(csc[:, self.columns == m].sum(axis=1) > 1)
        else:
            mask = csc[:, self.columns == markers].sum(axis=1) > 1
        return np.sum(mask)

    def cells_expressing(self, markers):
        """
        returns the number of cells expressing any marker in markers

        can be passed either a scid or gene name, but will return 0 if you give the
        wrong encoding, so pay attention!
        """
        csc = self._counts.tocsc()
        if isinstance(markers, Iterable) and not isinstance(markers, str):
            mask = np.zeros_like(self.index, dtype=np.bool)
            for m in markers:
                mask |= np.ravel(csc[:, self.columns == m].sum(axis=1) > 1)
        else:
            mask = csc[:, self.columns == markers].sum(axis=1) > 1
        return np.sum(mask)

    def threshold_cells(self, threshold):
        """
        return a new SparseCounts object with only cells that pass filters and features
        with non-zero observations, given the subset of cells"""
        cell_sums = np.ravel(self.counts.sum(axis=1) > threshold)
        csr = self._counts.tocsr()
        sub_cells = csr[cell_sums, :]
        csc = sub_cells.tocsc()
        gene_nonzero = np.ravel(csc.sum(axis=0) > 0)
        sub_counts = csc[:, gene_nonzero].tocoo()
        index = self._index[cell_sums]
        columns = self._columns[gene_nonzero]

        return SparseCounts(sub_counts, index, columns)

    def to_dense(self, threshold=None):
        """
        return a DenseCounts object with rows = cells and columns = scids or genes if
        convert_ids is True. for now, use threshold_cells to cut off "non-cell" barcodes

        #todo Improve Implementation once we have a better sense of how to threshold cells
        """

        # mask cells with fewer than threshold counts
        if threshold:
            subset = self.threshold_cells(threshold)
        else:
            subset = self

        # return a DenseCounts object
        dense_counts = np.asarray(subset.counts.todense())
        df = pd.DataFrame(dense_counts, subset.index, subset.columns)
        return DenseCounts(df)

    def convert_gene_ids(self, converter):
        """
        convert integer gene ids into gene names for downstream analysis

        args:
        -----
        converter: a seqc.convert_features.ConvertGeneCoordinates object. This is
         typically saved as output_stem + '_gene_id_map.p' during a SEQC run.

        returns:
        --------
        None

        modifies:
        ---------
        self.columns (int ids -> str ids)
        """
        ids = [converter.int2str_id(id_) for id_ in self.columns]
        self._columns = np.array(ids)

    def _convert_ids_deprecated(self, fgtf, scid_map):
        if not any([fgtf, scid_map]):
            raise ValueError('one of gtf or scid_map must be passed for conversion to'
                             'occur')

        # load the gene map (gmap)
        if not scid_map:
            r = seqc.gtf.Reader(fgtf)
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
        new_ids = np.zeros(len(self._columns), dtype=object)
        for i, val in enumerate(self._columns):
            try:
                new_ids[i] = gmap[val]
            except KeyError:
                new_ids[i] = None
        self._columns = np.array(new_ids)

    def plot_observations_per_cell(self, smooth=False, log=False, xlabel='', ylabel='',
                                   title=''):
        """
        Plot either a kde smoothed distribution or a histogram of the number of
        observations (either molecules or reads, depending on the construction of the
        SparseCounts object).

        returns:
        --------
        f: matplotlib.Figure object
        ax: matplotlib.Axes object
        """

        cellsums = self.counts.sum(axis=1)
        if log:
            cellsums = np.log(cellsums)
        f, ax = seqc.plot.get_axes()
        if not xlabel:
            xlabel = 'observations per cell'
        if not ylabel:
            ylabel = 'number of cells'
        if not title:
            title = 'distribution of observations per cell'
        if smooth:
            seqc.plot.kde(cellsums, fig=f, ax=ax, xlabel=xlabel, ylabel=ylabel,
                          title=title)
        else:
            seqc.plot.histogram(cellsums, fig=f, ax=ax, xlabel=xlabel, ylabel=ylabel,
                                title=title)
        return f, ax

    def plot_gc_bias(self, fig=None, ax=None, molecules=False, reads=False):
        """
        Scatter plot of GC content against library size. Pass either molecules or
        reads=True to correctly label the y-axis.
        """
        cell_sums = np.ravel(self.counts.sum(axis=1))
        return _plot_cell_gc_content_bias(cell_sums, self.index, fig=fig, ax=ax,
                                          molecules=molecules, reads=reads)

    def get_mitochondrial_fraction(self, return_sums=True):
        """return the fraction of reads that are mitochondrial"""
        cell_sums = self.counts.sum(axis=1)
        csc = self.counts.tocsc()
        mt_genes = np.array([True if g.startswith('MT-') else False
                             for g in self.columns],
                            dtype=np.bool)
        mt_sums = csc[:, mt_genes].tocsr().sum(axis=1)
        if np.all(mt_sums == 0):
            raise ValueError('No mitochondrial molecules detected, fractions are '
                             'undefined')
        mt_fraction = mt_sums / cell_sums
        if return_sums:
            return mt_fraction, cell_sums
        else:
            return mt_fraction

    def plot_fraction_mitochondrial(
            self, fig=None, ax=None, molecules=False, reads=False):

        # check that the ids are strings
        if not isinstance(self.columns[0], str):
            raise TypeError('Please convert gene ids to string identifiers using '
                            'self.convert_gene_ids() before calling this function.')
        cell_sums = self.counts.sum(axis=1)
        csc = self.counts.tocsc()
        mt_genes = np.array([True if g.startswith('MT-') else False
                             for g in self.columns],
                            dtype=np.bool)
        mt_sums = csc[:, mt_genes].tocsr().sum(axis=1)
        if np.all(mt_sums == 0):
            raise ValueError('No mitochondrial molecules detected, cannot plot.')

        f, ax = _plot_fraction_mitochondrial_rna(
            cell_sums, mt_sums, fig=fig, ax=ax, molecules=molecules, reads=reads)
        return f, ax

    def to_npz(self, npz_file):
        """
        Save an .npz archive containing all information necessary to reconstruct the
        SparseCounts object
        """
        data, row, col = self.counts.data, self.counts.row, self.counts.col
        np.savez(npz_file, data=data, row=row, col=col, index=self.index,
                 columns=self.columns)

    @classmethod
    def from_npz(cls, npz_file):
        """
        Load a SparseCounts object from a .npz file.
        """
        npz_data = np.load(npz_file)
        n = npz_data['index'].shape[0]
        m = npz_data['columns'].shape[0]
        data = coo_matrix((npz_data['data'], (npz_data['row'], npz_data['col'])),
                          shape=(n, m), dtype=npz_data['data'].dtype)
        index = npz_data['index']
        columns = npz_data['columns']
        return cls(data, index, columns)


class DenseCounts:
    """
    DenseCounts is a light wrapper for a pandas DataFrame that adds some relevant
    plotting and analysis methods

    """

    def __init__(self, df):
        """
        args:
        -----
        df: a Pandas DataFrame object whose rows are cells and columns are features
          (genes)
        """

        # todo preprocess the dataframe indices to enforce unique, sorted indices

        self._df = df
        self._tsne = None  # internal property to track if tsne has been run before

    def __repr__(self):
        return "<DenseCount wrapper of DataFrame with shape (%d, %d)" % self.df.shape

    @property
    def df(self):
        return self._df

    # For now, want the user to be able to manipulate the underlying DataFrame
    @df.setter
    def df(self, df):
        self._df = df

    @staticmethod
    def from_sparse(sparse_counts, threshold=None, fgtf=None, scid_map=None):
        """Convert a SparseCounts object into a DenseCounts object"""
        return sparse_counts.to_dense(threshold=threshold, fgtf=fgtf, scid_map=scid_map)

    def plot_gc_bias(self, fig=None, ax=None, molecules=False, reads=False):
        """
        Scatter plot of GC content against library size. Pass either molecules or
        reads=True to correctly label the y-axis.
        """
        cell_counts = self.df.sum(axis=1)
        return _plot_cell_gc_content_bias(cell_counts, self.df.index, fig=fig, ax=ax,
                                          molecules=molecules, reads=reads)

    def plot_tsne(self, markers=None, pca_d=None, new=False, fig=None, ax=None):
        """Generate a tsne plot of self.

        by default, plot_tsne() will color cells by density. If an optional marker or set
        of markers is passed, it will color the plot by the summed expression value of
        all markers.

        the tSNE coordinates will only be generated (1) if plot_tsne() has never been
        called before for this object, or (2) if new=True.

        tsne will automatically use PCA to reduce the dimensionality to 30. However,
        it is recommended that you plot the PCA results and visually inspect the
        explained variance to get the right number of components before setting this
        parameter.

        returns:
        --------
        f, ax: matplotlib Figure and Axes objects
        """

        # only calculate new tsne coordinates if requested
        if self._tsne is None or new:
            # sanitize data of nan/inf to avoid segfault
            sanitized = self.df.values.copy()
            sanitized[np.isnan(sanitized) | np.isinf(sanitized)] = 0
            sanitized = sanitized.astype(float)
            self._tsne = bh_sne(sanitized, pca_d=pca_d)

        # set title
        if isinstance(markers, Iterable) and not isinstance(markers, str):
            title = ('tSNE projection, colored by summed expression of:\n%s'
                     % ', '.join(markers))
        elif markers:
            title = 'tSNE projection, colored by %s' % markers
        else:
            title = 'tSNE projection, colored by cell density'
        xlabel = 'tSNE dimension 0'
        ylabel = 'tSNE dimension 1'

        # if no markers passed, seqc.plot density scatter
        if markers:
            x, y = self._tsne[:, 0], self._tsne[:, 1]
            sum_expression = self.df[markers].sum(axis=1)
            seqc.plot.scatter_colored_by_data(
                x, y, sum_expression, fig=fig, ax=ax, xlabel=xlabel, ylabel=ylabel,
                title=title)
        else:
            x, y = self._tsne[:, 0], self._tsne[:, 1]
            seqc.plot.scatter_density(
                x, y, fig=fig, ax=ax, xlabel=xlabel, ylabel=ylabel, title=title)
        return fig, ax

    def plot_fraction_mitochondrial(
            self, fig=None, ax=None, molecules=False, reads=False):
        cell_sums = self.df.sum(axis=1)
        mt_genes = np.array([True if g.startswith('MT-') else False
                             for g in self.df.columns],
                            dtype=np.bool)
        mt_sums = self.df[mt_genes].sum(axis=1)
        f, ax = _plot_fraction_mitochondrial_rna(
            cell_sums, mt_sums, fig=fig, ax=ax, molecules=molecules, reads=reads)
        return f, ax

    def to_npz(self, npz_file):
        """
        Save an .npz archive containing all information necessary to reconstruct the
        DenseCounts object
        """
        columns = np.asarray(self.df.columns)
        index = np.asarray(self.df.index)
        data = self.df.values
        np.savez(npz_file, data=data, index=index, columns=columns)

    @classmethod
    def from_npz(cls, npz_file):
        """load a DenseCounts object from a .npz archive"""
        npz_data = np.load(npz_file)
        return cls(npz_data['data'], index=npz_data['index'], columns=npz_data['columns'])

    def convert_gene_ids(self, converter):
        """
        Convert scids to gene identifiers either by parsing the gtf file (slow) or by
        directly mapping scids to genes (fast). In the latter case, the gtf_map must be
        pre-processed and saved to the index folder.

        see seqc.gtf.Reader.scid_to_gene(gtf, save=<output_file>) to save this gtf_map for
        repeat-use.
        """
        ids = [converter.int2str_id(id_) for id_ in self.df.columns]
        self.df.columns = ids

    def downsample(self, percentile: int=0):
        """downsample to equalize observations per row (cell)

        args:
        -----
        percentile: which percentile to downsample to (default 0). In some cases, it may
         be advantageous to set this higher (10-25) in order to reduce data loss from
         extremely low-value cells

        returns:
        seqc.analyze.DenseCounts object with downsampled observations
        """
        cell_sums = self.df.sum(axis=1)
        threshold = np.percentile(cell_sums, percentile)
        multinomial_probabilities = self.df.div(cell_sums, axis=0)
        downsampled = []
        for gene_probability_vector in multinomial_probabilities.values:
            downsampled.append(np.random.multinomial(threshold, gene_probability_vector))
        downsampled = pd.DataFrame(np.vstack(downsampled), index=self.df.index,
                                   columns=self.df.columns)
        return DenseCounts(downsampled)


class Experiment:
    """
    Some plotting methods require both molecule and cell counts; this class holds
    these values and implements related statistics and methods
    """

    def __init__(self, reads, molecules):
        for var in [reads, molecules]:
            if not any(isinstance(var, t) for t in [SparseCounts, DenseCounts]):
                raise TypeError('reads must be a SparseCounts or DenseCounts object')
        self._reads = reads
        self._molecules = molecules

    def __repr__(self):
        return ("<Experiment object:\nReads:%s\nMolecules:%s>" %
                (repr(self._reads), repr(self._molecules)))

    @property
    def reads(self):
        return self._reads

    @property
    def molecules(self):
        return self._molecules

    @staticmethod
    def from_read_array(ra, n, s):
        """create an experiment from a ReadArray object.

        args:
        -----
        ra: ReadArray object
        n: the number of T nulceotides that must follow the cell barcodes for the
         alignment to be considered valid
        s: the number of times a molecule must be observed to be included in the final
         read or molecule matrices

        """
        ua = ra.to_unique(n)
        return ua.to_experiment(s)

    @staticmethod
    def from_unique_read_array(ua):
        return ua.to_experiment()

    def to_npz(self, npz_file):
        """
        Save an .npz archive containing all information necessary to reconstruct the
        SparseCounts object
        """
        if isinstance(self.reads.index, np.ma.MaskedArray):
            self.reads._index = np.ma.compressed(self.reads.index)
        if isinstance(self.molecules.index, np.ma.MaskedArray):
            self.molecules._index = np.ma.compressed(self.molecules.index)

        data = {
            'rdata': self.reads.counts.data,
            'rrow': self.reads.counts.row,
            'rcol': self.reads.counts.col,
            'r_columns': self.reads.columns,
            'r_index': self.reads.index,
            'mdata': self.molecules.counts.data,
            'mrow': self.molecules.counts.row,
            'mcol': self.molecules.counts.col,
            'm_columns': self.molecules.columns,
            'm_index': self.molecules.index}
        np.savez(npz_file, **data)

    @classmethod
    def from_npz(cls, npz_file):
        """
        Load a SparseCounts object from a .npz file.
        """
        npz_data = np.load(npz_file)
        n = npz_data['r_index'].shape[0]
        m = npz_data['r_columns'].shape[0]
        data = coo_matrix((npz_data['rdata'], (npz_data['rrow'], npz_data['rcol'])),
                          shape=(n, m), dtype=npz_data['rdata'].dtype)
        index = npz_data['r_index']
        columns = npz_data['r_columns']
        reads = SparseCounts(data, index, columns)

        n = npz_data['m_index'].shape[0]
        m = npz_data['m_columns'].shape[0]
        data = coo_matrix((npz_data['mdata'], (npz_data['mrow'], npz_data['mcol'])),
                          shape=(n, m), dtype=npz_data['mdata'].dtype)
        index = npz_data['m_index']
        columns = npz_data['m_columns']
        molecules = SparseCounts(data, index, columns)
        return cls(reads, molecules)

    # todo add support for merging duplicate ids
    def convert_gene_ids(self, converter: str) -> None:
        """
        convert integer gene identifiers into string identifiers for downstream analysis

        args:
        -----
        converter: a filepath pointing to a serialized
         seqc.convert_features.ConvertGeneCoordinates object. or a loaded
         ConvertGeneCoordinates object. This object is usually saved as
         output_stem + '_gene_id_map.p' during a SEQC run.

        returns:
        --------
        None

        modifies:
        ---------
        self.reads.columns (int ids -> string ids)
        self.molecules.columns (int ids -> string ids)
        """
        if isinstance(converter, str):
            gmap = seqc.convert_features.ConvertGeneCoordinates.from_pickle(converter)
        elif isinstance(converter, seqc.convert_features.ConvertGeneCoordinates):
            gmap = converter
        else:
            raise TypeError('converter must either be a loaded ConvertGeneCoordinates '
                            'object or the filename of a serialized version')
        self._molecules.convert_gene_ids(gmap)
        self._reads.convert_gene_ids(gmap)

    def equalize_features(self):
        """
        Ensures columns of self.reads and self.molecules are equal

        Error correction removes any molecule the support < s (default 2). This can
        cause some features to drop out of the molecule matrix while still being present
        in the unfiltered reads. Many of the comparisons in this class are done across
        features, and therefore those features must first be equalized.

        returns:
        --------
        Experiment object
        """
        # get the read features that are in molecules; all molecules are in features
        if np.array_equal(self.reads.columns, self.molecules.columns):
            return self

        in_mols = np.in1d(self.reads.columns, self.molecules.columns)

        # sort both, index columns
        read_order = np.argsort(self.reads.columns)[in_mols]
        mol_order = np.argsort(self.molecules.columns)

        read_columns = self.reads.columns[read_order]
        mol_columns = self.molecules.columns[mol_order]

        read_data = self.reads.counts.tocsc()[:, read_order].tocoo()
        mol_data = self.molecules.counts.tocsc()[:, mol_order].tocoo()

        reads = SparseCounts(read_data, self.reads.index, read_columns)
        mols = SparseCounts(mol_data, self.molecules.index, mol_columns)
        return Experiment(reads, mols)

    def equalize_cells(self):
        """
        Ensures indices of self.reads and self.molecules are equal by eliminating any cell
        that does not appear in the molecule matrix. This happens when cells consist
        entirely of low-count (likely erroneous) reads
        """
        if np.array_equal(self.reads.index, self.molecules.index):
            return self

        in_mols = np.in1d(self.reads.index, self.molecules.index)

        # sort both, index columns
        read_order = np.argsort(self.reads.index)[in_mols]
        mol_order = np.argsort(self.molecules.index)

        read_index = self.reads.index[read_order]
        mol_index = self.molecules.index[mol_order]

        read_data = self.reads.counts.tocsr()[read_order, :].tocoo()
        mol_data = self.molecules.counts.tocsr()[mol_order, :].tocoo()

        reads = SparseCounts(read_data, read_index, self.reads.columns)
        mols = SparseCounts(mol_data, mol_index, self.molecules.columns)
        return Experiment(reads, mols)

    def plot_umi_correction(self, log=False):
        """
        Displays the impact of UMI correction on raw or log-transformed sc-RNA-seq data
        """

        # create figure
        gs = plt.GridSpec(nrows=2, ncols=2)
        fig = plt.figure(figsize=(6, 6))

        # plot # 1: Regression of Reads against Molecules
        # -----
        # get mean gene expression values for molecules and reads
        read_means = np.ravel(self.reads.mean(axis=0))
        mol_means = np.ravel(self.molecules.mean(axis=0))

        if log:
            read_means = np.log(read_means)
            mol_means = np.log(mol_means)

        # carry out polynomial regression
        x, y = mol_means, read_means  # used for reconstructing y_hat
        model = Pipeline([('poly', PolynomialFeatures(degree=3)),
                          ('linear', LinearRegression(fit_intercept=False))])
        model = model.fit(mol_means[:, np.newaxis], read_means)
        beta = model.named_steps['linear'].coef_

        # reconstruct y_hat and calculate the fit line.
        y_hat = beta[0] + beta[1] * x + beta[2] * x ** 2 + beta[3] * x ** 3
        xmin, xmax = np.min(mol_means), np.max(mol_means)
        f = np.linspace(xmin, xmax, 1000)
        fline = beta[0] + beta[1] * f + beta[2] * f ** 2 + beta[3] * f ** 3

        sx, sy, colors = seqc.plot.density_coloring(x, y)
        ax1 = plt.subplot(gs[0, 0])
        xlabel = 'Molecules  / Gene'
        ylabel = 'Reads / Gene'
        title = 'Polynomial Regression'
        seqc.plot.scatter_colored_by_data(sx, sy, colors, ax=ax1, xlabel=xlabel,
                                          ylabel=ylabel, title=title)

        # plot # 2: Regression Residuals
        # -----
        ax2 = plt.subplot(gs[0, 1])
        residuals = y - y_hat
        sx, s_residuals, colors = seqc.plot.density_coloring(x, residuals)
        ylabel = 'Residual'
        title = 'Residuals'
        seqc.plot.scatter_colored_by_data(sx, s_residuals, colors, ax=ax2, xlabel=xlabel,
                                          ylabel=ylabel, title=title)

        # plot # 3: Scaled Residuals
        # -----
        ax3 = plt.subplot(gs[1, 0])
        scaled_residuals = np.abs(y - y_hat) / y
        sx, s_scaled_residuals, colors = seqc.plot.density_coloring(x, scaled_residuals)
        ylabel = 'Residual / Reads per Gene'
        title = 'Normalized Residuals'
        seqc.plot.scatter_colored_by_data(sx, s_scaled_residuals, colors, ax=ax3,
                                          xlabel=xlabel, ylabel=ylabel, title=title)
        # adjust axes
        cmax = np.max(sx)
        ymax = np.max(scaled_residuals)
        plt.xlim((0, cmax))

        plt.ylim((0, ymax + min(1, ymax + 0.05)))

        # plot # 4: Fano Factors
        # -----
        ax4 = plt.subplot(gs[1, 1])
        idx = np.where(np.sum(self.molecules > 0, axis=0) > 2)
        read_fano = (np.var(self.reads[:, idx], axis=0) /
                     np.mean(self.reads[:, idx], axis=0))
        mols_fano = (np.var(self.molecules[:, idx], axis=0) /
                     np.mean(self.molecules[:, idx], axis=0))

        x, y, colors = seqc.plot.density_coloring(read_fano, mols_fano)
        xlabel = r'Fano ($\frac{\sigma^2}{\mu}$) Molecules / Gene'
        ylabel = r'Fano ($\frac{\sigma^2}{\mu}$) Reads / Gene'
        title = 'Normalied Variance of Reads \nvs. Molecules per Gene'
        seqc.plot.scatter_colored_by_data(x, y, colors, ax=ax4, xlabel=xlabel,
                                          ylabel=ylabel, title=title)
        # plot x = y
        cmin, cmax = np.min(np.hstack([x, y])), np.max(np.hstack([x, y]))
        x_y = np.linspace(cmin, cmax, 1000)
        ax4.plot(x_y, x_y, 'r--')
        plt.xlim((cmin, cmax))
        plt.ylim((cmin, cmax))

        seqc.plot.clean_figure(fig)

        return fig, dict(regression=ax1, residual=ax2, scaled_residual=ax3, fano=ax4)

    def plot_coverage_vs_rpm_ratio(self, min_ratio=2, min_mols=10):
        """plot cell selection"""

        # calculate cell sums
        mol_counts = self.molecules.counts.sum(axis=1)
        read_counts = self.reads.counts.sum(axis=1)

        # convert to series
        mol_counts = pd.Series(mol_counts.T.tolist()[0], index=self.molecules.index)
        read_counts = pd.Series(read_counts.T.tolist()[0], index=self.reads.index)

        # get ratios
        ratios = read_counts / mol_counts

        # filter out junk
        barcodes = mol_counts.index[(mol_counts > min_mols) & (ratios > min_ratio)]

        # set up plot
        fig = plt.figure(figsize=(8, 4))
        gs = plt.GridSpec(nrows=1, ncols=2)

        ax1 = plt.subplot(gs[0, 0])
        xlabel = 'log Molecule count per cell'
        ylabel = 'Average Reads per Molecule'
        title = 'Reads per Molecule vs. Library Size'
        fig, ax1 = seqc.plot.scatter_density(np.log(mol_counts[barcodes]), ratios[barcodes],
                                        fig=fig, ax=ax1, xlabel=xlabel, ylabel=ylabel,
                                        title=title)
        plt.ylim([0, 200])

        # fit GMM
        df = pd.DataFrame([np.log(mol_counts[barcodes]), ratios[barcodes]]).T
        n_clusters = 3
        gmm = GMM(n_components=n_clusters)
        gmm.fit(df)
        clusters = gmm.predict(df)

        # Plot clusters
        ax2 = plt.subplot(gs[0, 1])
        colors = seqc.plot.qualitative_colors[:3]
        # for i, c in enumerate(colors):
        #     ax2.scatter(df.ix[clusters == i, 0], df.ix[clusters == i, 1],
        #                 color=c, s=, edgecolor='none')
        xlabel = 'log Molecule count per cell'
        ylabel = 'Average Reads per Molecule'
        title = 'Cluster Classification'
        seqc.plot.scatter_colored_by_data(df.ix[:, 0], df.ix[:, 1], clusters,
                                          xlabel=xlabel, ylabel=ylabel, title=title)
        plt.ylim([0, 200])
        return fig, (ax1, ax2)

    def _select(self, mask, axis=1):
        """return a new array containing only objects passing mask"""
        mols = self.molecules._select(mask, axis=1)
        reads = self.reads._select(mask, axis=1)
        return Experiment(reads, mols)

    def select_cells(self, min_ratio=2, min_mols=10):

        # calculate cell sums
        mol_counts = np.ravel(self.molecules.sum(axis=1))
        read_counts = np.ravel(self.reads.sum(axis=1))

        # get ratios
        ratios = read_counts / mol_counts

        # filter out junk
        mask = (mol_counts > min_mols) & (ratios > min_ratio)

        # fit GMM
        train = np.vstack([np.log(mol_counts[mask]), ratios[mask]]).T
        n_clusters = 3
        gmm = GMM(n_components=n_clusters)
        gmm.fit(train)
        clusters = gmm.predict(train)

        # Choose the "real" cell cluster
        cell_cluster = np.where(gmm.means_[:, 0] == max(gmm.means_[:, 0]))[0][0]
        indices = np.where(mask)
        cell_indices = indices[clusters == cell_cluster]

        # get a mask for the original counts arrays
        counts_mask = np.zeros_like(mol_counts, dtype=np.bool)
        counts_mask[cell_indices] = 1

        # now, shrink the original array!
        return self._select(mask, axis=1)

    def plot_gc_content_vs_rpm_ratio(self, min_ratio=2, min_mols=10):

        # calculate cell sums
        mol_counts = np.ravel(self.molecules.counts.sum(axis=1))
        read_counts = np.ravel(self.reads.counts.sum(axis=1))

        # get ratios
        ratios = read_counts / mol_counts

        # filter out junk
        mask = (mol_counts > min_mols) & (ratios > min_ratio)

        cells = self.molecules.index[mask]
        gc_content = np.array([seqc.encodings.ThreeBit.gc_content(i) for i in cells])
        ratios = ratios[mask]

        xlabel = 'cell GC content'
        ylabel = 'Average Reads per Molecule'
        title = 'GC Content vs Cell Coverage'
        f, ax = seqc.plot.scatter_density(gc_content, ratios, xlabel=xlabel,
                                          ylabel=ylabel, title=title)
        plt.ylim(0, 300)
        return f, ax

    def summary(self, alignment_summary, fout=None):
        """print a basic summary of the run for technological debugging purposes

        args:
        -----
        alignment_summary:  name of alignment summary file. Found in output directory
         of SEQC run.
        fout: name of the file in which to save the summary.
        """
        read_counts = self.reads.counts.sum(axis=1)
        mol_counts = self.molecules.counts.sum(axis=1)
        stats = pd.Series()

        def get_val(frame, val):
            idx = np.array([True if val in l else False for l in frame.index])
            return frame.ix[idx][0]

        # parse alignment
        lines = pd.DataFrame.from_csv(alignment_summary, sep='\t').iloc[:, 0]
        stats['Reverse Read Length'] = get_val(lines, 'Average input read length')
        stats['No of reads'] = get_val(lines, 'Number of input reads')
        stats['Uniquely mapping reads'] = '%s (%s)' % (
            get_val(lines, 'Uniquely mapped reads number'),
            get_val(lines, 'Uniquely mapped reads %'))
        stats['Unmapped reads (includes phiX)'] = get_val(
            lines, '% of reads unmapped: too short')

        # Molecule counts
        stats['Total molecules'] = '%d' % mol_counts.sum()
        stats['Reads contributing to molecules'] = '%d (%.2f%%)' % (
            read_counts.sum(), read_counts.sum() / int(stats['No of reads']) * 100)
        stats['Reads per molecule'] = '%d' % (read_counts.sum() / mol_counts.sum())

        # Cell counts
        stats['Cells: > 500 mols.'] = str(int(sum(mol_counts > 500)[0]))
        stats['Cells: > 1k mols.'] = str(int(sum(mol_counts > 1000)[0]))
        stats['Cells: > 5k mols.'] = str(int(sum(mol_counts > 5000)[0]))
        stats['Cells: > 10k mols.'] = str(int(sum(mol_counts > 10000)[0]))

        if fout:
            with open(fout, 'w') as f:
                f.write(repr(stats))

        return stats


class CompareExperiments:
    """
    Some plotting methods require experiments themselves to be compareed; this class holds
    these values and implements related statistics and methods
    """
    pass
