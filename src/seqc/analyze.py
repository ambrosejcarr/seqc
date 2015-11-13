__author__ = 'ambrose'


import numpy as np
import pandas as pd
from collections import defaultdict, Mapping, Iterable
from scipy.sparse import coo_matrix
from seqc import plot, gtf, three_bit
from tsne import bh_sne
from sklearn.preprocessing import PolynomialFeatures
from sklearn.linear_model import LinearRegression
from sklearn.pipeline import Pipeline
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
    gc_content = np.array([three_bit.ThreeBit.gc_content(s) for s in sequences])
    fig, ax = plot.scatter_density(gc_content, cell_counts, fig=fig, ax=ax, xlabel=xlabel,
                                   ylabel=ylabel, title=title)
    return fig, ax


def _plot_fraction_mitochondrial_rna(mt_counts, cell_counts, fig=None, ax=None,
                                     molecules=False, reads=False):
    raise NotImplementedError


class SparseCounts:
    """
    SparseCounts is a thin wrapper around a scipy.sparse.coo_matrix object that provides
    additional functionality for the analysis of sc-RNA-seq data
    """

    def __init__(self, sparse_counts, index, columns):
        self._counts = sparse_counts
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

    @classmethod
    def from_read_array(cls, read_array, collapse_molecules=True, n_poly_t_required=0):
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

    def to_dense(self, threshold=None, fgtf=None, scid_map=None):
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

        # convert ids
        if any([fgtf, scid_map]):
            subset.convert_ids(fgtf=fgtf, scid_map=scid_map)

        # return a DenseCounts object
        dense_counts = np.asarray(subset.counts.todense())
        df = pd.DataFrame(dense_counts, subset.index, subset.columns)
        return DenseCounts(df)

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

    def plot_gc_bias(self, fig=None, ax=None, molecules=False, reads=False):
        """
        Scatter plot of GC content against library size. Pass either molecules or
        reads=True to correctly label the y-axis.
        """
        cell_sums = np.ravel(self.counts.sum(axis=1))
        return _plot_cell_gc_content_bias(cell_sums, self.index, fig=fig, ax=ax,
                                          molecules=molecules, reads=reads)

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

    def plot_tsne(self, markers=None, new=False, fig=None, ax=None):
        """Generate a tsne plot of self.

        by default, plot_tsne() will color cells by density. If an optional marker or set
        of markers is passed, it will color the plot by the summed expression value of
        all markers.

        the tSNE coordinates will only be generated (1) if plot_tsne() has never been
        called before for this object, or (2) if new=True.

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
            self._tsne = bh_sne(sanitized)

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

        # if no markers passed, plot density scatter
        if markers:
            x, y = self._tsne[:, 0], self._tsne[:, 1]
            sum_expression = self.df[markers].sum(axis=1)
            plot.scatter_colored_by_data(
                x, y, sum_expression, fig=fig, ax=ax, xlabel=xlabel, ylabel=ylabel,
                title=title)
        else:
            x, y = self._tsne[:, 0], self._tsne[:, 1]
            plot.scatter_density(
                x, y, fig=fig, ax=ax, xlabel=xlabel, ylabel=ylabel, title=title)
        return fig, ax

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

        sx, sy, colors = plot.density_coloring(x, y)
        ax1 = plt.subplot(gs[0, 0])
        xlabel = 'Molecules  / Gene'
        ylabel = 'Reads / Gene'
        title = 'Polynomial Regression'
        plot.scatter_colored_by_data(sx, sy, colors, ax=ax1, xlabel=xlabel, ylabel=ylabel,
                                     title=title)

        # plot # 2: Regression Residuals
        # -----
        ax2 = plt.subplot(gs[0, 1])
        residuals = y - y_hat
        sx, s_residuals, colors = plot.density_coloring(x, residuals)
        ylabel = 'Residual'
        title = 'Residuals'
        plot.scatter_colored_by_data(sx, s_residuals, colors, ax=ax2, xlabel=xlabel,
                                     ylabel=ylabel, title=title)

        # plot # 3: Scaled Residuals
        # -----
        ax3 = plt.subplot(gs[1, 0])
        scaled_residuals = np.abs(y - y_hat) / y
        sx, s_scaled_residuals, colors = plot.density_coloring(x, scaled_residuals)
        ylabel = 'Residual / Reads per Gene'
        title = 'Normalized Residuals'
        plot.scatter_colored_by_data(sx, s_scaled_residuals, colors, ax=ax3,
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

        x, y, colors = plot.density_coloring(read_fano, mols_fano)
        xlabel = r'Fano ($\frac{\sigma^2}{\mu}$) Molecules / Gene'
        ylabel = r'Fano ($\frac{\sigma^2}{\mu}$) Reads / Gene'
        title = 'Normalied Variance of Reads \nvs. Molecules per Gene'
        plot.scatter_colored_by_data(x, y, colors, ax=ax4, xlabel=xlabel, ylabel=ylabel,
                                     title=title)
        # plot x = y
        cmin, cmax = np.min(np.hstack([x, y])), np.max(np.hstack([x, y]))
        x_y = np.linspace(cmin, cmax, 1000)
        ax4.plot(x_y, x_y, 'r--')
        plt.xlim((cmin, cmax))
        plt.ylim((cmin, cmax))

        plot.clean_figure(fig)

        return fig, dict(regression=ax1, residual=ax2, scaled_residual=ax3, fano=ax4)


    def plot_coverage_vs_gc_content(self):
        raise NotImplementedError
        # TODO stopped here.


class CompareExperiments:
    """
    Some plotting methods require experiments themselves to be compareed; this class holds
    these values and implements related statistics and methods
    """


