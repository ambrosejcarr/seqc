import re
import os
import random
import pickle
import warnings
import multiprocessing
import shlex
import shutil
from itertools import combinations
from functools import partial
from copy import deepcopy
from collections import defaultdict
from subprocess import call, Popen, PIPE
import numpy as np
import pandas as pd
import tables as tb
from scipy.stats.mstats import kruskalwallis
from statsmodels.sandbox.stats.multicomp import multipletests
from scipy.sparse import coo_matrix
from scipy.stats import gaussian_kde, pearsonr, mannwhitneyu
import matplotlib
try:
    os.environ['DISPLAY']
except KeyError:
    matplotlib.use('Agg')
import matplotlib.pyplot as plt
with warnings.catch_warnings():
    warnings.simplefilter('ignore')  # catch experimental ipython widget warning
    import seaborn as sns
import phenograph
from tsne import bh_sne
import seqc


class SparseMatrixError(Exception):
    pass


# set plotting defaults
sns.set_style('ticks')
matplotlib.rc('font', **{'family': 'serif',
                         'serif': ['Computer Modern Roman'],
                         'monospace': ['Computer Modern Typewriter']
                         })

matplotlib.rc('figure', **{'figsize': (4, 4),
                           'dpi': 150})

matplotlib.rc('patch', **{'facecolor': 'royalblue',
                          'edgecolor': 'none'})

matplotlib.rc('lines', **{'color': 'royalblue',
                          'markersize': 7})

matplotlib.rc('savefig', **{'dpi': '150'})
cmap = matplotlib.cm.viridis
size = 8


def qualitative_colors(n):
    return sns.color_palette('husl', n)


def get_fig(fig=None, ax=None):
    """fills in any missing axis or figure with the currently active one
    :param ax: matplotlib Axis object
    :param fig: matplotlib Figure object
    """
    if not fig:
        fig = plt.gcf()
    if not ax:
        ax = plt.gca()
    return fig, ax


def density_2d(x, y):
    """return x and y and their density z, sorted by their density (smallest to largest)

    :param x:
    :param y:
    :return:
    """
    xy = np.vstack([np.ravel(x), np.ravel(y)])
    z = gaussian_kde(xy)(xy)
    i = np.argsort(z)
    return np.ravel(x)[i], np.ravel(y)[i], np.arcsinh(z[i])


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
        reader = seqc.alignment.sam.Reader(samfile)
        translator = seqc.sequence.gtf.GeneIntervals(gtf)

        data = np.recarray((len(reader),), cls._dtype)
        for i, alignment in enumerate(reader):
            gene = translator.translate(
                    alignment.rname, alignment.strand, alignment.pos)
            if gene is None:
                gene = 0
            pool = seqc.sequence.encodings.DNA3Bit.encode(alignment.pool)
            cell = seqc.sequence.encodings.DNA3Bit.encode(alignment.cell)
            rmt = seqc.sequence.encodings.DNA3Bit.encode(alignment.rmt)
            n_poly_t = alignment.poly_t.count(b'T')
            dust_score = alignment.dust_low_complexity_score
            data[i] = (pool, cell, rmt, n_poly_t, dust_score, gene, alignment.pos)

        return cls(data)

    def reads_passing_filters(self, min_poly_t: int, max_dust_score: int):
        """
        :param min_poly_t:
        :param max_dust_score:
        :return: ReadArray
        """
        data = self.data[((self.data['n_poly_t'] >= min_poly_t) &
                          (self.data['dust_score'] <= max_dust_score) &
                          (self.data['gene'] != 0) &
                          (self.data['cell'] != 0))]
        return ReadArray(data)

    def save(self, archive_name):
        """save a ReadArray in .h5 format

        :param archive_name:
        :return:
        """

        # create table
        blosc5 = tb.Filters(complevel=5, complib='blosc')
        f = tb.open_file(archive_name, mode='w', title='Data for seqc.ReadArray',
                         filters=blosc5)

        # store data
        f.create_table(f.root, 'data', self._data)
        f.close()

    @classmethod
    def load(cls, archive_name):
        f = tb.open_file(archive_name, mode='r')
        data = f.root.data.read()
        f.close()
        return cls(data)

    def __len__(self):
        return len(self.data)


class SparseFrame:

    def __init__(self, data, index, columns):

        if not isinstance(data, coo_matrix):
            raise TypeError('data must be type coo_matrix')
        if not isinstance(index, np.ndarray):
            raise TypeError('index must be type np.ndarray')
        if not isinstance(columns, np.ndarray):
            raise TypeError('columns must be type np.ndarray')

        self._data = data.astype(np.uint32)
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
        return self.data.sum(axis=axis)


class Experiment:

    def __init__(self, reads, molecules, metadata=None):
        """

        :param reads:  DataFrame or SparseFrame
        :param molecules: DataFrame or SparseFrame
        :param metadata: None or DataFrame
        :return:
        """
        if not (isinstance(reads, SparseFrame) or isinstance(reads, pd.DataFrame)):
            raise TypeError('reads must be of type SparseFrame or DataFrame')
        if not (isinstance(molecules, SparseFrame) or
                isinstance(molecules, pd.DataFrame)):
            raise TypeError('molecules must be of type SparseFrame or DataFrame')
        if metadata is None:
            metadata = pd.DataFrame(index=molecules.index, dtype='O')
        self._reads = reads
        self._molecules = molecules
        self._metadata = metadata
        self._normalized = False
        self._pca = None
        self._tsne = None
        self._diffusion_eigenvectors = None
        self._diffusion_eigenvalues = None
        self._diffusion_map_correlations = None
        self._cluster_assignments = None

    def save(self, fout: str) -> None:
        """
        :param fout: str, name of archive to store pickled Experiment data in. Should end
          in '.p'.
        :return: None
        """
        with open(fout, 'wb') as f:
            pickle.dump(vars(self), f)

    @classmethod
    def load(cls, fin):
        """

        :param fin: str, name of pickled archive containing Experiment data
        :return: Experiment
        """
        with open(fin, 'rb') as f:
            data = pickle.load(f)
        experiment = cls(data['_reads'], data['_molecules'], data['_metadata'])
        del data['_reads']
        del data['_molecules']
        del data['_metadata']
        for k, v in data.items():
            setattr(experiment, k[1:], v)
        return experiment

    def __repr__(self):
        c, g = self.molecules.shape
        _repr = ('Experiment: {c} cells x {g} genes\nsparse={s}'.format(
                s=self.is_sparse(), g=g, c=c))
        for k, v in sorted(vars(self).items()):
            if not (k == '_reads' or k == '_molecules}'):
                _repr += '\n{}={}'.format(k[1:], 'None' if v is None else 'True')
        return _repr

    @property
    def reads(self):
        return self._reads

    @reads.setter
    def reads(self, item):
        if not (isinstance(item, SparseFrame) or isinstance(item, pd.DataFrame)):
            raise TypeError('Experiment.reads must be of type SparseFrame or DataFrame')
        self._reads = item

    @property
    def molecules(self):
        return self._molecules

    @molecules.setter
    def molecules(self, item):
        if not (isinstance(item, SparseFrame) or isinstance(item, pd.DataFrame)):
            raise TypeError('Experiment.molecules must be of type SparseFrame or'
                            'DataFrame')
        self._molecules = item

    @property
    def metadata(self):
        return self._metadata

    @metadata.setter
    def metadata(self, item):
        if not isinstance(item, pd.DataFrame):
            raise TypeError('Experiment.metadata must be of type DataFrame')
        self._metadata = item

    @property
    def pca(self):
        return self._pca

    @pca.setter
    def pca(self, item):
        if not (isinstance(item, dict) or item is None):
            raise TypeError('self.pca must be a dictionary of pd.DataFrame object')
        self._pca = item

    @property
    def tsne(self):
        return self._tsne

    @tsne.setter
    def tsne(self, item):
        if not (isinstance(item, pd.DataFrame) or item is None):
            raise TypeError('self.tsne must be a pd.DataFrame object')
        self._tsne = item

    @property
    def diffusion_eigenvectors(self):
        return self._diffusion_eigenvectors

    @diffusion_eigenvectors.setter
    def diffusion_eigenvectors(self, item):
        if not (isinstance(item, pd.DataFrame) or item is None):
            raise TypeError('self.diffusion_eigenvectors must be a pd.DataFrame object')
        self._diffusion_eigenvectors = item

    @property
    def diffusion_eigenvalues(self):
        return self._diffusion_eigenvalues

    @diffusion_eigenvalues.setter
    def diffusion_eigenvalues(self, item):
        if not (isinstance(item, pd.DataFrame) or item is None):
            raise TypeError('self.diffusion_eigenvalues must be a pd.DataFrame object')
        self._diffusion_eigenvalues = item

    @property
    def diffusion_map_correlations(self):
        return self._diffusion_map_correlations

    @diffusion_map_correlations.setter
    def diffusion_map_correlations(self, item):
        if not (isinstance(item, pd.DataFrame) or item is None):
            raise TypeError('self.diffusion_map_correlations must be a pd.DataFrame'
                            'object')
        self._diffusion_map_correlations = item

    @property
    def cluster_assignments(self):
        return self._cluster_assignments

    @cluster_assignments.setter
    def cluster_assignments(self, item):
        if not (isinstance(item, pd.Series) or item is None):
            raise TypeError('self.cluster_assignments must be a pd.Series '
                            'object')
        self._cluster_assignments = item

    @classmethod
    def from_v012(cls, read_and_count_matrix: str):
        with open(read_and_count_matrix, 'rb') as f:
            d = pickle.load(f)
        reads = SparseFrame(d['reads']['matrix'], index=d['reads']['row_ids'],
                            columns=d['reads']['col_ids'])
        molecules = SparseFrame(d['molecules']['matrix'], index=d['molecules']['row_ids'],
                                columns=d['molecules']['col_ids'])
        return cls(reads, molecules)

    @classmethod
    def refresh(cls, experiment):
        """
        Generate a new Experiment object from an old instance of experiment. Convenience
        function for development, when class methods have changed.

        :param experiment:
        :return:
        """
        # todo not working because of copy vs. reference problems
        exp = cls(experiment.reads, experiment.molecules, experiment.metadata)
        properties = vars(experiment)
        del properties['_reads']
        del properties['_molecules']
        del properties['_metadata']
        for k, v in properties.items():
            setattr(exp, k[1:], v)

    def is_sparse(self):
        if all(isinstance(o, SparseFrame) for o in (self.reads, self.molecules)):
            return True
        else:
            return False

    def is_dense(self):
        return not self.is_sparse()

    # todo currently, this method has side-effects (mutates existing dataframes) -- change?
    # specifically, it mutates index (cell) ids
    @staticmethod
    def concatenate(experiments, metadata_labels=None):
        """
        Concatenate a set of Experiment objects. Each cell is considered to be unique,
        even if they share the same barcodes. If collisions are found, _%d will be added
        to each cell, where %d is equal to the position of the Experiment in the input
        experiments list.

        Concatenate should be run after filtering each dataset individually, since the
        datasets are considered to have distinct cells (even if they share the same
        barcodes)

        If metadata_labels are provided, a new entry in metadata will be included for
        each cell. This can be useful to mark each cell by which experiment it originated
        from in order to track batch effects, for example. metadata_labels should be a
        dictionary of lists, where the list has the same number of entries as the number
        of passed experiments. Thus, if one passed experiments=[e1, e2, e3] and
        metadata_labels={'batch': [1, 2, 3]}, batches 1, 2, and 3 would be propagated to
        each cell in the corresponding experiment

        :param experiments:
        :param metadata_labels: dict {label_name: [values], ...}
        :return:
        """

        if metadata_labels is None:
            metadata_labels = {}
        if not all(e.is_dense() for e in experiments):
            raise ValueError('merge may only be run on sparse inputs.')

        # mutate cell identifiers to ensure they are unique after merging
        for i, e in enumerate(experiments):
            suffix = '_{!s}'.format(i)
            new_cells = [str(c) + suffix for c in e.molecules.index]
            e.molecules.index = new_cells
            e.reads.index = new_cells
            e.metadata.index = new_cells

        # merge genes, create new cell list
        genes = set()
        cells = []
        meta_cols = set()
        for e in experiments:
            genes.update(list(e.molecules.columns))
            cells.extend(list(e.molecules.index))
            meta_cols.update(e.metadata.columns)

        # combine data
        empty_molc = np.zeros((len(cells), len(genes)), dtype=np.float64)
        empty_read = np.zeros((len(cells), len(genes)), dtype=np.float64)
        metadata = pd.DataFrame(index=cells, columns=meta_cols)
        m_combined = pd.DataFrame(empty_molc, index=cells, columns=genes)
        r_combined = pd.DataFrame(empty_read, index=cells, columns=genes)
        for e in experiments:
            m_combined.ix[e.molecules.index, e.molecules.columns] = e.molecules
            r_combined.ix[e.reads.index, e.reads.columns] = e.reads
            metadata.ix[e.metadata.index, e.metadata.columns] = e.metadata

        # add additional metadata
        for k, v in metadata_labels.items():
            metadata[k] = pd.Series(index=cells, dtype='O')
            for e in experiments:
                metadata[k].ix[e.metadata.index] = [v] * len(e.metadata.index)

        return Experiment(reads=r_combined, molecules=m_combined, metadata=metadata)

    @staticmethod
    def merge(experiments):
        """
        Merge a set of Experiment objects. Cells with duplicate cell barcodes will have
        their data summed across experiments.

        If metadata_labels are provided, a new entry in metadata will be included for
        each cell. This can be useful to mark each cell by which experiment it originated
        from in order to track batch effects, for example.

        :param experiments:
        :return:
        """

        if not all(e.is_sparse() for e in experiments):
            raise ValueError('merge may only be run on sparse inputs.')

        # get the intersection of all cell indices and gene columns
        cells = set()
        genes = set()
        for e in experiments:
            cells.update(e.molecules.index)
            genes.update(e.molecules.columns)

        # create new, global maps of genes and cells to indices
        cell_index = np.array(list(cells))
        gene_index = np.array(list(genes))
        cell2int = dict(zip(cell_index, range(len(cell_index))))
        gene2int = dict(zip(gene_index, range(len(gene_index))))

        # merge molecule counts
        new_molecule_data = defaultdict(int)
        new_read_data = defaultdict(int)
        for e in experiments:

            # map molecule and cell indices back to their respective cell and gene ids
            local_cell_map = dict(zip(range(len(e.molecules.index)), e.molecules.index))
            local_gene_map = dict(zip(range(len(e.molecules.columns)),
                                      e.molecules.columns))

            for mols, reads, cell, gene in zip(
                    e.molecules.data.data, e.reads.data.data, e.molecules.data.row,
                    e.molecules.data.col):

                # transform from index -> ids -> global index
                cell = cell2int[local_cell_map[cell]]
                gene = gene2int[local_gene_map[gene]]

                # add counts to global count matrices
                new_molecule_data[(cell, gene)] += mols
                new_read_data[(cell, gene)] += reads

        # get global row and col for coo_matrix construction
        row, col = (np.array(v) for v in zip(*new_molecule_data.keys()))

        # extract read and molecule data from dictionaries
        molecule_data = np.array(list(new_molecule_data.values()))
        read_data = np.array(list(new_read_data.values()))

        # get gene ids and cell ids for SparseFrame construction

        # make coo matrices
        shape = len(cell_index), len(gene_index)
        reads = coo_matrix((read_data, (row, col)), shape=shape, dtype=np.uint32)
        molecules = coo_matrix((molecule_data, (row, col)), shape=shape, dtype=np.uint32)

        sparse_reads = SparseFrame(reads, cell_index, gene_index)
        sparse_molecules = SparseFrame(molecules, cell_index, gene_index)

        return Experiment(reads=sparse_reads, molecules=sparse_molecules)

    def ensembl_gene_id_to_official_gene_symbol(self, gtf):
        """convert self.index containing scids into an index of gene names
        :param gtf:
        :return: Experiment
        """
        pattern = re.compile(
                r'(^.*?gene_id "[^0-9]*)([0-9]*)(\.?.*?gene_name ")(.*?)(".*?$)')

        gene_id_map = defaultdict(set)
        with open(gtf, 'r') as f:
            for line in f:
                match = re.match(pattern, line)
                if match:
                    gene_id_map[int(match.group(2))].add(match.group(4))

        reads = deepcopy(self.reads)
        reads.columns = ['-'.join(gene_id_map[i]) for i in self.reads.columns]
        molecules = deepcopy(self.molecules)
        molecules.columns = ['-'.join(gene_id_map[i]) for i in self.molecules.columns]

        return Experiment(reads=reads, molecules=molecules, metadata=self.metadata)

    def plot_molecules_vs_reads_per_molecule(self, fig=None, ax=None, min_molecules=10,
                                             ylim=(0, 150), title='Cell Coverage Plot'):
        """
        plots log10 molecules counts per barcode vs reads per molecule / molecules per
        barcode
        :param title: str, title for the plot (e.g. the sample name)
        :param ylim: tuple, indicates the min and max for the y-axis of the plot
        :param min_molecules: display only cells with this number of molecules or more
        :param ax: axis
        :param fig: figure

        :return: figure, axis
        """
        if self._normalized:
            raise RuntimeError('plot_molecules_vs_reads_per_molecule() should be run on '
                               'unnormalized data')

        fig, ax = get_fig(fig, ax)

        # get molecule and read counts per cell
        molecule_cell_sums = np.ravel(self.molecules.sum(axis=1))
        read_cell_sums = np.ravel(self.reads.sum(axis=1))

        # remove low expression reads
        mask = molecule_cell_sums >= min_molecules
        molecule_cell_sums = molecule_cell_sums[mask]
        read_cell_sums = read_cell_sums[mask]

        ratios = read_cell_sums / molecule_cell_sums
        x, y, z = density_2d(np.log10(molecule_cell_sums), ratios)
        ax.scatter(x, y, edgecolor='none', s=size, c=z, cmap=cmap)
        ax.set_xlabel('log10(Molecules per barcode)')
        ax.set_ylabel('Reads per barcode / Molecules per barcode')
        ax.set_title(title)
        ax.set_ylim(ylim)
        xlim = ax.get_xlim()
        ax.set_xlim((np.log10(min_molecules), xlim[1]))
        sns.despine(ax=ax)

        return fig, ax

    def remove_non_cell_barcodes(self, min_rpm=10, min_mpc=250):
        """removes low abundance cell barcodes

        remove barcodes with low molecules / barcode
        remove barcodes with low reads / molecule

        Note: defaults will often NOT work, must be set by visually inspecting
        plot_molecules_vs_reads_per_molecule()
        Note: changes data to dense_matrices

        :param min_mpc:
        :param min_rpm:
        :return: None
        """
        if self._normalized:
            raise RuntimeError('plot_molecules_vs_reads_per_molecule() should be run on '
                               'unnormalized data')

        # get molecule and read counts per cell
        molecule_cell_sums = pd.Series(np.ravel(self.molecules.sum(axis=1)),
                                       index=self.molecules.index, dtype=np.uint32)
        read_cell_sums = pd.Series(np.ravel(self.reads.sum(axis=1)),
                                   index=self.reads.index, dtype=np.uint32)

        ratios = read_cell_sums / molecule_cell_sums

        # select cells & convert to dense
        filter_ = (molecule_cell_sums >= min_mpc) & (ratios >= min_rpm)
        codes_passing_filter = molecule_cell_sums.index[filter_]

        row_inds = np.where(pd.Series(self.molecules.index).isin(codes_passing_filter))[0]
        read_data = self.reads.data.tocsr()[row_inds, :].todense()
        reads = pd.DataFrame(read_data, index=self.reads.index[row_inds],
                             columns=self.reads.columns)
        molecule_data = self.molecules.data.tocsr()[row_inds, :].todense()
        molecules = pd.DataFrame(molecule_data, index=self.molecules.index[row_inds],
                                 columns=self.molecules.columns)
        metadata = self.metadata.ix[row_inds, :]
        return Experiment(reads=reads, molecules=molecules, metadata=metadata)

    def plot_mitochondrial_molecule_fraction(self, fig=None, ax=None,
                                             title='Dead Cell Identification Plot'):
        """ plot the fraction

        :param title: title for the plot (e.g. the sample name)
        :param fig: figure
        :param ax: axis
        :return: fig, ax
        """
        if self._normalized:
            raise RuntimeError('plot_molecules_vs_reads_per_molecule() should be run on '
                               'unnormalized data')

        if self.is_sparse():
            raise SparseMatrixError('Must convert to dense matrix before calling this '
                                    'function. Use self.remove_non_cell_barcodes()')

        mt_genes = self.molecules.columns[self.molecules.columns.str.contains('MT-')]
        mt_counts = self.molecules[mt_genes].sum(axis=1)
        library_size = self.molecules.sum(axis=1)

        fig, ax = get_fig(fig=fig, ax=ax)
        x, y, z = density_2d(library_size, mt_counts / library_size)
        ax.scatter(x, y, s=size, edgecolors='none', c=z, cmap=cmap)
        ax.set_title('Mitochondrial Fraction')
        ax.set_xlabel('Molecules per cell')
        ax.set_ylabel('Mitochondrial molecules / Total molecules')
        ax.set_title(title)
        xlim = ax.get_xlim()
        ax.set_xlim((0, xlim[1]))
        ylim = ax.get_ylim()
        ax.set_ylim((0, ylim[1]))
        sns.despine(ax=ax)

        return fig, ax

    def exclude_dead_cells_with_high_mt_fraction(self, max_mt_fraction=0.2):
        """remove cells containing > max_mt_fraction mitochondrial molecules

        :param max_mt_fraction:
        :return: Experiment
        """
        if self._normalized:
            raise RuntimeError('plot_molecules_vs_reads_per_molecule() should be run on '
                               'unnormalized data')

        mt_genes = self.molecules.columns[self.molecules.columns.str.contains('MT-')]
        mt_counts = self.molecules[mt_genes].sum(axis=1)
        library_size = self.molecules.sum(axis=1)
        ratios = mt_counts / library_size
        pass_filter = ratios.index[ratios <= max_mt_fraction]

        molecules = self.molecules.ix[pass_filter]
        reads = self.molecules.ix[pass_filter]
        metadata = self.metadata.ix[pass_filter]

        return Experiment(reads=reads, molecules=molecules, metadata=metadata)

    def plot_molecules_vs_genes(self, fig=None, ax=None, title=None):
        """
        should be linear relationship
        :param ax:
        :param fig:
        :return:
        """
        if self._normalized:
            raise RuntimeError('plot_molecules_vs_reads_per_molecule() should be run on '
                               'unnormalized data')

        fig, ax = get_fig(fig=fig, ax=ax)
        molecule_counts = self.molecules.sum(axis=1)
        gene_counts = np.sum(self.molecules > 0, axis=1)
        x, y, z = density_2d(np.log10(molecule_counts), np.log10(gene_counts))
        ax.scatter(x, y, c=z, edgecolor='none', s=size, cmap=cmap)
        ax.set_xlabel('log10(Number of molecules)')
        ax.set_ylabel('log10(Number of genes)')
        ax.set_title('')
        sns.despine(ax=ax)

        return fig, ax

    def remove_low_complexity_cells(self):
        # todo implement, based on below-fit cells in the above plot
        if self._normalized:
            raise RuntimeError('plot_molecules_vs_reads_per_molecule() should be run on '
                               'unnormalized data')

        raise NotImplementedError

    def normalize_data(self):
        """
        uses AVO method of normalization
        :return: Experiment
        """

        if self.is_sparse():
            raise RuntimeError('Please convert to dense data before attempting to '
                               'normalize using Experiment.remove_non_cell_barcodes().')

        molecule_sums = self.molecules.sum(axis=1)
        molecules = self.molecules.div(molecule_sums, axis=0)\
            .mul(np.median(molecule_sums), axis=0)
        read_sums = self.reads.sum(axis=1)
        reads = self.reads.div(read_sums, axis=0)\
            .mul(np.median(read_sums), axis=0)
        exp = Experiment(reads=reads, molecules=molecules, metadata=self.metadata)
        exp._normalized = True

        # check that none of the genes are empty; if so remove them
        nonzero_genes = exp.molecules.sum(axis=0) != 0
        exp.molecules = exp.molecules.ix[:, nonzero_genes]
        exp.reads = exp.reads.ix[:, nonzero_genes]
        return exp

    def run_pca_python(self, n_components=100):

        # todo refactored no_comp -> n_components=100
        X = self.molecules.values  # todo added this
        # Make sure data is zero mean
        X = np.subtract(X, np.amin(X))
        X = np.divide(X, np.amax(X))

        # Compute covariance matrix
        if (X.shape[1] < X.shape[0]):
            C = np.cov(X, rowvar=0)
        # if N>D, we better use this matrix for the eigendecomposition
        else:
            C = np.multiply((1/X.shape[0]), np.dot(X, X.T))

        # Perform eigendecomposition of C
        C[np.where(np.isnan(C))] = 0
        C[np.where(np.isinf(C))] = 0
        l, M = np.linalg.eig(C)

        # Sort eigenvectors in descending order
        ind = np.argsort(l)[::-1]
        l = l[ind]
        if n_components < 1:
            n_components = np.where(np.cumsum(np.divide(l, np.sum(l)), axis=0) >= n_components)[0][0] + 1
            print('Embedding into ' + str(n_components) + ' dimensions.')
        if n_components > M.shape[1]:
            n_components = M.shape[1]
            print('Target dimensionality reduced to ' + str(n_components) + '.')

        M = M[:, ind[:n_components]]
        l = l[:n_components]

        # Apply mapping on the data
        if X.shape[1] >= X.shape[0]:
            M = np.multiply(np.dot(X.T, M), (1 / np.sqrt(X.shape[0] * l)).T)
        mappedX = np.dot(X, M)

        return mappedX, M, l

    def run_pca(self, no_components=100):
        """

        :param no_components:
        :return:
        """

        # Random tag to allow for multiple diffusion map runs
        rand_tag = random.random()

        # Write to csv file
        self.molecules.to_csv('/tmp/pc_data_%f.csv' % rand_tag)

        # Construct matlab script; Change directory to diffusion geometry
        matlab_cmd = "\" cd {};".format(os.path.expanduser('~/.seqc/tools'))
        # Set-up diffusion geometry
        matlab_cmd += "data = csvread('/tmp/pc_data_%f.csv', 1, 1);" % rand_tag
        # PCA
        matlab_cmd += "[mappedX, mapping]  = pca2(data, %d);" % no_components

        # Write results to file
        matlab_cmd += " csvwrite('/tmp/pc_mapped_x_%f.csv', mappedX);" % rand_tag
        matlab_cmd += " csvwrite('/tmp/pc_mapping_M_%f.csv', mapping.M);" % rand_tag
        matlab_cmd += (" csvwrite('/tmp/pc_mapping_lambda_%f.csv', mapping.lambda); exit;"
                       "\"" % rand_tag)

        # Run matlab command
        call(['matlab', '-nodesktop', '-nosplash', '-r %s' % matlab_cmd])

        # Read in results
        loadings = pd.DataFrame.from_csv('/tmp/pc_mapping_M_%f.csv' % rand_tag,
                                         header=None, index_col=None)
        loadings.index = self.molecules.columns
        eigenvalues = pd.DataFrame.from_csv('/tmp/pc_mapping_lambda_%f.csv' % rand_tag,
                                            header=None, index_col=None)
        self.pca = {'loadings': loadings, 'eigenvalues': eigenvalues}

        # Clean up
        os.remove('/tmp/pc_data_%f.csv' % rand_tag)
        os.remove('/tmp/pc_mapped_x_%f.csv' % rand_tag)
        os.remove('/tmp/pc_mapping_M_%f.csv' % rand_tag)
        os.remove('/tmp/pc_mapping_lambda_%f.csv' % rand_tag)

    def plot_pca_variance_explained(self, fig=None, ax=None, n_components=30,
                                    ylim=(0, 0.1)):
        # plot the eigenvalues
        fig, ax = get_fig(fig=fig, ax=ax)
        ax.plot(np.ravel(self.pca['eigenvalues'].values))
        plt.ylim(ylim)
        plt.xlim((0, n_components))
        sns.despine(ax=ax)
        return fig, ax

    def run_tsne(self, n_components=15):
        """
        normally run on PCA components; 1st component is normally excluded

        :param n_components:
        :return:
        """
        data = deepcopy(self.molecules)
        data -= np.min(np.ravel(data))
        data /= np.max(np.ravel(data))
        data = pd.DataFrame(np.inner(data, self.pca['loadings'].iloc[:, 0:n_components].T),
                            index=self.molecules.index)

        self.tsne = pd.DataFrame(bh_sne(data),
                                 index=self.molecules.index, columns=['x', 'y'])

    def plot_tsne(self, fig=None, ax=None, title='tSNE projection'):
        fig, ax = get_fig(fig=fig, ax=ax)
        x, y, z = density_2d(self.tsne['x'], self.tsne['y'])
        plt.scatter(x, y, c=z, s=size, cmap=cmap)
        ax.xaxis.set_major_locator(plt.NullLocator())
        ax.yaxis.set_major_locator(plt.NullLocator())
        ax.set_title(title)
        return fig, ax

    def run_phenograph(self, n_pca_components=15):
        """
        normally run on PCA components; 1st component is normally excluded
        :param n_pca_components:
        :return:
        """
        data = deepcopy(self.molecules)
        data -= np.min(np.ravel(data))
        data /= np.max(np.ravel(data))
        data = pd.DataFrame(np.inner(data, self.pca['loadings'].iloc[:, 0:n_pca_components].T),
                            index=self.molecules.index)

        communities, graph, Q = phenograph.cluster(data)
        self.cluster_assignments = pd.Series(communities, index=data.index)

    def plot_phenograph(self, fig=None, ax=None, labels=None):

        if self.tsne is None:
            raise RuntimeError('Cannot plot phenograph before generating tSNE'
                               'projection. Please call Experiment.run_tsne()')

        fig, ax = get_fig(fig=fig, ax=ax)
        clusters = sorted(set(self.cluster_assignments))
        colors = qualitative_colors(len(clusters))
        for i in range(len(clusters)):
            if labels:
                label=labels[i]
            else:
                label = clusters[i]
            data = self.tsne.ix[self.cluster_assignments == clusters[i], :]
            ax.plot(data['x'], data['y'], c=colors[i], linewidth=0, marker='o',
                    markersize=np.sqrt(size), label=label)
        plt.legend()
        ax.xaxis.set_major_locator(plt.NullLocator())
        ax.yaxis.set_major_locator(plt.NullLocator())
        return fig, ax

    # todo add plot for library sizes (see notebook ajc_worklogs/manuscript/MSKCC/TIL_processing.ipynb)
    # todo for this one, going to want access to original library size
    def remove_clusters_based_on_library_size_distribution(self, clusters: list):
        """
        throws out low library size clusters; normally there is only one, but worth
        visual inspection.

        note: retains cluster_assignments for clusters that are retained. Other analytical
        results are discarded, since they must be re-run on the reduced data.

        :param clusters: clusters to DISCARD
        :return: Experiment
        """
        retain = self.cluster_assignments.index[~self.cluster_assignments.isin(clusters)]
        molecules = self.molecules.ix[retain, :]
        reads = self.reads.ix[retain, :]
        metadata = self.metadata.ix[retain, :]
        experiment = Experiment(reads=reads, molecules=molecules, metadata=metadata)
        experiment._cluster_assignments = self.cluster_assignments[retain]

        # remove any genes that now have zero expression
        nonzero_genes = experiment.molecules.sum(axis=0) != 0
        experiment.molecules = experiment.molecules.ix[:, nonzero_genes]
        experiment.reads = experiment.reads.ix[:, nonzero_genes]
        return experiment

    def run_diffusion_map(self, knn=10, n_diffusion_components=20, n_pca_components=15,
                          params=pd.Series()):
        """
        :param knn:
        :param n_diffusion_components:
        :param params:
        :return:
        """

        if not self.pca:
            raise RuntimeError('Please run PCA before calculating diffusion components. '
                               'DiffusionGeometry is run on the PCA decomposition of '
                               'self.molecules.')

        # Random tag to allow for multiple diffusion map runs
        rand_tag = random.random()

        data = deepcopy(self.molecules)
        data -= np.min(np.ravel(data))
        data /= np.max(np.ravel(data))
        data = pd.DataFrame(np.inner(data, self.pca['loadings'].iloc[:, 0:n_pca_components].T),
                            index=self.molecules.index)

        # Write to csv file
        data.to_csv('/tmp/dm_data_%f.csv' % rand_tag)

        # Construct matlab script
        # Change directory to diffusion geometry
        matlab_cmd = "\" cd {};".format(os.path.expanduser(
                '~/.seqc/tools/DiffusionGeometry'))
        # Set up diffusion geometry
        matlab_cmd += "Startup_DiffusionGeometry;"
        matlab_cmd += "data = csvread('/tmp/dm_data_%f.csv', 1, 1);" % rand_tag

        # Diffusion map parameters
        dm_params = pd.Series()
        dm_params['Normalization'] = "'smarkov'"
        dm_params['Epsilon'] = str(1)
        dm_params['kNN'] = str(knn)
        dm_params['kEigenVecs'] = str(n_diffusion_components)
        dm_params['Symmetrization'] = "'W+Wt'"

        # Update additional parameters
        for i in params.index:
            dm_params[i] = str(params[i])

        # Add to matlab command
        print('Diffusion geometry parameters..')
        for i in dm_params.index:
            print("%s: %s" % (i, dm_params[i]))
            matlab_cmd += "GraphDiffOpts.%s = %s; " % (i, dm_params[i])

        # Run diffusion map
        matlab_cmd += "G = GraphDiffusion(data', 0, GraphDiffOpts);"

        # Write eigen vectors to file
        matlab_cmd += " csvwrite('/tmp/dm_eigs_%f.csv', G.EigenVecs);" % rand_tag
        matlab_cmd += " csvwrite('/tmp/dm_eig_vals_%f.csv', G.EigenVals); exit;\"" % (
            rand_tag)

        # Run matlab command
        call(['matlab', '-nodesktop', '-nosplash', '-r %s' % matlab_cmd])

        # Read in results
        eigs = pd.DataFrame.from_csv('/tmp/dm_eigs_%f.csv' % rand_tag, header=None,
                                     index_col=None)
        eigs.index = self.molecules.index
        vals = pd.DataFrame.from_csv('/tmp/dm_eig_vals_%f.csv' % rand_tag, header=None,
                                     index_col=None)

        self.diffusion_eigenvectors = eigs
        self.diffusion_eigenvalues = vals

        # Clean up
        os.remove('/tmp/dm_data_%f.csv' % rand_tag)
        os.remove('/tmp/dm_eigs_%f.csv' % rand_tag)
        os.remove('/tmp/dm_eig_vals_%f.csv' % rand_tag)

    def plot_diffusion_components(self, title='Diffusion Components'):
        if self.tsne is None:
            raise RuntimeError('Please run tSNE before plotting diffusion components.')

        height = int(2 * np.ceil(self.diffusion_eigenvalues.shape[0] / 5))
        width = 10
        fig = plt.figure(figsize=[width, height])
        n_rows = int(height / 2)
        n_cols = int(width / 2)
        gs = plt.GridSpec(n_rows, n_cols)

        for i in range(self.diffusion_eigenvectors.shape[1]):
            ax = plt.subplot(gs[i // n_cols, i % n_cols])
            plt.scatter(self.tsne['x'], self.tsne['y'], c=self.diffusion_eigenvectors[i],
                        cmap=cmap, edgecolors='none', s=size)
            ax.set_axis_off()

        fig.suptitle(title)
        return fig

    # todo extremely slow for large numbers of genes, components. parallelize
    # across components using multiprocessing pool similar to how the below GSEA function
    # was programmed.
    # NOTE: it will be more complex because larger arrays are involved in these
    # calculations. Need to use memmapped numpy arrays to avoid large-scale copying of the
    # input data.
    def determine_gene_diffusion_correlations(self, components=None, no_cells=10):
        """

        :param no_cells:
        :return:
        """

        if components is None:
            components = np.arange(self.diffusion_eigenvectors.shape[1])

        components = components[components != 0]

        # Empty container
        empty = np.zeros([self.molecules.shape[1], len(components)])
        self.diffusion_map_correlations = pd.DataFrame(
            empty, index=self.molecules.columns,
            columns=components)

        # Determine correlation for each component
        for component in components:

            # Cell order
            order = self.diffusion_eigenvectors.iloc[:, component]\
                .sort_values(inplace=False).index
            x = np.ravel(pd.rolling_mean(self.diffusion_eigenvectors.ix[order, component],
                                         no_cells))
            x = x[no_cells:]

            # For each gene
            for gene in self.molecules.columns:
                # Rolling mean
                vals = pd.rolling_mean(self.molecules.ix[order, gene],
                                       no_cells)[no_cells:]

                # Determine correlation
                self.diffusion_map_correlations.ix[gene, component] = pearsonr(x, vals)[0]

        # Reset NaNs
        self.diffusion_map_correlations = self.diffusion_map_correlations.fillna(0)

    def plot_gene_component_correlations(
            self, components=None, fig=None, ax=None,
            title='Gene vs. Diffusion Component Correlations'):
        """ plots gene-component correlations for a subset of components

        :param title: str, title for the plot
        :param components: Iterable of integer component numbers
        :param fig: Figure
        :param ax: Axis
        :return: fig, ax
        """
        fig, ax = get_fig(fig=fig, ax=ax)
        if self.diffusion_map_correlations is None:
            raise RuntimeError('Please run determine_gene_diffusion_correlations() '
                               'before attempting to visualize the correlations.')
        colors = qualitative_colors(self.diffusion_map_correlations.shape[1])

        if components is None:
            components = self.diffusion_map_correlations.columns

        for c in components:
            sns.kdeplot(self.diffusion_map_correlations.iloc[:, c].fillna(0), label=c,
                        ax=ax, color=colors[c])
        sns.despine(ax=ax)
        ax.set_title(title)
        ax.set_xlabel('correlation')
        ax.set_ylabel('gene density')
        plt.legend()
        return fig, ax

    @staticmethod
    def _gmt_options():
        mouse_options = os.listdir(os.path.expanduser('~/.seqc/tools/mouse'))
        human_options = os.listdir(os.path.expanduser('~/.seqc/tools/human'))
        print('Available GSEA .gmt files:\n\nmouse:\n{m}\n\nhuman:\n{h}\n'.format(
                m='\n'.join(mouse_options),
                h='\n'.join(human_options)))
        print('Please specify the gmt_file parameter as gmt_file=(organism, filename)')

    @staticmethod
    def _gsea_process(c, diffusion_map_correlations, output_stem, gmt_file):

        # save the .rnk file
        out_dir, out_prefix = os.path.split(output_stem)
        genes_file = '{stem}_cmpnt_{component}.rnk'.format(
                stem=output_stem, component=c)
        ranked_genes = diffusion_map_correlations.ix[:, c]\
            .sort_values(inplace=False, ascending=False)

        # set any NaN to 0
        ranked_genes = ranked_genes.fillna(0)

        # dump to file
        pd.DataFrame(ranked_genes).to_csv(genes_file, sep='\t', header=False)

        # Construct the GSEA call
        cmd = shlex.split(
            'java -cp {user}/.seqc/tools/gsea2-2.2.1.jar -Xmx1g '
            'xtools.gsea.GseaPreranked -collapse false -mode Max_probe -norm meandiv '
            '-nperm 1000 -include_only_symbols true -make_sets true -plot_top_x 0 '
            '-set_max 500 -set_min 50 -zip_report false -gui false -rnk {rnk} '
            '-rpt_label {out_prefix}_{component} -out {out_dir}/ -gmx {gmt_file}'
            ''.format(user=os.path.expanduser('~'), rnk=genes_file,
                      out_prefix=out_prefix, component=c, out_dir=out_dir,
                      gmt_file=gmt_file))

        # Call GSEA
        p = Popen(cmd, stderr=PIPE)
        _, err = p.communicate()

        # remove annoying suffix from GSEA
        if err:
            return err
        else:
            pattern = out_prefix + '_' + str(c) + '.GseaPreranked.[0-9]*'
            repl = out_prefix + '_' + str(c)
            files = os.listdir(out_dir)
            for f in files:
                mo = re.match(pattern, f)
                if mo:
                    curr_name = mo.group(0)
                    shutil.move('{}/{}'.format(out_dir, curr_name),
                                '{}/{}'.format(out_dir, repl))
                    return err

            # execute if file cannot be found
            return b'GSEA output pattern was not found, and could not be changed.'

    def run_gsea(self, output_stem, gmt_file=None, components=None):
        """

        :param output_stem:  the file location and prefix for the output of GSEA
        :param gmt_file:
        :param components:
        :return:
        """

        out_dir, out_prefix = os.path.split(output_stem)
        out_dir += '/'
        os.makedirs(out_dir, exist_ok=True)

        if self.diffusion_eigenvectors is None:
            raise RuntimeError('Please run self.calculate_diffusion_map_components() '
                               'before running GSEA to annotate those components.')

        if not gmt_file:
            self._gmt_options()
            return
        else:
            if not len(gmt_file) == 2:
                raise ValueError('gmt_file should be a tuple of (organism, filename).')
            gmt_file = os.path.expanduser('~/.seqc/tools/{}/{}').format(*gmt_file)

        if components is None:
            components = self.diffusion_map_correlations.columns

        n_workers = min((multiprocessing.cpu_count() - 1, len(components)))
        print('initializing {n_workers} worker processes and mapping GSEA to them.'
              ''.format(n_workers=n_workers))

        partial_gsea_process = partial(
                self._gsea_process,
                diffusion_map_correlations=self.diffusion_map_correlations,
                output_stem=output_stem,
                gmt_file=gmt_file)

        pool = multiprocessing.Pool(n_workers)
        errors = pool.map(partial_gsea_process, components)
        return errors

    def select_genes_from_diffusion_components(self, components, plot=False):
        """

        done based on GSEA enrichments

        :param components:
        :return:
        """
        if self.diffusion_map_correlations is None:
            raise RuntimeError('Please run self.determine_gene_diffusion_correlations() '
                               'before selecting genes based on those correlations.')

        # Plot the correlation distributions along the selected components
        if plot:
            _ = self.plot_gene_component_correlations(components)

        # Select the genes
        use_genes = list()
        for component in components:
            cutoff = (np.mean(self.diffusion_map_correlations.iloc[:, component]) +
                      2 * np.std(self.diffusion_map_correlations.iloc[:, component]))
            use_genes = use_genes + list(self.diffusion_map_correlations.index[abs(
                    self.diffusion_map_correlations.iloc[:, component]) > cutoff])

        # Unique genes
        use_genes = list(set(use_genes))

        # Create new scdata object
        subset_molecules = self.molecules.ix[:, use_genes]
        subset_reads = self.reads.ix[:, use_genes]
        metadata = self.metadata.ix[:, use_genes]

        # throws out analysis results; this makes sense
        return Experiment(subset_reads, subset_molecules, metadata)

    def pairwise_differential_expression(self, c1, c2, alpha=0.05):
        """
        carry out differential expression (post-hoc tests) between cells c1 and cells c2,
        using bh-FDR to correct for multiple tests
        :param alpha:
        :param g1:
        :param g2:
        :return:
        """
        g1 = self.molecules.loc[c1]
        g2 = self.molecules.loc[c2]

        res = pd.Series(index=self.molecules.columns, dtype=float)
        for g in self.molecules.columns:
            try:
                res[g] = mannwhitneyu(np.ravel(g1[g]), np.ravel(g2[g]))[1]
            except ValueError:
                res[g] = 1

        pval_corrected = multipletests(res.values, alpha=alpha, method='fdr_tsbh')[1]

        return pd.Series(pval_corrected, index=self.molecules.columns,
                         dtype=float).sort_values()

    # todo test
    def single_gene_differential_expression(self, gene, alpha):
        """
        carry out kruskal-wallis non-parametric (rank-wise) ANOVA with two-stage bh-FDR
        correction to determine if gene is differentially expressed in different clusters

        Explicitly, this tests the null hypothesis that all samples are drawn from the
        same distribution.

        if the KW-test is significant, Post-hoc rank-sum tests on the genes whose
        corrected p-values are below alpha determine the specific samples that are
        differentially expressed

        :param gene:
        :return: KW p-val, pairwise ranksums pvals
        """
        gene = self.molecules.ix[:, gene]

        clusters = np.unique(self.cluster_assignments)
        samples = [gene.ix[self.cluster_assignments == c].values for
                   c in clusters]

        pval = kruskalwallis(*samples)[1]

        if pval < alpha:
            pairs = list(combinations_with_replacement(np.arange(len(samples)), 2))
            pvals = pd.Series(index=pairs, dtype=float)
            for a, b in pairs:
                pvals.ix[(a, b)] = mannwhitneyu(gene[self.cluster_assignments == a],
                                                gene[self.cluster_assignments == b])[1]
            return pval, pvals.sort_values()
        else:
            return pval, None

    def differential_expression(self, alpha=0.05):
        """
        carry out kruskal-wallis non-parametric (rank-wise) ANOVA with two-stage bh-FDR
        correction to determine the genes that are differentially expressed in at least
        two populations, as defined by self._cluster_assignments

        Explicitly, this tests the null hypothesis that all samples are drawn from the
        same distribution.

        Post-hoc tests on the genes whose corrected p-values are below alpha determine
        the specific samples that are differentially expressed

        Need:
        (1) sense of an expression floor to call expressed vs. not expressed.

        Intuitively: molecules vs. population size? if differential expression exists
        then presumably there is a non-zero mean in at least one group.

        Does not deal with the circumstance where all means are non-zero. can call
        positive expression but not negative.

        :param experiments:
        :param gene:
        :return:
        """

        if not isinstance(self.cluster_assignments, pd.Series):
            raise RuntimeError('Please determine cluster assignments before carrying out '
                               'differential expression')

        # sort self.molecules by cluster assignment
        idx = np.argsort(self.cluster_assignments)
        data = self.molecules.values[idx, :]

        # get points to split the array
        split_indices = np.where(np.diff(self._cluster_assignments[idx]))[0] + 1

        # create the function for the kruskal test
        f = lambda v: kruskalwallis(*np.split(v, split_indices))[1]

        # get p-values
        pvals = np.apply_along_axis(f, 0, data)

        # correct the pvals
        alpha = 0.05
        reject, pval_corrected, _, _ = multipletests(pvals, alpha, method='fdr_tsbh')

        return pd.Series(pval_corrected, index=self.molecules.columns).sort_values(inplace=False)

    def plot_gene_expression(self, genes, suptitle='tSNE-projected Gene Expression'):

        not_in_dataframe = set(genes).difference(self.molecules.columns)
        if not_in_dataframe:
            if len(not_in_dataframe) < len(genes):
                print('The following genes were either not observed in the experiment, '
                      'or the wrong gene symbol was used: {!r}'.format(not_in_dataframe))
            else:
                print('None of the listed genes were observed in the experiment, or the '
                      'wrong symbols were used.')
                return

        # remove genes missing from experiment
        genes = set(genes).difference(not_in_dataframe)

        height = int(2 * np.ceil(len(genes) / 5))
        width = 10
        fig = plt.figure(figsize=[width, height])
        n_rows = int(height / 2)
        n_cols = int(width / 2)
        gs = plt.GridSpec(n_rows, n_cols)

        axes = []
        for i, g in enumerate(genes):
            ax = plt.subplot(gs[i // n_cols, i % n_cols])
            axes.append(ax)
            plt.scatter(self.tsne['x'], self.tsne['y'], c=np.arcsinh(self.molecules[g]),
                        cmap=cmap, edgecolors='none', s=size)
            ax.set_axis_off()
            ax.set_title(g)

        fig.suptitle(suptitle)

        return fig, axes

    def plot_aggregate_gene_expression(self, genes, fig=None, ax=None, title=None):

        # remove genes missing from experiment
        not_in_dataframe = set(genes).difference(self.molecules.columns)
        if not_in_dataframe:
            if len(not_in_dataframe) < len(genes):
                print('The following genes were either not observed in the experiment, '
                      'or the wrong gene symbol was used: {!r}'.format(not_in_dataframe))
            else:
                print('None of the listed genes were observed in the experiment, or the '
                      'wrong symbols were used.')
                return
        genes = set(genes).difference(not_in_dataframe)

        # create a generic title if not passed
        if not title:
            if len(repr(genes)) > 40:
                title = 'Aggregated gene expression of:\n{:.40s}...'.format(repr(genes))
            else:
                title = 'Aggregated gene expression of:\n{:.40s}'.format(repr(genes))

        # plot the data
        fig, ax = get_fig(fig=fig, ax=ax)
        ax.scatter(self.tsne['x'], self.tsne['y'],
                   c=np.arcsinh(self.molecules[genes].sum(axis=1)),
                   cmap=cmap, edgecolors='none', s=size)
        ax.xaxis.set_major_locator(plt.NullLocator())
        ax.yaxis.set_major_locator(plt.NullLocator())
        ax.set_title(title)
        return fig, ax

    def remove_housekeeping_genes(self, organism='human'):
        """
        Remove snoRNA and miRNAs because they are undetectable by the library construction
        procedure; their identification is a probable false-postive.

        Next, removes housekeeping genes based on GO annotations.
        :param organism: options: [human, mouse], default=human.
        """

        genes = self.molecules.columns
        # MT genes
        genes = genes[~genes.str.contains('MT-')]
        # "AL" miRNA
        genes = genes[~genes.str.contains('^AL[0-9]')]
        # "FP" miRNA
        genes = genes[~genes.str.contains('^FP[0-9]')]

        # housekeeping gene ontologies
        cats = ['GO_0016072|rRNA metabolic process', 'GO_0006412|translation',
                'GO_0006414|translational elongation', 'GO_0006364|rRNA processing']

        with open(os.path.expanduser('~/.seqc/tools/{organism}/gofat.bp.v1.0.gmt.txt'
                                     ''.format(organism=organism))) as f:
            data = {fields[0]: fields[1:] for fields in
                    [line[:-1].split('\t') for line in f.readlines()]}

        hk_genes = {'B2M', 'AC091053.2', 'MALAT1', 'TMSB4X', 'RPL13AP5-RPL13A',
                    'RP11-864I4.1-EEF1G', 'GAPDH', 'CTD-3035D6.1', 'CTC-575D19.1', 'ACTB',
                    'ACTG1', 'RP5-1056L3.3'}
        for c in cats:
            hk_genes.update(data[c])

        # remove genes
        genes = genes.difference(hk_genes)

        return Experiment(self.molecules[genes], self.reads[genes], self.metadata)

    def annotate_clusters_by_expression(self, alpha=0.05):
        """
        considering only genes which are differentially expressed across the
        population, label each cluster with the genes that are significantly expressed
        in one population relative to > 80% of other populations.

        :param alpha: allowable type-I error for ANOVA
        :return: boolean matrix (clusters x genes)
        """
        if not isinstance(self.cluster_assignments, pd.Series):
            raise RuntimeError('Please determine cluster assignments before carrying out '
                               'differential expression')

        # sort self.molecules by cluster assignment
        idx = np.argsort(self.cluster_assignments)
        data = self.molecules.values[idx, :]

        # get points to split the array
        split_indices = np.where(np.diff(self.cluster_assignments[idx]))[0] + 1

        # create the function for the kruskal test
        f = lambda v: kruskalwallis(*np.split(v, split_indices))[1]

        # get p-values
        pvals = np.apply_along_axis(f, 0, data)

        # correct the pvals
        reject, pval_corrected, _, _ = multipletests(pvals, alpha, method='fdr_tsbh')

        # restrict genes to ones with significant anova results
        data = data[:, reject]  # cells x genes

        # get cluster positions in data
        split_indices = np.concatenate(
                [[0], split_indices, [len(self.cluster_assignments)]])  # add end points
        cluster_slicers = [(i - 1, slice(split_indices[i - 1], split_indices[i]))
                           for i in range(1, len(split_indices))]

        # get pairs of clusters
        cluster_pairs = list(combinations(cluster_slicers, 2))
        n_clusters = len(set(self.cluster_assignments))
        n_genes = np.sum(reject)
        global_results = np.zeros((n_clusters, n_genes), dtype=np.bool)

        # need to know gene means in order to fill up the experiment matrix.
        mean_expression = np.zeros((n_clusters, n_genes), dtype=np.float)
        for id_, slice_ in cluster_slicers:
            mean_expression[id_, :] = np.mean(data[slice_, :], axis=0)

        # for each gene in data, carry out all vs all mann-whitney
        local_results_shape = (n_clusters, n_clusters)
        for i, gene in enumerate(data.T):
            gene_results = np.zeros(local_results_shape, dtype=np.bool)
            for (aid, a), (bid, b) in cluster_pairs:  # all pairwise tests
                try:
                    p = mannwhitneyu(gene[a], gene[b])[1]
                except ValueError:
                    continue  # neither cell expresses this gene
                if p < alpha:
                    if mean_expression[aid, i] > mean_expression[bid, i]:
                        gene_results[aid, bid] = 1
                    else:
                        gene_results[bid, aid] = 1

            # update global results
            res = gene_results.sum(axis=1) > gene_results.shape[0] * .8
            global_results[:, i] = res

        return pd.DataFrame(global_results, columns=self.molecules.columns[reject])
