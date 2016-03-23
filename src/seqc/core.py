import re
import os
import random
import pickle
import warnings
from copy import deepcopy
from collections import defaultdict
from subprocess import call
import numpy as np
import pandas as pd
import tables as tb
from scipy.sparse import coo_matrix
from scipy.stats import gaussian_kde, pearsonr
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

matplotlib.rc('savefig', **{'dpi': '300'})
cmap = matplotlib.cm.viridis
size = 6


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
            self._index = np.array(item)
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
        if not (isinstance(item, pd.DataFrame) or item is None):
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

    def is_sparse(self):
        if all(isinstance(o, SparseFrame) for o in (self.reads, self.molecules)):
            return True
        else:
            return False

    def is_dense(self):
        return not self.is_sparse()

    @staticmethod
    def concatenate(experiments, metadata_labels=None):
        """
        Concatenate a set of Experiment objects. Each cell is considered to be unique,
        even if they share the same barcodes. If collisions are found, _%d will be added
        to each cell, where %d is equal to the position of the Experiment in the input
        experiments list.

        If metadata_labels are provided, a new entry in metadata will be included for
        each cell. This can be useful to mark each cell by which experiment it originated
        from in order to track batch effects, for example.


        :param self:
        :param experiments:
        :param metadata_labels:
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
        cell2int = dict(zip(cells, range(len(cells))))
        int2cell = {v: k for k, v in cell2int.items()}
        gene2int = dict(zip(genes, range(len(genes))))
        int2gene = {v: k for k, v in gene2int.items()}

        # merge molecule counts
        new_molecule_data = defaultdict(int)
        new_read_data = defaultdict(int)
        for e in experiments:

            # map molecule and cell indices back to their respective cell and gene ids
            local_cell_map = dict(zip(e.index, range(len(e.index))))
            local_gene_map = dict(zip(e.columns, range(len(e.columns))))

            for mols, reads, cell, gene in zip(e.molecules.data, e.reads.data,
                                               e.molecules.row, e.molecules.col):
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

        # make coo matrices
        reads = coo_matrix((read_data, (row, col)), dtype=np.uint32)
        molecules = coo_matrix((molecule_data, (row, col)), dtype=np.uint32)

        # get gene ids and cell ids for SparseFrame construction
        row_index = [int2cell[v] for v in row]
        col_index = [int2gene[v] for v in col]

        sparse_reads = SparseFrame(reads, row_index, col_index)
        sparse_molecules = SparseFrame(molecules, row_index, col_index)

        return Experiment(reads=sparse_reads, molecules=sparse_molecules)

    @staticmethod
    def merge(experiments):
        """
        Merge a set of Experiment objects. Cells with duplicate cell barcodes will have
        their data summed across experiments.

        If metadata_labels are provided, a new entry in metadata will be included for
        each cell. This can be useful to mark each cell by which experiment it originated
        from in order to track batch effects, for example.

        :param self:
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
        cell_index = list(cells)
        gene_index = list(genes)
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
        sparse_molecules = SparseFrame(molecules, gene_index, cell_index)

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
                                             ylim=(0, 150)):
        """
        plots log10 molecules counts per barcode vs reads per molecule / molecules per
        barcode
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
        ax.set_ylim(ylim)
        xlim = ax.get_xlim()
        ax.set_xlim((np.log10(min_molecules), xlim[1]))  # todo if this fails later, this is the problem

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

    def plot_mitochondrial_molecule_fraction(self, fig=None, ax=None):
        """ plot the fraction

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
        xlim = ax.get_xlim()
        ax.set_xlim((0, xlim[1]))
        ylim = ax.get_ylim()
        ax.set_ylim((0, ylim[1]))

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

    def plot_molecules_vs_genes(self, fig=None, ax=None):
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
        return exp

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
        self.pca = {'loadings': loadings, eigenvalues: eigenvalues}

        # Clean up
        os.remove('/tmp/pc_data_%f.csv' % rand_tag)
        os.remove('/tmp/pc_mapped_x_%f.csv' % rand_tag)
        os.remove('/tmp/pc_mapping_M_%f.csv' % rand_tag)
        os.remove('/tmp/pc_mapping_lambda_%f.csv' % rand_tag)

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

    def plot_tsne(self, fig=None, ax=None):
        fig, ax = get_fig(fig=fig, ax=ax)
        x, y, z = density_2d(self.tsne['x'], self.tsne['y'])
        plt.scatter(x, y, c=z, s=size, cmap=cmap)
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
                    markersize=size, label=label)
        plt.legend()
        ax.xaxis.set_major_locator(plt.NullLocator())
        ax.yaxis.set_major_locator(plt.NullLocator())
        return fig, ax

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
        return experiment

    def calculate_diffusion_map_components(self, knn=10, n_diffusion_components=20,
                                           n_pca_components=15, params=pd.Series()):
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

    def plot_diffusion_components(self, fig=None, ax=None):
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

    def determine_gene_diffusion_correlations(self, no_cells=10):
        """

        :param no_cells:
        :return:
        """
        # Container
        empty = np.zeros([self.molecules.shape[1], self.diffusion_eigenvectors.shape[1]])
        self.diffusion_map_correlations = pd.DataFrame(
            empty, index=self.molecules.columns,
            columns=self.diffusion_eigenvectors.columns)

        # Determine correlation for each component
        for component in range(self.diffusion_eigenvectors.shape[1]):

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

    def plot_gene_component_correlations(self, components, fig=None, ax=None):
        """ plots gene-component correlations for a subset of components

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
        for c in components:
            sns.kdeplot(self.diffusion_map_correlations.iloc[:, c], label=c, ax=ax,
                        color=[colors[c]])
        sns.despine(ax=ax)
        return fig, ax

    @staticmethod
    def _gmt_options():
        mouse_options = os.listdir(os.path.expanduser('~/.seqc/tools/mouse'))
        human_options = os.listdir(os.path.expanduser('~/.seqc/tools/human'))
        print('Available GSEA .gmt files:\n\nmouse:\n{m}\n\nhuman:\n{h}\n'.format(
                m='\n'.join(mouse_options),
                h='\n'.join(human_options)))
        print('Please specify the gmt_file parameter as gmt_file=(organism, filename)')

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
            components = range(1, self.diffusion_eigenvectors.shape[1])

        # Run for each component
        for c in components:

            print(c)

            # Write genes to file
            genes_file = out_dir + out_prefix + '_cmpnt_%d.rnk' % c
            ranked_genes = self.diffusion_map_correlations.iloc[:, c]\
                .sort_values(inplace=False, ascending=False)
            pd.DataFrame(ranked_genes).to_csv(genes_file, sep='\t', header=False)

            # Construct the GSEA call
            cmd = list()
            cmd += ['java', '-cp', os.path.expanduser('~/.seqc/tools/gsea2-2.2.1.jar'),
                    '-Xmx1G', 'xtools.gsea.GseaPreranked']
            cmd += ['-collapse false', '-mode Max_probe', '-norm meandiv', '-nperm 1000']
            cmd += ['-include_only_symbols true', '-make_sets true', '-plot_top_x 20']
            cmd += ['-set_max 500', '-set_min 15', '-zip_report false -gui false']

            # Input arguments
            cmd += ['-rnk %s' % genes_file]
            cmd += ['-rpt_label %s' % out_prefix]
            cmd += ['-out %s' % out_dir]
            cmd += ['-gmx %s' % gmt_file]

            # Call GSEA
            call(cmd)

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
        subset_molecules = self.molecules.ix[use_genes]
        subset_reads = self.reads.ix[use_genes]
        metadata = self.metadata.ix[use_genes]

        # throws out analysis results; this makes sense
        return Experiment(subset_reads, subset_molecules, metadata)

    # todo implement
    def differential_expression(self, g1, g2):
        """
        carry out differential expression between cells g1 and cells g2
        :param g1:
        :param g2:
        :return:
        """
        raise NotImplementedError
