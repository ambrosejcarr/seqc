import multiprocessing
import os
import pickle
import random
import re
import shlex
import shutil
import warnings
from collections import defaultdict, OrderedDict
from copy import deepcopy, copy
from functools import partial
from itertools import combinations
from subprocess import call, Popen, PIPE
import matplotlib
import numpy as np
import pandas as pd
import tables as tb
from numpy.linalg import norm
from scipy.sparse import coo_matrix
from scipy.sparse import csr_matrix, find
from scipy.sparse.linalg import eigs
from scipy.stats import gaussian_kde, pearsonr, mannwhitneyu, rankdata
from scipy.stats.mstats import kruskalwallis
from sklearn import linear_model
from sklearn.manifold import TSNE
from sklearn.neighbors import NearestNeighbors
from statsmodels.sandbox.stats.multicomp import multipletests

try:
    os.environ['DISPLAY']
except KeyError:
    matplotlib.use('Agg')
import matplotlib.pyplot as plt
with warnings.catch_warnings():
    warnings.simplefilter('ignore')  # catch experimental ipython widget warning
    import seaborn as sns
import phenograph
import seqc
from seqc.exceptions import SparseMatrixError

# set plotting defaults
sns.set_style('ticks')
matplotlib.rc('font', **{
    'serif': ['Computer Modern Roman'],
    'monospace': ['Inconsolata'],
    'sans-serif': ['Lato']
    })

# todo change to kwargs
matplotlib.rc('figure', **{'figsize': (4, 4), 'dpi': 150})

matplotlib.rc('patch', **{'facecolor': 'royalblue', 'edgecolor': 'none'})

matplotlib.rc('lines', **{'color': 'royalblue', 'markersize': 7})

matplotlib.rc('savefig', **{'dpi': '150'})

matplotlib.rc('image', **{'cmap': 'viridis'})

# todo remove
size=7
cmap=plt.cm.viridis


def qualitative_colors(n):
    return sns.color_palette('rainbow', n)


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


class FigureGrid:

    def __init__(self, n: int, max_cols=3):
        self.n = n
        self.nrows = int(np.ceil(n / max_cols))
        self.ncols = int(min((max_cols, n)))
        figsize = self.ncols * 4, self.nrows * 4

        # create figure
        self.gs = plt.GridSpec(nrows=self.nrows, ncols=self.ncols)
        self.figure = plt.figure(figsize=figsize)

        # create axes
        self.axes = {}
        for i in range(n):
            row = int(i // self.ncols)
            col = int(i % self.ncols)
            self.axes[i] = plt.subplot(self.gs[row, col])

    def __getitem__(self, item):
        return self.axes[item]

    def __iter__(self):
        for i in range(self.n):
            yield self[i]

    def tight_layout(self):
        self.gs.tight_layout(self.figure)

    def despine(self, top=True, right=True, **kwargs):
        for i in range(self.n):
            sns.despine(ax=self[i], top=top, right=right, **kwargs)

    def savefig(self, *args, **kwargs):
        self.figure.savefig(*args, **kwargs)


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
        reader = seqc.alignment.sam.Reader(samfile)
        translator = seqc.sequence.gtf.GeneIntervals(gtf)
        num_records = 0

        data = np.recarray((len(reader),), cls._dtype)
        for i, alignment in enumerate(reader):
            num_records += 1
            gene = translator.translate(
                    alignment.rname, alignment.strand, alignment.pos)
            if gene is None:
                gene = 0
            pool = seqc.sequence.encodings.DNA3Bit.encode(alignment.pool)
            cell = seqc.sequence.encodings.DNA3Bit.encode(alignment.cell)
            rmt = seqc.sequence.encodings.DNA3Bit.encode(alignment.rmt)
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

        # filter out N's in cell barcode and rmt
        res = np.zeros(len(data), dtype=np.bool)
        cell = data['cell'].copy()
        rmt = data['rmt'].copy()
        while np.any(rmt):
            n_filter = rmt & 0b111 == 0b111
            res[n_filter] = True
            rmt >>= 3

        while np.any(cell):
            n_filter = cell & 0b111 == 0b111
            res[n_filter] = True
            cell >>= 3
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

        self._data = data  #.astype(np.uint32)  # was causing casting errors; molecules should be float32!
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


class Experiment:
    """
    Object to hold and process read and molecule information from a SEQC experiment.

    :property reads: a SparseFrame or DataFrame containing a cells x genes matrix of
      read counts.
    :property molecules: a SparseFrame or DataFrame containing a cells x genes matrix of
      molecule counts.
    :property metadata: a pd.DataFrame of cells x metadata; each column represents a
      categorical variable containing metadata pertaining to the cells of Experiment.
    :property pca: pd.DataFrame or None, a stored PCA decomposition of Experiment.
      Produced by calling #todo
    :property tsne: pd.DataFrame or None, a stored tSNE projection of Experiment.
      Produced by calling #todo
    :property diffusion_map_eigenvectors: pd.DataFrame or None, eigenvectors of a stored
      diffusion map decompsition of Experiment
    :property diffusion_map_eigenvalues: pd.DataFrame or None, eigenvalues of a stored
      diffusion map decompsition of Experiment
    :property diffusion_map_correlations: pd.DataFrame or None, stored correlations of
      each gene with each diffusion component
    :property cluster_assignments: pd.Series or None, stored phenograph cluster
      assignments.

    -- IO methods --
    :method from_count_matrices: Constructs an Experiment object from a pickled SEQC run
      output "*_read_and_count_matrices.p" file.
    :method save: save a serialized python pickle representation of the Experiment that
      can quickly be loaded at a later time.
    :method load: load a serialized python pickle representation of the Experiment.
    :method convert_ensemble: convert self.index containing scids into an index of gene
      names
    :method remove_non_cell_barcodes: removes low abundance cell barcodes
    :staticmethod merge: Merge multiple experiment objects together, merging duplicate
      cell values.
    :staticmethod concatenate: Merge multiple experiment objects together, treating
      duplicate cell ids as unique molecules and distinguishing them with numeric
      suffixes.
    :staticmethod ensembl_gene_id_to_official_gene_symbol: map integer ENSEMBL ids to
      official gene symbols, optionally using a pre-created map.

    -- Filtering methods --
    :method create_filtered_dense_count_matrix:  Applies several filters to an Experiment
      object to carry out an unsupervised conversion of a SparseFrame representation to a
      reduced dense pd.DataFrame representation. This is the final step of seqc and is
      carried out automatically.
    :method exclude_dead_cells_with_high_mt_fraction: remove cells with mitochondrial
      content suggestive of a dead or apoptotic cell
    :method remove_low_complexity_cells: Remove any cells for which the residual of the
      linear fit between molecules and genes detected is greater than .15.
    :method remove_clusters_based_on_library_size_distribution: # throws out low
      library size clusters; normally there is only one, but worth visual inspection.
    :method remove_housekeeping_genes: remove housekeeping genes from Experiment

    -- Statistical methods --
    :method run_pca: runs Van Der Maaten PCA on Experiment and stores the decomposition
      in self.pca
    :method run_tsne: runs tsne on a PCA-reduced projection of experiment
    :method run_phenograph: runs phenograph on a PCA-reduced projection of experiment
    :method run_diffusion_map: runs diffusion maps on a PCA-reduced projection of
      Experiment.
    :method run_diffusion_map_correlations: correlate a neighbor-based rolling mean for
      each gene with each diffusion map correlation.
    :method run_gsea: Run GSEA on diffusion components
    :method run_gsea_preranked_list: Helper function. Run GSEA on an already-ranked list
      of corrleations
    :method select_biological_components: Reports the top n enrichments with adjusted
      p-value < alpha for selected diffusion component GSEA enrichments.
    :method select_genes_from_diffusion_components: Carry out feature selection on genes,
      retaining only genes having significant correlation with a selected component.
    :method pairwise_differential_expression: carry out differential expression between
      two subsets of cells in self.
    :method single_gene_differential_expression: carry out non-parametric ANOVA across
      cluster identifiers for a single gene.
    :method differential_expression: carry out non-parametric ANOVA across cluster
      identifiers for all genes.
    :method annotate_clusters_by_expression: Return genes which are uniquely expressed
      in each cluster. Susceptible to differences in sampling and improper normalization.

    -- Plotting methods --
    :method plot_molecules_vs_reads_per_molecule: plot log10 molecule counts per barcode
      vs average cell coverage.
    :method plot_mitochondrial_molecule_fraction: plot the fraction of mRNA that are of
      mitochondrial origin for each cell.
    :method plot_molecules_vs_genes: display the relationship between molecules captured
      and number of genes captured.
    :method plot_pca_variance_explained: Plot the variance explained by each PCA
      component.
    :method plot_tsne: plot the tsne projection, coloring according to c. If c is not
      provided, default coloration is based on cell density.
    :method plot_tsne_by_metadata: plot tsne, coloring the projection by metadata
    :method plot_tsne_by_cell_sizes: plot tsne, coloring the projection by cell size
    :method plot_clusters: plot phenograph clusters
    :method plot_diffusion_maps: Visualize diffusion component loadings on the tSNE map
    :method plot_gene_component_correlations: Plots gene-component correlations for
      selected components
    :method plot_gene_expression: display gene expression values on the tSNE projection
      of experiment. One plot per provided gene.
    :method plot_aggregate_gene_expression: Display summed expression of a gene set on
      a tSNE projection.

    -- Miscellaneous methods --
    :method is_sparse: True if any data is in sparse format.
    :method is_dense: True if no data is in sparse format.
    :method normalize_data: Divide by cell size and multiply by median values. Sets
      self.normalized=True.


    """

    def __init__(self, molecules, reads=None, metadata=None):
        """
        Object to hold and process read and molecule information from a SEQC experiment.

        :param reads:  DataFrame or SparseFrame
        :param molecules: DataFrame or SparseFrame
        :param metadata: None or DataFrame
        :return:
        """
        if not (isinstance(molecules, SparseFrame) or
                isinstance(molecules, pd.DataFrame)):
            raise TypeError('molecules must be of type SparseFrame or DataFrame')
        if not (isinstance(reads, SparseFrame) or isinstance(reads, pd.DataFrame) or
                reads is None):
            raise TypeError('reads must be of type SparseFrame or DataFrame')
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

    def create_filtered_dense_count_matrix(self, max_mt_content=0.2):
        """
        filter cells with low molecule counts, low read coverage, high mitochondrial content,
        and low gene detection. Returns a dense pd.DataFrame of filtered counts, the total
        original number of molecules (int), the number of molecules lost with each filter
        (dict), and the number of cells lost with each filter (dict).

        :param max_mt_content: the maximum percentage of mitochondrial RNA that is
          considered to constitute a viable cell
        :return: (pd.DataFrame, int, dict, dict)
        """

        cells_lost = OrderedDict()
        molecules_lost = OrderedDict()

        if not self.molecules.columns.dtype.char == 'U':
            if self.molecules.sum().sum() == 0:
                raise seqc.exceptions.EmptyMatrixError(
                    'Matrix is empty, cannot create dense matrix')
            else:
                raise RuntimeError(
                    'non-string column names detected. Please convert column names into '
                    'string gene symbols before calling this function.')
        if not isinstance(max_mt_content, float):
            raise TypeError('Parameter max_mt_content must be of type float.')
        if not 0 <= max_mt_content <= 1:
            raise ValueError('Parameter max_mt_content must be in the interval [0, 1]')

        # set data structures and original molecule counts
        M = self.molecules.data
        R = self.reads.data
        G = self.molecules.columns
        is_invalid = np.zeros(M.shape[0], np.bool)
        total_molecules = np.sum(M.sum(axis=1))

        def additional_loss(new_filter, old_filter, D):
            new_cell_loss = np.sum(new_filter) - np.sum(old_filter)
            D = D.tocsr()
            total_molecule_loss = D[new_filter].sum().sum()
            old_molecule_loss = D[old_filter].sum().sum()
            new_molecule_loss = total_molecule_loss - old_molecule_loss
            return new_cell_loss, new_molecule_loss

        # filter low counts
        count_invalid = seqc.filter.low_count(M, is_invalid)
        cells_lost['low_count'], molecules_lost['low_count'] = additional_loss(
            count_invalid, is_invalid, M)

        # filter low coverage
        cov_invalid = seqc.filter.low_coverage(M, R, count_invalid)
        cells_lost['low_coverage'], molecules_lost['low_coverage'] = additional_loss(
            cov_invalid, count_invalid, M)

        # filter high_mt_content
        mt_invalid = seqc.filter.high_mitochondrial_rna(M, G, cov_invalid, max_mt_content)
        cells_lost['high_mt'], molecules_lost['high_mt'] = additional_loss(
            mt_invalid, cov_invalid, M)

        # filter low gene abundance
        gene_invalid = seqc.filter.low_gene_abundance(M, mt_invalid)
        cells_lost['low_gene_detection'], molecules_lost['low_gene_detection'] = additional_loss(
            gene_invalid, mt_invalid, M)

        # construct dense matrix
        dense = M.tocsr()[~gene_invalid, :].todense()
        nonzero_gene_count = np.ravel(dense.sum(axis=0) != 0)
        dense = dense[:, nonzero_gene_count]
        dense = Experiment(molecules=pd.DataFrame(
            dense,
            index=self.molecules.index[~gene_invalid],
            columns=self.molecules.columns[nonzero_gene_count]))

        # describe cells
        cell_description = dense.molecules.sum(axis=1).describe()

        return dense, total_molecules, molecules_lost, cells_lost, cell_description

    def save(self, fout: str) -> None:
        """
        save a serialized python pickle representation of the Experiment that
          can quickly be loaded at a later time.

        :param fout: str, name of archive to store pickled Experiment data in. Should end
          in '.p'.
        :return: None
        """
        with open(fout, 'wb') as f:
            pickle.dump(vars(self), f)

    @classmethod
    def load(cls, fin):
        """
        load a serialized python pickle representation of the Experiment.

        :param fin: str, name of pickled archive containing Experiment data
        :return: Experiment
        """
        with open(fin, 'rb') as f:
            data = pickle.load(f)
        experiment = cls(molecules=data['_molecules'], reads=data['_reads'],
                         metadata=data['_metadata'])
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
                if v is None:
                    _repr += '\n{}={}'.format(k[1:], 'None')
                elif v is False:
                    _repr += '\n{}={}'.format(k[1:], 'False')
                else:
                    _repr += '\n{}={}'.format(k[1:], 'True')
        return _repr

    @property
    def reads(self):
        return self._reads

    @reads.setter
    def reads(self, item):
        if not (isinstance(item, SparseFrame) or isinstance(item, pd.DataFrame)
                or item is None):
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

    @property
    def is_normalized(self):
        return self._normalized

    @classmethod
    def from_count_matrices(cls, read_and_count_matrix: str):
        """
        Construct an Experiment object from a SEQC reads_and_counts_matrix.p file

        :param read_and_count_matrix: str, filename of the pickled SEQC output
        :return:
        """

        with open(read_and_count_matrix, 'rb') as f:
            d = pickle.load(f)
        reads = SparseFrame(d['reads']['matrix'], index=d['reads']['row_ids'],
                            columns=d['reads']['col_ids'])
        molecules = SparseFrame(d['molecules']['matrix'], index=d['molecules']['row_ids'],
                                columns=d['molecules']['col_ids'])
        return cls(reads=reads, molecules=molecules)

    def is_sparse(self):
        """Return True if any of the Experiment's data is in sparse format"""
        if all(isinstance(o, SparseFrame) for o in (self.reads, self.molecules)):
            return True
        else:
            return False

    def is_dense(self):
        """Return True if none of the Experiment's data is in sparse format"""
        return not self.is_sparse()

    # todo may still be side-effecting
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
        molecules = coo_matrix((molecule_data, (row, col)), shape=shape, dtype=np.float32)

        sparse_reads = SparseFrame(reads, cell_index, gene_index)
        sparse_molecules = SparseFrame(molecules, cell_index, gene_index)

        return Experiment(reads=sparse_reads, molecules=sparse_molecules)

    # todo
    # currently, this method has side-effects (mutates existing row_ids, dataframes.
    # Need more copying, if want to avoid)
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

        # todo are these getting mutated twice??
        # mutate cell identifiers to ensure they are unique after merging
        old_index = []
        for i, e in enumerate(experiments):
            old_index.append(e.molecules.index)
            suffix = '_{!s}'.format(i)
            new_cells = [str(c) + suffix for c in e.molecules.index]
            e.molecules.index = new_cells
            e.metadata.index = new_cells

        # merge genes, create new cell list
        genes = set()
        cells = []
        meta_cols = set()
        for e in experiments:
            genes.update(list(e.molecules.columns))
            cells.extend(deepcopy(list(e.molecules.index)))
            meta_cols.update(list(e.metadata.columns))  # list not necessary?

        # combine data
        empty_molc = np.zeros((len(cells), len(genes)), dtype=np.float32)
        metadata = pd.DataFrame(index=cells, columns=meta_cols)
        m_combined = pd.DataFrame(empty_molc, index=cells, columns=genes)
        for e in experiments:
            m_combined.ix[e.molecules.index, e.molecules.columns] = e.molecules
            metadata.ix[e.metadata.index, e.metadata.columns] = e.metadata

        # add additional metadata
        for k, v in metadata_labels.items():
            metadata[k] = pd.Series(index=cells, dtype='O')
            for e, metadata_val in zip(experiments, v):
                metadata[k].ix[e.metadata.index] = [metadata_val] * len(e.metadata.index)

        # regenerate original columns in experiemnt
        for i, e in zip(old_index, experiments):
            e.molecules.index = i
            e.metadata.index = i

        return Experiment(molecules=m_combined, reads=None, metadata=metadata)

    @staticmethod
    def fast_concatenate(experiments: list, labels: tuple=None):
        """
        Concatenate a set of Experiment objects or DataFrames. Each cell is considered to
        be unique, even if they share the same barcodes. An hierarchical multi-index is
        constructed on the rows, where labels serve as the first level to distinguish
        cells from multiple experiments that may have gotten the same cell barcode.

        Concatenate should be run after filtering each dataset individually, since the
        datasets are considered to have distinct cells (even if they share the same
        barcodes)

        :param experiments: list or iterable of Experiment or DataFrame objects
        :param labels: tuple, labels for experiments, length must match number of
          experiments.
        :return merged: Experiment object concatenated along the row axis (0).
        """

        if labels is None:
            labels = np.arange(len(experiments))

        if all(isinstance(e, Experiment) for e in experiments):
            data = pd.concat([e.molecules for e in experiments], keys=labels)
            metadata = pd.concat([e.molecules for e in experiments], keys=labels)
        elif all(isinstance(e, pd.DataFrame) for e in experiments):
            data = pd.concat(experiments, keys=labels)
            metadata = None
        else:
            raise TypeError('experiments must be pandas DataFrame or Experiment objects.')

        return Experiment(molecules=data, reads=None, metadata=metadata)

    # todo @ambrosejcarr make this faster; it is way too slow.
    @staticmethod
    def create_gene_id_to_official_gene_symbol_map(gtf: str):
        """
        create_gene_id_to_official_gene_symbol_map: map integer ENSEMBL ids to
        official gene symbols.

        :param gtf: str, filename of gtf file from which to create the map.
        """
        pattern = re.compile(
            r'(^.*?gene_id "[^0-9]*)([0-9]*)(\.?.*?gene_name ")(.*?)(".*?$)')
        gene_id_map = defaultdict(set)
        with open(gtf, 'r') as f:
            for line in f:
                match = re.match(pattern, line)
                if match:
                    gene_id_map[int(match.group(2))].add(match.group(4).upper())
        return gene_id_map

    def ensembl_gene_id_to_official_gene_symbol(self, gtf=None, gene_id_map=None):
        """convert self.index containing scids into an index of gene names

        :param gtf: str, filename of a gtf file
        :param gene_id_map: gene_id_map constructed from
          Experiment.create_gene_id_to_official_gene_symbol_map. If converting multiple
          objects, it is much faster to only construct the map a single time.
        :return: Experiment
        """
        if gene_id_map is None:
            if gtf is None:
                raise ValueError('User must pass either GTF or a gene_id_map object')
            gene_id_map = self.create_gene_id_to_official_gene_symbol_map(gtf)
        if self.reads is not None:
            reads = deepcopy(self.reads)
            reads.columns = ['-'.join(gene_id_map[i]) for i in self.reads.columns]
        else:
            reads = None
        molecules = deepcopy(self.molecules)
        molecules.columns = ['-'.join(gene_id_map[i]) for i in self.molecules.columns]

        exp = Experiment(molecules=molecules, reads=reads, metadata=self.metadata)
        exp._normalized = self._normalized
        return exp

    def scid_to_official_gene_symbol(self, gtf: str):
        """
        Deprecated. convert scids to official gene symbols; useful for older v0.1.2 data.

        :param gtf: str gtf filename.
        :return:
        """
        warnings.warn('DeprecationWarning: This function is only useful for data created '
                      'with SEQC v0.1.2. You should re-process your data to take '
                      'advantage of new tools present in modern versions of SEQC')

        pattern = re.compile(
            r'(^.*?gene_name ")(.*?)(".*?scseq_id "SC)(.*?)(".*?$)')

        gene_id_map = defaultdict(set)
        with open(gtf, 'r') as f:
            for line in f:
                match = re.match(pattern, line)
                if match:
                    gene_id_map[int(match.group(4))].add(match.group(2))

        if self.reads is not None:
            reads = deepcopy(self.reads)
            reads.columns = ['-'.join(gene_id_map[i]) for i in self.reads.columns]
        else:
            reads = None
        molecules = deepcopy(self.molecules)
        molecules.columns = ['-'.join(gene_id_map[i]) for i in self.molecules.columns]

        exp = Experiment(molecules=molecules, reads=reads, metadata=self.metadata)
        exp._normalized = self._normalized
        return exp

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
            raise RuntimeError('remove_non_cell_barcodes() should be run on '
                               'unnormalized data')

        # get molecule and read counts per cell
        molecule_cell_sums = pd.Series(np.ravel(self.molecules.sum(axis=1)),
                                       index=self.molecules.index, dtype=np.float32)
        read_cell_sums = pd.Series(np.ravel(self.reads.sum(axis=1)),
                                   index=self.reads.index, dtype=np.uint32)

        ratios = read_cell_sums / molecule_cell_sums

        # select cells & convert to dense
        filter_ = (molecule_cell_sums >= min_mpc) & (ratios >= min_rpm)
        codes_passing_filter = molecule_cell_sums.index[filter_]

        row_inds = np.where(pd.Series(self.molecules.index).isin(codes_passing_filter))[0]

        molecule_data = self.molecules.data.tocsr()[row_inds, :].todense()
        molecules = pd.DataFrame(molecule_data, index=self.molecules.index[row_inds],
                                 columns=self.molecules.columns)

        metadata = self.metadata.ix[row_inds, :]
        exp = Experiment(molecules=molecules, reads=None, metadata=metadata)
        exp._normalized = self._normalized
        return exp

    def plot_mitochondrial_molecule_fraction(
            self, fig=None, ax=None, title='Dead Cell Identification Plot'):
        """ plot the fraction of mRNA that are of mitochondrial origin for each cell.

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

        :param max_mt_fraction: float, maximum percentage of mRNA that can come from the
          mitochondria before the cell is considered invalid.
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
        metadata = self.metadata.ix[pass_filter]

        exp = Experiment(molecules=molecules, reads=None, metadata=metadata)
        exp._normalized = self._normalized
        return exp

    def plot_molecules_vs_genes(self, fig=None, ax=None, title=''):
        """
        plot the (expected linear) relationship between number of molecules captured and
        the number of genes detected

        :param title: optional title to plot on axis
        :param ax: matplotlib axis
        :param fig: matplotlib figure
        :return: fig, ax
        """
        if self._normalized:
            raise RuntimeError('plot_molecules_vs_reads_per_molecule() should be run on '
                               'unnormalized data')

        fig, ax = get_fig(fig=fig, ax=ax)
        molecule_counts = self.molecules.sum(axis=1)
        gene_counts = np.sum(self.molecules > 0, axis=1)
        x, y, z = density_2d(np.log10(molecule_counts), np.log10(gene_counts))
        ax.scatter(x, y, c=z, edgecolor='none', s=size, cmap=cmap, alpha=0.65)
        ax.set_xlabel('log10(Number of molecules)')
        ax.set_ylabel('log10(Number of genes)')
        ax.set_title(title)
        sns.despine(ax=ax)

        # get line of best fit
        regr = linear_model.LinearRegression()
        regr.fit(x[:, np.newaxis], y)

        # predicted y
        yhat = regr.predict(x[:, np.newaxis])

        # plot fit line
        ax.plot(x, yhat, color='r', linewidth=1)

        # plot residuals that should be removed in red
        residuals = yhat - y
        remove = residuals > .15
        ax.scatter(x[remove], y[remove], color='r', s=12, alpha=0.65)

        return fig, ax

    def remove_low_complexity_cells(self):
        """Remove any cells for which the residual of the linear fit between molecules
        and genes detected is greater than .15.

        See plot_molecules_vs_genes for visualization of this fit.
        """
        if self._normalized:
            raise RuntimeError('plot_molecules_vs_reads_per_molecule() should be run on '
                               'unnormalized data')

        molecule_counts = self.molecules.sum(axis=1)
        gene_counts = np.sum(self.molecules > 0, axis=1)
        x = np.log10(molecule_counts)[:, np.newaxis]
        y = np.log10(gene_counts)

        # get line of best fit
        regr = linear_model.LinearRegression()
        regr.fit(x, y)

        # predicted y
        yhat = regr.predict(x)

        # plot residuals that should be removed in red
        residuals = yhat - y
        retain = residuals < .15

        self.molecules = self.molecules.ix[retain, :]
        self.metadata = self.metadata.ix[retain, :]

    def normalize_data(self):
        """
        uses AVO method of normalization; stores original library sizes. Sets
        self._normalized = True

        :return: Experiment
        """

        if self.is_sparse():
            raise RuntimeError('Please convert to dense data before attempting to '
                               'normalize using Experiment.remove_non_cell_barcodes().')

        molecule_sums = self.molecules.sum(axis=1)
        molecules = self.molecules.div(molecule_sums, axis=0)\
            .mul(np.median(molecule_sums), axis=0)
        exp = Experiment(molecules=molecules, reads=None, metadata=self.metadata)
        exp._normalized = True

        # check that none of the genes are empty; if so remove them
        nonzero_genes = exp.molecules.sum(axis=0) != 0
        exp.molecules = exp.molecules.ix[:, nonzero_genes]

        return exp

    def merge_duplicate_genes(self):
        """only necessary for v0.1.2 data where scids are duplicated

        :return:
        """
        self.molecules = self.molecules.groupby(axis=1).sum()
        if self.reads is not None:
            self.reads = self.reads.groupby(axis=1).sum()

    def run_pca(self, n_components=100):
        """
        run Van Der Maaten PCA on Experiment and stores the decomposition in self.pca

        :param n_components: number of PCA components to store.
        :return:
        """

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

        loadings = pd.DataFrame(data=M, index=self.molecules.columns)
        l = pd.DataFrame(l)

        self.pca = {'loadings': loadings, 'eigenvalues': l}

    def run_pca_matlab(self, no_components=100):
        """
        Deprecated. Runs a Matlab implementation of PCA. Requires Matlab.

        :param no_components: number of PCA components to retain.
        :return:
        """
        warnings.warn('DeprecationWarning: please use Experiment.run_pca()')

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

    def plot_pca_variance_explained(self, fig=None, ax=None, n_components=30):
        """
        Plot the variance explained by each PCA component. Requires that self.run_pca()
         has been execued.

        :param fig: matplotlib Figure
        :param ax: matplotlib Axis
        :param n_components: number of PCA components to plot
        """

        fig, ax = get_fig(fig=fig, ax=ax)
        ax.plot(np.ravel(self.pca['eigenvalues'].values))
        ax.set_ylim((0, float(np.max(self.pca['eigenvalues']))))
        ax.set_xlim((0, n_components))
        ax.set_title('PCA variance by component')
        ax.set_xlabel('Component')
        ax.set_ylabel('Loading')
        sns.despine(ax=ax)
        return fig, ax

    def run_tsne(self, n_components=15, **kwargs):
        """
        run tnse on a pca-reduced projection of experiment

        :param n_components: number of PCA components to use as input for tsne
        :return:
        """
        data = deepcopy(self.molecules)
        data -= np.min(np.ravel(data))
        data /= np.max(np.ravel(data))
        data = pd.DataFrame(np.dot(data, self.pca['loadings'].iloc[:, 0:n_components]),
                            index=self.molecules.index)
        tsne = TSNE(n_components=2, method='barnes_hut', n_iter=1000, perplexity=30,
                    angle=0.5, init='pca', random_state=0, **kwargs)
        res = tsne.fit_transform(data.values)
        self.tsne = pd.DataFrame(res, index=self.molecules.index, columns=['x', 'y'])

    def plot_tsne_by_metadata(self, label, fig=None, ax=None):
        """
        plot tsne, coloring the projection by metadata column "label"

        :param label: the metadata column to color the tsne projection with
        :param fig: Matplotlib Figure
        :param ax: Matplotlib Axis
        :return: fig, ax
        """

        cats = set(self.metadata[label])
        colors = qualitative_colors(len(cats))
        fig, ax = get_fig(fig=fig, ax=ax)
        for c, color in zip(cats, colors):
            rows = self.metadata[label] == c
            plt.plot(self.tsne.ix[rows, 'x'], self.tsne.ix[rows, 'y'], marker='o',
                     markersize=np.sqrt(7), linewidth=0, c=color, label=c)
        ax.xaxis.set_major_locator(plt.NullLocator())
        ax.yaxis.set_major_locator(plt.NullLocator())
        ax.set_title('tSNE projection colored by {}'.format(label))
        ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), markerscale=3)
        return fig, ax

    def plot_tsne(self, fig=None, ax=None, c=None, **kwargs):
        """
        plot the tsne projection, coloring according to c. If c is not provided, default
         coloration is based on cell density.

        :param fig: Matplotlib Figure
        :param ax: Matplotlib Axis
        :param c: numerical vector used to color the cells of experiment
        :param kwargs: additional key word arguments for plt.scatter
        :return: fig, ax
        """
        fig, ax = get_fig(fig=fig, ax=ax)

        if c is not None:  # get colormap
            if 'cmap' in kwargs.keys():
                cmap = getattr(plt.cm, kwargs['cmap'])
            else:
                cmap = plt.get_cmap()

        if c is None:
            x, y, z = density_2d(self.tsne['x'], self.tsne['y'])
        else:
            x, y, z = self.tsne['x'], self.tsne['y'], c

        plt.scatter(x, y, c=z, **kwargs)
        ax.xaxis.set_major_locator(plt.NullLocator())
        ax.yaxis.set_major_locator(plt.NullLocator())
        plt.colorbar()
        return fig, ax

    def plot_tsne_by_cell_sizes(self, fig=None, ax=None, vmin=None, vmax=None,
                                title='', sizes=None):
        """
        plot tsne, coloring the projection by cell size

        :param fig: Matplotlib Figure
        :param ax: Matplotlib Axis
        :param vmin: min value for min color
        :param vmax: max value for max color
        :param title: optional title to plot on axis
        :param sizes: library sizes (calculated from data if not provided)
        :return: (fig, ax)
        """
        fig, ax = get_fig(fig, ax)
        if self.tsne is None:
            raise RuntimeError('Please run self.run_tsne() before plotting.')
        if sizes is None:
            sizes = self.molecules.sum(axis=1)  # use current library sizes
        plt.scatter(self.tsne['x'], self.tsne['y'], s=7, c=sizes, cmap=cmap, vmax=vmax,
                    vmin=vmin, alpha=0.65)
        ax.xaxis.set_major_locator(plt.NullLocator())
        ax.yaxis.set_major_locator(plt.NullLocator())
        ax.set_title(title)
        plt.colorbar()
        return fig, ax

    def run_phenograph(self, n_pca_components=15, **kwargs):
        """
        Run phenograph clustering on PCA-reduced cells of experiment.

        one useful kwargs is 'k', which inversely correlates with the number of
        clusters that phenograph outputs.

        :param n_pca_components: number of PCA components to consider when projecting
          experiment data.
        :param **kwargs: keyword arguments to pass directly to phenograph

        :return:
        """
        data = deepcopy(self.molecules)
        data -= np.min(np.ravel(data))
        data /= np.max(np.ravel(data))
        data = pd.DataFrame(np.dot(data, self.pca['loadings'].iloc[:, 0:n_pca_components]),
                            index=self.molecules.index)

        communities, graph, Q = phenograph.cluster(data, **kwargs)
        self.cluster_assignments = pd.Series(communities, index=data.index)

    def plot_clusters(self, fig=None, ax=None, labels=None, **kwargs):
        """
        plot phenograph clusters

        :param fig: Matplotlib Figure
        :param ax: Matplotlib Axis
        :param labels: cluster labels
        :param kwargs: additional arguments to provide to scatter
        :return:
        """

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
                    markersize=np.sqrt(size), label=label, **kwargs)
        ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), markerscale=3)
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
        metadata = self.metadata.ix[retain, :]
        experiment = Experiment(molecules=deepcopy(molecules), reads=None,
                                metadata=deepcopy(metadata))
        experiment._cluster_assignments = deepcopy(self.cluster_assignments[retain])

        # remove any genes that now have zero expression
        nonzero_genes = experiment.molecules.sum(axis=0) != 0
        experiment.molecules = experiment.molecules.ix[:, nonzero_genes]
        experiment._normalized = self._normalized
        return experiment

    def run_diffusion_map(self, knn=20, n_diffusion_components=20,
                          n_pca_components=15, epsilon=1):
        """
        Decompose experiment using diffusion components

        :param knn: number of nearest neighbors. Default can be a bit low, use 60+ for
          large datasets.
        :param n_diffusion_components: number of componenets to calculate
        :param n_pca_components: The data is pre-reduced; this is the number of PCA
          components to calculate when doing that reduction
        :param epsilon: error term
        :return: None
        """


        data = deepcopy(self.molecules)
        data -= np.min(np.ravel(data))
        data /= np.max(np.ravel(data))
        data = pd.DataFrame(np.dot(data, self.pca['loadings'].iloc[:, 0:n_pca_components]),
                            index=self.molecules.index)

        # Nearest neighbors
        N = data.shape[0]
        nbrs = NearestNeighbors(n_neighbors=knn, algorithm='ball_tree').fit(data)
        distances, indices = nbrs.kneighbors(data)

        # Adjacency matrix
        rows = np.zeros(N * knn, dtype=np.int32)
        cols = np.zeros(N * knn, dtype=np.int32)
        dists = np.zeros(N * knn)
        location = 0
        for i in range(N):
            inds = range(location, location + knn)
            rows[inds] = indices[i, :]
            cols[inds] = i
            dists[inds] = distances[i, :]
            location += knn
        W = csr_matrix((dists, (rows, cols)), shape=[N, N])

        # Symmetrize W
        W = W + W.T

        # Convert to affinity (with selfloops)
        rows, cols, dists = find(W)
        rows = np.append(rows, range(N))
        cols = np.append(cols, range(N))
        dists = np.append(dists/(epsilon ** 2), np.zeros(N))
        W = csr_matrix((np.exp(-dists), (rows, cols)), shape=[N, N])

        # Create D
        D = np.ravel(W.sum(axis=1))
        D[D != 0] = 1 / D[D != 0]

        # Symmetric markov normalization
        D = csr_matrix((np.sqrt(D), (range(N), range(N))), shape=[N, N])
        P = D
        T = D.dot(W).dot(D)
        T = (T + T.T) / 2

        # Eigen value decomposition
        D, V = eigs(T, n_diffusion_components, tol=1e-4, maxiter=1000)
        D = np.real(D)
        V = np.real(V)
        inds = np.argsort(D)[::-1]
        D = D[inds]
        V = V[:, inds]
        V = P.dot(V)

        # Normalize
        for i in range(V.shape[1]):
            V[:, i] = V[:, i] / norm(V[:, i])
        V = np.round(V, 10)

        # Update object
        self.diffusion_eigenvectors = pd.DataFrame(V, index=self.molecules.index)
        self.diffusion_eigenvalues = pd.DataFrame(D)

    def plot_diffusion_components(self):
        """
        Visualize diffusion component loadings on the tSNE map

        :return:
        """

        if self.tsne is None:
            raise RuntimeError('Please run tSNE before plotting diffusion components.')

        height = int(2 * np.ceil(self.diffusion_eigenvalues.shape[0] / 5))
        width = 10
        fig = plt.figure(figsize=[width, height])
        n_rows = int(height / 2)
        n_cols = int(width / 2)
        gs = plt.GridSpec(n_rows, n_cols)
        axes = []

        for i in range(self.diffusion_eigenvectors.shape[1]):
            axes.append(plt.subplot(gs[i // n_cols, i % n_cols]))
            axes[i].scatter(self.tsne['x'], self.tsne['y'], c=self.diffusion_eigenvectors[i],
                        cmap=cmap, edgecolors='none', s=size, alpha=0.65)
            axes[i].set_axis_off()

            # get outliers, plot in red
            # mu = np.mean(self.diffusion_eigenvectors[i])
            # std = np.std(self.diffusion_eigenvectors[i])
            # emin = mu - 10 * std
            # emax = mu + 10 * std
            # is_outlier = (self.diffusion_eigenvectors[i] < emin) | (
            # self.diffusion_eigenvectors[i] > emax)
            # axes[i].plot(self.tsne.ix[is_outlier, 'x'], self.tsne.ix[is_outlier, 'y'],
            #              c='r', linewidth=0, marker='o',
            #              markersize=np.sqrt(size))

        return fig, axes

    def remove_diffusion_outliers(self):
        raise NotImplementedError

    @staticmethod
    def _correlation(x: np.array, vals: np.array):
        """
        private method, calculate pairwise correlation correlation between x and each
         column in vals

        :param x: vector to be correlated
        :param vals: matrix containing n columns to correlate with x
        :return: n-length vector of correlation coefficients.
        """

        x = x[:, np.newaxis]
        mu_x = x.mean()  # cells
        mu_vals = vals.mean(axis=0)  # cells by gene --> cells by genes
        sigma_x = x.std()
        sigma_vals = vals.std(axis=0)

        return ((vals * x).mean(axis=0) - mu_vals * mu_x) / (sigma_vals * sigma_x)

    def run_diffusion_map_correlations(self, components=None, no_cells=10):
        """
        correlate a neighbor-based rolling mean for each gene with each diffusion
        map correlation

        :param components: list or array of components to correlate with
        :param no_cells: number of cells over which to smooth expression values
        :return:
        """
        if components is None:
            components = np.arange(self.diffusion_eigenvectors.shape[1])
        else:
            components = np.array(components)
        components = components[components != 0]

        component_shape = self.diffusion_eigenvectors.shape[1]

        diffusion_map_correlations = np.empty((self.molecules.shape[1],
                                               self.diffusion_eigenvectors.shape[1]),
                                               dtype=np.float)

        for component_index in components:
            component_data = self.diffusion_eigenvectors.values[:, component_index]

            order = np.argsort(component_data)
            x = pd.rolling_mean(component_data[order], no_cells)[no_cells:]
            # assert x.shape == (cell_shape - no_cells,)

            # this fancy indexing will copy self.molecules
            vals = pd.rolling_mean(self.molecules.values[order, :], no_cells, axis=0)[no_cells:]
            # assert vals.shape == (cell_shape - no_cells, gene_shape)
            cor_res = self._correlation(x, vals)
            # assert cor_res.shape == (gene_shape,)
            diffusion_map_correlations[:, component_index] = self._correlation(x, vals)

        # this is sorted by order, need it in original order (reverse the sort)

        self.diffusion_map_correlations = pd.DataFrame(
                diffusion_map_correlations[:, components], index=self.molecules.columns,
                columns=components)

    def determine_gene_diffusion_correlations(self, components=None, no_cells=10):
        """Deprecated. Matlab-based run_diffusion_map_correlations. Requires Matlab."""
        warnings.warn('DeprecationWarning: please use '
                      'Experiment.run_diffusion_map_correlations(). It offers a 10x '
                      'speed-up')

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
        """ plots gene-component correlations for selected of components

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
        """
        Private method. identifies GMT files available for mouse or human genomes
        :return: str, file options
        """

        mouse_options = os.listdir(os.path.expanduser('~/.seqc/tools/mouse'))
        human_options = os.listdir(os.path.expanduser('~/.seqc/tools/human'))
        print('Available GSEA .gmt files:\n\nmouse:\n{m}\n\nhuman:\n{h}\n'.format(
                m='\n'.join(mouse_options),
                h='\n'.join(human_options)))
        print('Please specify the gmt_file parameter as gmt_file=(organism, filename)')

    @staticmethod
    def _gsea_process(c: int, diffusion_map_correlations: str, output_stem: str,
                      gmt_file: str):
        """
        Private method. function passed to multiprocessing.map() to carry out the
        GSEA enrichment for a single diffusion component.

        :param c: int, diffusion component to be tested
        :param diffusion_map_correlations: list of correlations
        :param output_stem: location for GSEA output
        :param gmt_file: location for .gmt file containing gene sets to be tested
        :return: None
        """


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

    @classmethod
    def run_gsea_preranked_list(cls, rnk_file: str, output_stem: str, gmt_file: str=None):
        """
        Helper function. Run GSEA on an already-ranked list of corrleations

        :param rnk_file: .rnk file containing correlations in ranked order according to
          an independent variable.
        :param output_stem: location for GSEA output
        :param gmt_file: location for .gmt file containing gene sets to be tested
        :return: None
        """

        if output_stem.find('-') > -1:
            raise ValueError('ouptput_stem cannot contain special characters such as -')

        out_dir, out_prefix = os.path.split(output_stem)
        os.makedirs(out_dir, exist_ok=True)

        if not gmt_file:
            cls._gmt_options()
            return
        else:
            if not len(gmt_file) == 2:
                raise ValueError('gmt_file should be a tuple of (organism, filename).')
            else:
                gmt_file = os.path.expanduser('~/.seqc/tools/{}/{}').format(*gmt_file)

        # Construct the GSEA call
        cmd = shlex.split(
            'java -cp {user}/.seqc/tools/gsea2-2.2.1.jar -Xmx1g '
            'xtools.gsea.GseaPreranked -collapse false -mode Max_probe -norm meandiv '
            '-nperm 1000 -include_only_symbols true -make_sets true -plot_top_x 0 '
            '-set_max 500 -set_min 50 -zip_report false -gui false -rnk {rnk} '
            '-rpt_label {out_prefix} -out {out_dir}/ -gmx {gmt_file}'
            ''.format(user=os.path.expanduser('~'), rnk=rnk_file, out_prefix=out_prefix,
                      out_dir=out_dir, gmt_file=gmt_file))

        # Call GSEA
        p = Popen(cmd, stderr=PIPE)
        _, err = p.communicate()

        # remove annoying suffix from GSEA
        if err:
            print(err.decode())
            return
        else:
            pattern = '{p}.GseaPreranked.[0-9]*'.format(p=out_prefix)
            repl = out_prefix
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

    # todo @ambrosjcarr review statistical assumptions
    def run_gsea_diffexpr(self, cluster_1, cluster_2, output_stem, gmt_file=None):
        """
        Warning: may contain bugs or improper statistical functions.
        Helper function to carry out GSEA based on differential expression between two
        clusters

        :param cluster_1, cluster_2:
        :param output_stem:  the file location and prefix for the output of GSEA
        :param gmt_file:
        :param components:
        :return:
        """

        if output_stem.find('-') > -1:
            raise ValueError('ouptput_stem cannot contain special characters such as -')

        out_dir, out_prefix = os.path.split(output_stem)
        out_dir += '/'
        os.makedirs(out_dir, exist_ok=True)

        if not gmt_file:
            self._gmt_options()
            return
        else:
            if not len(gmt_file) == 2:
                raise ValueError('gmt_file should be a tuple of (organism, filename).')
            else:
                gmt_file = os.path.expanduser('~/.seqc/tools/{}/{}').format(*gmt_file)

        # save the .rnk file
        out_dir, out_prefix = os.path.split(output_stem)
        genes_file = '{stem}_cluster_{c1}_vs_{c2}.rnk'.format(
                stem=output_stem, c1=cluster_1, c2=cluster_2)

        # get cells in each group and convert to ranks
        c1 = self.molecules.ix[self.cluster_assignments == cluster_1, :]
        c2 = self.molecules.ix[self.cluster_assignments == cluster_2, :]
        c1 = np.apply_along_axis(rankdata, 1, c1)
        c2 = np.apply_along_axis(rankdata, 1, c2)

        # get average ranks, compute fold change
        c1_mu = pd.Series(np.mean(c1, axis=0), index=self.molecules.columns)
        c2_mu = pd.Series(np.mean(c2, axis=0), index=self.molecules.columns)
        fold_change = (c2_mu - c1_mu) / c1_mu

        # dump to file
        ranked_genes = fold_change.sort_values(inplace=False, ascending=False)
        pd.DataFrame(ranked_genes).fillna(0).to_csv(genes_file, sep='\t', header=False)

        # Construct the GSEA call
        cmd = shlex.split(
            'java -cp {user}/.seqc/tools/gsea2-2.2.1.jar -Xmx1g '
            'xtools.gsea.GseaPreranked -collapse false -mode Max_probe -norm meandiv '
            '-nperm 1000 -include_only_symbols true -make_sets true -plot_top_x 0 '
            '-set_max 500 -set_min 50 -zip_report false -gui false -rnk {rnk} '
            '-rpt_label {out_prefix}_{c1}_vs_{c2} -out {out_dir}/ -gmx {gmt_file}'
            ''.format(user=os.path.expanduser('~'), rnk=genes_file,
                      out_prefix=out_prefix, c1=cluster_1, c2=cluster_2, out_dir=out_dir,
                      gmt_file=gmt_file))

        # Call GSEA
        p = Popen(cmd, stderr=PIPE)
        _, err = p.communicate()

        # remove annoying suffix from GSEA
        if err:
            print(err.decode())
            return
        else:
            pattern = '{p}_{c1}_vs_{c2}.GseaPreranked.[0-9]*'.format(
                    p=out_prefix, c1=cluster_1, c2=cluster_2)
            repl = '{p}_{c1}_vs_{c2}'.format(p=out_prefix, c1=cluster_1, c2=cluster_2)
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

    def run_gsea(self, output_stem: str, gmt_file: str=None, components=None):
        """
        Run Gene Set Enrichment Analysis on diffusion component correlations

        :param output_stem: location for GSEA output
        :param gmt_file: location for .gmt file containing gene sets to be tested
        :param components: list or array of components to be tested
        :return:
        """
        if output_stem.find('-') > -1:
            raise ValueError('ouptput_stem cannot contain special characters such as -')

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

    def select_biological_components(self, output_stem: str, n: int=10, alpha: float=0.5,
                                     components=None):
        """
        Reports the top n enrichments with adjusted p-value < alpha for selected
        diffusion component GSEA enrichments.

        :param output_stem: str, location of GSEA output
        :param n: int, number of significant enrichments to report
        :param alpha: float, significance threshold
        :param components: list or array of components to be reported
        :return:
        """
        if components is None:
            components = self.diffusion_map_correlations.columns

        # get positive enrichments for each component
        for c in components:
            pos_report_pattern = 'gsea_report_for_na_pos_[0-9]*?\.html'
            neg_report_pattern = 'gsea_report_for_na_neg_[0-9]*?\.html'
            try:
                files = os.listdir('{stem}_{c}'.format(stem=output_stem, c=c))
            except FileNotFoundError:
                print('No enrichments were calculated for component {c}'.format(c=c))
                continue
            print('component: {}'.format(c))
            for f in files:
                mo_pos = re.match(pos_report_pattern, f)
                mo_neg = re.match(neg_report_pattern, f)
                if mo_pos:
                    try:
                        report_data = pd.read_html('{stem}_{c}/{f}'.format(
                            stem=output_stem, c=c, f=f), header=0, index_col=0)[0]
                    except ValueError:
                        print('could not find table for component {}'.format(c))
                        break
                    # report top n significant enrichments
                    significant = np.where(report_data['FDR q-val'] < alpha)[0][:n]
                    print('positive enrichment:')
                    print(list(report_data.iloc[significant, 0]))
                elif mo_neg:
                    try:
                        report_data = pd.read_html('{stem}_{c}/{f}'.format(
                            stem=output_stem, c=c, f=f), header=0, index_col=0)[0]
                    except ValueError:
                        print('could not find table for component {}'.format(c))
                        break
                    # report top n significant enrichments
                    significant = np.where(report_data['FDR q-val'] < alpha)[0][:n]
                    print('negative enrichment:')
                    print(list(report_data.iloc[significant, 0]))
            print('')

    def select_genes_from_diffusion_components(self, components, plot=False):
        """
        Carry out feature selection on genes, retaining only genes having significant
        correlation with a selected component.

        :param components: list or array of selected components
        :param plot: if True, visualize the distribution of correlations
        :return: Experiment, contains selected features
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
            cutoff = (np.mean(self.diffusion_map_correlations.ix[:, component]) +
                      2 * np.std(self.diffusion_map_correlations.ix[:, component]))
            use_genes = use_genes + list(self.diffusion_map_correlations.index[abs(
                    self.diffusion_map_correlations.ix[:, component]) > cutoff])

        # Unique genes
        use_genes = list(set(use_genes))

        # Create new scdata object
        subset_molecules = deepcopy(self.molecules.ix[:, use_genes])
        metadata = deepcopy(self.metadata)

        # throws out analysis results; this makes sense
        exp = Experiment(molecules=subset_molecules, reads=None,
                         metadata=metadata)
        exp._normalized = self._normalized

        # remove duplicate genes
        exp.molecules = exp.molecules.groupby(axis=1, level=0).sum()

        return exp

    def pairwise_differential_expression(self, c1, c2, alpha=0.05):
        """
        carry out differential expression (post-hoc tests) between cells c1 and cells c2,
        using bh-FDR to correct for multiple tests
        :param alpha: type I error rate
        :param c1, c2: arrays of cell identifiers to be compared

        :return: pval_sorted, fold_change
        """
        # get genes expressed in either set
        expressed = self.molecules.loc[c1 | c2].sum(axis=0) != 0

        g1 = self.molecules.loc[c1].ix[:, expressed]
        g2 = self.molecules.loc[c2].ix[:, expressed]

        g1mu = g1.mean(axis=0)
        g2mu = g2.mean(axis=0)
        fc = (g2mu - g1mu) / g1mu

        res = pd.Series(index=self.molecules.columns[expressed], dtype=float)
        for g in self.molecules.columns[expressed]:
            try:
                res[g] = mannwhitneyu(np.ravel(g1[g]), np.ravel(g2[g]))[1]
            except ValueError:
                res[g] = 1

        pval_corrected = multipletests(res.values, alpha=alpha, method='fdr_tsbh')[1]

        pval_sorted = pd.Series(pval_corrected, index=self.molecules.columns[expressed],
                                dtype=float).sort_values()
        fc = fc.ix[pval_sorted.index]
        return pval_sorted, fc

    def single_gene_differential_expression(self, gene, alpha):
        """
        carry out kruskal-wallis non-parametric (rank-wise) ANOVA with two-stage bh-FDR
        correction to determine if gene is differentially expressed in different clusters

        Explicitly, this tests the null hypothesis that all samples are drawn from the
        same distribution.

        if the KW-test is significant, Post-hoc rank-sum tests on the genes whose
        corrected p-values are below alpha determine the specific samples that are
        differentially expressed

        :param gene: str: gene name
        :param alpha: type I error.
        :return: KW p-val, pairwise ranksums pvals
        """
        gene = self.molecules.ix[:, gene]

        clusters = np.unique(self.cluster_assignments)
        samples = [gene.ix[self.cluster_assignments == c].values for
                   c in clusters]

        pval = kruskalwallis(*samples)[1]

        if pval < alpha:
            pairs = list(combinations(np.arange(len(samples)), 2))
            pvals = pd.Series(index=pairs, dtype=float)
            for a, b in pairs:
                pvals.ix[(a, b)] = mannwhitneyu(gene[self.cluster_assignments == a],
                                                gene[self.cluster_assignments == b])[1]
            return pval, pvals.sort_values()
        else:
            return pval, None

    def differential_expression(self, alpha: float=0.05):
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

        :param alpha: float, type I error rate (default=0.05)
        :return: pval_corrected
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
        reject, pval_corrected, _, _ = multipletests(pvals, alpha, method='fdr_tsbh')

        return pd.Series(pval_corrected, index=self.molecules.columns).sort_values(inplace=False)

    def plot_gene_expression(self, genes, log=False):
        """
        Display gene expression on the tSNE projection. One plot is made per gene.

        :param genes: array of string gene_ids to be displayed
        :param log: values should be log transformed (default False)
        :return: fig, axes
        """

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
        if log:
            molecules = self.molecules
        else:
            molecules = np.arcsinh(self.molecules)

        axes = []
        for i, g in enumerate(genes):
            ax = plt.subplot(gs[i // n_cols, i % n_cols])
            axes.append(ax)
            plt.scatter(self.tsne['x'], self.tsne['y'], c=molecules[g],
                        cmap=cmap, edgecolors='none', s=size, alpha=0.65)
            ax.set_axis_off()
            ax.set_title(g)

        return fig, axes

    def plot_aggregate_gene_expression(self, genes, fig=None, ax=None, title=None):
        """
        Display summed expression of a gene set on a tSNE projection.

        :param genes: genes to be summed over
        :param fig: Matplotlib Figure
        :param ax: Matplotlib Axis
        :param title: Plot this title on the figure
        :return: fig, ax
        """

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
                   cmap=cmap, edgecolors='none', s=size, alpha=0.65)
        ax.xaxis.set_major_locator(plt.NullLocator())
        ax.yaxis.set_major_locator(plt.NullLocator())
        ax.set_title(title)
        return fig, ax

    def plot_molecules_vs_reads_per_molecule(
            self, fig=None, ax=None, min_molecules=10, ylim=(0, 150),
            title='Cell Coverage Plot'):
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

    def remove_housekeeping_genes(self, organism='human', additional_hk_genes=set()):
        """
        Remove snoRNA and miRNAs because they are undetectable by the library construction
        procedure; their identification is a probable false-postive. Then remove
        housekeeping genes based on GO annotations.

        :param additional_hk_genes: set of additional string gene ids to remove
        :param organism: options: [human, mouse], default=human.
        """
        genes = self.molecules.columns
        # MT genes
        genes = genes[~genes.str.contains('MT-')]
        # "AL" miRNA
        genes = genes[~genes.str.contains('^AL[0-9]')]
        # "FP" miRNA
        genes = genes[~genes.str.contains('^FP[0-9]')]
        # SNO RNAs
        genes = genes[~genes.str.contains('RNU|SNO|RNVU')]
        # Ribosomal RNA
        genes = genes[~genes.str.contains('RNA5S')]

        # housekeeping gene ontologies
        cats = ['GO_0016072|rRNA metabolic process', 'GO_0006412|translation',
                'GO_0006414|translational elongation', 'GO_0006364|rRNA processing']

        with open(os.path.expanduser('~/.seqc/tools/{organism}/gofat.bp.v1.0.gmt.txt'
                                     ''.format(organism=organism))) as f:
            data = {fields[0]: fields[1:] for fields in
                    [line[:-1].split('\t') for line in f.readlines()]}

        hk_genes = {'AC091053.2', 'RPL13AP5-RPL13A', 'RP11-864I4.1-EEF1G', 'GAPDH',
                    'CTD-3035D6.1', 'CTC-575D19.1', 'ACTB', 'ACTG1', 'RP5-1056L3.3',
                    'U1', 'U2', 'U3', 'U4', 'U6', 'U7'}
        hk_genes.update(additional_hk_genes)
        for c in cats:
            hk_genes.update(data[c])

        # remove genes
        genes = genes.difference(hk_genes)

        exp = Experiment(molecules=deepcopy(self.molecules[genes]), reads=None,
                         metadata=deepcopy(self.metadata))
        exp._normalized = self._normalized
        return exp

    def annotate_clusters_by_expression(self, alpha=0.05):
        """
        Annotate clusters with genes that are uniquely expressed in them.

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
        split_indices = np.where(np.diff(self.cluster_assignments.iloc[idx]))[0] + 1

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
            res = gene_results.sum(axis=1) >= gene_results.shape[0] - 1
            global_results[:, i] = res

        return pd.DataFrame(global_results, columns=self.molecules.columns[reject],
                            index=np.unique(self.cluster_assignments))
