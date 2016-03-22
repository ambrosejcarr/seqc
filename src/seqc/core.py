import re
import os
import random
import pickle
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
import seaborn as sns
import phenograph
from tsne import bh_sne
import seqc

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

matplotlib.rc('lines', **{'color': 'royalblue'})

cmap = matplotlib.cm.viridis


class SparseMatrixError(Exception):
    pass


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
    xy = np.vstack([x, y])
    z = gaussian_kde(xy)(xy)
    i = np.argsort(z)
    return x[i], y[i], np.arcsinh(z[i])


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
        if not isinstance(index, np.array):
            raise TypeError('index must be type np.array')
        if not isinstance(columns, np.array):
            raise TypeError('columns must be type np.array')

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
            self._index = np.array(item)
        except:
            raise TypeError('self.columns must be convertible into a np.array object')

    def sum(self, axis=0):
        return self.data.sum(axis=axis)


class Experiment:

    def __init__(self, reads, molecules):
        if not isinstance(reads, SparseFrame) or isinstance(reads, pd.DataFrame):
            raise TypeError('reads must be of type SparseFrame or DataFrame')
        if not isinstance(molecules, SparseFrame) or isinstance(molecules, pd.DataFrame):
            raise TypeError('molecules must be of type SparseFrame or DataFrame')
        self._reads = reads
        self._molecules = molecules
        self._normalized_data = None  # filled by normalize_data; will be DataFrame
        self._pca = None
        self._tsne = None
        self._diffusion_eigenvectors = None
        self._diffusion_eigenvalues = None
        self._diffusion_map_correlations = None

    @property
    def reads(self):
        return self._reads

    @reads.setter
    def reads(self, item):
        if not isinstance(item, SparseFrame) or isinstance(item, pd.DataFrame):
            raise TypeError('Experiment.reads must be of type SparseFrame or DataFrame')
        self._reads = item

    @property
    def molecules(self):
        return self._molecules

    @molecules.setter
    def molecules(self, item):
        if not isinstance(item, SparseFrame) or isinstance(item, pd.DataFrame):
            raise TypeError('Experiment.molecules must be of type SparseFrame or'
                            'DataFrame')
        self._molecules = item

    @property
    def pca(self):
        return self._pca

    @pca.setter
    def pca(self, item):
        if not isinstance(item, dict):
            raise TypeError('self.pca must be a dictionary of pd.DataFrame object')
        self._pca = item

    @property
    def tsne(self):
        return self._tsne

    @tsne.setter
    def tsne(self, item):
        if not isinstance(item, pd.DataFrame):
            raise TypeError('self.tsne must be a pd.DataFrame object')
        self._tsne = item

    @property
    def diffusion_eigenvectors(self):
        return self._diffusion_eigenvectors

    @diffusion_eigenvectors.setter
    def diffusion_eigenvectors(self, item):
        if not isinstance(item, pd.DataFrame):
            raise TypeError('self.diffusion_eigenvectors must be a pd.DataFrame object')
        self._diffusion_eigenvectors = item

    @property
    def diffusion_eigenvalues(self):
        return self._diffusion_eigenvalues

    @diffusion_eigenvalues.setter
    def diffusion_eigenvalues(self, item):
        if not isinstance(item, pd.DataFrame):
            raise TypeError('self.diffusion_eigenvalues must be a pd.DataFrame object')
        self._diffusion_eigenvalues = item

    @property
    def diffusion_map_correlations(self):
        return self._diffusion_map_correlations

    @diffusion_map_correlations.setter
    def diffusion_map_correlations(self, item):
        if not isinstance(item, pd.DataFrame):
            raise TypeError('self.diffusion_map_correlations must be a pd.DataFrame'
                            'object')
        self._diffusion_map_correlations = item

    @classmethod
    def from_v012(cls, read_and_count_matrix: dict):
        d = read_and_count_matrix
        reads = SparseFrame(d['reads']['matrix'], index=d['reads']['row_ids'],
                            columns=d['reads']['col_ids'])
        molecules = SparseFrame(d['molecules']['matrix'], index=d['molecules']['row_ids'],
                                columns=d['molecules']['col_ids'])
        return cls(reads, molecules)

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
        experiment = cls(data['_reads'], data['_molecules'])
        experiment.normalized_data = data['_normalized_data']
        experiment.pca = data['_pca']
        experiment.tsne = data['_tsne']
        experiment.diffusion_eigenvectors = data['_diffusion_eigenvectors']
        experiment.diffusion_eigenvalues = data['_diffusion_eigenvalues']
        experiment.diffusion_map_correlations = ['_diffusion_map_correlations']
        return experiment

    def is_sparse(self):
        if all(isinstance(o, SparseFrame) for o in (self.reads, self.molecules)):
            return True
        else:
            return False

    def is_dense(self):
        return not self.is_sparse()

    def ensembl_gene_id_to_official_gene_symbol(self, gtf) -> None:
        """convert self.index containing scids into an index of gene names
        :param gtf:
        """
        pattern = re.compile(
                r'(^.*?gene_id "[^0-9]*)([0-9]*)(\..*?gene_name ")(.*?)(".*?$)')

        gene_id_map = defaultdict(set)
        with open(gtf, 'r') as f:
            for line in f:
                match = re.match(pattern, line)
                if match:
                    gene_id_map[int(match.group(2))].add(match.group(4))

        self.reads.index = ['-'.join(gene_id_map[i]) for i in self.reads.index]
        self.molecules.index = ['-'.join(gene_id_map[i]) for i in self.molecules.index]

    def plot_molecules_vs_reads_per_molecule(self, fig=None, ax=None, min_molecules=10):

        """
        plots log10 molecules counts per barcode vs reads per molecule / molecules per
        barcode
        :param min_molecules:
        :param ax:
        :param fig:

        """
        fig, ax = get_fig(fig, ax)

        # get molecule and read counts per cell
        molecule_cell_sums = self.molecules.sum(axis=1)
        read_cell_sums = self.reads.sum(axis=1)

        # remove low expression reads
        molecule_cell_sums = molecule_cell_sums[molecule_cell_sums >= min_molecules]
        read_cell_sums = read_cell_sums[molecule_cell_sums >= min_molecules]

        ratios = read_cell_sums / molecule_cell_sums
        x, y, z = density_2d(np.log10(molecule_cell_sums), ratios)
        ax.scatter(x, y, edgecolor='none', s=10, c=z, cmap=cmap)
        ax.set_xlabel('log10(Molecules per barcode)')
        ax.set_ylabel('Reads per barcode / Molecules per barcode')

        return fig, ax

    # todo does this raise errors if carried out on a DataFrame?
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
        self.reads = pd.DataFrame(read_data, index=self.reads.index[row_inds],
                                  columns=self.reads.columns)
        molecule_data = self.molecules.data.tocsr()[row_inds, :].todense()
        self.molecules = pd.DataFrame(molecule_data, index=self.molecules.index[row_inds],
                                      columns=self.molecules.columns)

    # todo check if works on both sparse and dense matrices
    def normalize_data(self):
        """
        use AVO method
        :return:
        """
        self.molecules = (self.molecules / self.molecules.sum(axis=1) /
                          np.median(self.molecules.sum(axis=1)))
        self.reads = (self.reads / self.reads.sum(axis=1) /
                      np.median(self.reads.sum(axis=1)))

    def plot_mitochondrial_molecule_fraction(self, fig=None, ax=None):
        """

        :param fig:
        :param ax:
        :return: fig, ax
        """

        if self.is_sparse():
            raise SparseMatrixError('Must convert to dense matrix before calling this '
                                    'function. Use self.remove_non_cell_barcodes()')

        mt_genes = self.molecules.index[self.molecules.columns.str.contains('MT-')]
        mt_counts = self.molecules[mt_genes].sum(axis=1)
        library_size = self.molecules.sum(axis=1)

        fig, ax = get_fig(fig=fig, ax=ax)
        x, y, z = density_2d(library_size, mt_counts / library_size)
        ax.scatter(x, y, s=5, edgecolors='none', c=z)
        ax.set_title('Mitochondrial Fraction')
        ax.set_xlabel('Molecules per cell')
        ax.set_ylabel('Mitochondrial molecules / Total molecules')

        return fig, ax

    def exclude_dead_cells_with_high_mt_fraction(self, max_mt_fraction=0.2):
        """

        :param max_mt_fraction:
        :return: Experiment
        """
        mt_genes = self.molecules.index[self.molecules.columns.str.contains('MT-')]
        mt_counts = self.molecules[mt_genes].sum(axis=1)
        library_size = self.molecules.sum(axis=1)
        ratios = mt_counts / library_size
        pass_filter = ratios.index[ratios <= max_mt_fraction]

        molecules = self.molecules.ix[pass_filter]
        reads = self.molecules.ix[pass_filter]

        return Experiment(reads=reads, molecules=molecules)

    def plot_molecules_vs_genes(self, fig=None, ax=None):
        """
        should be linear relationship
        :param ax:
        :param fig:
        :return:
        """
        fig, ax = get_fig(fig=fig, ax=ax)
        molecule_counts = self.molecules.sum(axis=1)
        gene_counts = np.sum(self.molecules > 0, axis=1)
        x, y, z = density_2d(np.log10(molecule_counts), np.log10(gene_counts))
        ax.scatter(x, y, c=z, edgecolor='none', s=5)
        ax.set_xlabel('log10(Number of molecules)')
        ax.set_ylabel('log10(Number of genes)')

        return fig, ax

    def remove_low_complexity_cells(self):
        # todo implement, based on below-fit cells in the above plot
        raise NotImplementedError

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
        matlab_cmd += "data = csvread('/tmp/pc_data_%d.csv', 1, 1);" % rand_tag
        # PCA
        matlab_cmd += "[mappedX, mapping]  = pca2(data, %d);" % no_components

        # Write results to file
        matlab_cmd += " csvwrite('/tmp/pc_mapped_x_%f.csv', mappedX);" % rand_tag
        matlab_cmd += " csvwrite('/tmp/pc_mapping_M_%f.csv', mapping.M);" % rand_tag
        matlab_cmd += (" csvwrite('/tmp/pc_mapping_lambda_%f.csv', mapping.lambda); exit;"
                       "\"" % rand_tag)

        # Run matlab command
        call(['matlab', '-nodesktop', '-nosplash', '-r %s' % matlab_cmd])

        # todo | ask manu about these files & their meanings; if possible, make similar to
        # todo | sklearn's implementation
        # Read in results
        self.pca = pd.Series()
        self.pca['mappedX'] = pd.DataFrame.from_csv(
                '/tmp/pc_mapped_x_%f.csv' % rand_tag, header=None, index_col=None)
        self.pca['mappedX'].index = self.molecules.index
        self.pca['mappingM'] = pd.DataFrame.from_csv(
                '/tmp/pc_mapping_M_%f.csv' % rand_tag, header=None, index_col=None)
        self.pca['mappingM'].index = self.molecules.columns
        self.pca['mappinglambda'] = pd.DataFrame.from_csv(
                '/tmp/pc_mapping_lambda_%f.csv' % rand_tag, header=None, index_col=None)

        # Clean up
        os.remove('/tmp/pc_data_%f.csv' % rand_tag)
        os.remove('/tmp/pc_mapped_x_%f.csv' % rand_tag)
        os.remove('/tmp/pc_mapping_M_%f.csv' % rand_tag)
        os.remove('/tmp/pc_mapping_lambda_%f.csv' % rand_tag)

    def run_tsne(self, n_components):
        """
        normally run on PCA components; 1st component is normally excluded

        :param n_components:
        :return:
        """
        data = deepcopy(self.molecules)
        data -= np.min(np.ravel(data))
        data /= np.max(np.ravel(data))  # todo I changed the below 0 -> 1 (exclude c#1)
        data = pd.DataFrame(np.inner(data, self.pca['mappingM'].iloc[:, 1:n_components]),
                            index=self.molecules.columns)

        self.tsne = pd.DataFrame(bh_sne(data),
                                 index=self.molecules.index, columns=['x', 'y'])

    # todo stopped here
    def phenograph(self):
        """
        normally run on PCA components; 1st component is normally excluded
        :return:
        """
        raise NotImplementedError

    def select_retained_clusters_based_on_library_size_distribution(self):
        """
        throws out low library size clusters; normally there is only one, but worth
        visual inspection.

        :return:
        """
        raise NotImplementedError

    def calculate_diffusion_map_components(self, knn=10, no_eigs=20, params=pd.Series()):
        """
        :param knn:
        :param no_eigs:
        :param params:
        :return:
        """

        # Random tag to allow for multiple diffusion map runs
        rand_tag = random.random()

        # Write to csv file
        self.molecules.to_csv('/tmp/dm_data_%f.csv' % rand_tag)

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
        dm_params['kEigenVecs'] = str(no_eigs)
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
        if not self.tsne:
            raise RuntimeError('Please run tSNE before plotting diffusion components.')

        fig, ax = get_fig(fig=fig, ax=ax)

        plt.figure(figsize=[12, 6])  # todo set adaptively
        no_cols = int(self.diffusion_eigenvectors.shape[1]/2)
        gs = plt.GridSpec(2, no_cols)

        for i in range(self.diffusion_eigenvectors.shape[1]):
            plt.subplot(gs[int(np.floor(i/no_cols)), i % no_cols])
            plt.scatter(self.tsne['x'], self.tsne['y'], c=self.diffusion_eigenvectors[i],
                        cmap=cmap, edgecolors='none', s=10)

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
                .sort(inplace=False).index
            x = range(no_cells, len(order))

            # For each gene
            for gene in self.molecules.columns:
                # Rolling mean
                vals = pd.rolling_mean(self.molecules.ix[order, gene],
                                       no_cells)[no_cells:]

                # Determine correlation
                self.diffusion_map_correlations.ix[gene, component] = pearsonr(x, vals)[0]

        # Reset NaNs
        self.diffusion_map_correlations.fillna(0)

    def run_gsea(self, out_dir, out_prefix, gmt_file, components=None,
                 gsea_jar='~/tools/gsea2-2.2.1.jar'):
        if not self.diffusion_eigenvectors:
            raise RuntimeError('Please run self.calculate_diffusion_map_components() '
                               'before running GSEA to annotate those components.')

        if components is None:
            components = range(1, self.diffusion_eigenvectors.shape[1])

        # Run for each component
        for c in components:

            print(c)

            # Write genes to file
            genes_file = out_dir + out_prefix + '_cmpnt_%d.rnk' % c
            ranked_genes = self.diffusion_map_correlations.iloc[:, c]\
                .sort(inplace=False, ascending=False)
            pd.DataFrame(ranked_genes).to_csv(genes_file, sep='\t', header=False)

            # Construct the GSEA call
            cmd = list()
            cmd += ['java', '-cp', gsea_jar,  '-Xmx1G', 'xtools.gsea.GseaPreranked']
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

    def select_genes_from_diffusion_components(self, components):
        """

        done based on GSEA enrichments

        :param components:
        :return:
        """
        if not self.diffusion_map_correlations:
            raise RuntimeError('Please run self.determine_gene_diffusion_correlations() '
                               'before selecting genes based on those correlations.')

        # Plot the correlation distributions along the selected components
        plt.figure(figsize=[5, 5])
        for c in components:
            sns.kdeplot(self.diffusion_map_correlations.iloc[:, c], label=c)

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
        return Experiment(subset_reads, subset_molecules)

    # todo add these later
    def differential_expression(self, g1, g2):
        """
        carry out differential expression between cells g1 and cells g2
        :param g1:
        :param g2:
        :return:
        """
        raise NotImplementedError
