import numpy as np
from numpy import ma
import pandas as pd
import multiprocessing
import fastcluster
from sklearn.neighbors import NearestNeighbors
from sklearn.manifold import TSNE
import warnings
import time
from scipy.sparse import csr_matrix, find
from scipy.sparse.linalg import eigs
from scipy.spatial.distance import pdist
from scipy.cluster.hierarchy import dendrogram
from numpy.linalg import norm
from statsmodels.sandbox.stats.multicomp import multipletests
from scipy.stats.mstats import kruskalwallis, rankdata
from scipy.stats import t
from tinydb import TinyDB
from functools import partial
from collections import namedtuple


class GraphDiffusion:
    def __init__(self, knn=10, normalization='smarkov', epsilon=1,
                 n_diffusion_components=10):
        """
        Run diffusion maps on the data. This implementation is based on the
        diffusion geometry library in Matlab:
        https://services.math.duke.edu/~mauro/code.html#DiffusionGeom and was implemented
        by Pooja Kathail

        :param knn: Number of neighbors for graph construction to determine distances between cells
        :param normalization: method for normalizing the matrix of weights
             'bimarkov'            force row and column sums to be 1
             'markov'              force row sums to be 1
             'smarkov'             symmetric conjugate to markov
             'beltrami'            Laplace-Beltrami normalization ala Coifman-Lafon
             'sbeltrami'           symmetric conjugate to beltrami
             'FokkerPlanck'        Fokker-Planck normalization
             'sFokkerPlanck'       symmetric conjugate to Fokker-Planck normalization
        :param epsilon: Gaussian standard deviation for converting distances to affinities
        :param n_diffusion_components: Number of diffusion components to generate
        """
        if normalization not in ['bimarkov', 'smarkov', 'markov', 'sbeltrami', 'beltrami',
                                 'FokkerPlanck', 'sFokkerPlanck']:
            raise ValueError(
                'Unsupported normalization. Please refer to the docstring for the '
                'supported methods')

        self.knn = knn
        self.normalization = normalization
        self.epsilon = epsilon
        self.n_diffusion_components = n_diffusion_components
        self.eigenvectors = None
        self.eigenvalues = None
        self.diffusion_operator = None
        self.weights = None

    @staticmethod
    def keigs(T, k, P, take_diagonal=0):
        """ return k largest magnitude eigenvalues for the matrix T.
        :param T: Matrix to find eigen values/vectors of
        :param k: number of eigen values/vectors to return
        :param P: in the case of symmetric normalizations,
                  this is the NxN diagonal matrix which relates the nonsymmetric
                  version to the symmetric form via conjugation
        :param take_diagonal: if 1, returns the eigenvalues as a vector rather than as a
                              diagonal matrix.
        """
        D, V = eigs(T, k, tol=1e-4, maxiter=1000)
        D = np.real(D)
        V = np.real(V)
        inds = np.argsort(D)[::-1]
        D = D[inds]
        V = V[:, inds]
        if P is not None:
            V = P.dot(V)

        # Normalize
        for i in range(V.shape[1]):
            V[:, i] = V[:, i] / norm(V[:, i])
        V = np.round(V, 10)

        if take_diagonal == 0:
            D = np.diag(D)

        return V, D

    @staticmethod  # todo fix; what is S?
    def bimarkov(W, max_iters=100, abs_error=0.00001, verbose=False, **kwargs):
        """normalization method for GraphDiffusion"""

        if W.size == 0:
            return

        # process input
        if W.shape[0] != W.shape[1]:
            raise ValueError('Bimarkov.py: kernel must be NxN\n')

        N = W.shape[0]

        # initialize
        p = np.ones(N)

        # iterative
        for i in range(max_iters):

            S = np.ravel(W.sum(axis=1))
            err = np.max(np.absolute(1.0 - np.max(S)), np.absolute(1.0 - np.min(S)))

            if err < abs_error:
                break

            D = csr_matrix((np.divide(1, np.sqrt(S)), (range(N), range(N))), shape=[N, N])
            p *= S
            W = D.dot(W).dot(D)

        # iron out numerical errors
        T = (W + W.T) / 2
        return T, p

    @staticmethod
    def smarkov(D, N, W):
        """normalization method for GraphDiffusion"""
        D = csr_matrix((np.sqrt(D), (range(N), range(N))), shape=[N, N])
        P = D
        T = D.dot(W).dot(D)
        T = (T + T.T) / 2
        return T, P

    @staticmethod
    def markov(D, N, W):
        """normalization method for GraphDiffusion"""
        T = csr_matrix((D, (range(N), range(N))), shape=[N, N]).dot(W)
        return T, None

    @staticmethod
    def sbeltrami(D, N, W):
        """normalization method for GraphDiffusion"""
        P = csr_matrix((D, (range(N), range(N))), shape=[N, N])
        K = P.dot(W).dot(P)

        D = np.ravel(K.sum(axis=1))
        D[D != 0] = 1 / D[D != 0]

        D = csr_matrix((D, (range(N), range(N))), shape=[N, N])
        P = D
        T = D.dot(K).dot(D)

        T = (T + T.T) / 2
        return T, P

    @staticmethod
    def beltrami(D, N, W):
        """normalization method for GraphDiffusion"""
        D = csr_matrix((D, (range(N), range(N))), shape=[N, N])
        K = D.dot(W).dot(D)

        D = np.ravel(K.sum(axis=1))
        D[D != 0] = 1 / D[D != 0]

        V = csr_matrix((D, (range(N), range(N))), shape=[N, N])
        T = V.dot(K)
        return T, None

    @staticmethod
    def FokkerPlanck(D, N, W):
        """normalization method for GraphDiffusion"""
        D = csr_matrix((np.sqrt(D), (range(N), range(N))), shape=[N, N])
        K = D.dot(W).dot(D)

        D = np.ravel(K.sum(axis=1))
        D[D != 0] = 1 / D[D != 0]

        D = csr_matrix((D, (range(N), range(N))), shape=[N, N])
        T = D.dot(K)
        return T, None

    def sFokkerPlanck(self, D, N, W):
        """normalization method for GraphDiffusion"""
        print('(sFokkerPlanck) ... ')

        D = csr_matrix((np.sqrt(D), (range(N), range(N))), shape=[N, N])
        K = D.dot(W).dot(D)

        D = np.ravel(K.sum(axis=1))
        D[D != 0] = 1 / D[D != 0]

        D = csr_matrix((np.sqrt(D), (range(N), range(N))), shape=[N, N])
        P = D
        T = D.dot(K).dot(D)

        T = (T + T.T) / 2
        return T, P

    def fit(self, data, verbose=True):
        """
        :param data: Data matrix of samples X features
        :param verbose: print progress report

        :return: Dictionary containing diffusion operator, weight matrix,
                 diffusion eigen vectors, and diffusion eigen values
        """
        if verbose:
            print('Running Diffusion maps with the following parameters:')
            print('Normalization: %s' % self.normalization)
            print('Number of nearest neighbors k: %d' % self.knn)
            print('Epsilon: %.4f' % self.epsilon)

        # Nearest neighbors
        start = time.process_time()
        N = data.shape[0]
        nbrs = NearestNeighbors(n_neighbors=self.knn).fit(data)
        distances, indices = nbrs.kneighbors(data)

        # Adjacency matrix
        rows = np.zeros(N * self.knn, dtype=np.int32)
        cols = np.zeros(N * self.knn, dtype=np.int32)
        dists = np.zeros(N * self.knn)
        location = 0
        for i in range(N):
            inds = range(location, location + self.knn)
            rows[inds] = indices[i, :]
            cols[inds] = i
            dists[inds] = distances[i, :]
            location += self.knn
        W = csr_matrix((dists, (rows, cols)), shape=[N, N])

        # Symmetrize W
        W = W + W.T

        # Convert to affinity (with selfloops)
        rows, cols, dists = find(W)
        rows = np.append(rows, range(N))
        cols = np.append(cols, range(N))
        dists = np.append(dists / (self.epsilon ** 2), np.zeros(N))
        W = csr_matrix((np.exp(-dists), (rows, cols)), shape=[N, N])

        # Create D
        D = np.ravel(W.sum(axis=1))
        D[D != 0] = 1 / D[D != 0]

        # Go through the various normalizations
        fnorm = getattr(self, self.normalization)
        T, P = fnorm(D=D, N=N, W=W)

        if self.normalization != 'bimarkov' and verbose:
            print('%.2f seconds' % (time.process_time() - start))

        # Eigen value decomposition
        V, D = GraphDiffusion.keigs(T, self.n_diffusion_components, P, take_diagonal=1)
        self.eigenvectors = V
        self.eigenvalues = D
        self.diffusion_operator = T
        self.weights = W


class NormalizeCells:
    """
    Contains row normalization techniques for equalizing cell sampling.
    """

    @staticmethod
    def size(data: np.ndarray, percentile=None):
        """
        Normalize data such that each cell has an identical library size

        :param data: n observation x k feature array
        :param log: if True, log transform the data
        :param percentile: if None, multiply each library by the total number
         of nonzero genes observed across the experiment. If float in [0, 1),
         multiply by the cell size observed at
         np.percentile(cellsums, percentile)
        :return: normalized n observation x k feature array
        """
        if percentile is None:
            data = data.div(data.sum(axis=1), axis=0) * data.shape[1]
        else:
            cell_sums = data.sum(axis=1)
            data = data.div(cell_sums, axis=0).mul(np.percentile(cell_sums, percentile), axis=0)
        return data


class ScaleFeatures:
    """
    Contains column normalization techniques for equalizing feature sampling
    """

    @staticmethod
    def unit_size(data: np.ndarray, copy=True):
        """
        scales features in data between 0 and 1

        :param data: n observation x k feature array
        :param copy: bool, if True, creates a copy of data. Otherwise carries out
          standardization in-place
        :return: n x k scaled array
        """
        if copy:
            data = data.copy().astype(np.float32)
        else:
            data = data.astype(np.float32)
        d_min, d_max = np.min(data), np.max(data)
        data -= d_min
        data /= d_max - d_min
        return data

    @staticmethod
    def standardize(data: np.ndarray, copy=True):
        """
        scales data such that each feature (column) has mean=0 and standard deviation=1

        :param data: n observation x k feature array
        :param copy: bool, if True, creates a copy of data. Otherwise carries out
          standardization in-place
        :return: n x k scaled array
        """
        if copy:
            data = data.copy().astype(np.float32)
        else:
            data = data.astype(np.float32)
        data -= data.mean(axis=0)
        data /= data.std(axis=0)

    @staticmethod
    def unit_length(data: np.ndarray, copy=True):
        """
        scales each feature so that its euclidean length is 1

        :param data: n observation x k feature array
        :param copy: bool, if True, creates a copy of data. Otherwise carries out
          standardization in-place
        :return: n x k scaled array
        """
        if copy:
            data = data.copy().astype(np.float32)
        else:
            data = data.astype(np.float32)
        data /= np.sqrt(np.sum(data ** 2, axis=0))
        return data


class tSNE:

    def __init__(self, n_components: int=2, scale: bool=True, run_pca: bool=True,
                 n_pca_components: int=10, **kwargs):
        """
        Carry out t-stochastic neighbor embedding

        :param n_components: number of components to which data should be projected
        :param normalize: if True, scales features to unit size
        :param run_pca: if True, runs PCA on the input data and runs tSNE on the
          components retained by PCA.
        :param n_pca_components:  Number of components retained by PCA
        :param kwargs:  additional keyword arguments to pass sklearn.manifold.TSNE
        """
        self.scale = scale
        self.run_pca = run_pca
        self.n_components = n_components
        self.n_pca_components = n_pca_components
        self.kwargs = kwargs
        self.tsne = None
        self.pca = None

    def fit_transform(self, data: np.ndarray or pd.DataFrame) -> None:
            """
            fit the tSNE model to data given the parameters provided during
             initialization and transform the output

            :param data: n observation x k feature data array
            :return: None
            """
            if self.scale:
                data = ScaleFeatures.unit_size(data, copy=True)
            if self.run_pca:
                self.pca = PCA(n_components=self.n_pca_components)
                data = self.pca.fit_transform(data)

            if not self.kwargs:
                self.kwargs = dict(angle=0.4, init='pca', random_state=0, n_iter=500,
                                   perplexity=30)

            tsne = TSNE(n_components=self.n_components, method='barnes_hut',
                        **self.kwargs)
            res = tsne.fit_transform(data.values)
            if isinstance(data, pd.DataFrame):
                self.tsne = pd.DataFrame(res, index=data.index)
            else:
                self.tsne = res


class PCA:

    def __init__(self, n_components=100):
        """
        construct a model for Principle Component Analysis

        :param n_components: number of principle components to retain
        """
        self.n_components = n_components
        self.loadings = None
        self.eigenvalues = None

    def fit(self, data):
        """
        Fit the model to data

        :param data: n observation x k feature data array
        :return: None
        """

        if isinstance(data, pd.DataFrame):
            X = data.values
        elif isinstance(data, np.ndarray):
            X = data
        else:
            raise TypeError('data must be a pd.DataFrame or np.ndarray')

        # Make sure data is zero mean
        X = ScaleFeatures.unit_size(X)
        # X = np.subtract(X, np.amin(X))
        # X = np.divide(X, np.amax(X))

        # Compute covariance matrix
        if X.shape[1] < X.shape[0]:
            C = np.cov(X, rowvar=0)
        # if N > D, we better use this matrix for the eigendecomposition
        else:
            C = np.multiply((1 / X.shape[0]), np.dot(X, X.T))

        # Perform eigendecomposition of C
        C[np.where(np.isnan(C))] = 0
        C[np.where(np.isinf(C))] = 0
        l, M = np.linalg.eig(C)

        # Sort eigenvectors in descending order
        ind = np.argsort(l)[::-1]
        l = l[ind]
        if self.n_components < 1:
            self.n_components = (
                np.where(np.cumsum(np.divide(l, np.sum(l)), axis=0) >=
                         self.n_components)[0][0] + 1)
            print('Embedding into ' + str(self.n_components) + ' dimensions.')
        elif self.n_components > M.shape[1]:
            self.n_components = M.shape[1]
            print('Target dimensionality reduced to ' + str(self.n_components) + '.')

        M = M[:, ind[:self.n_components]]
        l = l[:self.n_components]

        # Apply mapping on the data
        if X.shape[1] >= X.shape[0]:
            M = np.multiply(np.dot(X.T, M), (1 / np.sqrt(X.shape[0] * l)).T)

        self.loadings = M
        self.eigenvalues = l

    def transform(self, data, components=None) -> np.ndarray or pd.DataFrame:
        """
        Transform data using the fit PCA model.

        :param data:  n observation x k feature data array
        :param components:  components to retain when transforming
          data, if None, uses all components
        :return: np.ndarray containing transformed data
        """

        if components is None:
            components = np.arange(self.n_components)
        projected = np.dot(data, self.loadings[:, components])
        if isinstance(data, pd.DataFrame):
            return pd.DataFrame(projected, index=data.index)
        else:
            return projected

    def fit_transform(self, data, n_components=None) -> np.ndarray or pd.DataFrame:
        """
        Fit the model to data and transform the data using the fit model

        :param data:  n observation x k feature data array
        :param n_components:  number of components to retain when transforming
          data
        :return: np.ndarray containing transformed data
        """

        self.fit(data)
        return self.transform(data, components=n_components)


class correlation:

    @staticmethod
    def vector(x: np.array, y: np.array):
        """
        Correlate each column in y with a vector x

        :param x: np.ndarray vector of length n
        :param y: np.ndarray matrix of shape (n, k)
        :returns: vector of length n
        """
        # x = x[:, np.newaxis]  # for working with matrices
        mu_x = x.mean()  # cells
        mu_y = y.mean(axis=0)  # cells by gene --> cells by genes
        sigma_x = x.std()
        sigma_y = y.std(axis=0)

        return ((y * x).mean(axis=0) - mu_y * mu_x) / (sigma_y * sigma_x)


    @staticmethod
    def map(x, y):
        """Correlate each column of x with each column of y

        :param x: np.array; shape N x T.
        :param y: np.array; shape M x T.
        :returns: np.array; shape N x M in which each element is a correlation
                            coefficient.
        """
        N = y.shape[1]
        sum_x = x.sum(axis=1)
        sum_y = y.sum(axis=1)
        p1 = N * np.dot(y, x.T)
        p2 = sum_x * sum_y[:, np.newaxis]
        p3 = N * ((y ** 2).sum(axis=1)) - (sum_y ** 2)
        p4 = N * ((x ** 2).sum(axis=1)) - (sum_x ** 2)
        return (p1 - p2) / np.sqrt(p4 * p3[:, None])

    @staticmethod
    def eigv(evec, data, components=tuple(), knn=10):
        """
        get pairwise correlations of eigenvectors with columns in data

        :param evec: eigenvectors
        :param data: np.ndarray genes x cells data matrix
        :param components: which eigenvectors to select
        :param knn: number of neighbors to smooth gene expression values over
        :return:
        """
        if isinstance(data, pd.DataFrame):
            D = data.values
            df = True
        elif isinstance(data, np.ndarray):
            D = data
            df = False
        else:
            raise TypeError('data must be a pd.DataFrame or np.ndarray')

        # set components, remove zero if it was specified
        if not components:
            components = np.arange(evec.shape[1])
        else:
            components = np.array(components)
        components = components[components != 0]

        eigv_corr = np.empty((D.shape[1], evec.shape[1]), dtype=np.float)

        for component_index in components:
            component_data = evec[:, component_index]

            order = np.argsort(component_data)
            x = pd.DataFrame(component_data[order]).rolling(
                window=knn, center=False).mean()[knn:].values
            # this fancy indexing will copy self.molecules
            vals = pd.DataFrame(D[order, :]).rolling(
                window=knn, center=False, axis=0).mean()[knn:].values
            eigv_corr[:, component_index] = correlation.vector(x, vals)

        # this is sorted by order, need it in original order (reverse the sort)
        eigv_corr = eigv_corr[:, components]
        if df:
            eigv_corr = pd.DataFrame(eigv_corr, index=data.columns, columns=components)
        return eigv_corr


class ExperimentalYield:

    output = (
        '{divide}\nINPUT\n{divide}\n'
        'Total input reads:\t{n_fastq}\n'
        '{divide}\nALIGNMENT (% FROM INPUT)\n{divide}\n'
        'Total reads aligned:\t{n_sam} ({prop_al}%)\n'
        ' - Genomic alignments:\t{genomic} ({prop_gen}%)\n'
        ' - PhiX alignments:\t{phi_x} ({prop_phix}%)\n'
        ' - Transcriptome alignments:\t{trans} ({prop_trans}%)\n'
        '{divide}\nFILTERING (% FROM ALIGNMENT)\n{divide}\n'
        'Genomic alignments:\t{genomic} ({bad_gen}%)\n'
        'PhiX alignments:\t{phi_x} ({bad_phi}%)\n'
        'Incorrect barcodes:\t{wrong_cb} ({bad_cb}%)\n'
        'Missing cell barcodes/RMT:\t{no_cell} ({bad_cell}%)\n'
        'N present in RMT:\t{rmt_N} ({bad_rmtN}%)\n'
        'Insufficient poly(T):\t{poly_t} ({bad_polyt}%)\n'
        '{divide}\nCELL/MOLECULE COUNT DISTRIBUTION\n{divide}\n'
        'Total molecules:\t\t{tot_mc}\n'
        'Molecules lost:\t{mols_lost}\n'
        'Cells lost:\t{cells_lost}\n'
        'Cell description:\n{cell_desc}\n'
        '{divide}\nSUMMARY\n{divide}\n'
        'Total retained reads:\t{n_good} ({prop_good}%)\n'
        'Total reads unaligned:\t{lost_al} ({prop_un}%)\n'
        'Total reads filtered:\t{n_bad} ({prop_bad}%)\n'
        '{divide}\n')

    @classmethod
    def construct_run_summary(cls, summary: dict):
        """
        calculates basic loss statistics and constructs a summary
        that will be sent to the user after the SEQC run has completed.

        :param summary: dictionary constructed during error correction
        :return: output of basic summary statistics
        """
        if not summary:
            return

        # obtain values from summary
        n_fastq = summary['n_fastq']
        n_sam = summary['n_sam']
        genomic = summary['gene_0']
        phix = summary['phi_x']
        no_cell = summary['cell_0']
        # no_rmt = summary['rmt_0']
        rmt_N = summary['rmt_N']
        poly_t = summary['poly_t']
        tot_mc = summary['total_mc']
        mols_lost = list(summary['mols_lost'].items())
        cells_lost = list(summary['cells_lost'].items())
        cell_desc = summary['cell_desc'].to_string()
        divide = '-' * 40

        # run summary will not be calculated if user started SEQC midway
        if n_fastq == 'NA' or n_sam == 'NA':
            return

        # calculate summary statistics
        trans = n_sam - genomic - phix
        prop_al = round((n_sam/n_fastq) * 100, 1)
        prop_gen = round((genomic/n_sam) * 100, 1)
        prop_phix = round((phix/n_sam) * 100, 1)
        prop_trans = round((trans/n_sam) * 100, 1)
        lost_al = n_fastq - n_sam
        prop_un = round(100 - prop_al, 1)
        n_bad = genomic + phix + no_cell + rmt_N + poly_t
        # n_bad = genomic + phix + no_cell + no_rmt + rmt_N + poly_t
        # wrong_cb does not apply to drop-seq
        try:
            wrong_cb = summary['cb_wrong']
            n_bad += wrong_cb
            bad_cb = round((wrong_cb/n_bad) * 100, 1)
        except KeyError:
            wrong_cb = 0
            bad_cb = 0
        # continue with calculations
        n_good = n_sam - n_bad
        bad_gen = round((genomic/n_bad) * 100, 1)
        bad_phi = round((phix/n_bad) * 100, 1)
        bad_cell = round((no_cell/n_bad) * 100, 1)
        # bad_rmt = round((no_rmt/n_bad) * 100, 1)
        bad_rmtN = round((rmt_N/n_bad) * 100, 1)
        bad_polyt = round((poly_t/n_bad) * 100, 1)
        prop_bad = round((n_bad/n_fastq) * 100, 1)
        prop_good = round((n_good/n_fastq) * 100, 1)

        # format output
        output = cls.output.format(
            n_fastq=n_fastq, n_sam=n_sam, genomic=genomic, phi_x=phix, no_cell=no_cell,
            wrong_cb=wrong_cb, rmt_N=rmt_N, poly_t=poly_t, divide=divide,
            prop_al=prop_al, prop_gen=prop_gen, prop_phix=prop_phix, lost_al=lost_al,
            n_bad=n_bad, n_good=n_good, prop_good=prop_good, prop_bad=prop_bad,
            prop_un=prop_un, bad_gen=bad_gen, bad_phi=bad_phi, bad_cb=bad_cb,
            bad_cell=bad_cell, bad_rmtN=bad_rmtN, bad_polyt=bad_polyt, trans=trans,
            prop_trans=prop_trans, tot_mc=tot_mc, mols_lost=mols_lost,
            cells_lost=cells_lost, cell_desc=cell_desc)
        return output


class smoothing:

    @staticmethod
    def kneighbors(data, n_neighbors=50):
        """
        Smooth gene expression values by setting the expression of each gene in each
        cell equal to the mean value of itself and its n_neighbors

        :param data: np.ndarray | pd.DataFrame; genes x cells array
        :param n_neighbors: int; number of neighbors to smooth over
        :return: np.ndarray | pd.DataFrame; same as input
        """

        if isinstance(data, pd.DataFrame):
            df = True
            data_ = data.values
        elif isinstance(data, np.ndarray):
            df = False
            data_ = data
        else:
            raise TypeError("data must be a pd.DataFrame or np.ndarray")

        knn = NearestNeighbors(
            n_neighbors=n_neighbors,
            n_jobs=multiprocessing.cpu_count() - 1)

        knn.fit(data_)
        dist, inds = knn.kneighbors(data_)

        # set values equal to their means across neighbors
        res = data_[inds, :].mean(axis=1)

        if df:
            res = pd.DataFrame(res, index=data.index,
                               columns=data.columns)
        return res


class GSEA:

    def __init__(self, correlations, labels=None):
        """
        Create a Gene Set Enrichment Object

        :param correlations: pd.Series or iterable of pd.Series;
          objects containing correlations. If a single series is passed, the main process
          will compute gene set enrichments. If an iterable is passed, processes will be
          spawned to process each row in the correlation array. If the pd.Series contain
          name fields, these labels will be retained in the output.
        :param labels: iterable; alternatively, labels can be passed directly to GSEA.
          these will overwrite any names given in provided Series, and must be unique.

        Producing correlation vectors:
        (1) Correlation is a natural measure of similarity for diffusion components or
            PCA components. One can correlate the expression value of each gene with its
            loading for that vector.
        (2) Correlation can be used to measure association with a phenotype label (e.g.
            a cluster). For multi-cluster association, a one vs other binary labeling is
            appropriate in many cases. In circumstances where clusters have known
            progression, correlation can be run against ordered clusters.

        """
        if isinstance(correlations, pd.Series):
            self.correlations = [correlations,]
        elif all(isinstance(obj, pd.Series) for obj in correlations):
            self.correlations = list(correlations)
            # self.nproc = min((multiprocessing.cpu_count() - 1, len(correlations)))
        else:
            raise ValueError('correlations must be passed as pd.Series objects.')

        # set correlation labels
        if labels:
            if not len(labels) == len(self.correlations):
                raise ValueError(
                    'number of passed labels ({!s}) does not match number of passed '
                    'correlations ({!s})'.format(
                        len(labels), len(self.correlations)))
            if not len(labels) == len(set(labels)):
                raise ValueError('labels cannot contain duplicates.')
            for c, l in zip(self.correlations, labels):
                c.name = l
        else:
            for c in self.correlations:
                try:
                    getattr(c, 'name')
                except AttributeError:  # at least one series lacks a name; reset all
                    for i, c in enumerate(self.correlations):
                        self.correlations[i] = c.copy()  # don't mutate original vector
                        self.correlations[i].name = i
                    break
            if not len(set(c.name for c in self.correlations)) == len(self.correlations):
                raise ValueError('correlation Series names must be unique')

    class Database:

        @staticmethod
        def create(gene_set_files, file_labels, db_name):
            """
            Create a searchable database of Gene Sets from .gmt files provided by GSEA

            :param gene_set_files: [str, ...]; GSEA .gmt files to add to the database
            :param file_labels: [str, ...]; labels corresponding to each of the
              gene_set_files. Must be the same length as gene_set_files
            :param db_name: str; path and name of .json file in which to store database.
              ".json" will be appended to to the filename if db_name does not end with
              that suffix.
            :return: TinyDB database (stored in .json)
            """
            records = []
            for file_, label in zip(gene_set_files, file_labels):
                with open(file_, 'r') as f:
                    sets = [l.strip().split('\t') for l in f.readlines()]
                records.extend(
                    {'type': label.upper(), 'name': fields[0], 'genes': fields[2:]} for
                    fields in sets)

            if not db_name.endswith('.json'):
                db_name += '.json'
            db = TinyDB(db_name)
            db.insert_multiple(records)
            return db

        @staticmethod
        def load(db_name):
            """
            wrapper for TinyDB loading; loads an existing database from a .json file

            :param db_name:  database to load
            :return: TinyDB database
            """
            return TinyDB(db_name)

        @staticmethod
        def select(data, field):
            """
            return field from each data entry extracted from a database

            e.g. to get genes from a database query,
            $> data = Query(<query>)
            $> genes = GSEA.Database.select(data, 'genes')

            :param data: [dict, ...]; output of a database Query
            :param field: str; key for field to be extracted
            :return: list; field for each entry in data
            """
            return list(map(lambda x: x[field], data))

        @staticmethod
        def keys(db):
            """Return the available fields (keys) in a homogeneous database object"""
            return db.table('_default').all()[0].keys()

    class SetMasks:

        @staticmethod
        def from_db(sets, correlations):
            """
            Convenience function for GSEA. Map sets of genes to the index of
            correlatations for masking unused genes.

            :param sets: k sets of genes to map to the index of correlations
            :param correlations: correlations across n genes
            :return: np.ndarray of n genes x k sets
            """
            set_genes = GSEA.Database.select(sets, 'genes')

            # map genes to index; should the map consider the union?
            gene_map = dict(zip(
                correlations.index, np.arange(correlations.shape[0])))

            masks = np.zeros((len(sets), len(gene_map)), dtype=np.bool)
            for i, set_ in enumerate(set_genes):
                masks[i, [gene_map[g] for g in
                          set(set_).intersection(gene_map.keys())]] = 1
            return masks

    @staticmethod
    def _gsea_process(correlations, sets, n_perm, alpha):
        """
        Initialize a GSEA process where sets are tested against correlations.

        :param correlations: pd.Series; an ordered list of correlations
        :param sets: gene sets extracted from GSEA database
        :param n_perm: int; number of permutations to use in constructing the null
          distribution for significance testing
        :alpha: float; percentage of type I error to use when computing False Discovery
          Rate correction.
        """
        set_names = GSEA.Database.select(sets, 'name')
        masks = GSEA.SetMasks.from_db(sets, correlations)
        pos, neg = GSEA.construct_normalized_correlation_matrices(
            correlations, masks)
        p, es = GSEA.calculate_enrichment_significance(
            pos, neg, n_perm=n_perm, alpha=alpha)
        return pd.DataFrame({'p': p, 'es': es}, index=set_names)

    def test(self, sets, n_perm=500, alpha=0.05):
        """
        Test sets of genes against the correlations used to initialize the GSEA object.

        :param sets: gene sets to test against ordered correlations; these can be
          extracted from a Database using Query() syntax. e.g.:
            # retrieve a database
            db = seqc.stats.GSEA.Database.load('MSigDB.json')
            # search for sets of type 'HALLMARK' that have > 50 genes
            sets = db.search(
                (Query().type == 'HALLMARK') &
                (Query().genes.test(lambda x: len(x) > 50))
            )
            # search for sets with 'IMMUNE' in name:
            sets = db.search(Query().name.matches('.*?IMMUNE.*'))
        :param n_perm: int; number of permutations to use in constructing the null
          distribution for significance testing; note that n_perm sets the granularity of
          significance values; if n_perm=500, significance bins are of size 1/500.
        :param alpha: float; percentage of type I error to use when computing False
          Discovery Rate correction.
        :return: dict
        """

        def split(a, n):
            k, m = len(a) // n, len(a) % n
            return (a[i * k + min(i, m):(i + 1) * k + min(i + 1, m)] for i in range(n))

        def group(a, n):
            args = [iter(a)] * n
            return zip(*args)

        def join(a, n):
            return [pd.concat(s, axis=0) for s in group(a, n)]

        nproc = multiprocessing.cpu_count()
        ncorr = len(self.correlations)
        nsets = len(sets)

        assert(nsets != 0)
        assert(ncorr != 0)
        assert(nproc != 0)

        if ncorr < nproc:
            p = min(nproc // ncorr, nsets)
            correlations = self.correlations * p
            sets = list(split(sets, p))
        else:
            correlations = self.correlations
            sets = [sets]
            p = 1
        args = list(zip(correlations, sets * ncorr * p, [n_perm] * ncorr * p,
                        [alpha] * ncorr * p))

        # adjust nproc in case n * p != nproc
        nproc = min(ncorr * p, nproc)
        pool = multiprocessing.Pool(nproc)
        try:
            res = pool.starmap(GSEA._gsea_process, args)
        except:
            pool.close()
            raise
        pool.close()

        # merge sets
        res = join(res, p)

        res_df = {c.name: df for (c, df) in zip(self.correlations, res)}

        for df in res_df.values():
            df['adj_p'] = multipletests(df['p'], alpha=alpha / len(res_df),
                                        method='fdr_tsbh')[1]

        return res_df

    @staticmethod
    def construct_normalized_correlation_matrices(correlations, set_mask):
        """
        construct normalized matrices of genes in set_masks and not in set_masks

        :param correlations: pd.Series of correlations
        :param set_mask: single set mask
        :return: pos, neg: normalized sets
        """

        N = np.tile(correlations, (set_mask.shape[0], 1))  # tile one row per set mask

        # create arrays masked based on whether or not genes were in the gene set
        pos = ma.masked_array(N, mask=set_mask)
        # copy for neg, else can't use in-place operations
        neg = ma.masked_array(N.copy(), mask=~set_mask.copy())
        neg.fill(1)  # negative decrements are binary, not based on correlations

        # normalize positive and negative arrays
        pos /= ma.sum(pos, axis=1)[:, np.newaxis]
        neg /= ma.sum(neg, axis=1)[:, np.newaxis]

        return pos, neg

    @staticmethod
    def calculate_enrichment_score(pos, neg):
        """
        calculate a gene set enrichment score from normalized correlation matrices.

        :param pos: n x k matrix of n gene_sets x k normalized gene correlations where
          genes were in sets, else 0.
        :param neg: n x k matrix of n gene_sets x k genes where entries are equal to
          1 / k where genes were not in sets, else 0.
        :return: vector of maximum enrichment scores
        """
        cumpos = pos.filled(0).cumsum(axis=1)
        cumneg = neg.filled(0).cumsum(axis=1)
        es = cumpos - cumneg
        return es[np.arange(es.shape[0]), np.argmax(np.abs(es), axis=1)]

    @staticmethod
    def permute_correlation_ranks(pos_corr_matrix, neg_corr_matrix):
        """wrapper for np.permutation that catches warnings about a desired feature"""
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            np.random.shuffle(pos_corr_matrix)
            np.random.shuffle(neg_corr_matrix)

    @staticmethod
    def calculate_enrichment_significance(pos, neg, n_perm=1000, alpha=0.05):
        """
        calculate enrichment scores and significance of the enrichment

        :param pos: np.ndarray, normalized correlations of genes in set
        :param neg: np.ndarray, normalized correlations of genes not in set
        :param n_perm: int, number of permutations
        :param alpha: float (0, 1]
        :return: pvals, es
        """
        es = GSEA.calculate_enrichment_score(pos, neg)  # score

        # todo here's the slow part. Any way to parallelize?
        # calculate permutations
        perms = np.zeros((n_perm, pos.shape[0]), dtype=np.float)  # n_perm x num sets
        for i in np.arange(n_perm):
            GSEA.permute_correlation_ranks(pos, neg)
            perms[i, :] = GSEA.calculate_enrichment_score(pos, neg)

        # calculate p-vals
        pvals = 1 - np.sum(perms < np.abs(es[np.newaxis, :]), axis=0) / perms.shape[0]

        return pvals, es


def _sampling_function(n_iter, n_molecules, theta, n_cells):

    def online_mean_var(nb, mu_b, var_b, na, mu_a, var_a):
        nx = na + nb
        delta = mu_b - mu_a
        mu_x = mu_a + delta * nb / nx
        var_x = (na * (var_a + mu_a ** 2) + nb * (var_b + mu_b ** 2)) / nx - mu_x ** 2
        return nx, mu_x, var_x

    res_mu = np.zeros((n_iter, theta.shape[0]), dtype=np.float32)
    res_var = np.zeros((n_iter, theta.shape[0]), dtype=np.float32)
    n_cells //= 10
    for i in np.arange(n_iter):
        # break sampling (n_cells) into 10 pieces
        obs = np.random.multinomial(n_molecules, theta, n_cells)
        mu_x = np.mean(obs, axis=0)
        var_x = np.mean(obs, axis=0)
        n_x = obs.shape[0]
        for _ in np.arange(9):
            obs = np.random.multinomial(n_molecules, theta, n_cells)
            mu = np.mean(obs, axis=0)
            var = np.mean(obs, axis=0)
            n = obs.shape[0]
            n_x, mu_x, var_x = online_mean_var(n, mu, var, n_x, mu_x, var_x)
        res_mu[i, :] = mu_x
        res_var[i, :] = var_x / n_x
    return res_mu, res_var


class ExpressionTree:

    def __init__(self, data, group_assignments, alpha=0.05):
        """
        Create an hierarchical tree and calculate differential expression at each branch
        point, agglomerating over the leaves of the branch. Uses multinomial resampling
        of cells to eliminate non-identical sampling across cell types

        :param data: n cells x k genes 2d array
        :param group_assignments: n cells 1d vector
        :param alpha: float (0, 1], acceptable type I error
        """

        self._tree = None
        self._results = None

        # make sure group_assignments and data have the same length
        if not data.shape[0] == group_assignments.shape[0]:
            raise ValueError(
                'Group assignments shape ({!s}) must equal the number of rows in data '
                '({!s}).'.format(group_assignments.shape[0], data.shape[0]))

        # store index if both data and group_assignments are pandas objects
        if isinstance(data, pd.DataFrame) and isinstance(group_assignments, pd.Series):
            # ensure assignments and data indices are aligned
            try:
                ordered_assignments = group_assignments[data.index]
                if not len(ordered_assignments) == data.shape[0]:
                    raise ValueError(
                        'Index mismatch between data and group_assignments detected when '
                        'aligning indices. check for duplicates.')
            except:
                raise ValueError('Index mismatch between data and group_assignments.')

            # sort data by cluster assignment
            idx = np.argsort(ordered_assignments.values)
            if not np.array_equal(idx, np.arange(idx.shape[0])):
                self._data = data.iloc[idx, :].values
                ordered_assignments = ordered_assignments.iloc[idx]
                self._group_assignments = ordered_assignments.values
            else:
                self._data = data.values
                self._group_assignments = group_assignments.values
            self._index = data.columns

        else:  # get arrays from input values
            self._index = None  # inputs were not all indexed pandas objects

            if isinstance(data, np.ndarray):
                self._data = data
            elif isinstance(data, pd.DataFrame):
                self._data = data.values
            else:
                raise ValueError('data must be a pd.DataFrame or np.ndarray')

            try:
                group_assignments = np.array(group_assignments)
            except:
                raise ValueError('group_assignments must be convertible to a np.ndarray')

            idx = np.argsort(group_assignments)
            if not np.array_equal(idx, np.arange(idx.shape[0])):
                self._data = data[idx, :]
                self._group_assignments = group_assignments[idx]
            else:
                self._data = data
                self._group_assignments = group_assignments

        # set nan values = 0
        self._data[np.isnan(self._data) | np.isinf(self._data)] = 0

        # get points to split the array, create slicers for each group
        self._split_indices = np.where(np.diff(self._group_assignments))[0] + 1

        if not 0 < alpha <= 1:
            raise ValueError('Parameter alpha must fall within the interval (0, 1].')
        self._alpha = alpha

    @property
    def data(self):
        return self._data

    @property
    def group_assignments(self):
        return self._group_assignments

    @property
    def alpha(self):
        return self._alpha

    @property
    def tree(self):
        return self._tree

    @property
    def results(self):
        return self._results

    class Node:
        """Simple hashable node class that for tree leaves"""
        __slots__ = ['node_id', 'left', 'right', 'distance', 'n_daughters']

        def __init__(self, node_id, left=None, right=None, distance=None,
                     n_daughters=None):
            self.node_id = node_id
            self.left = left
            self.right = right
            self.distance = distance
            self.n_daughters = n_daughters

        def __repr__(self):
            return 'Node id: {!s}, left: {!s}, right: {!s}, distance: {!s}'.format(
                self.node_id, self.left, self.right, self.distance)

        def __hash__(self):
            return hash(self.node_id)

    class Tree:

        def __init__(self, linkage, data=None):
            """
            Simple tree class to represent the hierarchical clustering of data

            :param linkage: result of hierarchical clustering. Used to build the tree
            :param data: np.ndarray, of cells x genes, and basis of the clustering.
              Examined for differential expression.

            """


            # store expression data
            self.linkage = linkage
            self.data = data

            # set up tree
            root_id = np.max(linkage[:, :2]) + 1  # node names
            current_id = root_id
            self.tree = {}
            self.n_nodes = linkage.shape[0]

            # construct tree
            for row in linkage[::-1, ]:
                self.tree[current_id] = ExpressionTree.Node(current_id, *row)
                current_id -= 1

            # add leaf nodes
            for i in np.arange(self.n_nodes + 1, dtype=float):
                self.tree[i] = ExpressionTree.Node(i)

            self.root = self.tree[root_id]

        def __getitem__(self, item):
            if isinstance(item, ExpressionTree.Node):
                return self.tree[item.node_id]
            elif isinstance(item, (float, int)):
                return self.tree[item]
            elif item is None:
                return None
            else:
                raise TypeError(
                    'Tree objects support indexing by Node, Node.node_id, or None, not '
                    '{!s}'.format(type(item)))

        def __repr__(self):
            return 'Hierarchical Tree: {:.0f} nodes'.format(self.n_nodes)

        def dfs(self, node=None, visited=None):
            """depth first ordered search of Tree"""
            if visited is None:
                visited = []
            if node is None:
                node = self.root
            if isinstance(node, (int, float)):
                node = self[node]
            visited.append(node)
            for next_ in (self[node.left], self[node.right]):
                if next_ is not None:
                    self.dfs(next_, visited)
            return visited

        def bfs(self, node=None):
            """breadth-first ordered search of tree"""

            if node is None:
                node = self.root
            elif isinstance(node, (int, float)):
                node = self[node]
            elif not isinstance(node, ExpressionTree.Node):
                raise ValueError('node must be a node object, a node_id, or None.')

            visited, queue = [], [node]
            while queue:
                node = queue.pop(0)
                visited.append(node)
                for next_ in [self[node.left], self[node.right]]:
                    if next_ is not None:
                        queue.append(next_)
            return visited

        def agglomerate(self, node):
            """returns all leaves of node"""
            return set(n.node_id for n in self.dfs(node))

        def plot(self, ax, **kwargs):
            """plot a dendrogram on ax"""
            dendrogram(self.linkage, ax=ax, **kwargs)

    @staticmethod
    def resampled_welchs_t(
            a, b, n_iter=100, n_cells=1000, downsample_value_function='median',
            alpha=0.01):
        """
        Return t-statistic from resampled populations a and b

        :param a: n cells x k genes np.ndarray
        :param b: m cells x k genes np.ndarray
        :param n_iter: number of times to resample a and b
        :param downsample_value_function: options: 'mean', 'median' or float percentile
          in the interval (0, 1]. Function to apply to find the value that a & b should
          be sampled to. Default: median. A and B are subjected to these functions, and
          the minimum value between the two is selected.
        """

        # get downsample function
        if isinstance(downsample_value_function, float):
            if 0 < downsample_value_function <= 1:
                downsample_value_function = partial(np.percentile, q=downsample_value_function * 100)
            else:
                raise ValueError('If using a float value for downsample_value_function, the value must fall in the '
                                 'interval (0, 1].')
        elif downsample_value_function is 'mean':
            downsample_value_function = np.mean
        elif downsample_value_function is 'median' or downsample_value_function is None:
            downsample_value_function = np.median

        # get sampling value
        n = min(downsample_value_function(a.sum(axis=1)), downsample_value_function(b.sum(axis=1))).astype(int)

        # normalize cells and obtain probability vectors
        a_prob = (a / a.sum(axis=1)[:, np.newaxis]).mean(axis=0)
        b_prob = (b / b.sum(axis=1)[:, np.newaxis]).mean(axis=0)

        # determine iterations per process
        n_cpu = multiprocessing.cpu_count() - 1
        if n_iter > n_cpu:
            size = np.array([n_iter // n_cpu] * n_cpu)
            size[:n_iter % n_cpu] += 1
        else:
            size = np.ones((n_iter,))

        # map to a process pool
        a_sampler = partial(_sampling_function, n_molecules=n, theta=a_prob, n_cells=n_cells)
        b_sampler = partial(_sampling_function, n_molecules=n, theta=b_prob, n_cells=n_cells)
        pool = multiprocessing.Pool(n_cpu)

        a_res = pool.map(a_sampler, size)
        b_res = pool.map(b_sampler, size)
        a_mu, a_var = (np.vstack(mats) for mats in zip(*a_res))
        b_mu, b_var = (np.vstack(mats) for mats in zip(*b_res))

        # no nans should be produced by the mean
        assert np.sum(np.isnan(a_mu)) == 0
        assert np.sum(np.isnan(b_mu)) == 0

        # in cases where variance is np.nan, we can safely set the variance to zero since
        # the mean will also be zero; this will reduce singularities caused by one tissue
        # never expressing a protein.
        a_var[np.isnan(a_var)] = 0
        b_var[np.isnan(b_var)] = 0

        numerator = a_mu - b_mu  # (samples, genes)
        denominator = np.sqrt(a_var + b_var)  # (samples, genes)
        statistic = numerator / denominator  # (samples, genes)

        # statistic has NaNs where there are no observations of a or b (DivideByZeroError)
        statistic[np.isnan(statistic)] = 0  # calculate df

        df = 2

        statistic = statistic.mean(axis=0)
        p = t.cdf(np.abs(statistic), df)

        q = multipletests(p, alpha=alpha, method='fdr_tsbh')[1]

        return statistic, q

    def create_hierarchy(self):
        """
        hierarchically cluster cluster centroids and store the results as self.tree
        """
        fmean = partial(np.mean, axis=0)

        cluster_array_views = np.array_split(self.data, self._split_indices, axis=0)
        means = np.vstack(list(map(fmean, cluster_array_views)))
        distance = pdist(means)
        linkage = fastcluster.linkage(distance, method='single')

        self._tree = self.Tree(linkage)

    def compare_hierarchy(self, n_iter=100, central_value_func=None, n_cells=1000):
        """
        carry out differential expression at each stage of the Hierarchy.

        :param n_iter: int, number of iterations
        :param central_value_func: function,
          Indicates that clusters should be downsampled to the median of all populations
          prior to comparison to control for sampling bias. This can be done to check
          clustering stability and to get reliable gene expression differences.
        :param n_cells: int, number of cells to sample at each iteration
        :return: dictionary of pd.DataFrame objects containing results
        """

        # get expression means
        if self.tree is None:
            self.create_hierarchy()

        self._results = {}
        for node in self.tree.bfs():
            l_clusters = np.array(self.tree.agglomerate(node.left))
            r_clusters = np.array(self.tree.agglomerate(node.right))
            l_cells = np.in1d(self.group_assignments, l_clusters)
            r_cells = np.in1d(self.group_assignments, r_clusters)
            alpha_adj = self.alpha / (len(set(self.group_assignments)) - 1)
            res = self.resampled_welchs_t(
                self.data[l_cells],
                self.data[r_cells],
                n_iter=n_iter,
                downsample_value_function=central_value_func,
                n_cells=n_cells,
                alpha=alpha_adj)
            if self._index is not None:
                res = pd.Series(res, index=self._index).sort_values(inplace=False)
            self._results[(node.left.node_id, node.right.node_id)] = res
        return self._results

    def plot_hierarchy(self, ax=None, alternate_labels=None):
        """Plot the dendrogram corresponding to hierarchical clustering of self."""
        if ax is None:
            import matplotlib.pyplot as plt
            ax = plt.gca()
        if self.tree is None:
            self.create_hierarchy()
        dendrogram(self.tree.linkage, ax=ax)
        if alternate_labels and isinstance(alternate_labels, dict):
            xtl = ax.get_xticklabels()
            xtl = list(map(lambda x: alternate_labels[x.get_text()], xtl))
            ax.set_xticklabels(xtl, rotation='vertical')


class DifferentialExpression:

    def __init__(self, data, group_assignments, alpha=0.05):
        """
        Carry out ANOVA between the groups of data

        :param data: n cells x k genes 2d array
        :param group_assignments: n cells 1d vector
        :param alpha: float (0, 1], acceptable type I error
        """
        # make sure group_assignments and data have the same length
        if not data.shape[0] == group_assignments.shape[0]:
            raise ValueError(
                'Group assignments shape ({!s}) must equal the number of rows in data '
                '({!s}).'.format(group_assignments.shape[0], data.shape[0]))

        # todo
        # may want to verify that each group has at least two observations
        # (else variance won't work)

        # store index if both data and group_assignments are pandas objects
        if isinstance(data, pd.DataFrame) and isinstance(group_assignments, pd.Series):
            # ensure assignments and data indices are aligned
            try:
                ordered_assignments = group_assignments[data.index]
                if not len(ordered_assignments) == data.shape[0]:
                    raise ValueError(
                        'Index mismatch between data and group_assignments detected when '
                        'aligning indices. check for duplicates.')
            except:
                raise ValueError('Index mismatch between data and group_assignments.')

            # sort data by cluster assignment
            idx = np.argsort(ordered_assignments.values)
            self.data = data.iloc[idx, :].values
            ordered_assignments = ordered_assignments.iloc[idx]
            self.group_assignments = ordered_assignments.values
            self.index = data.columns

        else:  # get arrays from input values
            self.index = None  # inputs were not all indexed pandas objects

            try:
                data = np.array(data)
            except:
                raise ValueError('data must be convertible to a np.ndarray')

            try:
                group_assignments = np.array(group_assignments)
            except:
                raise ValueError('group_assignments must be convertible to a np.ndarray')

            idx = np.argsort(group_assignments)
            self.data = data[idx, :]
            self.group_assignments = group_assignments[idx]

        self.post_hoc = None
        self.groups = np.unique(group_assignments)

        # get points to split the array, create slicers for each group
        self.split_indices = np.where(np.diff(self.group_assignments))[0] + 1
        # todo is this a faster way of calculating the below anova?
        # self.array_views = np.array_split(self.data, self.split_indices, axis=0)

        if not 0 < alpha <= 1:
            raise ValueError('Parameter alpha must fall within the interval (0, 1].')
        self.alpha = alpha

        self._anova = None

    def anova(self, min_mean_expr=None):
        """
        carry out non-parametric ANOVA across the groups of self.

        :param min_mean_expr: minimum average gene expression value that must be reached
          in at least one cluster for the gene to be considered
        :return:
        """
        if self._anova is not None:
            return self._anova

        # run anova
        f = lambda v: kruskalwallis(*np.split(v, self.split_indices))[1]
        pvals = np.apply_along_axis(f, 0, self.data)  # todo could shunt to a multiprocessing pool

        # correct the pvals
        _, pval_corrected, _, _ = multipletests(pvals, self.alpha, method='fdr_tsbh')

        # store data & return
        if self.index is not None:
            self._anova = pd.Series(pval_corrected, index=self.index)
        else:
            self._anova = pval_corrected
        return self._anova

    def post_hoc_tests(self):
        """
        carries out post-hoc tests between genes with significant ANOVA results using
        Welch's U-test on ranked data.
        """
        if self._anova is None:
            self.anova()

        anova_significant = np.array(self._anova) < 1  # call array in case it is a Series

        # limit to significant data, convert to column-wise ranks.
        data = self.data[:, anova_significant]
        rank_data = np.apply_along_axis(rankdata, 0, data)
        # assignments = self.group_assignments[anova_significant]

        split_indices = np.where(np.diff(self.group_assignments))[0] + 1
        array_views = np.array_split(rank_data, split_indices, axis=0)

        # get mean and standard deviations of each
        fmean = partial(np.mean, axis=0)
        fvar = partial(np.var, axis=0)
        mu = np.vstack(list(map(fmean, array_views))).T  # transpose to get gene rows
        n = np.array(list(map(lambda x: x.shape[0], array_views)))
        s = np.vstack(list(map(fvar, array_views))).T
        s_norm = s / n  # transpose to get gene rows

        # calculate T
        numerator = mu[:, np.newaxis, :] - mu[:, :, np.newaxis]
        denominator = np.sqrt(s_norm[:, np.newaxis, :] + s_norm[:, :, np.newaxis])
        statistic = numerator / denominator

        # calculate df
        s_norm2 = s**2 / (n**2 * n-1)
        numerator = (s_norm[:, np.newaxis, :] + s_norm[:, :, np.newaxis]) ** 2
        denominator = (s_norm2[:, np.newaxis, :] + s_norm2[:, :, np.newaxis])
        df = np.floor(numerator / denominator)

        # get significance
        p = t.cdf(np.abs(statistic), df)  # note, two tailed test

        # calculate fdr correction; because above uses 2-tails, alpha here is halved
        # because each test is evaluated twice due to the symmetry of vectorization.
        p_adj = multipletests(np.ravel(p), alpha=self.alpha, method='fdr_tsbh')[1]
        p_adj = p_adj.reshape(*p.shape)

        phr = namedtuple('PostHocResults', ['p_adj', 'statistic', 'mu'])
        self.post_hoc = phr(p_adj, statistic, mu)

        if self.index is not None:
            p_adj = pd.Panel(
                p_adj, items=self.index[anova_significant], major_axis=self.groups,
                minor_axis=self.groups)
            statistic = pd.Panel(
                statistic, items=self.index[anova_significant], major_axis=self.groups,
                minor_axis=self.groups)
            mu = pd.DataFrame(mu, self.index[anova_significant], columns=self.groups)

        return p_adj, statistic, mu

    def population_markers(self, p_crit=0.0):
        """
        Return markers that are significantly differentially expressed in one
        population vs all others

        :param p_crit: float, fraction populations that may be indistinguishable from the
          highest expressing population for each gene. If zero, each marker gene is
          significantly higher expressed in one population relative to all others.
          If 0.1, 10% of populations may share high expression of a gene, and those
          populations will be marked as expressing that gene.

        """
        if self.post_hoc is None:
            self.post_hoc_tests()

        # get highest mean for each gene
        top_gene_idx = np.argmax(self.post_hoc.mu, axis=1)

        # index p_adj first dimension with each sample, will reduce to 2d genes x samples
        top_gene_sig = self.post_hoc.p_adj[:, top_gene_idx, :]

        # for each gene, count the number of non-significant DE results.
        sig = np.array(top_gene_sig < self.alpha)
        num_sig = np.sum(sig, axis=2)

        # if this is greater than N - 1 * p_crit, discard the gene.
        n = self.post_hoc.p_adj.shape[2] - 1  # number of genes, sub 1 for self
        idx_marker_genes = np.where(num_sig < n * (1 - p_crit))
        marker_genes = sig[idx_marker_genes, :]

        # correctly index these genes
        if self.index:
            pass  # todo fix this

        return marker_genes
