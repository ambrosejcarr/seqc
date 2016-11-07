import os
import re
import glob
import warnings
import shlex
import numpy as np
import pandas as pd
import multiprocessing
from sklearn.neighbors import NearestNeighbors
import time
import bhtsne
from scipy.special import expit
from scipy.sparse import csr_matrix, find
from scipy.sparse.linalg import eigs
from numpy.linalg import norm
from statsmodels.sandbox.stats.multicomp import multipletests
from scipy.stats.mstats import kruskalwallis, rankdata
from scipy.stats import t
from functools import partial
from collections import namedtuple
import subprocess

# add some functions to the namespace
from seqc.ttest import bootstrap_t

# safe bug in external dependency, suppress warning.
warnings.filterwarnings(action="ignore", module="scipy", message="^internal gelsd")

class GraphDiffusion:
    def __init__(self, knn=10, normalization='smarkov', epsilon=1,
                 n_diffusion_components=10):
        """
        Run diffusion maps on the data. This implementation is based on the
        diffusion geometry library in Matlab:
        https://services.math.duke.edu/~mauro/code.html#DiffusionGeom and was implemented
        by Pooja Kathail

        :param knn: Number of neighbors for graph construction to determine distances
          between cells
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
    def bimarkov(W, max_iters=100, abs_error=0.00001, **kwargs):
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

    @staticmethod
    def sFokkerPlanck(D, N, W):
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
        return {'operator': T, 'eigval': D, 'eigvec': V, 'weights': W}


class NormalizeCells:
    """Contains row normalization techniques for equalizing cell sampling.

    :method size: divides by library size
    """

    @staticmethod
    def size(data: np.ndarray, percentile=None):
        """
        Normalize data such that each cell has an identical library size

        :param data: n observation x k feature array
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
            data = data.div(cell_sums, axis=0).mul(np.percentile(cell_sums, percentile),
                                                   axis=0)
        return data


class ScaleFeatures:
    """
    Contains column normalization techniques for equalizing feature sampling

    :method unit_size: scale each feature between 0 and 1
    :method standardize: center the mean of each feature at zero and scale by standard
      deviation
    :method unit_length: scale each feature so that its euclidean length is 1
    """

    @staticmethod
    def unit_size(data: np.ndarray, copy=True, dtype=np.float32):
        """
        scales features in data between 0 and 1


        :param data: n observation x k feature array
        :param copy: bool, if True, creates a copy of data. Otherwise carries out
          standardization in-place
        :param dtype: optional datatype parameter, options [np.float32, np.float64,
          float], defaults to np.float32 to save space.
        :return: n x k scaled array
        """

        if dtype not in [float, np.float32, np.float64]:
            raise TypeError('dtype must be one of np.float32, np.float64 or float')
        if copy:
            data = data.copy().astype(dtype)
        else:
            data = data.astype(dtype)
        d_min = np.amin(data, axis=0)
        data -= d_min
        d_max = np.amax(data, axis=0)

        # divide by max unless max is zero, then divide by 1
        data /= np.amax(np.array([d_max, np.ones_like(d_max)]), axis=0)
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


class TSNE:

    def __init__(self, n_components: int=2, run_pca: bool=False,
                 n_pca_components: int=20, fillna: float=None, **kwargs):
        """
        t-stochastic neighbor embedding


        :param normalize: if True, scales features to unit size
        :param run_pca: if True, runs PCA on the input data and runs tSNE on the
          components retained by PCA.
        :param n_components: number of tSNE components to return
        :param n_pca_components: number of components to which data should be projected,
          if run_pca is True
        :param fillna: fills np.nan values with this float value
        :param kwargs:  additional keyword arguments to pass tsne

        :method fit_transform: fits the tSNE model to data and returns the transformed
          result

        """

        self.run_pca = run_pca
        self.n_components = n_components
        self.n_pca_components = n_pca_components
        self.kwargs = kwargs
        self.tsne = None
        self.pca = None
        self.fillna = fillna

    def fit_transform(self, data: np.ndarray or pd.DataFrame) -> None:
        """
        fit the tSNE model to data given the parameters provided during
         initialization and transform the output

        :param data: n observation x k feature data array
        :return np.ndarray or pd.DataFrame: tsne results
        """
        if isinstance(data, pd.DataFrame):
            data_ = data.values
        else:
            data_ = data

        if self.fillna is not None:
            data_[np.where(np.isnan(data_))] = self.fillna
            data_[np.where(np.isinf(data_))] = self.fillna
        if self.run_pca:
            self.pca = PCA(n_components=self.n_pca_components)
            data_ = self.pca.fit_transform(data_)

        res = bhtsne.tsne(data_.astype(float), dimensions=self.n_components, **self.kwargs)

        if isinstance(data, pd.DataFrame):
            self.tsne = pd.DataFrame(res, index=data.index)
        else:
            self.tsne = res
        return self.tsne


class PCA:

    def __init__(self, n_components=30):
        """
        construct a model for Principle Component Analysis

        :param n_components: number of principle components to retain

        :property eigenvalues: stores the eigenvalues computed by fit()
        :property loadings: stores the eigenvectors of the pca decomposition computed by
          fit()
        :method fit: fit the model to the data
        :method transform: project the data onto a subset of the principle components
          (default: all components other than the first)
        :method fit_transform: fit and transform the data, returning the projected result
        """
        self.n_components = n_components
        self.loadings = None
        self.eigenvalues = None

    def fit(self, data: np.ndarray or pd.DataFrame, fillna=0):
        """
        Fit the model to data

        :param data: n observation x k feature data array
        :param fillna: fill np.NaN values with this value. If None, will not fill.
        :return: None
        """

        if isinstance(data, pd.DataFrame):
            X = data.values
        elif isinstance(data, np.ndarray):
            X = data
        else:
            raise TypeError('data must be a pd.DataFrame or np.ndarray')

        if fillna is not None:
            X[np.where(np.isnan(X))] = fillna
            X[np.where(np.isinf(X))] = fillna

        # Compute covariance matrix
        if X.shape[1] < X.shape[0]:
            C = np.cov(X, rowvar=False)
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
          data, if None, uses all components except for the first
        :return: np.ndarray containing transformed data
        """

        if components is None:
            components = np.arange(1, self.n_components)

        projected = np.dot(data, self.loadings[:, components])
        if isinstance(data, pd.DataFrame):
            return pd.DataFrame(projected, index=data.index, columns=components)
        else:
            return projected

    def fit_transform(self, data: np.ndarray or pd.DataFrame, n_components=None) -> \
            np.ndarray or pd.DataFrame:
        """
        Fit the model to data and transform the data using the fit model

        :param data:  n observation x k feature data array
        :param n_components:  number of components to retain when transforming
          data
        :return np.ndarray or pd.DataFrame: transformed data
        """

        self.fit(data)
        return self.transform(data, components=n_components)


class correlation:
    """Fast vectorized correlation methods

    :method vector(x, y): correlate each column in y with a vector in x
    :method map(x, y): correlate each column of x with each column in y
    :method eigv(evec, data): get pairwise correlations of eigenvectors
      with columns of data
    """

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
    def map(x: np.ndarray, y: np.ndarray):
        """Correlate each column of x with each column of y

        :param x: np.array; shape N x T.
        :param y: np.array; shape M x T.
        :returns: np.array; shape N x M in which each element is a correlation
                            coefficient.
        """
        assert(x.shape[1] == y.shape[1])
        n = x.shape[1]
        x_diff = x - x.mean(axis=-1)[:, None]
        y_diff = y - y.mean(axis=-1)[:, None]
        x_std = x.std(axis=-1)
        y_std = y.std(axis=-1)
        return np.dot(x_diff, y_diff.T) / (n * x_std[:, np.newaxis] * y_std)

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
        elif isinstance(data, np.ndarray):
            D = data
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
        if isinstance(data, pd.DataFrame):
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
        'N present in CB:\t{cell_N} ({bad_cellN}%)\n'
        'Insufficient poly(T):\t{poly_t} ({bad_polyt}%)\n'
        'High dust score:\t{dust} ({bad_dust}%)\n'
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
        cell_N = summary['cell_N']
        dust = summary['dust']
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
        n_bad = genomic + phix + no_cell + rmt_N + cell_N + poly_t + dust
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
        bad_cellN = round((cell_N/n_bad) * 100, 1)
        bad_polyt = round((poly_t/n_bad) * 100, 1)
        bad_dust = round((dust/n_bad) * 100, 1)
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
            cell_N=cell_N, bad_cellN=bad_cellN, dust=dust, bad_dust=bad_dust,
            prop_trans=prop_trans, tot_mc=tot_mc, mols_lost=mols_lost,
            cells_lost=cells_lost, cell_desc=cell_desc)
        return output


class smoothing:
    """Data smoothing kernels

    :method kneighbors: transforms each observation (row) of data by setting it
      equal to the average of its k-nearest neighbors
    """

    @staticmethod
    def kneighbors(data: np.array or pd.DataFrame, n_neighbors=50, pca=None, **kwargs):
        """
        Smooth gene expression values by setting the expression of each gene in each
        cell equal to the mean value of itself and its n_neighbors

        :param data: np.ndarray | pd.DataFrame; genes x cells array
        :param n_neighbors: int; number of neighbors to smooth over
        :param pca: dimensionality reduced matrix, knn will be run on this and applied
          to data (runs much faster)
        :param kwargs: keyword arguments to pass sklearn.NearestNeighbors
        :return: np.ndarray | pd.DataFrame; same as input
        """

        if isinstance(data, pd.DataFrame):
            data_ = data.values
        elif isinstance(data, np.ndarray):
            data_ = data
        else:
            raise TypeError("data must be a pd.DataFrame or np.ndarray")

        knn = NearestNeighbors(
            n_neighbors=n_neighbors,
            n_jobs=multiprocessing.cpu_count() - 1,
            **kwargs)

        if pca is not None:
            knn.fit(pca)
            inds = knn.kneighbors(pca, return_distance=False)
        else:
            knn.fit(data_)
            inds = knn.kneighbors(data_, return_distance=False)

        # smoothing creates large intermediates; break up to avoid memory errors
        pieces = []
        num_partitions = np.round(data_.shape[0] / 2000) + 1
        if num_partitions > 2:  # 2 partitions produces start + end, need a third to split
            sep = np.linspace(0, data_.shape[0] + 1, num_partitions, dtype=int)
            for start, end in zip(sep, sep[1:]):
                pieces.append(data_[inds[start:end, :], :].mean(axis=1))
            res = np.vstack(pieces)
        else:
            res = data_[inds, :].mean(axis=1)

        if isinstance(data, pd.DataFrame):
            res = pd.DataFrame(res, index=data.index, columns=data.columns)

        return res


class JavaGSEA:

    def __init__(self, correlations):
        """initialize a gsea object
        :param correlations: a pandas series of correlations in the range of [-1, 1] whose
          index contains gene names

        :method linear_scale: method to linearly scale a vector to lie on the interval
          [-1, 1]
        :method logisitc_scale: method to scale a vector by the logistic function to lie
          on the interval [-1, 1]
        :method run: run GSEA
        """
        if not isinstance(correlations, pd.Series):
            raise TypeError('correlations must be a pandas series')
        if not ((np.min(correlations) >= -1) & (np.max(correlations) <= 1)):
            raise RuntimeError(
                'input correlations were not contained within the interval [-1, 1]. '
                'Please use JavaGSEA.linear_scale() or JavaGSEA.logistic_scale() to '
                'scale values to this interval before running.')
        self._correlations = correlations.sort_values()
        self._rnk = None
        self._output_stem = os.environ['TMPDIR'] + 'gsea_corr_{!s}'.format(
            np.random.randint(0, 1000000))
        self._results = {}

    @property
    def correlations(self):
        return self._correlations

    @correlations.setter
    def correlations(self):
        raise RuntimeError('Please create a new object to compare different correlations')

    @property
    def results(self):
        return self._results

    @staticmethod
    def linear_scale(data: pd.Series) -> pd.Series:
        """scale input vector to interval [-1, 1] using a linear scaling
        :return correlations: pd.Series, data scaled to the interval [-1, 1]
        """
        data = data.copy()
        data -= np.min(data, axis=0)
        data /= np.max(data, axis=0) / 2
        data -= 1
        return data

    @staticmethod
    def logistic_scale(data: pd.Series) -> pd.Series:
        """scale input vector to interval [-1, 1] using a sigmoid scaling
        :return correlations: pd.Series, data scaled to the interval [-1, 1]
        """
        return pd.Series((expit(data.values) * 2) - 1, index=data.index)

    def _save_rank_file(self) -> None:
        """save the correlations to a .rnk file in TMPDIR

        :return: None
        """
        self._rnk = self._output_stem + '.rnk'
        df = pd.DataFrame(self._correlations).fillna(0)
        df.to_csv(self._rnk, sep='\t', header=False)

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

    def run(self, gmt_file: (str, str) = None, output_stem: str = None) -> (
            pd.DataFrame, pd.DataFrame):
        """
        Helper function. Run GSEA on an already-ranked list of corrleations. To see
        available files, leave gmt_file parameter empty

        :param output_stem: location for GSEA output
        :param gmt_file: (organism: str, filename: str) organism and filename of the gmt
          file to use
        :return: positive, negative: pd.DataFrames containing GSEA enrichment results
        """

        if output_stem is None:
            output_stem = self._output_stem

        if output_stem.find('-') > -1:
            raise ValueError('ouptput_stem cannot contain special characters such as -')

        out_dir, out_prefix = os.path.split(output_stem)
        os.makedirs(out_dir, exist_ok=True)

        if self._rnk is None:
            self._save_rank_file()

        if not gmt_file:
            self._gmt_options()
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
            ''.format(user=os.path.expanduser('~'), rnk=self._rnk, out_prefix=out_prefix,
                      out_dir=out_dir, gmt_file=gmt_file))

        # Call GSEA
        p = subprocess.Popen(cmd, stderr=subprocess.PIPE)
        _, err = p.communicate()

        # find the file that GSEA created
        if err:
            print(err.decode())
            return
        else:
            pattern = '{p}.GseaPreranked.[0-9]*'.format(p=out_prefix)
            files = os.listdir(out_dir)
            folder = None
            for f in files:
                mo = re.match(pattern, f)
                if mo:
                    folder = os.environ['TMPDIR'] + mo.group(0)
        if folder is None:
            raise RuntimeError(
                'seqc.JavaGSEA was not able to recover the output of the Java '
                'executable. This likely represents a bug.')

        # recover information from run
        names = ['size', 'es', 'nes', 'p', 'fdr_q', 'fwer_p', 'rank_at_max',
                 'leading_edge']
        pos = pd.DataFrame.from_csv(glob.glob(folder + '/gsea*pos*xls')[0],
                                    sep='\t').iloc[:, :-1]
        pos.drop(['GS<br> follow link to MSigDB', 'GS DETAILS'], axis=1, inplace=True)
        neg = pd.DataFrame.from_csv(glob.glob(folder + '/gsea*neg*xls')[0],
                                    sep='\t').iloc[:, :-1]
        neg.drop(['GS<br> follow link to MSigDB', 'GS DETAILS'], axis=1, inplace=True)
        pos.columns, neg.columns = names, names
        self._results[gmt_file] = {'positive': pos, 'negative': neg}
        return self._results[gmt_file].values()


class DifferentialExpression:

    def __init__(self, data, group_assignments, alpha=0.05):
        """
        Carry out ANOVA between the groups of data

        :param data: n cells x k genes 2d array
        :param group_assignments: n cells 1d vector
        :param alpha: float (0, 1], acceptable type I error
        """
        # make sure group_assignments and data have the same length
        warnings.warn('DeprecationWarning: This function is deprecated.')
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
