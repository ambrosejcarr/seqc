import time
import numpy as np
from scipy.sparse.linalg import eigs
from numpy.linalg import norm
from scipy.sparse import csr_matrix, find
from sklearn.neighbors import NearestNeighbors


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
