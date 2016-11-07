import numpy as np
import pandas as pd


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
