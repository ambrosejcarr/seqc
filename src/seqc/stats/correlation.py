import numpy as np
import pandas as pd


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
