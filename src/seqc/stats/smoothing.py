import numpy as np
import pandas as pd
import multiprocessing
from sklearn.neighbors import NearestNeighbors


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
