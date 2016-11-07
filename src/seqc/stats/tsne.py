import numpy as np
import pandas as pd
import bhtsne
from seqc.stats.pca import PCA

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
