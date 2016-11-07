import numpy as np
import pandas as pd
from contextlib import closing
from multiprocessing import Pool
from sklearn.cluster import KMeans


def _assign(d):
    """

    :param np.ndarray d: 1d vector of scaled differences
    :return np.ndarray: 1d boolean gene-enrichment assignment vector
    """
    km = KMeans(n_clusters=2)
    km.fit(d[:, np.newaxis])
    assignments = km.labels_.astype(bool)
    if np.argmax(km.cluster_centers_) == 0:
        return assignments
    else:
        return ~assignments


def g_test(data, labels, log=False):
    """

    :param pd.DataFrame data:
    :param labels:
    :param log:
    :return:
    """

    if log:
        data = np.log(data + 1)

    data = pd.DataFrame(data.values / data.values.sum(axis=1)[:, np.newaxis],
                        index=labels, columns=data.columns)

    # calculate data that are useful for determining observed and expected values
    gene_sums = data.sum(axis=0)
    grouped = data.groupby(axis=0, level=0)  # group only once
    category_sizes = grouped.size()
    category_fractions = category_sizes / category_sizes.sum()  # normalize

    # get observed, expected
    expected = pd.DataFrame(
        data=np.dot(category_fractions.values[:, np.newaxis],
                    gene_sums.values[np.newaxis, :]),
        index=category_sizes.index,
        columns=gene_sums.index)
    observed = grouped.sum()

    # scaled ratios are used in both g-test, and partitioning of expressed vs. not
    logratio = np.log(observed / expected)
    logratio.values[~np.isfinite(logratio.values)] = 0
    scaled_diff = observed * logratio

    g = 2 * np.sum(scaled_diff, axis=0)  # g-test

    # todo only assign significant values
    # todo calculate significance
    with closing(Pool()) as pool:
        assignments = pool.map(_assign, scaled_diff.values.T)

    assignments = pd.DataFrame(
        data=np.vstack(assignments).T,
        index=category_sizes.index,
        columns=data.columns
    )

    return g, assignments
