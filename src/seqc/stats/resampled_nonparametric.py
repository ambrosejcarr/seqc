import os
from functools import partial
from multiprocessing import Pool
from contextlib import closing
from itertools import repeat
import numpy as np
import numpy.ma as ma
import pandas as pd
from scipy.stats.mstats import count_tied_groups, rankdata
from scipy.stats.mstats import kruskalwallis as _kruskalwallis
from scipy.special import erfc
from statsmodels.sandbox.stats.multicomp import multipletests


def get_memory():
    """
    """
    return os.sysconf('SC_PAGE_SIZE') * os.sysconf('SC_PHYS_PAGES') / (1024 ** 3)


def _mannwhitneyu(x, y, use_continuity=True):
    """
    Computes the Mann-Whitney statistic
    Missing values in `x` and/or `y` are discarded.
    Parameters
    ----------
    x : ndarray,
        Input, vector or observations x features matrix
    y : ndarray,
        Input, vector or observations x features matrix. If matrix, must have
        same number of features as x
    use_continuity : {True, False}, optional
        Whether a continuity correction (1/2.) should be taken into account.
    Returns
    -------
    statistic : float
        The Mann-Whitney statistic
    approx z : float
        The normal-approximated z-score for U.
    pvalue : float
        Approximate p-value assuming a normal distribution.
    """
    if x.ndim == 1 and y.ndim == 1:
        x, y = x[:, np.newaxis], y[:, np.newaxis]
    ranks = rankdata(np.concatenate([x, y]), axis=0)
    nx, ny = x.shape[0], y.shape[0]
    nt = nx + ny
    U = ranks[:nx].sum(0) - nx * (nx + 1) / 2.

    mu = (nx * ny) / 2.
    u = np.amin([U, nx*ny - U], axis=0)  # get smaller U by convention

    sigsq = np.ones(ranks.shape[1]) * (nt ** 3 - nt) / 12.

    for i in np.arange(len(sigsq)):
        ties = count_tied_groups(ranks[:, i])
        sigsq[i] -= np.sum(v * (k ** 3 - k) for (k, v) in ties.items()) / 12.
    sigsq *= nx * ny / float(nt * (nt - 1))

    if use_continuity:
        z = (U - 1 / 2. - mu) / np.sqrt(sigsq)
    else:
        z = (U - mu) / np.sqrt(sigsq)

    prob = erfc(abs(z) / np.sqrt(2))
    return np.vstack([u, z, prob]).T


def find_sampling_value(group_data, percentile):
    """

    :param group_data:
    :param int percentile:
    :return:
    """
    return min(np.percentile(g.sum(axis=1), percentile) for g in group_data)


def normalize(data, downsample_value, upsample=False, labels=None):
    """
    :param data:
    :param downsample_value: value to normalize cell counts to. In current implementation,
        a small number of cells (10%) are upsampled to this value.
    :param upsample: if False, all observations with size < downsample_value are excluded.
        if True, those cells are upsampled to downsample_value.
    :return:
    """
    obs_size = data.sum(axis=1)
    if not upsample:
        keep = obs_size >= downsample_value
        data = data[keep, :]
        if labels is not None:
            labels = labels[keep]
    norm = (data * downsample_value) / data.sum(axis=1)[:, np.newaxis]
    if labels is not None:
        return norm, labels
    else:
        return norm


def _draw_sample(normalized_data, n):
    """
    :param normalized_data:
    :param n:
    """
    np.random.seed()
    idx = np.random.randint(0, normalized_data.shape[0], n)
    sample = normalized_data[idx, :]
    p = np.random.sample(sample.shape)  # round samples probabilistically

    return np.floor(sample) + (sample % 1 > p).astype(int)


def _mw_sampling_function(norm_data, n_cell):
    """
    :param norm_data:
    :param n_cell:
    :return:
    """
    a, b = (_draw_sample(d, n_cell) for d in norm_data)
    return _mannwhitneyu(a, b)  # dim = (n_genes, 3)


def confidence_interval(z):
    """

    :param z:
    :return:
    """
    return np.percentile(z, [2.5, 97.5], axis=0).T


def mannwhitneyu(
        x, y, n_iter=50, sampling_percentile=10, alpha=0.05, verbose=False,
        upsample=False):
    """
    :param x: observations by features array or DataFrame (ndim must be 2, although there
        needn't be more than one feature)
    :param y: observations by features array or DataFrama. Features must be the same as x
    :param n_iter: number of times to sample x and y
    :param sampling_percentile: percentile to downsample to. observations with row sums
        lower than this value will be excluded
    :param alpha: significance threshold for FDR correction
    :param verbose: if True, report number of cells sampled in each iteration and the
        integer value to which cells are downsampled
    :param upsample: if False, cells with size lower than sampling_percentile are
        discarded. If True, those cells are upsampled.
    :return pd.DataFrame: DataFrame with columns:
        U: median u-statistic over the n_iter iterations of the test
        z_approx: median approximate tie-corrected z-score for the mann-whitney U-test
        z_lo: lower bound, 95% confidence interval over z
        z_hi: upper bound, 95% confidence interval over z
        p: p-value for z_approx
        q: FDR-corrected q-value over all tests in output, using two-stage BH-FDR.
    """

    # do some sanity checks on input data
    if isinstance(x, pd.DataFrame) and isinstance(y, pd.DataFrame):
        assert np.array_equal(x.columns, y.columns)
        labels = x.columns
        x = x.values
        y = y.values
    elif x.ndim > 1:
        assert x.shape[1] == y.shape[1]
        labels = None
    else:
        labels = None

    # calculate sampling values
    v = find_sampling_value([x, y], sampling_percentile)
    norm_data = [normalize(d, v, upsample) for d in [x, y]]
    n_cell = min(d.shape[0] for d in norm_data)
    sampling_function = partial(_mw_sampling_function, n_cell=n_cell)

    if verbose:  # report sampling values
        print('sampling %d cells (with replacement) per iteration' % n_cell)
        print('sampling %d molecules per cell' % v)

    with closing(Pool()) as pool:
        results = pool.map(sampling_function, repeat(norm_data, n_iter))

    results = np.stack(results)  # u, z, p

    ci = confidence_interval(results[:, :, 1])
    results = pd.DataFrame(
        data=np.concatenate([np.median(results, axis=0), ci], axis=1),
        index=labels,
        columns=['U', 'z_approx', 'p', 'z_lo', 'z_hi'])

    # add multiple-testing correction
    results['q'] = multipletests(results['p'], alpha=alpha, method='fdr_tsbh')[1]

    # remove low-value genes whose median sampling value is -inf
    neginf = np.isneginf(results['z_approx'])
    results.ix[neginf, 'z_lo'] = np.nan
    results.ix[neginf, 'z_approx'] = 0
    results.ix[neginf, ['p', 'q']] = 1.

    results = results[['U', 'z_approx', 'z_lo', 'z_hi', 'p', 'q']].sort_values('q')
    results.iloc[:, 1:4] = np.round(results.iloc[:, 1:4], 2)

    return results


def _kw_sampling_function(data, splits, n_cell):
    data = [_draw_sample(d, n_cell) for d in np.split(data, splits)]
    return _kruskal(data)


def _kruskal(data):
    """
    Compute the Kruskal-Wallis H-test for independent samples
    Parameters
    ----------
    sample1, sample2, ... : array_like
       Two or more arrays with the sample measurements can be given as
       arguments.
    Returns
    -------
    statistic : float
       The Kruskal-Wallis H statistic, corrected for ties
    pvalue : float
       The p-value for the test using the assumption that H has a chi
       square distribution
    Notes
    -----
    For more details on `kruskal`, see `stats.kruskal`.
    """
    results = []
    for i in np.arange(data[0].shape[1]):
        args = [d[:, i] for d in data]
        try:
            results.append(_kruskalwallis(*args))
        except ValueError:
            results.append([0, 1.])
    return results


def category_to_numeric(labels):
    """transform categorical labels to a numeric array"""
    labels = np.array(labels)
    if np.issubdtype(labels.dtype, np.integer):
        return labels
    else:
        cats = np.unique(labels)
        map_ = dict(zip(cats, np.arange(cats.shape[0])))
        return np.array([map_[i] for i in labels])


def kruskalwallis(
        data, labels, n_iter=50, sampling_percentile=10, alpha=0.05, verbose=False,
        upsample=False):
    """
    :param data: np.ndarray or pd.DataFrame of observations x features
    :param labels: observation labels for categories to be compared
    :param n_iter: number of times to sample x and y
    :param sampling_percentile: percentile to downsample to. observations with row sums
        lower than this value will be excluded
    :param alpha: significance threshold for FDR correction
    :param verbose: if True, report number of cells sampled in each iteration and the
        integer value to which cells are downsampled
    :param upsample: if False, cells with size lower than sampling_percentile are
        discarded. If True, those cells are upsampled.
    :return pd.DataFrame: DataFrame with columns:
        H: median u-statistic over the n_iter iterations of the test
        z_approx: median approximate tie-corrected z-score for the mann-whitney U-test
        z_lo: lower bound, 95% confidence interval over z
        z_hi: upper bound, 95% confidence interval over z
        p: p-value for z_approx
        q: FDR-corrected q-value over all tests in output, using two-stage BH-FDR.
    """

    if isinstance(data, pd.DataFrame):
        features = data.columns
        data = data.values
    elif isinstance(data, np.ndarray):
        features = None
    else:
        raise ValueError('data must be a np.ndarray or pd.DataFrame, not %s' %
                         repr(type(data)))

    # if labels are not numeric, transform to numeric categories
    labels = category_to_numeric(labels)
    if not labels.shape[0] == data.shape[0]:
        raise ValueError('labels (shape=%s) must match dimension 0 of data (shape=%s)' %
                         (repr(labels.shape), repr(labels.data)))

    idx = np.argsort(labels)
    data = data[idx, :]  # will copy
    labels = labels[idx]

    splits = np.where(np.diff(labels))[0] + 1

    # calculate sampling values and downsample data
    v = find_sampling_value(np.split(data, splits), sampling_percentile)
    norm_data, labels = normalize(data, v, upsample, labels)

    splits = np.where(np.diff(labels))[0] + 1  # rediff, norm_data causes loss

    n_cell = min(d.shape[0] for d in np.split(norm_data, splits))
    sampling_function = partial(_kw_sampling_function, n_cell=n_cell, splits=splits)

    if verbose:  # report sampling values
        print('sampling %d cells (with replacement) per iteration' % n_cell)
        print('sampling %d molecules per cell' % v)

    with closing(Pool()) as pool:
        results = pool.map(sampling_function, repeat(norm_data, n_iter))

    results = np.stack(results)  # H, p

    ci = confidence_interval(results[:, :, 0])  # around H
    results = pd.DataFrame(
        data=np.concatenate([np.median(results, axis=0), ci], axis=1),
        index=features,
        columns=['H', 'p', 'H_lo', 'H_hi'])

    results['q'] = multipletests(results['p'], alpha=alpha, method='fdr_tsbh')[1]
    results = results[['H', 'H_lo', 'H_hi', 'p', 'q']]
    return results

