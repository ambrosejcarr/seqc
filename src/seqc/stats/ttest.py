import numpy as np
from collections.abc import Callable
from multiprocessing import Pool, cpu_count
from functools import partial
from contextlib import closing
from scipy.stats import t
import pandas as pd
from statsmodels.sandbox.stats.multicomp import multipletests


def estimate_multinomial(x):
    """estimate empirical multinomial expectation for a set of cells with each cell
     normalized to contribute equally to the expectation.

    :param np.ndarray x: cell x gene array containing expression data
    :return np.ndarray: multinomial expectation over genes of x
    """
    return (x / x.sum(axis=1)[:, np.newaxis]).mean(axis=0)


def assert_input_non_negative(*args):
    """
    :param [np.ndarray] args: input numpy arrays
    :return None:
    """
    if any(np.any(np.less(a, 0)) for a in args):
        raise ValueError('input data must be non-negative')


def _sampling_function(n_iter, n_molecules, theta, n_cells):
    """

    :param n_iter:
    :param n_molecules:
    :param theta:
    :param n_cells:
    :return:
    """

    def online_mean_var(nb, mu_b, var_b, na, mu_a, var_a):
        nx = na + nb
        delta = mu_b - mu_a
        mu_x_ = mu_a + delta * nb / nx
        var_x_ = (na * (var_a + mu_a ** 2) + nb * (var_b + mu_b ** 2)) / nx - mu_x_ ** 2
        return nx, mu_x_, var_x_

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


def sample_moments(mult_probability, n_samples, n_cells, n_molecules):
    """sample mean and variance of n_cells, each containing n_molecules. n_samples mean/
    variance pairs are sampled on each call.

    :param mult_probability:
    :param n_samples:
    :param n_cells:
    :param n_molecules:
    :return:
    """

    # parition iterations among available compute cores
    ncpu = cpu_count()
    if n_samples > ncpu:
        samples_per_process = np.array([n_samples // ncpu] * ncpu)
        samples_per_process[:n_samples % ncpu] += 1
    else:
        samples_per_process = np.ones((n_samples,))

    # map iterations across compute cores
    sampler = partial(
        _sampling_function, n_molecules=n_molecules, theta=mult_probability,
        n_cells=n_cells)
    with closing(Pool(ncpu)) as pool:
        results = pool.map(sampler, samples_per_process)
        mu, var = (np.vstack(mats) for mats in zip(*results))

    # all means should be finite
    assert np.sum(np.isnan(mu)) == 0

    # in cases where variance is np.nan, we can safely set the variance to zero since the
    # mean for that tissue will also be zero; this will eliminate singularities caused by
    # one tissue never expressing a protein.
    var[np.isnan(var)] = 0

    return mu, var


def whelch_satterthwaite_df(a_var, b_var, a_n, b_n):
    t1 = a_var.mean(axis=0)
    t2 = b_var.mean(axis=0)
    numerator = (t1 / a_n + t2 / b_n) ** 2
    denominator = t1 ** 2 / (a_n ** 2 * (a_n - 1)) + t2 ** 2 / (b_n ** 2 * (b_n - 1))
    df = numerator / denominator
    return df


def whelchs_t(a_mu, a_var, b_mu, b_var, a_n, b_n):
    """

    :param np.ndarray a_mu:
    :param np.ndarray a_var:
    :param np.ndarray b_mu:
    :param np.ndarray b_var:
    :param int a_n:
    :param int b_n:
    :return float, float: statistic and p-value
    """
    df = whelch_satterthwaite_df(a_var, b_var, a_n, b_n)
    numerator = a_mu - b_mu  # (samples, genes)
    denominator = np.sqrt(a_var + b_var)  # (samples, genes)
    statistic = numerator / denominator  # (samples, genes)

    # statistic has NaNs where there are no observations of a or b (DivideByZeroError)
    statistic[np.isnan(statistic)] = 0
    median_statistic = np.median(np.abs(statistic), axis=0)
    p = (1 - t.cdf(median_statistic, df)) * 2  # p-value
    ci_95 = np.percentile(np.abs(statistic), [2.5, 97.5], axis=0).T

    return median_statistic, p, ci_95


def bootstrap_t(a, b, n_samples=100, n_cells=None, alpha=0.05,
                downsample_value_function=np.median, labels=None):
    """

    :param np.ndarray a:
    :param np.ndarray b:
    :param int n_samples:
    :param int n_cells:
    :param float alpha: acceptable type-I error (default = 0.05)
    :param Callable downsample_value_function: function that identifies the number of
      molecules n to sample from a and b. the sampling number will be the minimum of the
      result across a and b. default = np.median. Other values include np.mean and np.max.
    :param labels: feature labels for columns of a & b
    :return (int, int) statistic, q_val:
    """
    assert_input_non_negative(a, b)
    mult_a = estimate_multinomial(a)
    mult_b = estimate_multinomial(b)

    # get number of molecules to sample
    a_sizes = a.sum(axis=1)
    b_sizes = b.sum(axis=1)
    n_molecules = min(
        map(lambda x: downsample_value_function(x).astype(int), [a_sizes, b_sizes]))

    # set n_cells to the smaller of the two passed samples (e.g. if comparing two sets,
    # one with 130 cells, and one with 1902 cells, n_cells = 130).
    if n_cells is None:
        n_cells = min(a.shape[0], b.shape[0])

    a_mu, a_var = sample_moments(mult_a, n_samples, n_cells, n_molecules)
    b_mu, b_var = sample_moments(mult_b, n_samples, n_cells, n_molecules)

    statistic, p, ci_95 = whelchs_t(a_mu, a_var, b_mu, b_var, a.shape[0], b.shape[0])

    q = multipletests(p, alpha=alpha, method='fdr_tsbh')[1]

    results = pd.DataFrame(
        data=np.vstack([statistic, ci_95.T, p, q]).T,
        index=labels,
        columns=['t', 't_ci95_low', 't_ci95_high', 'p', 'q'])

    return results
