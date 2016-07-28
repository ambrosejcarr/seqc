import numpy as np


def jsd(p, q) -> float:
    """Jensen Shannon distance of two variables normalized variables p and q

    Note that if p and q are not normalized, this function will not return a proper distance,
    so matrices should be normalized prior to use

    use with sklearn.NearestNeighbors:

    >>> from sklearn.neighbors import NearestNeighbors
    >>> # set some dummy variables
    >>> data = np.random.random((100, 100))
    >>> data = data / data.sum(axis=1)[:, np.newaxis]  # norm rows
    >>> assert(np.all(np.array(data.sum(axis=1) == 1)))3
    >>> k = 10
    >>>
    >>> nn = NearestNeighbors(k=k, metric='pyfunc', algorithm='ball_tree',
    >>>          metric_params={'func': jsd})
    >>> nn.fit(data)

    Parameters
    ----------
    p, q : np.array

    Returns
    -------
    float : kl divergence between p and q
    """
    idx = np.logical_or(p != 0, q != 0)
    p = p[idx]
    q = q[idx]
    m = (p + q) / 2
    return np.sqrt((.5 * kldiv(p, m)) + (.5 * kldiv(q, m)))


def kldiv(x: np.ndarray, m: np.ndarray) -> float:
    """Modified Kullback-Liebler divergence of two variables x and m.

    depends upon normalization done by jsd parent function, namely that (1) there are no zero-
    valued entries in m, and (2) both x and m are probability distributions that sum to 1

    Parameters
    ----------
    x, m : normalized probability vectors

    Returns
    -------
    float : kl divergence between p and q
    """
    return np.nansum(x * np.log2(x / m))