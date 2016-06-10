from seqc.sequence.fastq import Reader
from math import floor
import numpy as np
import pandas as pd
from sklearn.mixture import GMM
from sklearn.linear_model import LinearRegression

_primer_lengths = dict(
    in_drop=47,
    in_drop_v2=49,
    drop_seq=20,
    mars1_seq=None,  # mars-seq has not been provided to us as pre-demultiplexed data
    mars2_seq=None
)


def estimate_min_poly_t(fastq_files: list, platform: str) -> int:
    """
    estimate the minimum size of poly-t tail that should be present on a properly captured molecule's
     forward read. If multiple fastq files are passed, the minimum value across all files will be
     returned

    :param fastq_files: list of fastq filenames
    :param platform: the platform used to generate this library
    :return: int minimum number of poly-t expected from a valid capture primer
    """
    min_vals = []
    try:
        primer_length = _primer_lengths[platform]
    except KeyError:
        raise ValueError('provided value {} is not a valid argument for platform.'.format(platform))
    if primer_length is None:
        raise RuntimeError('provided platform does not have a defined primer length, and thus the '
                           'min_poly_t parameter cannot be estimated. Please provide --min-poly-t '
                           'explicitly in process_experiment.py.')
    for f in fastq_files:
        mean = Reader(f).estimate_sequence_length()[0]
        available_nucleotides = max(0, mean - primer_length)
        min_vals.append(floor(min(available_nucleotides * .8, 20)))
    return min(min_vals)


def low_count(molecules, is_invalid):
    """
    :param molecules:
    :param is_invalid:
    :return: np.array(dtype=bool), listing whether each cell is valid
    """

    # copy, sort, and normalize molecule sums
    ms = np.ravel(molecules[~is_invalid, :].sum(axis=1))
    ms = ms[ms > 0]  # throw out zero values
    idx = np.argsort(ms)
    norm_ms = ms[idx] / ms[idx].sum()  # sorted, normalized array

    # identify inflection point form second derivative
    cms = np.cumsum(norm_ms)
    d1 = np.diff(pd.Series(cms).rolling(10).mean()[10:])
    d2 = np.diff(pd.Series(d1).rolling(10).mean()[10:])
    inflection_pt = np.min(np.where(np.abs(d2) == 0)[0])
    vcrit = ms[idx][inflection_pt]
    print(vcrit)

    is_invalid[ms < vcrit] = True

    return is_invalid


def low_coverage(molecules, reads, is_invalid):
    """
    For best results, should be run after filter.low_count()

    :param molecules:
    :param reads:
    :param is_valid:
    :return:
    """
    ms = molecules[~is_invalid, :].sum(axis=1)
    rs = reads[~is_invalid, :].sum(axis=1)

    # get read / cell ratio, filter out low coverage cells
    ratio = rs / ms

    # fit two GMMs on one and two modes
    col_ratio = ratio[:, np.newaxis]
    gmm1 = GMM(n_components=1)
    gmm2 = GMM(n_components=2)
    gmm1.fit(col_ratio)
    gmm2.fit(col_ratio)

    # check if adding a second component is necessary; if not, filter is pass-through
    if gmm2.bic(col_ratio) / gmm1.bic(col_ratio) < 0.85:
        res = gmm2.fit_predict(col_ratio)
        failing = np.where(res == np.argmin(gmm2.means_))[0]

        # set smaller mean as invalid
        is_invalid[np.where(~is_invalid)[0][failing]] = True

    return is_invalid


def high_mitochondrial_rna(molecules, is_invalid, max_mt_content=0.2):
    """
    :param molecules:
    :param is_invalid:
    :param max_mt_content:
    :return:
    """
    # identify % genes that are mitochondrial
    mt_genes = molecules.columns[molecules.columns.str.contains('MT-')]
    mt_molecules = molecules[~is_invalid, mt_genes].sum(axis=1)
    ms = molecules[~is_invalid, :].sum(axis=1)
    ratios = mt_molecules / ms

    failing = ratios.index[ratios > max_mt_content]
    is_invalid[np.where(~is_invalid)[0][failing]] = True

    return is_invalid


def low_gene_abundance(molecules, is_invalid):
    """
    :param molecules:
    :param is_invalid:
    :return:
    """

    ms = molecules[~is_invalid, :].sum(axis=1)
    genes = np.array(molecules[~is_invalid, :] > 0).sum(axis=1)
    x = np.log10(ms)[:, np.newaxis]
    y = np.log10(genes)

    # get line of best fit
    regr = LinearRegression()
    regr.fit(x, y)

    # mark large residuals as failing
    yhat = regr.predict(x)
    residuals = yhat - y
    failing = residuals > .15

    is_invalid[np.where(~is_invalid)[0][failing]] = True

    return is_invalid
