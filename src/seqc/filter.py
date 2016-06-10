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
    updates is_invalid to reflect cells whose molecule counts are below the inflection
    point of an ecdf constructed from cell molecule counts. Typically this reflects cells
    whose molecule counts are approximately <= 100.

    :param molecules: scipy.stats.coo_matrix, molecule count matrix
    :param is_invalid:  np.ndarray(dtype=bool), declares valid and invalid cells
    :return: is_invalid, np.ndarray(dtype=bool), updated valid and invalid cells
    """

    # copy, sort, and normalize molecule sums
    ms = np.ravel(molecules.tocsr()[~is_invalid, :].sum(axis=1))
    idx = np.argsort(ms)[::-1]  # largest cells first
    norm_ms = ms[idx] / ms[idx].sum()  # sorted, normalized array

    # identify inflection point form second derivative
    cms = np.cumsum(norm_ms)
    d1 = np.diff(pd.Series(cms).rolling(10).mean()[10:])
    d2 = np.diff(pd.Series(d1).rolling(10).mean()[10:])
    inflection_pt = np.min(np.where(np.abs(d2) == 0)[0])
    vcrit = ms[idx][inflection_pt]

    is_invalid = is_invalid.copy()
    is_invalid[ms < vcrit] = True

    return is_invalid


def low_coverage(molecules, reads, is_invalid):
    """
    Fits a two-component gaussian mixture model to the data. If a component is found
    to fit a low-coverage fraction of the data, this fraction is set as invalid. Not
    all datasets contain this fraction.

    For best results, should be run after filter.low_count()

    :param molecules: scipy.stats.coo_matrix, molecule count matrix
    :param reads: scipy.stats.coo_matrix, read count matrix
    :param is_invalid:  np.ndarray(dtype=bool), declares valid and invalid cells
    :return: is_invalid, np.ndarray(dtype=bool), updated valid and invalid cells
    """
    ms = np.ravel(molecules.tocsr()[~is_invalid, :].sum(axis=1))
    rs = np.ravel(reads.tocsr()[~is_invalid, :].sum(axis=1))

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
        is_invalid = is_invalid.copy()
        is_invalid[np.where(~is_invalid)[0][failing]] = True

    return is_invalid


def high_mitochondrial_rna(molecules, gene_ids, is_invalid, max_mt_content=0.2):
    """
    Sets any cell with a fraction of mitochondrial mRNA greater than max_mt_content to
    invalid.

    :param molecules: scipy.stats.coo_matrix, molecule count matrix
    :param gene_ids: np.ndarray(dtype=str) containing string gene identifiers
    :param is_invalid:  np.ndarray(dtype=bool), declares valid and invalid cells
    :param max_mt_content: float, maximum percentage of reads that can come from
      mitochondria in a valid cell
    :return: is_invalid, np.ndarray(dtype=bool), updated valid and invalid cells
    """
    # identify % genes that are mitochondrial
    mt_genes = np.fromiter(map(lambda x: x.startswith('MT-'), gene_ids), dtype=np.bool)
    mt_molecules = np.ravel(molecules.tocsr()[~is_invalid, :].tocsc()[:, mt_genes].sum(axis=1))
    ms = np.ravel(molecules.tocsr()[~is_invalid, :].sum(axis=1))
    ratios = mt_molecules / ms

    failing = ratios > max_mt_content
    is_invalid = is_invalid.copy()
    is_invalid[np.where(~is_invalid)[0][failing]] = True

    return is_invalid


def low_gene_abundance(molecules, is_invalid):
    """
    Fits a linear model to the relationship between number of genes detected and number
    of molecules detected. Cells with a lower than expected number of detected genes
    are set as invalid.

    :param molecules: scipy.stats.coo_matrix, molecule count matrix
    :param is_invalid:  np.ndarray(dtype=bool), declares valid and invalid cells
    :return: is_invalid, np.ndarray(dtype=bool), updated valid and invalid cells
    """

    ms = np.ravel(molecules.tocsr()[~is_invalid, :].sum(axis=1))
    genes = np.ravel(molecules.tocsr()[~is_invalid, :].getnnz(axis=1))
    x = np.log10(ms)[:, np.newaxis]
    y = np.log10(genes)

    # get line of best fit
    regr = LinearRegression()
    regr.fit(x, y)

    # mark large residuals as failing
    yhat = regr.predict(x)
    residuals = yhat - y
    failing = residuals > .15

    is_invalid = is_invalid.copy()
    is_invalid[np.where(~is_invalid)[0][failing]] = True

    return is_invalid
