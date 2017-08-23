import warnings
from collections import OrderedDict
from seqc.sequence.fastq import Reader
from math import floor
import numpy as np
import pandas as pd
from sklearn.mixture import GaussianMixture
from sklearn.linear_model import LinearRegression
from seqc.exceptions import EmptyMatrixError
from seqc.sparse_frame import SparseFrame
from numpy.linalg import LinAlgError
import seqc.plot
from seqc import log


def estimate_min_poly_t(fastq_files: list, platform) -> int:
    """
    estimate the minimum size of poly-t tail that should be present on a properly captured
    molecule's forward read. If multiple fastq files are passed, the minimum value across
    all files will be returned

    :param fastq_files: list of fastq filenames
    :param platform: the platform used to generate this library
    :return: int minimum number of poly-t expected from a valid capture primer
    """
    min_vals = []
    primer_length = platform.primer_length()
    if primer_length is None:
        raise RuntimeError(
            'provided platform does not have a defined primer length, and thus the '
            'min_poly_t parameter cannot be estimated. Please provide --min-poly-t '
            'explicitly in process_experiment.py.')
    for f in fastq_files:
        mean = Reader(f).estimate_sequence_length()[0]
        available_nucleotides = max(0, mean - primer_length)
        min_vals.append(floor(min(available_nucleotides * .8, 20)))
    return min(min_vals)


def low_count(molecules, is_invalid, plot=False, ax=None):
    """
    updates is_invalid to reflect cells whose molecule counts are below the inflection
    point of an ecdf constructed from cell molecule counts. Typically this reflects cells
    whose molecule counts are approximately <= 100.

    :param molecules: scipy.stats.coo_matrix, molecule count matrix
    :param is_invalid:  np.ndarray(dtype=bool), declares valid and invalid cells
    :param bool plot: if True, plot a summary of the filter
    :param ax: Must be passed if plot is True. Indicates the axis on which to plot the
      summary.
    :return: is_invalid, np.ndarray(dtype=bool), updated valid and invalid cells
    """

    # copy, sort, and normalize molecule sums
    ms = np.ravel(molecules.tocsr()[~is_invalid, :].sum(axis=1))
    idx = np.argsort(ms)[::-1]  # largest cells first
    norm_ms = ms[idx] / ms[idx].sum()  # sorted, normalized array

    # identify inflection point from second derivative
    cms = np.cumsum(norm_ms)
    d1 = np.diff(pd.Series(cms).rolling(10).mean()[10:])
    d2 = np.diff(pd.Series(d1).rolling(10).mean()[10:])
    try:
        # throw out an extra 5% of cells from where the inflection point is found.
        # these cells are empirically determined to have "transition" library sizes
        # that confound downstream analysis
        inflection_pt = np.min(np.where(np.abs(d2) == 0)[0])
        inflection_pt = int(inflection_pt * .9)
    except ValueError as e:
        if e.args[0] == ('zero-size array to reduction operation minimum which has no '
                         'identity'):
            log.notify('Low count filter passed-through; too few cells to estimate '
                       'inflection point.')
            return is_invalid  # can't estimate validity
        else:
            raise

    vcrit = ms[idx][inflection_pt]

    is_invalid = is_invalid.copy()
    is_invalid[ms < vcrit] = True

    if plot and ax:
        cms /= np.max(cms)  # normalize to one
        ax.plot(np.arange(len(cms))[:inflection_pt], cms[:inflection_pt])
        ax.plot(np.arange(len(cms))[inflection_pt:], cms[inflection_pt:], c='indianred')
        ax.hlines(cms[inflection_pt], *ax.get_xlim(), linestyle='--')
        ax.vlines(inflection_pt, *ax.get_ylim(), linestyle='--')
        ax.set_xticklabels([])
        ax.set_xlabel('putative cell')
        ax.set_ylabel('ECDF (Cell Size)')
        ax.set_title('Cell Size')
        ax.set_ylim((0, 1))
        ax.set_xlim((0, len(cms)))

    return is_invalid


def low_coverage(molecules, reads, is_invalid, plot=False, ax=None, filter_on=True):
    """
    Fits a two-component gaussian mixture model to the data. If a component is found
    to fit a low-coverage fraction of the data, this fraction is set as invalid. Not
    all datasets contain this fraction.

    For best results, should be run after filter.low_count()

    :param molecules: scipy.stats.coo_matrix, molecule count matrix
    :param reads: scipy.stats.coo_matrix, read count matrix
    :param is_invalid:  np.ndarray(dtype=bool), declares valid and invalid cells
    :param bool plot: if True, plot a summary of the filter
    :param ax: Must be passed if plot is True. Indicates the axis on which to plot the
      summary.
    :param filter_on: indicate whether low coverage filter is on
    :return: is_invalid, np.ndarray(dtype=bool), updated valid and invalid cells
    """
    ms = np.ravel(molecules.tocsr()[~is_invalid, :].sum(axis=1))
    rs = np.ravel(reads.tocsr()[~is_invalid, :].sum(axis=1))

    if ms.shape[0] < 10 or rs.shape[0] < 10:
        log.notify(
            'Low coverage filter passed-through; too few cells to calculate '
            'mixture model.')
        return is_invalid

    # get read / cell ratio, filter out low coverage cells
    ratio = rs / ms

    # fit two GMMs on one and two modes
    col_ratio = ratio[:, np.newaxis]
    gmm1 = GaussianMixture(n_components=1)
    gmm2 = GaussianMixture(n_components=2)
    gmm1.fit(col_ratio)
    gmm2.fit(col_ratio)

    if filter_on:
    # check if adding a second component is necessary; if not, filter is pass-through
        filter_on = gmm2.bic(col_ratio) / gmm1.bic(col_ratio) < 0.95

    if filter_on:
        res = gmm2.predict(col_ratio)

        # Molecule sum means
        means = [np.mean(ms[res == 0]), np.mean(ms[res == 1])]
        failing = np.where(res == np.argmin(means))[0]

        # set smaller mean as invalid
        is_invalid = is_invalid.copy()
        is_invalid[np.where(~is_invalid)[0][failing]] = True
    else:
        res, means = None, None

    if plot and ax:
        logms = np.log10(ms)
        try:
            seqc.plot.scatter.continuous(logms, ratio, colorbar=False, ax=ax, s=3)
        except LinAlgError:
            warnings.warn('SEQC: Insufficient number of cells to calculate density for '
                          'coverage plot')
            ax.scatter(logms, ratio, s=3)
        ax.set_xlabel('log10(molecules)')
        ax.set_ylabel('reads / molecule')
        ax.set_title('Coverage')
        xmin, xmax = np.min(logms), np.max(logms)
        ymax = np.max(ratio)
        ax.set_xlim((xmin, xmax))
        ax.set_ylim((0, ymax))
        seqc.plot.xtick_vertical(ax=ax)

        # plot 1d conditional densities of two-component model
        # todo figure out how to do this!!

        # plot the discarded cells in red, like other filters
        if filter_on:
            ax.scatter(
                logms[res == np.argmin(means)], ratio[res == np.argmin(means)],
                s=4, c='indianred')

    return is_invalid


def high_mitochondrial_rna(molecules, gene_ids, is_invalid, max_mt_content=0.2,
                           plot=False, ax=None):
    """
    Sets any cell with a fraction of mitochondrial mRNA greater than max_mt_content to
    invalid.

    :param molecules: scipy.stats.coo_matrix, molecule count matrix
    :param gene_ids: np.ndarray(dtype=str) containing string gene identifiers
    :param is_invalid:  np.ndarray(dtype=bool), declares valid and invalid cells
    :param max_mt_content: float, maximum percentage of reads that can come from
      mitochondria in a valid cell
    :param bool plot: if True, plot a summary of the filter
    :param ax: Must be passed if plot is True. Indicates the axis on which to plot the
      summary.
    :return: is_invalid, np.ndarray(dtype=bool), updated valid and invalid cells
    """
    # identify % genes that are mitochondrial
    mt_genes = np.fromiter(map(lambda x: x.startswith('MT-'), gene_ids), dtype=np.bool)
    mt_molecules = np.ravel(molecules.tocsr()[~is_invalid, :].tocsc()[:, mt_genes].sum(
        axis=1))
    ms = np.ravel(molecules.tocsr()[~is_invalid, :].sum(axis=1))
    ratios = mt_molecules / ms

    failing = ratios > max_mt_content
    is_invalid = is_invalid.copy()
    is_invalid[np.where(~is_invalid)[0][failing]] = True

    if plot and ax:
        if ms.shape[0] and ratios.shape[0]:
            try:
                seqc.plot.scatter.continuous(ms, ratios, colorbar=False, ax=ax, s=3)
            except LinAlgError:
                log.notify('Inadequate number of cells or MT gene abundance to plot MT '
                           'filter, no visual will be produced, but filter has been '
                           'applied.')
                return is_invalid
        else:
            return is_invalid  # nothing else to do here
        if np.sum(failing) != 0:
            ax.scatter(ms[failing], ratios[failing], c='indianred', s=3)  # failing cells
        xmax = np.max(ms)
        ymax = np.max(ratios)
        ax.set_xlim((0, xmax))
        ax.set_ylim((0, ymax))
        ax.hlines(max_mt_content, *ax.get_xlim(), linestyle='--', colors='indianred')
        ax.set_xlabel('total molecules')
        ax.set_ylabel('fraction mitochondrial\nmolecules')
        ax.set_title(
            'MT-RNA Fraction: {:.2}%'.format(np.sum(failing) / len(failing) * 100))
        seqc.plot.xtick_vertical(ax=ax)

    return is_invalid


def low_gene_abundance(molecules, is_invalid, plot=False, ax=None):
    """
    Fits a linear model to the relationship between number of genes detected and number
    of molecules detected. Cells with a lower than expected number of detected genes
    are set as invalid.

    :param molecules: scipy.stats.coo_matrix, molecule count matrix
    :param is_invalid:  np.ndarray(dtype=bool), declares valid and invalid cells
    :param bool plot: if True, plot a summary of the filter
    :param ax: Must be passed if plot is True. Indicates the axis on which to plot the
      summary.
    :return: is_invalid, np.ndarray(dtype=bool), updated valid and invalid cells
    """

    ms = np.ravel(molecules.tocsr()[~is_invalid, :].sum(axis=1))
    genes = np.ravel(molecules.tocsr()[~is_invalid, :].getnnz(axis=1))
    x = np.log10(ms)[:, np.newaxis]
    y = np.log10(genes)

    if not (x.shape[0] or y.shape[0]):
        return is_invalid

    # get line of best fit
    with warnings.catch_warnings():  # ignore scipy LinAlg warning about LAPACK bug.
        warnings.simplefilter('ignore')
        regr = LinearRegression()
        regr.fit(x, y)

    # mark large residuals as failing
    yhat = regr.predict(x)
    residuals = yhat - y
    failing = residuals > .15

    is_invalid = is_invalid.copy()
    is_invalid[np.where(~is_invalid)[0][failing]] = True

    if plot and ax:
        m, b = regr.coef_, regr.intercept_
        try:
            seqc.plot.scatter.continuous(x, y, ax=ax, colorbar=False, s=3)
        except LinAlgError:
            log.notify('Inadequate number of cells to plot low coverage filter no visual '
                       'will be produced, but filter has been applied.')
            return is_invalid
        xmin, xmax = np.min(x), np.max(x)
        ymin, ymax = np.min(y), np.max(y)
        lx = np.linspace(xmin, xmax, 200)
        ly = m * lx + b
        ax.plot(lx, np.ravel(ly), linestyle='--', c='indianred')
        ax.scatter(x[failing], y[failing], c='indianred', s=3)
        ax.set_ylim((ymin, ymax))
        ax.set_xlim((xmin, xmax))
        ax.set_xlabel('molecules (cell)')
        ax.set_ylabel('genes (cell)')
        ax.set_title('Low Complexity')
        seqc.plot.xtick_vertical(ax=ax)

    return is_invalid


def create_filtered_dense_count_matrix(
        molecules: SparseFrame, reads: SparseFrame, max_mt_content=0.2, plot=False,
        figname=None, filter_mitochondrial_rna: bool=True, filter_low_count: bool=True, filter_low_coverage: bool=True):
    """
    filter cells with low molecule counts, low read coverage, high mitochondrial content,
    and low gene detection. Returns a dense pd.DataFrame of filtered counts, the total
    original number of molecules (int), the number of molecules lost with each filter
    (dict), and the number of cells lost with each filter (dict).

    :param filter_mitochondrial_rna: if True, run the mitochondrial RNA filter.
    :param filter_low_count: if True, run the low count cell filter.
    :param filter_low_coverage: if True, run the low coverage filter.
    :param molecules: SparseFrame
    :param reads: SparseFrame
    :param max_mt_content: the maximum percentage of mitochondrial RNA that is
    :param plot: if True, plot filtering summaries.
    :param figname: if plot is True, name of the figure to save.
    :return: (pd.DataFrame, int, dict, dict)
    """

    cells_lost = OrderedDict()
    molecules_lost = OrderedDict()

    if not molecules.columns.dtype.char == 'U':
        if molecules.sum().sum() == 0:
            raise EmptyMatrixError('Matrix is empty, cannot create dense matrix')
        else:
            raise RuntimeError(
                'non-string column names detected. Please convert column names into '
                'string gene symbols before calling this function.')
    if not isinstance(max_mt_content, float):
        raise TypeError('Parameter max_mt_content must be of type float.')
    if not 0 <= max_mt_content <= 1:
        raise ValueError('Parameter max_mt_content must be in the interval [0, 1]')

    # set data structures and original molecule counts
    molecules_data = molecules.data
    reads_data = reads.data
    molecules_columns = molecules.columns
    is_invalid = np.zeros(molecules_data.shape[0], np.bool)
    total_molecules = np.sum(molecules_data.sum(axis=1))

    def additional_loss(new_filter, old_filter, data_matrix):
        new_cell_loss = np.sum(new_filter) - np.sum(old_filter)
        data_matrix = data_matrix.tocsr()
        total_molecule_loss = data_matrix[new_filter].sum().sum()
        old_molecule_loss = data_matrix[old_filter].sum().sum()
        new_molecule_loss = total_molecule_loss - old_molecule_loss
        return new_cell_loss, new_molecule_loss

    if plot:
        fig = seqc.plot.FigureGrid(4, max_cols=2)
        ax_count, ax_cov, ax_mt, ax_gene = iter(fig)  # get axes
    else:
        fig, ax_count, ax_cov, ax_mt, ax_gene = [None] * 5  # dummy figure

    # filter low counts
    if filter_low_count:
        count_invalid = low_count(molecules_data, is_invalid, plot, ax_count)
        cells_lost['low_count'], molecules_lost['low_count'] = additional_loss(
            count_invalid, is_invalid, molecules_data)
    else:
        count_invalid = is_invalid


    # filter low coverage
    cov_invalid = low_coverage(molecules_data, reads_data, count_invalid, plot, ax_cov,filter_low_coverage)
    cells_lost['low_coverage'], molecules_lost['low_coverage'] = additional_loss(
        cov_invalid, count_invalid, molecules_data)

    # filter high_mt_content if requested
    if filter_mitochondrial_rna:
        mt_invalid = high_mitochondrial_rna(
            molecules_data, molecules_columns, cov_invalid, max_mt_content, plot, ax_mt)
        cells_lost['high_mt'], molecules_lost['high_mt'] = additional_loss(
            mt_invalid, cov_invalid, molecules_data)
    else:
        mt_invalid = cov_invalid

    # filter low gene abundance
    gene_invalid = low_gene_abundance(molecules_data, mt_invalid, plot, ax_gene)
    cells_lost['low_gene_detection'], molecules_lost[
        'low_gene_detection'] = additional_loss(
        gene_invalid, mt_invalid, molecules_data)

    # construct dense matrix
    dense = molecules_data.tocsr()[~gene_invalid, :].todense()
    nonzero_gene_count = np.ravel(np.array(dense.sum(axis=0) != 0))
    dense = dense[:, nonzero_gene_count]
    dense = pd.DataFrame(
        dense,
        index=molecules.index[~gene_invalid],
        columns=molecules.columns[nonzero_gene_count])

    # describe cells
    cell_description = dense.sum(axis=1).describe()

    if plot:
        fig.tight_layout()
        fig.savefig(figname, dpi=300, transparent=True)

    return dense, total_molecules, molecules_lost, cells_lost, cell_description
