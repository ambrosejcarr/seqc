__author__ = 'ambrose'

import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.use('AGG')
from scipy.stats import gaussian_kde
from sklearn.preprocessing import PolynomialFeatures
from sklearn.linear_model import LinearRegression
from sklearn.pipeline import Pipeline
import matplotlib.pyplot as plt
import seaborn as sns
from copy import copy
from collections import defaultdict
sns.set_style('ticks')


def _savefig(f, fname):
    if fname.endswith('.jpg') or fname.endswith('.png'):
        f.savefig(fname, dpi=450)
    else:
        f.savefig(fname)


def deobfuscate(df):
    """exchange ObfuscatedTuple classes for regular tuples"""
    if isinstance(df, np.ndarray):
        df = pd.DataFrame(df)
    features = df['features'].apply(lambda x: x.to_tuple())
    positions = df['positions'].apply(lambda x: x.to_tuple())
    df['positions'] = positions
    df['features'] = features
    return df


def visualize_filter_correlation(arr, fname, experiment_name=None):

    # df = _load_array(arr)
    df = arr
    corr_df = df.copy()
    corr_df = corr_df.drop(['features', 'positions', 'cell'], axis=1)
    corr_df['trimmed_bases'][corr_df['trimmed_bases'] > 0] = 1
    corr_df['rmt'][corr_df['rmt'] > 0] = 1
    corr_df['disambiguation_results'][corr_df['disambiguation_results'] < 5] = 1
    corr_df['disambiguation_results'][corr_df['disambiguation_results'] == 5] = 0
    corr_df['n_poly_t'] = [0 if n < 4 else 1 for n in corr_df['n_poly_t']]

    corr = corr_df.corr()

    corr.values[np.triu_indices_from(corr.values)] = 0
    corr.values[np.diag_indices_from(corr.values)] = 0

    labels = [l.replace('_', ' ') for l in corr_df.columns]
    labels[-1] = 'is ambiguous'

    f = plt.figure(figsize=(5, 4))
    with sns.axes_style('whitegrid'):
        sns.heatmap(corr, linewidths=0.5, xticklabels=labels, yticklabels=labels)
        title = 'read characteristic correlation'
        if experiment_name:
            title = experiment_name + ' ' + title
        plt.title(title)
        plt.tight_layout()
        _savefig(f, fname)


def _load_array(arr):
    # open the file
    if isinstance(arr, str):
        df = pd.DataFrame(np.load(arr))
    elif isinstance(arr, np.ndarray):
        df = pd.DataFrame(arr)
    elif isinstance(arr, pd.DataFrame):
        df = arr
    else:
        raise ValueError('parameter "arr" must be the filename of a serialized array,'
                         'a structured array, a recarray or a pandas dataframe')
    return df


def _density_plot(data, x, y, ax, xlab, ylab, arcsinh=True):
    xd = data[x].copy()
    yd = data[y].copy()  # todo | these are a bit of a hack, want a copying slice
    xd[np.isnan(xd)] = 0
    yd[np.isnan(yd)] = 0
    xd[np.isinf(xd)] = 0
    yd[np.isinf(yd)] = 0
    h, xe, ye = np.histogram2d(xd, yd, bins=(50, 50))
    h = h.T  # image flips axes
    if arcsinh:
        h = np.arcsinh(h)
    plt.imshow(h, interpolation='nearest', origin='lower', cmap='viridis')
    if xlab:
        xmin, xmax = ax.get_xlim()
        xtick_positions = np.linspace(xmin, xmax, 5)
        xdmin, xdmax = np.min(xe), np.max(xe)

        # find amount to round
        n = 2
        i = copy(xdmax)
        while i > 1:
            i /= 10
            n -= 1

        xtick_labels = np.round(np.linspace(xdmin, xdmax, 5), n)
        ax.set_xticks(xtick_positions)
        ax.set_xticklabels(xtick_labels, rotation='vertical')

    if ylab:
        ymin, ymax = ax.get_ylim()
        ytick_positions = np.linspace(ymin, ymax, 5)
        ydmin, ydmax = np.min(ye), np.max(ye)

        # find amount to round
        n = 2
        i = copy(ydmax)
        while i > 1:
            i /= 10
            n -= 1

        ytick_labels = np.round(np.linspace(ydmin, ydmax, 5), n)
        ax.set_yticks(ytick_positions)
        ax.set_yticklabels(ytick_labels)


def barcode_summary(arr, aggregation_col_name, save):

    df = _load_array(arr)

    # mask reads failing filters

    aggregation_functions = {
        ('cell', lambda a: np.unique(a).shape[0]),  # number of unique cells
        ('rmt', lambda a: np.unique(a).shape[0]),  # number of unique rmts
        ('n_poly_t', np.mean),
        ('valid_cell', np.mean),  # will be between 0 and 1
        ('trimmed_bases', np.mean),
        ('rev_quality', np.mean),
        ('fwd_quality', np.mean),
        ('features', lambda a: np.unique(a).shape[0]),
        ('positions', lambda a: np.unique(a).shape[0]),
        ('is_aligned', np.sum),  # this is counts
        ('alignment_score', np.mean)
    }

    # group by barcode_type
    grouped = df.groupby(aggregation_col_name)
    barcode_data = grouped.agg(aggregation_functions)

    # drop the column that was aggregated over
    barcode_data = barcode_data.drop(aggregation_col_name, axis=1)

    # create an appropriately sized figure
    n = len(barcode_data.columns)
    width, height = [1.75 * n] * 2
    fig = plt.figure(figsize=(width, height))
    gs = plt.GridSpec(nrows=n, ncols=n)
    ax = []
    sns.set_style('ticks')

    # plot everything using subplot/gridspec
    for row, y in enumerate(barcode_data.columns):
        for col, x in enumerate(barcode_data.columns):
            if col >= row:  # don't fill lower triangular
                continue

            cur_ax = plt.subplot(gs[row, col])
            ax.append(cur_ax)
            if row == n - 1:
                plt.xlabel(x.replace('_', ' ').replace('pct', '%').upper())
                xticks = True
            else:
                xticks = False
                plt.xticks([])
            if col == 0:
                plt.ylabel(y.replace('_', ' ').replace('pct', '%').upper())
                yticks = True
            else:
                yticks = False
                plt.yticks([])

            _density_plot(barcode_data, x, y, cur_ax, xticks, yticks)

            sns.despine()

    plt.suptitle(save.replace('.png', '').replace('_', ' '), fontsize=12)
    gs.tight_layout(fig, rect=[0, 0.03, 1, 0.95])
    plt.savefig(save, dpi=300)


def umi_correction(reads, molecules, save):
    """
    show impact of umi correction. reads and molecules should come from the same dataset
    """

    # set-up the plot
    gs = plt.GridSpec(nrows=2, ncols=2)
    fig = plt.figure(figsize=(6, 6))

    # calculate 1d vectors of mean values for molecules and reads
    if isinstance(reads, pd.DataFrame):
        reads = reads.values
    if isinstance(molecules, pd.DataFrame):
        molecules = molecules.values
    mol_means = np.ravel(np.mean(reads, axis=0))
    read_means = np.ravel(np.mean(molecules, axis=0))

    # carry out polynomial regression
    model = Pipeline([('poly', PolynomialFeatures(degree=3)),
                      ('linear', LinearRegression(fit_intercept=False))])
    model = model.fit(mol_means[:, np.newaxis], read_means)

    # plot non-linear relationship
    ax1 = plt.subplot(gs[0, 0])
    beta = model.named_steps['linear'].coef_
    x = mol_means
    y = read_means
    y_hat = beta[0] + beta[1] * x + beta[2] * x ** 2 + beta[3] * x ** 3
    xmin, xmax = np.min(mol_means), np.max(mol_means)
    f = np.linspace(xmin, xmax, 1000)
    fline = beta[0] + beta[1] * f + beta[2] * f ** 2 + beta[3] * f ** 3

    # get density
    xy = np.vstack([x, y])
    z = gaussian_kde(xy)(xy)
    i = np.argsort(z)
    sort_x, sort_y, z = x[i], y[i], np.arcsinh(z[i])

    ax1.scatter(sort_x, sort_y, edgecolor='none', c=z, cmap='viridis')
    ax1.plot(f, fline, linestyle='--', color='indianred')
    plt.xlabel('molecules per gene')
    plt.ylabel('reads per gene')
    plt.title('Polynomial Regression')
    sns.despine()

    # ----
    # plot residuals from the fit
    ax2 = plt.subplot(gs[0, 1])
    residuals = y - y_hat
    xy = np.vstack([x, residuals])
    z = gaussian_kde(xy)(xy)
    i = np.argsort(z)
    sort_x, residuals, z = x[i], residuals[i], np.arcsinh(z[i])

    # plot residuals
    ax2.scatter(sort_x, residuals, edgecolor='none', c=z, cmap='viridis')
    cmax = np.max(sort_x)
    ax2.hlines(0, 0, cmax, color='red', linestyle='dashed')
    plt.xlim((0, cmax))
    plt.xlabel('Molecules per Gene')
    plt.ylabel('Residual')
    plt.title('Residual Plot')
    sns.despine()

    # ----
    # show residuals as a function of data magnitude
    ax3 = plt.subplot(gs[1, 0])
    pct_residuals = np.abs(y - y_hat) / y
    xy = np.vstack([x, residuals])
    z = gaussian_kde(xy)(xy)
    i = np.argsort(z)
    sort_x, pct_residuals, z = x[i], pct_residuals[i], np.arcsinh(z[i])

    # plot residuals
    ax3.scatter(sort_x, pct_residuals, edgecolor='none', c=z, cmap='viridis')
    cmax = np.max(sort_x)
    ymax = np.max(pct_residuals)
    plt.xlim((0, cmax))
    plt.ylim((0, ymax + min(1, ymax + 0.05)))
    plt.xlabel('Molecules per Gene')
    plt.ylabel('Residual / Reads per gene)')
    plt.title('Normalized Residual Plot')
    sns.despine()

    # ----
    # show comparison of fano factors
    ax4 = plt.subplot(gs[1, 1])
    idx = np.where(np.sum(molecules > 0, axis=0) > 2)
    read_fano = np.var(reads[:, idx], axis=0) / np.mean(reads[:, idx], axis=0)
    mols_fano = np.var(molecules[:, idx], axis=0) / np.mean(molecules[:, idx], axis=0)

    # get density
    xy = np.vstack([np.ravel(read_fano), np.ravel(mols_fano)])
    z = gaussian_kde(xy)(xy)
    i = np.argsort(z)
    x, y, z = mols_fano[:, i], read_fano[:, i], z[i]

    # plot stuff
    ax4.scatter(x, y, edgecolor='none', color=z, cmap='viridis')
    cmin, cmax = np.min(np.hstack([x, y])), np.max(np.hstack([x, y]))
    x_y = np.linspace(cmin, cmax, 1000)
    ax4.plot(x_y, x_y, 'r--')
    plt.xlim((cmin, cmax))
    plt.ylim((cmin, cmax))
    plt.xlabel(r'Fano ($\frac{\sigma^2}{\mu}$) molecules / gene', fontsize=12)
    plt.ylabel(r'Fano ($\frac{\sigma^2}{\mu}$) reads / gene', fontsize=12)
    plt.title('Normalied Variance of Reads \nvs. Molecules per Gene', fontsize=12)
    sns.despine()

    gs.tight_layout(fig)
    plt.savefig(save, dpi=300)


def _vertical_marginal(data, ax):
    marginals = data.sum(axis=0)
    with sns.axes_style('ticks'):
        idx = np.arange(len(marginals))
        width = 1
        height = marginals
        ax.bar(idx, height, width, color=plt.cm.viridis(0), edgecolor='w')
        ylim = ax.get_ylim()
        ax.set_yticks(ylim)
        ax.set_xticks([])
        sns.despine(bottom=False, right=True, top=True, left=False, ax=ax)


def _horizontal_marginal(data, ax):
    # get data
    marginals = data.sum(axis=1)
    with sns.axes_style('ticks'):
        idx = np.arange(len(marginals))
        width = 1
        height = marginals
        ax.barh(idx, height, width, color=plt.cm.viridis(0), edgecolor='w')
        xlim = ax.get_xlim()
        # ax.xaxis.tick_top()
        ax.set_xticks(xlim)
        ax.set_yticks([])
        sns.despine(bottom=False, right=True, left=False, top=True, ax=ax)


def error_rate(arr, fname, experiment_name=None):
    """plot error rates with marginal distributions for each experiment type"""

    df = _load_array(arr)

    df[np.isnan(df)] = 0
    df[np.isinf(df)] = 0

    # set up figure
    width, height = 7, 6
    f = plt.figure(figsize=(width, height))
    nrows, ncols = 20, 20
    gs = plt.GridSpec(nrows=nrows, ncols=ncols)

    # set main axis
    hm_ax = plt.subplot(gs[5:20, 0:15])
    pc = plt.pcolor(df * 100, cmap='viridis', edgecolors='w', vmax=4.0)

    # set ylabels
    labels = ['%s to %s' % (v[0], v[1]) for v in list(df.index)]
    hm_ax.set_yticks(np.arange(df.shape[0]) + 0.5, minor=False)
    hm_ax.set_yticklabels(labels, minor=False)
    plt.ylabel('Error')
    for t in hm_ax.yaxis.get_major_ticks():
        t.tick10n = False
        t.tick20n = False

    # set xlabels
    xstep = 4 if df.shape[1] >= 15 else 2
    hm_ax.set_xticks(np.arange(0, df.shape[1] + 1, xstep) + 0.5)
    hm_ax.set_xticklabels(np.arange(0, df.shape[1] + 1, xstep))
    plt.xlabel('Position')
    for t in hm_ax.xaxis.get_major_ticks():
        t.tick10n = False
        t.tick20n = False

    # turn off frame & grid
    hm_ax.set_frame_on(False)
    hm_ax.grid(False)

    # plot vertical marginals
    with sns.axes_style('ticks'):
        vm_ax = plt.subplot(gs[1:4, 0:15])
        _vertical_marginal(df, vm_ax)

    # plot horizontal marginals
    with sns.axes_style('ticks'):
        hm_ax = plt.subplot(gs[5:20, 15:18])
        _horizontal_marginal(df, hm_ax)

    # set cmap axis
    cm_ax = plt.subplot(gs[5:20, 18:20])
    plt.colorbar(pc, cax=cm_ax)  # fraction=0.045, pad=0.1)

    # save figure
    # gs.tight_layout(f)
    if experiment_name:
        plt.suptitle('Cell Barcode Error Rates for %s' % experiment_name)
    else:
        plt.suptitle('Cell Barcode Error Rates')
    plt.savefig(fname, dpi=300)


def multimap_rate_vs_read_length(h5db):
    """display how multimapping rate decreases with read length"""
    # assume that h5db is structured as follows:
    _ = h5db  # dummy
    raise NotImplementedError


def _pie(ax, sizes, explode, labels, colors=None):
    """general pie chart formatting"""

    # set colors
    if colors is not None:
        colors = plt.cm.Set2(np.linspace(0, 1, len(sizes)))

    # reformat labels (change '_' to ' ', first letter capital.)

    p, t, a = ax.pie(sizes, explode=explode, labels=labels, colors=colors,
                     autopct='%1.1f%%', shadow=False, startangle=90)
    ax.axis('equal')
    for i in range(len(a)):
        if float(a[i].get_text()[:-1]) > 1.5:
            t[i].set_size(7)
            a[i].set_size(7)
        else:
            t[i].set_text('')
            a[i].set_text('')

    sns.despine(top=True, bottom=True, left=True, right=True)
    return ax


def molecule_yeild(arr, metadata, fname, experiment_name=None):
    """get molecular yield information for exp_name"""

    df = _load_array(arr)

    raise NotImplementedError

    with open(metadata, 'r') as f:
        short_filtered_reads = int(f.read().strip().split('\t')[-1])

    # set up figure
    f = plt.figure()
    gs = plt.GridSpec(nrows=2, ncols=2)

    ################################ INITIAL FILTERING ##################################

    trimmed = df['trimmed']
    too_short = df['too_short']
    full_length = df['total_reads'] - df['trimmed']

    ax1 = plt.subplot(gs[0, 0])
    colors = plt.cm.RdYlGn([0.1, 0.45, 0.9])
    _pie(ax1, [trimmed, too_short, full_length], [0, 0, 0.05],
         ['Trimmed', 'Too Short', 'Full Length'], colors)

    #################################### ALIGNMENT ######################################
    not_aligned = df['unmapped_reads']
    genomic = df['aligned_genomic']
    no_cell = df['no_cell']  # todo | should filter this during error correction?
    no_rmt = df['no_rmt']
    multimapped = df['mmapped_reads']
    useful_reads = df['total_reads']

    ax2 = plt.subplot(gs[0, 1])
    greens = plt.cm.Greens([0.5, 0.7])
    reds = plt.cm.Reds([0.2, 0.4, 0.6, 0.8])
    colors = np.vstack([reds, greens])
    values = [not_aligned, genomic, no_cell, no_rmt, multimapped, useful_reads]
    labels = ['Unaligned', 'Genomic', 'No Cell Barcode', 'No RMT', 'Multi-aligned',
              'Unique-aligned']
    _pie(ax2, values, [0, 0, 0, 0, 0.05, 0.05], labels, colors)

    ######################### DISAMBIGUATION & ERROR CORRECTION #########################
    no_model = df['no_model']
    errors_corrected = df['errors_corrected']
    not_disambiguated = df['not_disambiguated']
    partially_disambiguated = df['partially_disambiguated']
    disambiguated = df['disambiguated']
    unambiguous = df['total_molecules'] - df['disambiguated']

    ax3 = plt.subplot(gs[1, 0])
    greens = plt.cm.Greens([0.5, 0.7])
    yellows = plt.cm.YlOrBr([0.1, 0.2])
    reds = plt.cm.Reds([0.3, 0.5])
    colors = np.vstack([reds, yellows, greens])
    sizes = [no_model, errors_corrected, not_disambiguated, partially_disambiguated,
             disambiguated, unambiguous]
    labels = ['no_model', 'errors_corrected', 'not_disambiguated',
              'partially_disambiguated', 'disambiguated', 'unambiguous']
    _pie(ax3, sizes, [0, 0, 0, 0, 0.05, 0.05], labels, colors)

    ################################## CELL SELECTION ###################################
    ax4 = plt.subplot(gs[1, 1])
    labels = ['molecules_in_subthreshold_cells', 'molecules_in_cells']
    values = [df[k] for k in labels]
    explode = [0, 0.05]
    colors = ['r', 'b']
    _pie(ax4, values, explode, labels, colors)

    # mark the final graph with a yield percentage
    fastq_filter_yeild = df['total_reads'] / (df['total_reads'] + df['trimmed'])
    alignment_yield = (df['mmapped_reads'] + df['unique_reads']) / df['total_reads']
    filter_yeild = (df['disambiguated'] + unambiguous) / (
        df['partially_disambiguated'] + df['not_disambiguated'] + df['no_model'] +
        df['errors_corrected'] + df['disambiguated'] + unambiguous)
    cell_yield = df['molecules_in_cells'] / (df['molecules_in_subthreshold_cells'] +
                                               df['molecules_in_cells'])
    yeilds = [fastq_filter_yeild, alignment_yield, filter_yeild, cell_yield]

    gs.tight_layout(f)

    if fname.endswith('.jpg') or fname.endswith('.png'):
        plt.savefig(fname, dpi=450)
    else:
        plt.savefig(fname)

    return f, yeilds


def comparative_alignment(data, experiment_labels, fname):

    height, width = 4, 8
    sns.set_style('ticks')
    f = plt.figure(figsize=(width, height))
    gs = plt.GridSpec(nrows=1, ncols=20)
    ax = list()
    bars = defaultdict(list)

    # get a colormap to color each experiment
    colors = plt.cm.Set2(np.linspace(0, 1, len(data)))

    def plot_bar(category_labels, axes, label, bar_heights):
        n = len(category_labels)
        ind = np.arange(n)  # locations of bars
        _width = 1 / (n + 1)
        experiments = [[e[c] for c in category_labels] for e in data]
        for i, e in enumerate(experiments):
            bar_heights[label].append(axes.bar(
                ind + (_width * i), e, width, color=colors[i]
            ))
        axes.set_ylabel(label)
        axes.set_xticks(ind + ((width * len(experiments)) / 2))
        axes.set_xticklabels([c.replace('_', ' ') for c in categories],
                             rotation='vertical')
        sns.despine()
        return bars

    # plot alignment rates
    ax.append(plt.subplot(gs[0, 0:5]))
    categories = ['mmap_rate', 'uniq_rate', 'unmapped_rate']
    bars = plot_bar(categories, ax[0], 'Alignment Rate (%)', bars)

    # plot alignment numbers
    ax.append(plt.subplot(gs[0, 5:11]))
    categories = ['total_reads', 'unique_reads', 'mmapped_reads', 'unmapped_reads']
    bars = plot_bar(categories, ax[1], 'Alignments (counts)', bars)
    ax[1].yaxis.get_major_formatter().set_powerlimits((0, 1))

    # plot alignment errors
    ax.append(plt.subplot(gs[0, 11:16]))
    categories = ['mismatch_rate', 'insertion_rate', 'deletion_rate']
    bars = plot_bar(categories, ax[2], 'Error Rates (% per base)', bars)

    # plot alignment errors
    ax.append(plt.subplot(gs[0, 16:20]))
    categories = ['deletion_size', 'insertion_size']
    bars = plot_bar(categories, ax[3], 'Error sizes (bases)', bars)

    bars_for_legend = [b[0] for b in bars['Alignment Rate']]
    plt.legend(handles=bars_for_legend, labels=experiment_labels, loc=2,
               bbox_to_anchor=(1.05, 1))

    plt.suptitle('Library Construction Alignment Information')
    gs.tight_layout(f, rect=[0, 0, 0.85, 0.95])
    plt.savefig(fname, dpi=300)
    return f
