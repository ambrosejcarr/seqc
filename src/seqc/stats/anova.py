import warnings
from collections import namedtuple
import numpy as np
import pandas as pd
from functools import partial
from scipy.stats.mstats import kruskalwallis, rankdata
from scipy.stats import t
from statsmodels.sandbox.stats.multicomp import multipletests

class ANOVA:

    def __init__(self, data, group_assignments, alpha=0.05):
        """
        Carry out ANOVA between the groups of data

        :param data: n cells x k genes 2d array
        :param group_assignments: n cells 1d vector
        :param alpha: float (0, 1], acceptable type I error
        """
        # make sure group_assignments and data have the same length
        warnings.warn('DeprecationWarning: This function is deprecated.')
        if not data.shape[0] == group_assignments.shape[0]:
            raise ValueError(
                'Group assignments shape ({!s}) must equal the number of rows in data '
                '({!s}).'.format(group_assignments.shape[0], data.shape[0]))

        # todo
        # may want to verify that each group has at least two observations
        # (else variance won't work)

        # store index if both data and group_assignments are pandas objects
        if isinstance(data, pd.DataFrame) and isinstance(group_assignments, pd.Series):
            # ensure assignments and data indices are aligned
            try:
                ordered_assignments = group_assignments[data.index]
                if not len(ordered_assignments) == data.shape[0]:
                    raise ValueError(
                        'Index mismatch between data and group_assignments detected when '
                        'aligning indices. check for duplicates.')
            except:
                raise ValueError('Index mismatch between data and group_assignments.')

            # sort data by cluster assignment
            idx = np.argsort(ordered_assignments.values)
            self.data = data.iloc[idx, :].values
            ordered_assignments = ordered_assignments.iloc[idx]
            self.group_assignments = ordered_assignments.values
            self.index = data.columns

        else:  # get arrays from input values
            self.index = None  # inputs were not all indexed pandas objects

            try:
                data = np.array(data)
            except:
                raise ValueError('data must be convertible to a np.ndarray')

            try:
                group_assignments = np.array(group_assignments)
            except:
                raise ValueError('group_assignments must be convertible to a np.ndarray')

            idx = np.argsort(group_assignments)
            self.data = data[idx, :]
            self.group_assignments = group_assignments[idx]

        self.post_hoc = None
        self.groups = np.unique(group_assignments)

        # get points to split the array, create slicers for each group
        self.split_indices = np.where(np.diff(self.group_assignments))[0] + 1
        # todo is this a faster way of calculating the below anova?
        # self.array_views = np.array_split(self.data, self.split_indices, axis=0)

        if not 0 < alpha <= 1:
            raise ValueError('Parameter alpha must fall within the interval (0, 1].')
        self.alpha = alpha

        self._anova = None

    def anova(self, min_mean_expr=None):
        """
        carry out non-parametric ANOVA across the groups of self.

        :param min_mean_expr: minimum average gene expression value that must be reached
          in at least one cluster for the gene to be considered
        :return:
        """
        if self._anova is not None:
            return self._anova

        # run anova
        f = lambda v: kruskalwallis(*np.split(v, self.split_indices))[1]
        pvals = np.apply_along_axis(f, 0, self.data)  # todo could shunt to a multiprocessing pool

        # correct the pvals
        _, pval_corrected, _, _ = multipletests(pvals, self.alpha, method='fdr_tsbh')

        # store data & return
        if self.index is not None:
            self._anova = pd.Series(pval_corrected, index=self.index)
        else:
            self._anova = pval_corrected
        return self._anova

    def post_hoc_tests(self):
        """
        carries out post-hoc tests between genes with significant ANOVA results using
        Welch's U-test on ranked data.
        """
        if self._anova is None:
            self.anova()

        anova_significant = np.array(self._anova) < 1  # call array in case it is a Series

        # limit to significant data, convert to column-wise ranks.
        data = self.data[:, anova_significant]
        rank_data = np.apply_along_axis(rankdata, 0, data)
        # assignments = self.group_assignments[anova_significant]

        split_indices = np.where(np.diff(self.group_assignments))[0] + 1
        array_views = np.array_split(rank_data, split_indices, axis=0)

        # get mean and standard deviations of each
        fmean = partial(np.mean, axis=0)
        fvar = partial(np.var, axis=0)
        mu = np.vstack(list(map(fmean, array_views))).T  # transpose to get gene rows
        n = np.array(list(map(lambda x: x.shape[0], array_views)))
        s = np.vstack(list(map(fvar, array_views))).T
        s_norm = s / n  # transpose to get gene rows

        # calculate T
        numerator = mu[:, np.newaxis, :] - mu[:, :, np.newaxis]
        denominator = np.sqrt(s_norm[:, np.newaxis, :] + s_norm[:, :, np.newaxis])
        statistic = numerator / denominator

        # calculate df
        s_norm2 = s**2 / (n**2 * n-1)
        numerator = (s_norm[:, np.newaxis, :] + s_norm[:, :, np.newaxis]) ** 2
        denominator = (s_norm2[:, np.newaxis, :] + s_norm2[:, :, np.newaxis])
        df = np.floor(numerator / denominator)

        # get significance
        p = t.cdf(np.abs(statistic), df)  # note, two tailed test

        # calculate fdr correction; because above uses 2-tails, alpha here is halved
        # because each test is evaluated twice due to the symmetry of vectorization.
        p_adj = multipletests(np.ravel(p), alpha=self.alpha, method='fdr_tsbh')[1]
        p_adj = p_adj.reshape(*p.shape)

        phr = namedtuple('PostHocResults', ['p_adj', 'statistic', 'mu'])
        self.post_hoc = phr(p_adj, statistic, mu)

        if self.index is not None:
            p_adj = pd.Panel(
                p_adj, items=self.index[anova_significant], major_axis=self.groups,
                minor_axis=self.groups)
            statistic = pd.Panel(
                statistic, items=self.index[anova_significant], major_axis=self.groups,
                minor_axis=self.groups)
            mu = pd.DataFrame(mu, self.index[anova_significant], columns=self.groups)

        return p_adj, statistic, mu

    def population_markers(self, p_crit=0.0):
        """
        Return markers that are significantly differentially expressed in one
        population vs all others

        :param p_crit: float, fraction populations that may be indistinguishable from the
          highest expressing population for each gene. If zero, each marker gene is
          significantly higher expressed in one population relative to all others.
          If 0.1, 10% of populations may share high expression of a gene, and those
          populations will be marked as expressing that gene.

        """
        if self.post_hoc is None:
            self.post_hoc_tests()

        # get highest mean for each gene
        top_gene_idx = np.argmax(self.post_hoc.mu, axis=1)

        # index p_adj first dimension with each sample, will reduce to 2d genes x samples
        top_gene_sig = self.post_hoc.p_adj[:, top_gene_idx, :]

        # for each gene, count the number of non-significant DE results.
        sig = np.array(top_gene_sig < self.alpha)
        num_sig = np.sum(sig, axis=2)

        # if this is greater than N - 1 * p_crit, discard the gene.
        n = self.post_hoc.p_adj.shape[2] - 1  # number of genes, sub 1 for self
        idx_marker_genes = np.where(num_sig < n * (1 - p_crit))
        marker_genes = sig[idx_marker_genes, :]

        # correctly index these genes
        if self.index:
            pass  # todo fix this

        return marker_genes
