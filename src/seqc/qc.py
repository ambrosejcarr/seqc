__author__ = 'ambrose'

import numpy as np
import numpy.lib.recfunctions as rf
import pandas as pd
from collections import defaultdict
import pickle
from seqc import three_bit
from scipy.special import gammaln, gammaincinv
from scipy.stats import chi2
from itertools import chain
import sys
from itertools import permutations


def multinomial_loglikelihood(x, probs):
    """return the multinomial log-likelihood of probs given x log(L) = Mult(probs | x)"""
    return (gammaln(np.sum(x) + 1) - np.sum(gammaln(x + 1))) + np.sum(x * np.log(probs))


def likelihood_ratio_test(current_model, best_model, df):
    """likelihood ratio test, evaluated with a chi2 distribution with df=df

    note: assumes input parameters are log-likelihoods"""
    ratio = -2 * current_model + 2 * best_model
    return chi2.cdf(ratio, df)


class UnionFind:
    """Union-find data structure.

    Each unionFind instance X maintains a family of disjoint sets of
    hashable objects, supporting the following two methods:

    - X[item] returns a name for the set containing the given item.
      Each set is named by an arbitrarily-chosen one of its members; as
      long as the set remains unchanged it will keep the same name. If
      the item is not yet part of a set in X, a new singleton set is
      created for it.

    - X.union(item1, item2, ...) merges the sets containing each item
      into a single larger set.  If any item is not yet part of a set
      in X, it is added to X as one of the members of the merged set.
    """

    def __init__(self):
        """Create a new empty union-find structure."""
        self.weights = {}
        self.parents = {}

    def __getitem__(self, obj):
        """Find and return the name of the set containing the object."""

        # check for previously unknown object
        if obj not in self.parents:
            self.parents[obj] = obj
            self.weights[obj] = 1
            return obj

        # find path of objects leading to the root
        path = [obj]
        root = self.parents[obj]
        while root != path[-1]:
            path.append(root)
            root = self.parents[root]

        # compress the path and return
        for ancestor in path:
            self.parents[ancestor] = root
        return root

    def __iter__(self):
        """Iterate through all items ever found or unioned by this structure."""
        return iter(self.parents)

    def union(self, *objects):
        """Find the sets containing the objects and merge them all."""
        roots = [self[x] for x in objects]
        heaviest = max([(self.weights[r], r) for r in roots])[1]
        for r in roots:
            if r != heaviest:
                self.weights[heaviest] += self.weights[r]
                self.parents[r] = heaviest

    def union_all(self, iterable):
        self.union(i for i in iterable)

    def find_all(self, vals):
        vals = [self.find_component(v) for v in vals]
        unique = set(vals)
        reindex = dict(zip(unique, range(len(unique))))
        return np.array([reindex[v] for v in vals]), np.array(list(reindex.values()))

    def find_component(self, iterable):
        """Return the set that obj belongs to

        If the iterable contains items that have been unioned, then any entry in the
         iterable will be sufficient to identify the set that obj belongs to. Use the
         first entry, and return the set associated with iterable.

        If the iterable has not been entered into the structure, this method can yield
         incorrect results
        """
        return self[next(iter(iterable))]


def disambiguate(data, expectations, alpha=0.1):
    """Resolve molecules that are ambiguously aligned. Input is an ndarray

    Need a way of grouping the ndarray such that I can assign the original positions
    their resolved positions + a column that indicates the type of resolution that
    was achieved.

    Need: groups of gene - (rmt / cell)

    results:
    0: unambiguous, 1: disambiguated by disjoint separation,
    2: partial disambiguation by model, 3: complete disambiguation by model
    4: ambiguous, 5: no model, likely contaminant

    args:
    -----

    returns:
    --------
    """
    # load the expectations
    with open(expectations, 'rb') as f:
        expectations = pickle.load(f)

    # determine some constants
    eps = 10e-10

    # create array to hold disambiguation results and a counter for total molecules
    results = np.zeros(data.shape, dtype=np.int8)
    # total_molecules = 0

    # merge cell & rmt
    seq_data = np.vstack([data['cell'], data['rmt'].astype(np.int64)]).T
    seq = np.apply_along_axis(three_bit.ThreeBit.ints2int, axis=0,
                              arr=seq_data)

    # todo mask filtered reads

    # get indices of reads associated with each putative molecule (rmt/cell pair)
    molecules = defaultdict(list)
    for i, s in enumerate(seq):
        molecules[s].append(i)
    for k, v in molecules.items():
        molecules[k] = np.array(v)  # covert all lists to arrays for indexing

    # set structured array dtype for counting reads per molecule occurences
    obs_dtype = [('features', np.object), ('nobs', np.int)]
    exp_dtype = [('features', np.object), ('pexp', np.float)]

    # loop over potential molecules
    for read_indices in molecules.values():

        # get a view of the structarray
        group = data[read_indices]

        # check if all reads have a single, identical feature
        check_trivial = np.unique(group['features'])
        if check_trivial.shape[0] == 1 and len(check_trivial[0]) == 1:
            continue  # no changes necessary, result == 0 (unambiguous)

        # get disjoint set membership
        uf = UnionFind()
        uf.union_all(group['features'])
        set_membership, sets = uf.find_all(group['features'])

        # loop over disjoint molecules
        for s in sets:

            disjoint_group = group[set_membership == s]
            disjoint_group_idx = read_indices[set_membership == s]  # track index

            # get observed counts
            obs_features, obs_counts, = np.unique(group['features'], return_counts=True)

            # check that creating disjoint sets haven't made this trivial
            if obs_features.shape[0] == 1 and len(obs_features[0]) == 1:
                results[disjoint_group_idx] = 1
                continue  # no need to calculate probabilities for single model

            # convert observed counts to a structured array
            obs = np.core.records.fromarrays([obs_features, obs_counts], dtype=obs_dtype)

            # get all potential molecule identities, create container for likelihoods, df
            possible_models = np.array(list(set(chain(*disjoint_group['features']))))
            model_likelihoods = np.empty_like(possible_models, dtype=np.float)
            df = np.empty_like(possible_models, dtype=np.int)

            for i, m in enumerate(possible_models):

                # get model probabilities todo this should be pre-created in pickled index
                exp_features = np.array(list(expectations[m].keys()))
                exp_probs = np.array(list(expectations[m].values()))
                exp = np.core.records.fromarrays(
                    [exp_features, exp_probs], dtype=exp_dtype)

                # join on features
                ma = rf.join_by('features', obs, exp, jointype='outer')
                ma.set_fill_value = ['?', 0, eps]
                ma = ma.filled()

                # calculate model probability
                model_likelihoods[i] = multinomial_loglikelihood(ma['nobs'], ma['pexp'])
                df[i] = ma.shape[0]

            likelihood_ordering = np.argsort(model_likelihoods)[::-1]
            models = possible_models[likelihood_ordering]
            model_likelihoods = model_likelihoods[likelihood_ordering]

            # get top model
            passing_models, top_likelihood = [models[0]], model_likelihoods[0]

            # gauge relative likelihoods
            for i in range(1, len(model_likelihoods)):
                model = models[i]
                lh = model_likelihoods[i]
                df = df[i]
                p = likelihood_ratio_test(lh, top_likelihood, df)
                if p < alpha:
                    passing_models.append(model)
                else:
                    break  # no models with lower likelihoods will pass, terminate loop

            # adjust models, record results
            if len(passing_models) == 1:
                res = 3
            elif 1 < len(passing_models) < len(models):
                res = 2
            else:  # len(passing_models == len(models); no improvement was made.
                res = 4

            # set results
            results[disjoint_group_idx] = res

            # change features
            new_features = tuple(passing_models)

            data['features'][disjoint_group_idx] = new_features

    return results, data


def disambiguate_old(data, expectations, alpha=0.1):
    """Resolve molecules that are ambiguously aligned. Input is an ndarray

    Need a way of grouping the ndarray such that I can assign the original positions
    their resolved positions + a column that indicates the type of resolution that
    was achieved.

    Need: groups of gene - (rmt / cell)

    results:
    0: unambiguous, 1: disambiguated by disjoint separation,
    2: partial disambiguation by model, 3: complete disambiguation by model
    4: ambiguous, 5: no model, likely contaminant

    args:
    -----

    returns:
    --------
    """
    # load the expectations
    with open(expectations, 'rb') as f:
        expectations = pickle.load(f)

    # determine some constants
    eps = 10e-10

    # transform data type to DataFrame
    data = pd.DataFrame(data)
    data.index.name = 'idx'
    data.reset_index(inplace=True)  # get a column of indices for tracking results

    # array for tracking disambiguation results and total molecules
    results = np.zeros_like(data.index, dtype=np.int8)
    total_molecules = 0

    # groupby is only interested in reads passing filters; create a filter mask here

    # set aggregation functions
    aggregation_functions = {
        'features': np.size,
        'idx': lambda x: tuple(x),
    }

    # groupby cell / rmt
    # todo: this copies data, should exchange for a more efficient structure.
    ambiguous_groups = data.groupby(['cell', 'rmt'])

    for _, g in ambiguous_groups:

        assert(isinstance(g, pd.DataFrame))

        feature_counts = g.groupby('features').agg(aggregation_functions)
        feature_counts.columns = ['idx', 'counts']

        if feature_counts.shape[0] == 1:
            total_molecules += 1
            continue  # unambiguous; enter 0 in results, but start as zeros so no action.

        # todo: there is likely a way to streamline unionset usage here so that it
        # todo: naturally uses the np.array / dataframe data structure
        # separate into disjoint subsets
        disjoint = UnionFind()
        for feature in feature_counts.index:
            disjoint.union(*feature)

        # assign to sets
        observed_subsets = np.empty_like(feature_counts.index)
        for i, feature in enumerate(feature_counts.index):
            observed_subsets[i] = disjoint.find_component(feature)
        all_subsets = np.unique(observed_subsets)

        for s in all_subsets:
            observed = feature_counts.loc[observed_subsets == s]

            total_molecules += 1  # increment molecule counter.

            # check if still ambiguous
            if observed.shape[0] == 1:
                # todo might be a bug here; include a trivial case in the toy data
                for i in observed['idx']:
                    results[i] = 1  # disambiguated by disjoint separation
                continue

            # get the molecule models that need to be checked
            models = np.array(list(set(chain(*observed.index))))

            # set container for likelihoods, expectations and observations
            model_likelihoods = np.zeros_like(models, dtype=float)
            model_data = {}

            for i, m in enumerate(models):
                expected = expectations[m]
                eidx = pd.Index(list(expected.keys()), tupleize_cols=False)
                evals = list(expected.values())
                expected = pd.Series(evals, index=eidx)
                expected.name = 'expected'
                filled = pd.DataFrame(observed['counts']).join(expected)

                # fill any nan's
                filled[np.isnan(filled.values)] = eps
                model_data[m] = filled  # can remove if not comparing models

                lh = multinomial_loglikelihood(filled['counts'], filled['expected'])

                # check for malformed results
                if np.isnan(lh) or np.isinf(lh):
                    lh = -np.inf
                model_likelihoods[i] = lh

            # sort likelihoods by best model
            idx = np.argsort(model_likelihoods)[::-1]
            models = models[idx]
            model_likelihoods = model_likelihoods[idx]

            # get top model
            passing_models, top_likelihood = [models[0]], model_likelihoods[0]

            # gauge relative likelihoods
            for i in range(1, len(model_likelihoods)):
                model = models[i]
                lh = model_likelihoods[i]
                df = model_data[model].shape[0]
                p = likelihood_ratio_test(lh, top_likelihood, df)
                if p < alpha:
                    passing_models.append(model)
                else:
                    break  # no other models will pass either, have lower likelihoods

            # store results of disambiguation
            # need to track both the result of disambiguation and the new model (features)
            if len(passing_models) == 1:
                res = 3
            elif len(passing_models) >= 1:
                if len(passing_models) < len(models):
                    res = 2
                else:
                    res = 4

            # get all of the indices
            indices = list(chain(*observed['idx']))

            # set results
            results[indices] = res

            # change features
            new_features = tuple(passing_models)

            data.ix[indices, 'features'] = new_features

    return results, data


high_value = sys.maxsize


def hamming_dist_bin(c1, c2):
    if three_bit.ThreeBit.seq_len(c1) != three_bit.ThreeBit.seq_len(c2):
        return high_value
    d = 0
    while c1 > 0:
        if c1 & 0b111 != c2 & 0b111:
            d += 1
        c1 >>= 3
        c2 >>= 3
    return d


def generate_close_seq(seq):
    """Return a list of all sequences that are up to 2 hamm distance from seq"""
    res = []
    l = three_bit.ThreeBit.seq_len(seq)

    # genereta all sequences that are dist 1
    for i in range(l):
        mask = 0b111 << (i * 3)
        cur_chr = (seq & mask) >> (i * 3)
        res += [seq & (~mask) | (new_chr << (i * 3)) for new_chr in
                three_bit.ThreeBit.bin_nums if new_chr != cur_chr]

    # genereta all sequences that are dist 2
    for i in range(l):
        mask_i = 0b111 << (i * 3)
        chr_i = (seq & mask_i) >> (i * 3)
        for j in range(i + 1, l):
            mask_j = 0b111 << (j * 3)
            chr_j = (seq & mask_j) >> (j * 3)
            mask = mask_i | mask_j
            res += [seq & (~mask) | (new_chr_i << (i * 3)) | (new_chr_j << (j * 3)) for
                    new_chr_i in three_bit.ThreeBit.bin_nums if new_chr_i != chr_i for
                    new_chr_j in three_bit.ThreeBit.bin_nums if new_chr_j != chr_j]

    return res


def prob_d_to_r_bin(d_seq, r_seq, err_rate):
    """
    Return the probability of d_seq turning into r_seq based on the err_rate table
    (all binary)
    """

    if three_bit.ThreeBit.seq_len(d_seq) != three_bit.ThreeBit.seq_len(r_seq):
        return 1

    p = 1.0
    while d_seq > 0:
        if d_seq & 0b111 != r_seq & 0b111:
            p *= err_rate[three_bit.ThreeBit.ints2int([d_seq & 0b111, r_seq & 0b111])]
        d_seq >>= 3
        r_seq >>= 3
    return p


def estimate_error_rate(cell_barcodes, cell_barcode_data):

    """
    Estimate the error rate based on the barcodes in the data and the correct barcodes in
    the barcode file. Return an error_rate table.
    """

    if isinstance(cell_barcodes, str):
        with open(cell_barcodes, 'rb') as f:
            cell_barcodes = pickle.load(f)

    # create table of error types to hold number of occurrences
    errors = list(
        three_bit.ThreeBit.ints2int([p[0], p[1]]) for p in
        permutations(three_bit.ThreeBit.bin_nums, r=2))
    error_table = dict(zip(errors, [0] * len(errors)))

    correct_instances = 0
    for code in cell_barcode_data:
        # try:
        errors = cell_barcodes.map_errors(code)
        # except TypeError:
        #     print(code, type(code))
        #     raise
        if errors:
            for err_type in errors:
                error_table[err_type] += 1
            # some of the barcodes were correct
            correct_instances += three_bit.ThreeBit.seq_len(code) - len(errors)
        else:
            correct_instances += three_bit.ThreeBit.seq_len(code)

    # convert to error rates
    default_error_rate = 0.02
    err_rate = dict(zip(errors, [0.0] * len(errors)))
    if sum(error_table.values()) == 0:
        print('No errors were detected, using %f uniform error chance.' % (
            default_error_rate))
        err_rate = dict(zip(errors, [default_error_rate] * len(errors)))
    for k, v in error_table.items():
        try:
            err_rate[k] = v / (sum(n for err_type, n in error_table.items() if
                                   err_type & 0b111000 == k & 0b111000) +
                               correct_instances)
        except ZeroDivisionError:
            print('Warning: too few reads to estimate error rate for %r '
                  'setting default rate of %f' % (k, default_error_rate))
            err_rate[k] = default_error_rate
    return err_rate


# def list_errors(code, correct_barcodes):
#     """
#     For a given code and a list of correct barcodes - find the correct barcode
#     that is closest to code and return the list of errors that turned it into
#     code. An error is a six bit int representing a two chr string of type "AG","CT", etc.
#     """
#     # find the closest correct barcode
#     min_dist = high_value
#     donor = 0
#     for cor_code in correct_barcodes:
#         hamm_d = hamming_dist_bin(code, cor_code)
#
#         if hamm_d < min_dist:
#             min_dist = hamm_d
#             donor = cor_code
#             if hamm_d == 1:
#                 break
#
#     if donor == 0:
#         print('Error: no donor code was found to be closest. code = ',
#               three_bit.ThreeBit.bin2str(code))
#     # return the actual error
#     err_list = []
#     while code > 0:
#         if code & 0b111 != donor & 0b111:
#             err_list.append(three_bit.ThreeBit.ints2int([code & 0b111, donor & 0b111]))
#         code >>= 3
#
#     return err_list


def correct_errors(data, cell_barcodes, donor_cutoff=1, p_val=0.1):
    """calculate and correct errors in barcode sets"""

    # todo mask (1) invalid barcodes and (2) reads failing filters

    is_error = np.zeros(data.shape[0], dtype=np.bool)
    err_rate = estimate_error_rate(cell_barcodes, data['cell'])

    N = three_bit.ThreeBit._str2bindict['N']

    # merge sequences
    seq_data = np.vstack([data['cell'], data['rmt'].astype(np.int64)]).T
    seq = np.apply_along_axis(three_bit.ThreeBit.ints2int, axis=1,
                              arr=seq_data)
    # assert(seq.shape[0] == data.shape[0])

    # group by features
    genes = defaultdict(list)
    for i, f in enumerate(data['features']):
        genes[f].append(i)
    for k, v in genes.items():
        genes[k] = np.array(v)  # covert all lists to arrays for indexing todo might be unnecessary

    for gene_indices in genes.values():

        # mask gene_indices for reads failing filters to get potential recipient sequences
        gene_data = data[gene_indices]
        gene_seqs = seq[gene_indices]

        # get donor sequences (any with valid cell barcodes)
        d_mask = gene_data['valid_cell'] == 1
        # d_indices = gene_indices[d_mask]  # todo I commented a bunch of stuff that I don't think is necessary; np.unique shouls suffice for donors.

        # group donor reads by indices
        donor_counts = dict(zip(*np.unique(gene_seqs[d_mask], return_counts=True)))
        # donor_counts = defaultdict(list)
        # for i in d_indices:
        #     donor_counts[seq[i]].append(i)

        # get recipient sequences passing all filters todo include all filters
        r_mask = (gene_data['valid_cell'] == 1) & (gene_data['alignment_score'] > 30)
        r_indices = gene_indices[r_mask]

        recipient_data = defaultdict(list)
        for i in r_indices:
            recipient_data[seq[i]].append(i)

        for r_seq, r_indices in recipient_data.items():
            r_num_occurences = len(r_indices)

            """
            We want to calculate the error rate (L) required such that we expect
            to see k observations of r_seq or more 5% of the time.
            This is:
            Pois(K >= k) = 0.1; P(r_seq == error)
            = Pois(K <= k) = 0.9
            Pois(K <= k) is the poisson cdf function, which is numerically
            approximated by the regularized incomplete gamma function:
            Pois(K <= k | L) = gammainc(k + 1, L), x >= 0
                                                      = 0, otherwise
            thus L = gammaincinv(k + 1, 0.9)

            Rami: after checking it should be L = gammaincinv(k + 1, 0.1)
            summary - threshold is the mean of the dist for which the chance of getting
            r_num_occurences is P_VALUE if the expected errors is bigger than that,
            than the chance of getting r_num_occurences by chance of error is higher
            than P_VALUE a and so we need to regect r_seq as an error.
            lower P_VALUE will increase the number of errors (we require higher level of
            certainty)
            """

            threshold = gammaincinv(r_num_occurences, p_val)

            expected_errors = 0
            for d_seq in generate_close_seq(r_seq):
                try:
                    d_num_occurences = donor_counts[d_seq]
                except KeyError:
                    continue
                if d_num_occurences <= donor_cutoff:
                    continue  # don't bother calculating probability if counts < cutoff

                # get probability of conversion
                p_dtr = prob_d_to_r_bin(d_seq, r_seq, err_rate)

                expected_errors += d_num_occurences * p_dtr
                if expected_errors > threshold:
                    is_error[r_indices] = 1
                    break

    return is_error, err_rate
