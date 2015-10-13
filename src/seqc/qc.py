__author__ = 'ambrose'

import numpy as np
import numpy.lib.recfunctions as rf
import pandas as pd
from collections import defaultdict, Counter
import pickle
from seqc import three_bit
from seqc.sam import ObfuscatedTuple
from seqc.convert_features import GeneTable
from scipy.special import gammaln, gammaincinv
from scipy.stats import chi2
from scipy.sparse import coo_matrix
from itertools import chain
import sys
from itertools import permutations
import subprocess


def deobfuscate(df):
    """exchange ObfuscatedTuple classes for regular tuples"""
    if (isinstance(df.iloc[0, :]['features'], tuple) and
            isinstance(df.iloc[0, :]['positions'], tuple)):
        return df
    if isinstance(df, np.ndarray):
        df = pd.DataFrame(df)
    features = df['features'].apply(lambda x: x.to_tuple())
    positions = df['positions'].apply(lambda x: x.to_tuple())
    df['positions'] = positions
    df['features'] = features
    return df


def mask_failing_cells(df_or_array, n_poly_t_required):
    return ((df_or_array['cell'] != 0) &
            (df_or_array['rmt'] != 0) &
            (df_or_array['n_poly_t'] >= n_poly_t_required) &
            (df_or_array['is_aligned'])
            )


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
        for i in iterable:
            self.union(*i)

    def find_all(self, vals):
        vals = [self.find_component(v) for v in vals]
        unique = set(vals)
        reindex = dict(zip(unique, range(len(unique))))
        set_membership = np.array([reindex[v] for v in vals])
        sets = np.array(list(reindex.values()))
        return set_membership, sets

    def find_component(self, iterable):
        """Return the set that obj belongs to

        If the iterable contains items that have been unioned, then any entry in the
         iterable will be sufficient to identify the set that obj belongs to. Use the
         first entry, and return the set associated with iterable.

        If the iterable has not been entered into the structure, this method can yield
         incorrect results
        """
        return self[next(iter(iterable))]

# determine some constants
_eps = 10e-10


def merge_observations_and_expectations(o, e):
    """may need to modify some stuff to deal with the weird feature object"""
    all_keys = set(o.keys()).union(e.keys())
    n = len(all_keys)
    obs = np.zeros(n, dtype=np.float)
    exp = np.zeros(n, dtype=np.float)
    for i, k in enumerate(all_keys):
        try:
            obs[i] = o[k]
        except KeyError:
            obs[i] = _eps
        try:
            exp[i] = e[k]
        except KeyError:
            exp[i] = _eps
    return obs, exp


def disambiguate(data, expectations, n_poly_t_required=0, alpha=0.1):
    """Resolve molecules that are ambiguously aligned. Input is an ndarray

    Need a way of grouping the ndarray such that I can assign the original positions
    their resolved positions + a column that indicates the type of resolution that
    was achieved.

    Need: groups of gene - (rmt / cell)

    results:
    1: unambiguous, 2: disambiguated by disjoint separation,
    3: partial disambiguation by model, 4: complete disambiguation by model
    5: ambiguous, 6: no model, likely contaminant

    args:
    -----

    returns:
    --------
    """
    # load the expectations
    if isinstance(expectations, str):
        with open(expectations, 'rb') as f:
            expectations = pickle.load(f)

    # create array to hold disambiguation results and a counter for total molecules
    indices = np.arange(data.shape[0])
    results = np.zeros(data.shape[0], dtype=np.int8)

    # mask = ((data['cell'] != 0) & (data['rmt'] != 0) &
    #         (data['n_poly_t'] >= n_poly_t_required) & data['is_aligned'])

    mask = mask_failing_cells(data, n_poly_t_required=n_poly_t_required)

    # filter, then merge cell & rmt
    seq_data = np.vstack([data['cell'][mask], data['rmt'][mask].astype(np.int64)]).T
    seq = np.apply_along_axis(three_bit.ThreeBit.ints2int, axis=1,
                              arr=seq_data)
    indices = indices[mask]

    # get indices of reads associated with each putative molecule (rmt/cell pair)
    molecules = defaultdict(list)
    for i, s in zip(indices, seq):
        molecules[s].append(i)
    for k, v in molecules.items():
        molecules[k] = np.array(v)  # covert all lists to arrays for indexing

    # set structured array dtype for counting reads per molecule occurences
    # obs_dtype = [('features', np.object), ('nobs', np.int)]
    # exp_dtype = [('features', np.object), ('pexp', np.float)]

    # loop over potential molecules
    for read_indices in molecules.values():

        # when indices is a single molecule, numpy does very weird things when converting
        # to arrays, so we need to keep track of this and structure array construction
        # accordingly

        # get a view of the structarray
        putative_molecule = data[read_indices]

        # check if all reads have a single, identical feature
        check_trivial = Counter(putative_molecule['features'])
        if len(check_trivial) == 1 and len(next(iter(check_trivial))) == 1:
            results[read_indices] = 1  # trivial
            continue

        # get disjoint set membership
        uf = UnionFind()
        uf.union_all(putative_molecule['features'])
        set_membership, sets = uf.find_all(putative_molecule['features'])

        # Loop Over Disjoint molecules
        for s in sets:

            disjoint_group = putative_molecule[set_membership == s]
            disjoint_group_idx = read_indices[set_membership == s]  # track index

            # get observed counts
            obs_counter = Counter(f.to_tuple() for f in putative_molecule['features'])
            # obs_features = np.array(list(obs_counter.keys()))
            # obs_counts = np.array(list(obs_counter.values()))

            # check that creating disjoint sets haven't made this trivial
            if len(obs_counter) == 1 and len(next(iter(obs_counter))) == 1:
                results[disjoint_group_idx] = 2
                continue  # no need to calculate probabilities for single model

            # convert observed counts to a structured array
            # obs = np.core.records.fromarrays([obs_features, obs_counts], dtype=obs_dtype)

            # get all potential molecule identities, create container for likelihoods, df
            possible_models = np.array(
                list(set(chain(*(f.to_tuple() for f in disjoint_group['features']))))
            )

            model_likelihoods = np.empty_like(possible_models, dtype=np.float)
            df = np.empty_like(possible_models, dtype=np.int)

            for i, m in enumerate(possible_models):

                # get model probabilities todo this should be pre-created in pickled index
                # todo this does the wrong thing with several types of input:
                # (1, 2): 1. --> array([1, 2]) instead of array([(1, 2)])
                # (1,) : 1. --> array([[1]]) (shape (1, 1) instead of shape (1,))

                exp_dict = expectations[m]
                obs, exp = merge_observations_and_expectations(obs_counter, exp_dict)

                # calculate model probability
                model_likelihoods[i] = multinomial_loglikelihood(obs, exp)
                df[i] = len(obs) - 1

            likelihood_ordering = np.argsort(model_likelihoods)[::-1]
            models = possible_models[likelihood_ordering]
            model_likelihoods = model_likelihoods[likelihood_ordering]

            # get top model
            passing_models, top_likelihood = [models[0]], model_likelihoods[0]

            # gauge relative likelihoods
            for i in range(1, len(model_likelihoods)):
                model = models[i]
                lh = model_likelihoods[i]
                degrees_freedom = df[i]
                p = likelihood_ratio_test(lh, top_likelihood, degrees_freedom)
                if p < alpha:
                    passing_models.append(model)
                else:
                    break  # no models with lower likelihoods will pass, terminate loop

            # adjust models, record results
            if len(passing_models) == 1:
                res = 4
            elif 1 < len(passing_models) < len(models):
                res = 3
            else:  # len(passing_models == len(models); no improvement was made.
                res = 5

            # set results
            results[disjoint_group_idx] = res

            # change features
            new_features = ObfuscatedTuple(passing_models)

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
            new_features = ObfuscatedTuple(passing_models)

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


def length_bias(arr, gtf):
    """calculate what, if any, length bias exists in the experiment found in arr"""
    raise NotImplementedError


def counts_matrix(arr, collapse_molecules, n_poly_t_required):

    # mask array for failed molecules
    arr = arr[mask_failing_cells(arr, n_poly_t_required)]

    # deobfuscate the array
    df = deobfuscate(arr)

    # group data by cell, rmt, feature to get molecules/counts
    counts = df.groupby(['cell', 'features', 'rmt']).size()

    # process data into d[feature][cell] : molecule count
    data = defaultdict(dict)
    if collapse_molecules:
        for (cell, features, _), _ in counts.iteritems():
                try:
                    data[features][cell] += 1
                except KeyError:
                    data[features][cell] = 1
    else:
        for (cell, features, _), count in counts.iteritems():
                try:
                    data[features][cell] += count
                except KeyError:
                    data[features][cell] = count

    # convert to values, row, col form for scipy.coo
    # pre-allocate arrays
    size = sum(len(c) for c in data.values())
    values = np.empty(size, dtype=int)
    row = np.empty(size, dtype=int)
    col = np.empty(size, dtype=int)
    i = 0
    for gene in data:
        for cell, count in data[gene].items():
            values[i] = count
            row[i] = cell
            col[i] = gene[0]  # todo this is super arbitrary, I'm selecting the first gene of the multi-alignment.
            i += 1

    # get max count to shrink dtype if possible
    maxcount = np.max(values)

    # set dtype
    if 0 < maxcount < 2 ** 8:
        dtype = np.uint8
    elif maxcount < 2 ** 16:
        dtype = np.uint16
    elif maxcount < 2 ** 32:
        dtype = np.uint32
    elif maxcount < 2 ** 64:
        dtype = np.uint64
    elif maxcount < 0:
        raise ValueError('Negative count value encountered. These values are not'
                         'defined and indicate a probable upstream bug')
    else:
        raise ValueError('Count values too large to fit in int64. This is very '
                         'unlikely, and will often cause Memory errors. Please check '
                         'input data.')

    # map row and cell to integer values for indexing
    unq_row = np.unique(row)  # these are the ids for the new rows / cols of the array
    unq_col = np.unique(col)
    row_map = dict(zip(unq_row, np.arange(unq_row.shape[0])))
    col_map = dict(zip(unq_col, np.arange(unq_col.shape[0])))
    row_ind = np.array([row_map[i] for i in row])
    col_ind = np.array([col_map[i] for i in col])

    # change dtype, set shape
    values = values.astype(dtype)
    shape = (unq_row.shape[0], unq_col.shape[0])

    # return a sparse array
    coo = coo_matrix((values, (row_ind, col_ind)), shape=shape, dtype=dtype)
    return coo, unq_row, unq_col


def set_filter_thresholds():
    """identify intelligent thresholds for each filter

    for each filter, determine a good threshold for eliminating a read based on its
    probability of contributing to a valid alignment

    which leads to the question: what is a valid alignment? Intuitively, it is one with
    (1) good sequence quality
    (2) good alignment score
    (3) minimal homopolymer trimming
    (4) associated with a cell that has enough reads
    (5) has a cell barcode and an RMT
    (6) is a part of an RMT that we are confident is real

    so I guess we could regress each of these predictors against a final gold standard
    if we could find one.

    What would that gold standard be? We could start with "is in a good cell" but most
    reads probably are in good cells, so we would want to correct for this.

    Perhaps we go back to the correlation plot? Could we add some additional columns?
    """
    raise NotImplementedError


def sam_to_count_multiple_files(sam_files, gtf_file, read_length):
    """count genes in each cell"""
    gt = GeneTable(gtf_file)
    all_genes = gt.all_genes()

    # map genes to ids
    n_genes = len(all_genes)
    gene_to_int_id = dict(zip(sorted(all_genes), range(n_genes)))
    cell_number = 1
    read_count = defaultdict(int)

    # add metadata fields to mimic htseq output; remember to remove these in the final
    # analysis
    gene_to_int_id['ambiguous'] = n_genes
    gene_to_int_id['no_feature'] = n_genes + 1
    gene_to_int_id['not_aligned'] = n_genes + 2

    for sam_file in sam_files:
    # pile up counts
        with open(sam_file) as f:
            for record in f:

                # discard headers
                if record.startswith('@'):
                    continue
                record = record.strip().split('\t')

                # get start, end, chrom, strand
                flag = int(record[1])
                if flag & 4:
                    int_gene_id = n_genes + 2  # not aligned
                else:
                    chromosome = record[2]
                    if flag & 16:
                        strand = '-'
                        end = int(record[3])
                        start = end - read_length
                    else:
                        strand = '+'
                        start = int(record[3])
                        end = start + read_length

                    genes = gt.coordinates_to_gene_ids(chromosome, start, end, strand)
                    if len(genes) == 1:
                        int_gene_id = gene_to_int_id[genes[0]]
                    if len(genes) == 0:
                        int_gene_id = n_genes + 1
                    if len(genes) > 1:
                        int_gene_id = n_genes
                read_count[(cell_number, int_gene_id)] += 1
        cell_number += 1

    # create sparse matrix
    cell_row, gene_col = zip(*read_count.keys())
    data = list(read_count.values())
    m = cell_number
    n = n_genes + 3

    coo = coo_matrix((data, (cell_row, gene_col)), shape=(m, n), dtype=np.int32)
    gene_index = sorted(all_genes) + ['ambiguous', 'no_feature', 'not_aligned']
    cell_index = ['no_cell'] + list(range(1, cell_number))

    return coo, gene_index, cell_index

    # out = subprocess.check_output(['which', 'htseq-count'])
    # if not out:
    #     raise RuntimeError('htseq-count not found in PATH')
    #
    # for sam_file in sam_files:
    #     htseq_cmd = ['htseq-count', sam_file, gtf_file]
    #
    #     p = subprocess.Popen(htseq_cmd, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
    #     out, err = p.communicate()
    #     if err:
    #         raise ChildProcessError('htseq error: %s' % err.decode())
    #
    #     # process the count file
    #
    # return fout


def sam_to_count_single_file(sam_file, gtf_file, read_length):
    """cannot separate file due to operating system limitations. Instead, implement
    a mimic of htseq-count that uses the default 'union' approach to counting, given the
    same gtf file"""

    # get conversion table, all possible genes for the count matrix
    gt = GeneTable(gtf_file)
    all_genes = gt.all_genes()

    # map genes to ids
    n_genes = len(all_genes)
    gene_to_int_id = dict(zip(sorted(all_genes), range(n_genes)))
    cell_to_int_id = {'no_cell': 0}
    cell_number = 1
    read_count = defaultdict(int)

    # add metadata fields to mimic htseq output; remember to remove these in the final
    # analysis
    gene_to_int_id['ambiguous'] = n_genes
    gene_to_int_id['no_feature'] = n_genes + 1
    gene_to_int_id['not_aligned'] = n_genes + 2

    # pile up counts
    with open(sam_file) as f:
        for record in f:

            # discard headers
            if record.startswith('@'):
                continue
            record = record.strip().split('\t')

            # get cell id, discard if no cell, add new cell if not found.
            cell = record[0].split(':')[0]
            if cell == 0:
                int_cell_id = 0
            else:
                try:
                    int_cell_id = cell_to_int_id[cell]
                except KeyError:
                    cell_to_int_id[cell] = cell_number
                    int_cell_id = cell_number
                    cell_number += 1

            # get start, end, chrom, strand
            flag = int(record[1])
            if flag & 4:
                int_gene_id = n_genes + 2  # not aligned
            else:
                chromosome = record[2]
                if flag & 16:
                    strand = '-'
                    end = int(record[3])
                    start = end - read_length
                else:
                    strand = '+'
                    start = int(record[3])
                    end = start + read_length

                genes = gt.coordinates_to_gene_ids(chromosome, start, end, strand)
                if len(genes) == 1:
                    int_gene_id = gene_to_int_id[genes[0]]
                if len(genes) == 0:
                    int_gene_id = n_genes + 1
                if len(genes) > 1:
                    int_gene_id = n_genes
            read_count[(int_cell_id, int_gene_id)] += 1

    # create sparse matrix
    cell_row, gene_col = zip(*read_count.keys())
    data = list(read_count.values())
    m = cell_number
    n = n_genes + 3

    coo = coo_matrix((data, (cell_row, gene_col)), shape=(m, n), dtype=np.int32)
    gene_index = sorted(all_genes) + ['ambiguous', 'no_feature', 'not_aligned']
    cell_index = ['no_cell'] + list(range(1, cell_number))

    return coo, gene_index, cell_index
