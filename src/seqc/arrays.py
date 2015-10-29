__author__ = 'Ambrose J. Carr'

from operator import itemgetter
import numpy as np
from copy import copy
from numpy.lib.recfunctions import append_fields
from itertools import islice, chain
from seqc.three_bit import ThreeBit
from seqc.qc import UnionFind, multinomial_loglikelihood, likelihood_ratio_test
from seqc.sam import average_quality
from collections import deque, defaultdict
from seqc import h5
from scipy.sparse import coo_matrix
import tables as tb
import pickle
import os


# todo | it would save more memory to skip "empty" features and give them equal indices
# todo | to sparsify the structure. e.g. empty index 5 would index into data[7:7],
# todo | returning array([])
class JaggedArray:

    def __init__(self, data, index):
        diff = index[:, 0] - index[:, 1]
        if np.any(diff == 0):
            raise ValueError('Invalid index: each index pair must have a non-zero size')
        self._data = data
        self._index = index

    @property
    def data(self):
        return self._data

    @property
    def index(self):
        return self._index

    def __len__(self):
        return self._index.shape[0]

    def __getitem__(self, item):
        """
        returns a selection of the array. There are three valid types with different
        return types:

        (1) int: returns a np.array object containing the values corresponding to the
        index
        (3) np.array: assumes the array contains multiple indices. Returns an array of
        np.array objects corresponding to each index
        (2) slice: returns a JaggedArray object containing the subsection of indices
        indicated by the slice
        """
        if isinstance(item, slice):
            min_index = self.index[slice.start, 0]
            max_index = self.index[slice.stop, 1]
            new_index = self.index[slice, :]
            adj_index = new_index - min_index  # zero the index
            new_data = self.data[min_index:max_index]
            return self.__init__(new_data, adj_index)
        elif isinstance(item, np.ndarray):
            indices = self._index[item, :]
            return np.array([self.data[i[0]:i[1]].view(HashableArray) for i in indices])
        else:
            i = self._index[item, :]
            try:
                return self._data[i[0]:i[1]].view(HashableArray)
            except IndexError:
                raise IndexError('invalid slice: index %d: %s for index of length %d' %
                                 (item, repr(i), self.__len__()))

    def _setitem(self, i, value):
        old = self[i]
        len_old = len(old)
        if len(value) > len_old:
            raise ValueError
        elif len(value) == len_old:
            idx = self.index[i]
            self._data[idx[0]:idx[1]] = value
        else:
            idx_start = self.index[i, 0]
            idx_end = idx_start + len(value)
            self._data[idx_start:idx_end] = value
            self._index[i, :] = (idx_start, idx_end)

    def __setitem__(self, index, value):
        """
        allows a feature or set of features to be re-set by adjusting the index and
        values. Note that increasing the number of features is not supported, and will
        raise a ValueError. Features can only be changed or decreased in number. Note
        that setting features in this way leaves unused space in the array.
        """
        if isinstance(index, slice):
            for i in self._index[index.start:index.stop]:
                self._setitem(i, value)
        elif isinstance(index, np.ndarray):
            for i in index:
                self._setitem(i, value)
        else:
            self._setitem(index, value)

    def __repr__(self):
        if len(self.index) < 10:
            return ('<JaggedArray:\n[' + '\n '.join(str(self[i])
                    for i in range(len(self))) + ']>')
        else:
            return ('<JaggedArray:\n[' + '\n '.join(str(self[i]) for i in range(5)) +
                    '\n ...\n ' +
                    '\n '.join(str(self[i]) for i in range(len(self) - 5, len(self))) +
                    ']>'
                    )

    def extend(self, other):
        new_data = np.concatenate([self.data, other.data])
        l = len(self.index)
        other_index = other.index + l
        new_index = np.vstack([self.index, other_index])
        return JaggedArray(new_data, new_index)

    @staticmethod
    def size_to_dtype(n):
        if n > 2 * (2 ** 64):
            raise ValueError('Too many input reads for int64, n must be <= 2^64')
        elif n > 2 * (2 ** 32):
            return np.uint64
        elif n > 2 * (2 ** 16):
            return np.uint32
        elif n > 2 * (2 ** 8):
            return np.uint16
        else:
            return np.uint8

    @classmethod
    def from_iterable(cls, nested_iterable, dtype=np.int32):
        if not nested_iterable:
            raise ValueError('Cannot initialize empty array')
        k = len(nested_iterable)
        n = sum(len(i) for i in nested_iterable)
        data = np.zeros(n, dtype=dtype)
        index = np.zeros((k, 2), dtype=cls.size_to_dtype(n))
        i = 0
        for j, val in enumerate(nested_iterable):
            size = len(val)
            index[j, :] = (i, i + size)
            data[i:i + size] = val
            i += size
        return cls(data, index)


class ArraySet:

    def __init__(self, iterable_of_arrays):
        """identifies a set of unique arrays"""
        self._seen = set()
        self._arrays = []
        for array in iterable_of_arrays:
            hashed = hash(array.tobytes())
            if not hashed in self._seen:
                self._seen.add(hashed)
                self._arrays.append(array)

    def __iter__(self):
        return iter(self._arrays)

    def __len__(self):
        return len(self._arrays)

    def __repr__(self):
        if len(self) > 10:
            return '<ArraySet {%s ...}>' % ', '.join(repr(a) for a in self._arrays)
        else:
            return '<ArraySet {%s}>' % ', '.join(repr(a) for a in self._arrays)
            # repr(np.array(self._arrays, dtype=object))

    def __contains__(self, array):
        return True if hash(array.tobytes()) in self._seen else False

    def add(self, array):
        try:
            hashed = hash(array.tobytes())
        except AttributeError:
            raise ValueError('array must be a np.ndarray object, not %s' % type(array))
        if array.tobytes() not in self._seen:
            self._arrays.append(array)
            self._seen.add(hashed)

    def union(self, other):
        union_ = copy(self)
        for v in other:
            union_.add(v)
        return union_

    def tolist(self):
        return self._arrays


class HashableArray(np.ndarray):
    """custom class that simply adds a __hash__ method to numpy arrays"""

    def __new__(cls, input_array_or_iterable):
        """returns a view of the input_array with hash capability"""
        try:
            return input_array_or_iterable.view(cls)
        except AttributeError:  # input_array_or_iterable was not an ndarray
            return np.array(input_array_or_iterable).view(cls)

    def __hash__(self):
        """this is fast for small arrays, but slow for large arrays. It is guaranteed
        to be unique provided that all hashed values are arrays. It may collide when
        hashed with mixed types"""
        return hash(self.tobytes())


class ArrayCounter:

    def __init__(self, iterable_of_arrays):
        self.counts = {}
        self._arrays = []
        self._keys = set()
        for a in iterable_of_arrays:
            hashed = hash(a.tobytes())
            if hashed not in self._keys:
                self._keys.add(hashed)
                self._arrays.append(a)
                self.counts[hashed] = 1
            else:
                self.counts[hashed] += 1

    def __getitem__(self, item):
        try:
            return self.counts[hash(item.tobytes())]
        except AttributeError:
            if not isinstance(item, np.ndarray):
                raise ValueError('item must be a np.array object')
            else:
                raise

    def __iter__(self):
        return iter(self._arrays)

    def __repr__(self):
        sorted_data = sorted(self.items(), key=itemgetter(1), reverse=True)
        top10 = ', '.join('%s: %d' % (repr(i[0]), i[1]) for i in sorted_data[:10])
        if len(sorted_data) > 10:
            return '<ArrayCounter {%s ...}>' % top10
        else:
            return '<ArrayCounter {%s}>' % top10

    def __len__(self):
        return len(self._arrays)

    def keys(self):
        for k in self._arrays:
            yield k

    def values(self):
        for k in self._keys:
            yield self[k]

    def items(self):
        for k in self._arrays:
            yield k, self[k]


class ReadArray:

    _dtype = [
        ('cell', np.int64),
        ('rmt', np.int32),
        ('n_poly_t', np.uint8),
        ('valid_cell', np.bool),
        ('trimmed_bases', np.uint8),
        ('rev_quality', np.uint8),
        ('fwd_quality', np.uint8),
        ('is_aligned', np.bool),
        ('alignment_score', np.uint8)]

    def __init__(self, data, features, positions):
        """
        n is the total number of reads. The array will contain m multialignments, where
        m == n when all reads are uniquely aligned, and m < n otherwise.
        """

        # check that arrays are aligned
        if not data.shape[0] == len(features) == len(positions):
            raise ValueError('Input arrays are not aligned. data (%d), features (%d), '
                             'and positions (%d) must all have the same number of rows' %
                             (data.shape[0], len(features), len(positions)))

        self._data = data
        self._features = features
        self._positions = positions

    def __getitem__(self, item):
        return self._data[item], self._features[item], self._positions[item]

    def __len__(self):
        return self._data.shape[0]

    def __repr__(self):
        return '<ReadArray object with %d items>' % len(self)

    @property
    def data(self):
        return self._data

    @property
    def features(self):
        return self._features

    @property
    def positions(self):
        return self._positions

    @property
    def nbytes(self):
        return (self._data.nbytes + self.features._index.nbytes +
                self.features._data.nbytes + self.positions._data.nbytes +
                self.positions._index.nbytes)

    @classmethod
    # todo this is losing either the first or the last read
    def from_samfile(cls, samfile, feature_table, feature_positions):

        # first pass, get size of file
        # second pass, group up the multi-alignments


        # first pass through file: get the total number of alignments and the number of
        # multialignments for array construction; store indices for faster second pass
        with open(samfile, 'r') as f:

            # first, get header size and get rid of the headers
            line = f.readline()
            header_size = 0
            while line.startswith('@'):
                header_size += 1
                line = f.readline()

            # next, count alignments and multi-alignments, keep track of where the MAs
            # switch
            n_reads = 0
            n_mult = 0
            ma_sizes = []

            current_name = line.split('\t')[0]
            current_size = 1
            for line in f:
                n_reads += 1
                new_name = line.split('\t')[0]
                if current_name != new_name:
                    n_mult += 1
                    current_name = new_name
                    ma_sizes.append(current_size)
                    current_size = 1
                else:
                    current_size += 1

        # construct the arrays in contiguous memory
        read_data = np.zeros((n_mult,), dtype=cls._dtype)

        # both features and positions can reference the same index object
        index_dtype = JaggedArray.size_to_dtype(n_reads)
        index = np.zeros((n_mult, 2), dtype=index_dtype)

        # not all reads will have valid features; this size estimate is an upper bound
        # that will shrink when the array is filled
        features = np.zeros(n_reads, dtype=np.uint32)
        positions = np.zeros(n_reads, dtype=np.uint32)

        # 2nd pass through the file: populate the arrays
        current_index_position = 0
        with open(samfile, 'r') as f:
            records = iter(f)
            _ = list(islice(records, header_size))  # get rid of header
            for i, s in enumerate(ma_sizes):
                # get all records associated with a multialignment
                multi_alignment = list(islice(records, s))
                if not multi_alignment:  # iterator is exhausted
                    break
                recd, feat, pos = cls._process_multialignment(
                    multi_alignment, feature_table, feature_positions)
                read_data[i] = recd
                features[current_index_position:current_index_position + len(feat)] = feat
                positions[current_index_position:current_index_position + len(pos)] = pos
                index[i, :] = (current_index_position, current_index_position + s)
                current_index_position += s

        # eliminate any extra size from the features & positions arrays
        features = features[:current_index_position]
        positions = positions[:current_index_position]

        jagged_features = JaggedArray(features, index)
        jagged_positions = JaggedArray(positions, index)

        return cls(read_data, jagged_features, jagged_positions)

    @classmethod
    def _process_multialignment(cls, alignments, feature_table, feature_positions):
        """translate a sam multi-alignment into a ReadArray row"""

        # all fields are identical except feature; get from first alignment
        first = alignments[0].strip().split('\t')

        rev_quality = average_quality(first[10])
        alignment_score = int(first[13].split(':')[-1])

        # parse data from name field, previously extracted from forward read
        forward_metadata = (int(f) for f in first[0].strip().split(';')[0].split(':'))
        cell, rmt, n_poly_t, valid_cell, trimmed_bases, fwd_quality = forward_metadata

        # get all features and positions
        if len(alignments) == 1:

            # check if read is unaligned
            flag = int(first[1])
            if flag & 4:
                is_aligned = False
                rec = (cell, rmt, n_poly_t, valid_cell, trimmed_bases, rev_quality,
                       fwd_quality, is_aligned, alignment_score)
                features, positions = [0], [0]
                return rec, features, positions
            else:
                true_position = int(first[3])
                strand = '-' if (flag & 16) else '+'
                reference_name = first[2]
                features = [cls.translate_feature(reference_name, strand, true_position,
                                                  feature_table, feature_positions)]
                positions = [true_position]
        else:
            features = []
            positions = []
            for alignment in alignments:
                alignment = alignment.strip().split('\t')
                true_position = int(alignment[3])
                flag = int(alignment[1])
                strand = '-' if (flag & 16) else '+'
                reference_name = alignment[2]
                features.append(cls.translate_feature(
                    reference_name, strand, true_position, feature_table,
                    feature_positions))

                positions.append(true_position)

        delete = deque()
        for i in range(len(features)):
            if features[i] == 0:
                delete.append(i)

        features, positions = cls.multi_delete(delete, features, positions)

        is_aligned = True if features else False

        if not features:
            features, positions = [0], [0]

        rec = (cell, rmt, n_poly_t, valid_cell, trimmed_bases, rev_quality,
               fwd_quality, is_aligned, alignment_score)

        return rec, features, positions

    def stack(self, other):
        """returns a stacked copy of self and other"""
        data = np.hstack(self.data, other.data)
        features = self.features.extend(other.features)
        positions = self.positions.extend(other.positions)
        return ReadArray(data, features, positions)

    def group_for_disambiguation(self, required_poly_t):
        """
        generate an heirarchical dictionary that groups the indices of reads according to
        their gene and sequence information

        any indices associated with a failing (masked) read is not included in the
        index dictionary, and so is implicity given a zero value in the results vector
        """
        indices = np.arange(self.data.shape[0])

        mask = self.mask_failing_cells(n_poly_t_required=required_poly_t)

        # filter, then merge cell & rmt; this creates a copy of cell + rmt
        # memory overhead = n * 16b
        seq_data = np.vstack([self.data['cell'][mask],
                              self.data['rmt'][mask].astype(np.int64)]).T
        seq = np.apply_along_axis(ThreeBit.ints2int, axis=1,
                                  arr=seq_data)
        indices = indices[mask]

        # get indices of reads associated with each putative molecule (rmt/cell pair)
        molecules = defaultdict(list)
        for i, s in zip(indices, seq):
            molecules[s].append(i)
        for k, v in molecules.items():
            molecules[k] = np.array(v)  # covert all lists to arrays for indexing

        return molecules, seq

    def group_for_error_correction(self, required_poly_t=0):
        """
        returns a two-level dictionary that maps features to their sequences and sequences
        to all of the positions in self that feature-sequence positions are found.

        The dictionary is created from post-filtered sequences. required_poly_t is one
        such filter. Setting it equal to zero (default value) will result in a
        pass-through filter

        In detail, the returned dictionary has the structure:
        molecules (dict)
            |
            v
            hashed byte-representation of feature np.ndarray
                |
                v
                concatenated 3-bit representation of cb1-cb2-rmt, where the least
                significant bits of the sequence contain the rmt
                    |
                    v
                    np.ndarray containing the indices where the feature-seq combinations
                    are found in the original array

        >> group = molecules[features][seq]
        >> group
        [1] np.array([1, 5, 11, 224, 123412])
        >> group.shape[0]
        [2] 5  # number of counts associated with seq, given feature.
        """
        indices = np.arange(self.data.shape[0])

        mask = self.mask_failing_cells(n_poly_t_required=required_poly_t)

        # filter and merge cell/rmt
        seq_data = np.vstack([self.data['cell'][mask],
                              self.data['rmt'][mask].astype(np.int64)]).T
        # todo may not need to build this
        seq = np.apply_along_axis(ThreeBit.ints2int, axis=1,
                                  arr=seq_data)

        indices = indices[mask]

        # get indices of reads associated with each putative molecule (gene/rmt/cell)
        # cannot hash numpy arrays; instead hash the bytes representation of the array
        # this means that the feature (f) is giberish, but that's not important for
        # the experiment.

        molecules = defaultdict(dict)
        for i, s, f in zip(indices, seq, self.features):
            try:
                molecules[hash(f.tobytes())][s].append(i)
            except AttributeError:  # f is not an array
                molecules[hash(np.array(f).tobytes())][s].append(i)
            except KeyError:  # s has not yet been seen, cannot append.
                try:
                    molecules[hash(f.tobytes())][s] = [i]
                except AttributeError:
                    molecules[hash(np.array(f).tobytes())][s].append(i)
        for f in molecules.keys():
            for s, v in molecules[f].items():
                molecules[f][s] = np.array(v)

        return dict(molecules)

    def resolve_alignments(self, expectations, required_poly_t=0, alpha=0.1):
        molecules, seqs = self.group_for_disambiguation(required_poly_t)
        results = np.zeros(self.data.shape[0], dtype=np.int8)

        # load the expectations
        if os.path.isfile(expectations):
            with open(expectations, 'rb') as f:
                expectations = pickle.load(f)

        # loop over potential molecules
        for read_indices in molecules.values():

            # get a view of the structarray
            putative_molecule = self.data[read_indices]
            putative_features = self.features[read_indices]

            # check if all reads have a single, identical feature
            # can't hash numpy arrays, but can hash .tobytes() of the array
            check_trivial = ArrayCounter(putative_features)
            if len(check_trivial) == 1 and len(next(iter(check_trivial))) == 1:
                results[read_indices] = 1  # trivial
                continue

            # get disjoint set membership
            uf = UnionFind()
            uf.union_all(putative_features)
            set_membership, sets = uf.find_all(putative_features)

            # Loop Over Disjoint molecules
            for s in sets:

                disjoint_group = putative_molecule[set_membership == s]
                disjoint_features = putative_features[set_membership == s]  # todo doesn't like this because it's a tuple!
                disjoint_group_idx = read_indices[set_membership == s]  # track index

                # get observed counts
                obs_counter = ArrayCounter(putative_features)

                # check that creating disjoint sets haven't made this trivial
                if len(obs_counter) == 1 and len(next(iter(obs_counter))) == 1:
                    results[disjoint_group_idx] = 2
                    continue  # no need to calculate probabilities for single model

                # get all potential molecule identities, create container for lh, df
                possible_models = np.array(
                    list(set(chain(*disjoint_features)))
                )

                model_likelihoods = np.empty_like(possible_models, dtype=np.float)
                df = np.empty_like(possible_models, dtype=np.int)

                for i, m in enumerate(possible_models):

                    # get model probabilities todo this should be pre-created in pickled index
                    # todo this does the wrong thing with several types of input:
                    # (1, 2): 1. --> array([1, 2]) instead of array([(1, 2)])
                    # (1,) : 1. --> array([[1]]) (shape (1, 1) instead of shape (1,))

                    try:
                        exp_dict = expectations[m]
                    except KeyError:
                        continue
                    except IndexError:
                        print(m)
                        print(type(m))
                        raise
                    obs, exp = outer_join(obs_counter, exp_dict)

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
                self.features[disjoint_group_idx] = passing_models

        self._data = append_fields(self.data, 'disambiguation_results', results)

    @staticmethod
    def multi_delete(sorted_deque, *lists):
        while sorted_deque:
            i = sorted_deque.pop()
            for l in lists:
                del l[i]
        return lists

    def mask_failing_cells(self, n_poly_t_required):
        return ((self.data['cell'] != 0) &
                (self.data['rmt'] != 0) &
                (self.data['n_poly_t'] >= n_poly_t_required) &
                (self.data['is_aligned'])
                )

    @staticmethod
    def translate_feature(reference_name, strand, true_position, feature_table,
                          feature_positions):
        rounded_position = true_position // 100 * 100
        try:
            potential_scids = feature_table[(strand, reference_name, rounded_position)]
        except KeyError:
            return 0

        # only one feature in stranded libraries
        for scid in potential_scids:
            if any(s < true_position < e for (s, e) in feature_positions[scid]):
                return scid  # todo | test if ever > 1 feature
            else:
                return 0

    def save_h5(self, archive_name):
        """
        efficiently save the ReadArray to a compressed h5 archive; note that this
        will overwrite an existing archive with the same name
        """

        def store_carray(h5f, array, where, name):
            atom = tb.Atom.from_dtype(array.dtype)
            store = h5f.createCArray(where, name, atom, array.shape)
            store[:] = array

        blosc5 = tb.Filters(complevel=5, complib='blosc')
        f = tb.open_file(archive_name, mode='w', title='Data for seqc.ReadArray',
                         filters=blosc5)

        # store data
        f.create_table(f.root, 'data', self._data)
        # except tb.exceptions.NodeError:
        #     f.remove_node(f.root, 'data')
        #     f.create_table(f.root, 'data', self._data)

        # create group for feature-related data and store.
        feature_group = h5.create_group(
            f, 'features', '/', 'Data for ReadArray._features')
        store_carray(f, self._features._data, feature_group, 'data')
        store_carray(f, self._features._index, feature_group, 'index')

        # store position data
        positions_group = h5.create_group(
            f, 'positions', '/', 'Data for ReadArray._features')
        store_carray(f, self._positions._data, positions_group, 'data')
        store_carray(f, self._positions._index, positions_group, 'index')

        f.close()

    @classmethod
    def from_h5(cls, archive_name):
        f = tb.open_file(archive_name, mode='r')

        data = f.root.data.read()

        fdata = f.root.features.data.read()
        findex = f.root.features.index.read()

        pdata = f.root.positions.data.read()
        pindex = f.root.positions.index.read()

        features = JaggedArray(fdata, findex)
        positions = JaggedArray(pdata, pindex)

        return cls(data, features, positions)

    # todo write tests
    def to_sparse_counts(self, collapse_molecules, n_poly_t_required):

        mask = self.mask_failing_cells(n_poly_t_required)
        unmasked_inds = np.arange(self.data.shape[0])[mask]
        molecule_counts = defaultdict(dict)

        if collapse_molecules:
            # get molecule counts

            for i in unmasked_inds:
                feature = self._features[i]
                cell = self.data['cell'][i]
                rmt = self.data['rmt'][i]
                try:
                    molecule_counts[feature][cell].add(rmt)
                except KeyError:
                    molecule_counts[feature][cell] = {rmt}

            # convert to molecule counts
            for f in molecule_counts.keys():
                for c, rmts in molecule_counts[f].items():
                    molecule_counts[f][c] = len(rmts)
        else:
            for i in unmasked_inds:
                feature = self._features[i]
                cell = self.data['cell'][i]
                try:
                    molecule_counts[feature][cell] += 1
                except KeyError:
                    molecule_counts[feature][cell] = 1

        # convert to values, row, col form for scipy.coo
        # pre-allocate arrays
        size = sum(len(c) for c in molecule_counts.values())
        values = np.empty(size, dtype=int)
        row = np.empty(size, dtype=int)
        col = np.empty(size, dtype=int)
        i = 0
        for gene in molecule_counts:
            for cell, count in molecule_counts[gene].items():
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


def outer_join(left, right):
    """
    outer join two dictionaries of arrays on the union of their keys and return
    the expanded arrays in same order"""

    # get all keys
    all_keys = ArraySet(left.keys()).union(right.keys())

    # allocate empty arrays
    n = len(all_keys)
    left_array = np.zeros(n, dtype=np.float)
    right_array = np.zeros(n, dtype=np.float)

    # fill arrays
    for i, k in enumerate(all_keys):
        try:
            left_array[i] = left[k]
        except KeyError:
            left_array[i] = 0.
    for i, k in enumerate(all_keys):
        try:
            right_array[i] = left[k]
        except KeyError:
            right_array[i] = 0.

    return left_array, right_array
