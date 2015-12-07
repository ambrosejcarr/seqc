__author__ = 'Ambrose J. Carr'

# import warnings
# with warnings.catch_warnings():
from operator import itemgetter
import numpy as np
from copy import copy
from numpy.lib.recfunctions import append_fields
from itertools import islice
from collections import defaultdict
from scipy.sparse import coo_matrix
import numbers
import tables as tb
import pickle
from types import *
import os
from scipy.special import gammaln
from scipy.stats import chi2
import seqc

# register numpy integers as Integrals
numbers.Integral.register(np.integer)


def _multinomial_loglikelihood(x, probs):
    """return the multinomial log-likelihood of probs given x log(L) = Mult(probs | x)"""
    return (gammaln(np.sum(x) + 1) - np.sum(gammaln(x + 1))) + np.sum(x * np.log(probs))


def _likelihood_ratio_test(current_model, best_model, df):
    """likelihood ratio test, evaluated with a chi2 distribution with df=df

    note: assumes input parameters are log-likelihoods"""
    ratio = -2 * current_model + 2 * best_model
    return chi2.cdf(ratio, df)


class _UnionFind:
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


# todo | change jagged array slicing to return another jagged array
# todo | change from_iterable to be able to construct from a jagged array slice output
# todo test if id() is faster than .tobytes()
# todo | it would save more memory to skip "empty" features and give them equal indices
class JaggedArray:

    def __init__(self, data, index):
        self._data = data
        self._index = index

    @property
    def shape(self):
        return tuple(len(self))

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
        # if isinstance(item, numbers.Integral):
        i = self._index[item, :]
        try:
            return self._data[i[0]:i[1]]
        except IndexError:
            try:
                indices = self._index[item, :]
                return np.array([self.data[i[0]:i[1]].view(HashableArray)
                                 for i in indices])
            except:
                raise TypeError('__getitem__ supports int or arraytypes, not %s. If a '
                                'slice is desired us self.get_slice()'
                                % type(item))

    def __iter__(self):
        for i in range(len(self)):
            yield self[i]

    def get_slice(self):
        min_index = self.index[slice.start, 0]
        max_index = self.index[slice.stop, 1]
        new_index = self.index[slice, :]
        adj_index = new_index - min_index  # zero the index
        new_data = self.data[min_index:max_index]
        return self.__init__(new_data, adj_index)

    def _setitem(self, i, value):
        old = self[i]
        len_old = len(old)
        if len(value) > len_old:
            raise ValueError('Attempt to set longer sequence (%s; len=%d) than '
                             'allocated: %s; len=%d' %
                             (repr(value), len(value), repr(old), len(old)))
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

    def shrink(self):
        """eliminate any unused space in self._data

        if a method, such as ReadArray.resolve_alignments() has reduced the size of
        internal arrays via __setitem__(), there will be unused space in the array. This
        method will generate a new array object that uses the minimum amount of memory

        Note that this will generate a copy of self with memory usage <= self. The user
        is resposible for deleting the original object so that it can be garbage
        collected.
        """
        raise NotImplementedError  # todo implement
        # build a new index; identify any empty space
        index = np.zeros(self.index.shape, dtype=self.index.dtype)
        empty_size = 0
        current_position = 0
        for i in self._index:
            if i[0] > current_position:
                empty_size += i[0] - current_position
            index[i] = [current_position, len(self._data[i])]

        # check that there is no excess in the final position
        last_index = self._index[-1][1]
        if last_index > current_position:
            empty_size += last_index - current_position

        # if no empty space, return self
        if not empty_size:
            return self

        # otherwise, build new data
        new_size = self.data.shape[0] - empty_size
        data = np.zeros(new_size, dtype=self.data.dtype)
        for i, v in enumerate(self):
            data[index[i][0]:index[i][1]] = v

        return JaggedArray(data, index)

    def has_feature(self):
        """
        return a boolean vector indicating whether each entry in index is associated
        with a non-zero number of features
        """
        return self.index[:, 0] != self.index[:, 1]

    def is_unique(self):
        """
        return a boolean vector indicating whether each entry in index is assocaited
        with a single, unique feature.
        """
        return self.index[:, 0] + 1 == self.index[:, 1]

    def to_unique(self, bool_index=None, dtype=None):
        """
        return unique array features as a new np.ndarray.

        args:
        -----
        bool_index: If None, self.is_unique() is called to generate this index, otherwise
         it is assumed that this has been precalculated
        dtype: data type of new array. If None, uses the datatype of the old array.

        returns:
        --------
        np.array
        """
        if dtype is None:
            dtype = self.data.dtype

        if bool_index is None:
            index = self.index[self.is_unique(), 0]
        else:
            index = self.index[bool_index, 0]

        return self._data[index].astype(dtype)

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
    def from_iterable(cls, nested_iterable, dtype=np.uint32):
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

    @staticmethod
    def build_index():
        pass


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

    def __delitem__(self, key):
        # if key not in self, do nothing.
        if key in self:
            del self.counts[key]
            self._keys.remove(key)
            for i, v in enumerate(self._arrays):
                if v == key:
                    del self._arrays[i]

    def __contains__(self, key):
        return key in self._keys

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
        ('cell', np.uint64),
        ('rmt', np.uint32),
        ('n_poly_t', np.uint8),
        ('valid_cell', np.bool),
        ('dust_score', np.uint8),
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
    def shape(self):
        return self._data.shape

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
        return (self._data.nbytes + self.features.index.nbytes +
                self.features.data.nbytes + self.positions.data.nbytes +
                self.positions.index.nbytes)

    @classmethod
    # todo this is losing either the first or the last read
    # todo if necessary, this can be built out-of-memory using an h5 table
    def from_samfile(cls, samfile, feature_converter):

        # first pass through file: get the total number of alignments and the number of
        # multialignments for array construction; store indices for faster second pass
        with open(samfile) as f:

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
        idx = 0  # current index in jagged arrays
        with open(samfile) as f:
            records = iter(f)
            _ = list(islice(records, header_size))  # get rid of header
            for i, s in enumerate(ma_sizes):
                # get all records associated with a multialignment
                multi_alignment = list(islice(records, s))
                if not multi_alignment:  # iterator is exhausted
                    break
                recd, feat, pos = cls._process_multialignment(
                    multi_alignment, feature_converter)
                read_data[i] = recd
                if feat:
                    nfeatures = len(feat)
                    features[idx:idx + nfeatures] = feat
                    positions[idx:idx + nfeatures] = pos
                    index[i, :] = (idx, idx + nfeatures)
                    idx += nfeatures
                else:
                    index[i, :] = (idx, idx)

        # eliminate any extra size from the features & positions arrays
        features = features[:idx]
        positions = positions[:idx]

        jagged_features = JaggedArray(features, index)
        jagged_positions = JaggedArray(positions, index)

        return cls(read_data, jagged_features, jagged_positions)

    @staticmethod
    def _average_quality(quality_string):
        """calculate the average quality of a sequencing read from ASCII quality string"""
        n_bases = len(quality_string)
        return (sum(ord(q) for q in quality_string) - n_bases * 33) // n_bases

    @classmethod
    def _process_multialignment(cls, alignments, feature_converter):
        """translate a sam multi-alignment into a ReadArray row"""

        # all fields are identical except feature; get from first alignment
        first = alignments[0].strip().split('\t')

        rev_quality = cls._average_quality(first[10])
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
                features, positions = [], []
                return rec, features, positions
            else:  # single alignment
                pos = int(first[3])
                strand = '-' if (flag & 16) else '+'
                reference_name = first[2]
                features = feature_converter.translate(strand, reference_name, pos)
                if features:
                    positions = [pos]
                else:
                    positions = []
        else:
            features = []
            positions = []
            for alignment in alignments:
                alignment = alignment.strip().split('\t')
                pos = int(alignment[3])
                flag = int(alignment[1])
                strand = '-' if (flag & 16) else '+'
                reference_name = alignment[2]
                feature = feature_converter.translate(strand, reference_name, pos)
                if feature:
                    features += feature
                    positions.append(pos)

        # delete = deque()
        # for i in range(len(features)):
        #     if features[i] == 0:
        #         delete.append(i)
        #
        # features, positions = cls._multi_delete(delete, features, positions)

        is_aligned = True if features else False

        # if not features:
        #     features, positions = [0], [0]

        rec = (cell, rmt, n_poly_t, valid_cell, trimmed_bases, rev_quality,
               fwd_quality, is_aligned, alignment_score)

        # todo remove these checks
        assert(isinstance(features, list))
        assert(isinstance(positions, list))
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

        mask = self.mask_failing_records(n_poly_t_required=required_poly_t)

        # filter, then merge cell & rmt; this creates a copy of cell + rmt
        # memory overhead = n * 16b
        cells = self.data['cell'][mask]
        rmts = self.data['rmt'][mask]

        seq = [seqc.three_bit.ThreeBit.ints2int([int(c), int(r)])
               for c, r in zip(cells, rmts)]
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

        mask = self.mask_failing_records(n_poly_t_required=required_poly_t)

        # filter and merge cell/rmt
        # seq_data = np.vstack([self.data['cell'][mask],
        #                       self.data['rmt'][mask].astype(np.uint64)]).T
        cells = self.data['cell'][mask]
        rmts = self.data['rmt'][mask]

        seq = [seqc.three_bit.ThreeBit.ints2int([int(c), int(r)])
               for c, r in zip(cells, rmts)]
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
        """
        Resolve ambiguously aligned molecules and edit the ReadArray data structures
        in-place to reflect the more specific gene assignments.

        args:
        expectations: serialized coalignment expectation file. Often found in
          index/p_coalignment_array.p. Accepts either the loaded object or the filesystem
          location of the serialized object
        required_poly_t (default: 0): the number of poly_t required for a read to be
          considered a valid input for resolution. Reads lacking a poly_t tail are likely
          to randomly prime into non-mRNA sequences, so requiring the presence of poly_t
          in the 5' end of the forward read can reduce unwanted contamination
        alpha (default: 0.1): significance threshold for differentiating potential donor
          genes. lowering the alpha value causes the algorithm to require that models
          are better differentiated from inferior alternatives to declare the model
          differentiated. If the significance threshold is not reached, then the inferior
          model cannot be distinguished from the best model, and resolve_alignments fails
          to disambiguate the molecule in question.

        actions:
        edits self.features and self.data in-place to reflect resolved alignments. adds
          a column to self.data called 'disambiguation_results', this is an indicator
          column that tracks the result of this method.
          0: no model is fully supported by the data. This represents a molecule with
             higher-order complexity. These umi-cell combinations have two molecules
             associated with overlapping alignments and cannot be differentiated without
             considering multiple molecule contributions. This problem can be mitigated
             by increasing the UMI length or decreasing the concentration of input RNA.
          1: trivial result, molecule was unambiguous
          2: separating molecules into disjoint sets and limiting the algorithm to
          supported gene fully disambiguated the molecule
          3: molecule was partially disambiguated using the multinomial expectation model
          4: molecule was fully disambiguated using the multinomial expectation model
          5: model failed to remove any ambiguity
        """

        molecules, seqs = self.group_for_disambiguation(required_poly_t)
        results = np.zeros(self.data.shape[0], dtype=np.int8)

        # load the expectations
        if os.path.isfile(expectations):
            with open(expectations, 'rb') as f:
                expectations = pickle.load(f)
        elif isinstance(expectations, str):
            raise FileNotFoundError('could not locate serialized expectations object: %s'
                                    % expectations)
        elif isinstance(expectations, dict):
            pass  # expectations already loaded
        else:
            raise TypeError('invalid expecatation object type, must be a dict expectation'
                            ' object or a string filepath')

        # loop over potential molecules
        for read_indices in molecules.values():

            # get a view of the structarray
            # putative_molecule = self.data[read_indices]
            putative_features = self.features[read_indices]

            # check if all reads have a single, identical feature
            # can't hash numpy arrays, but can hash .tobytes() of the array
            check_trivial = ArrayCounter(putative_features)
            if len(check_trivial) == 1 and len(next(iter(check_trivial))) == 1:
                results[read_indices] = 1  # trivial
                continue

            # get disjoint set membership
            uf = _UnionFind()
            uf.union_all(putative_features)
            set_membership, sets = uf.find_all(putative_features)

            # Loop Over Disjoint molecules
            for s in sets:

                # disjoint_group = putative_molecule[set_membership == s]
                disjoint_features = putative_features[set_membership == s]
                disjoint_group_idx = read_indices[set_membership == s]  # track index

                # get observed counts
                obs_counter = ArrayCounter(putative_features)

                # get the subset of features that are supported by the data
                # these features should be present in all multi-alignments
                possible_models = set(disjoint_features[0])
                for molecule_group in disjoint_features:
                    possible_models.intersection_update(molecule_group)
                    if not possible_models:  # no molecule is supported by data.
                        continue  # could not resolve BUT a subset could be unique!

                # check that creating disjoint sets or eliminating unsupported molecules
                # haven't made this trivial
                if len(possible_models) == 1:
                    results[disjoint_group_idx] = 2
                    continue
                # if we no longer have any supported models, we cannot disambiguate this
                # molecule
                elif len(possible_models) == 0:
                    results[disjoint_group_idx] = 5
                    continue

                # convert possible models to np.array for downstream indexing
                possible_models = np.array(list(possible_models))

                model_likelihoods = np.empty_like(possible_models, dtype=np.float)
                df = np.empty_like(possible_models, dtype=np.int)

                for i, m in enumerate(possible_models):
                    try:
                        exp_dict = expectations[m]
                    except KeyError:
                        continue
                    except IndexError:
                        raise
                    obs, exp = outer_join(obs_counter, exp_dict)

                    # calculate model probability
                    model_likelihoods[i] = _multinomial_loglikelihood(obs, exp)
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
                    p = _likelihood_ratio_test(lh, top_likelihood, degrees_freedom)
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

        # todo can I avoid this? append_fields is very slow
        self._data = append_fields(self.data, 'disambiguation_results', results)
        # self._features = self.features.shrink()

    @staticmethod
    def _multi_delete(sorted_deque, *lists):
        while sorted_deque:
            i = sorted_deque.pop()
            for l in lists:
                del l[i]
        return lists

    # todo could be a bug here.
    def to_unique(self, n_poly_t_required):
        """Create a UniqueReadArray containing copying only unique reads from self"""

        fbool = ((self.data['cell'] != 0) &
                 (self.data['rmt'] != 0) &
                 (self.data['n_poly_t'] >= n_poly_t_required) &
                 (self.data['is_aligned']) &
                 (self.features.is_unique()))

        if not np.sum(fbool):
            raise ValueError('Cannot create UniqueReadArray; no records pass filters')

        data = self._data[fbool]
        features = self.features.to_unique(fbool)
        positions = self.positions.to_unique(fbool)

        return UniqueReadArray(data, features, positions)

    def mask_failing_records(self, n_poly_t_required):
        """Return a mask for failing records"""
        return ((self.data['cell'] != 0) &
                (self.data['rmt'] != 0) &
                (self.data['n_poly_t'] >= n_poly_t_required) &
                (self.data['is_aligned']) &
                (self.features.has_feature()))

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

        # check that there is data in array
        if not self.positions.data.shape[0]:
            raise ValueError('No feature positions were assigned to this object. Cannot '
                             'save an empty array.')
        if not self.features.data.shape[0]:
            raise ValueError('No features were assigned to this object. Cannot '
                             'save an empty array.')


        def store_carray(h5f, array, where, name):
            atom = tb.Atom.from_dtype(array.dtype)
            store = h5f.createCArray(where, name, atom, array.shape)
            store[:] = array

        blosc5 = tb.Filters(complevel=5, complib='blosc')
        f = tb.open_file(archive_name, mode='w', title='Data for seqc.ReadArray',
                         filters=blosc5)

        # store data
        f.create_table(f.root, 'data', self._data)

        # create group for feature-related data and store.
        store_carray(f, self.features.data, '/', 'features')
        store_carray(f, self.features.index, '/', 'index')
        store_carray(f, self.positions.data, '/', 'positions')

        f.close()

    @classmethod
    def from_h5(cls, archive_name):
        f = tb.open_file(archive_name, mode='r')

        data = f.root.data.read()

        fdata = f.root.features.read()
        findex = f.root.index.read()

        pdata = f.root.positions.read()
        pindex = f.root.index.read()

        features = JaggedArray(fdata, findex)
        positions = JaggedArray(pdata, pindex)

        f.close()

        return cls(data, features, positions)

    def plot_filter_correlations(self):
        """plot filter correlations"""
        raise NotImplementedError

    def plot_forward_quality_dist(self):
        data = self.data['fwd_quality']
        f, ax = seqc.plot.histogram(
            data, bins=40, fig=None, ax=None, xlabel='Quality', ylabel='Number of Reads',
            title='Forward Read Quality')
        return f, ax

    def plot_reverse_quality_dist(self):
        data = self.data['rev_quality']
        f, ax = seqc.plot.histogram(
            data, bins=40, fig=None, ax=None, xlabel='Quality', ylabel='Number of Reads',
            title='Reverse Read Quality')
        return f, ax

    def plot_number_aligned_features(self):
        """excludes single alignments"""
        data = np.array([len(f) for f in self.features])
        non_unique = data[data > 1]  # exclude unique alignments
        n_unique = non_unique.shape[0] / data.shape[0] * 100
        f, ax = seqc.plot.histogram(
            data, bins=40, fig=None, ax=None, xlabel='Quality', ylabel='Number of Reads',
            title='Reverse Read Quality')
        ax.text(0.9, 0.9, '%d Unique' % int(n_unique), horizontalalignment='right',
                verticalalignment='center', transform=ax.transAxes)
        return f, ax

    def plot_alignment_score(self):
        data = self.data['alignment_score']
        f, ax = seqc.plot.histogram(
            data, bins=40, fig=None, ax=None, xlabel='Quality', ylabel='Number of Reads',
            title='Reverse Read Quality')
        return f, ax

    def summarize(self, gtf, alignment_metadata=None, save_plots=False):
        """return descriptive statistics for this dataset

        Note that this function can take a very long time to run, particularly if
        save_plots is not False.

        summarize calculates several statistics:
        1. average reads per molecule, and the read per molecule distribution
        2. average molecules per cell, and the molecule per cell distribution
        3. percentage of records that:
           a. align to the transcriptome
           b. are phiX
           c. are mitochondrial
           d. that result from technical errors
           e. have a valid cell barcode
           f. contribute to molecules (pass all filters)
           g. aligned ambiguously
           h. aligned ambiguously but were corrected
           i. aligned to unique features
        4. nucleotide error rates for:
           a. each error type
           b. overall average error rate

        # todo implement
        if a folder location is provided to save_plots, summarize also generates several
        plots:
        1. Cell Yield across a moving threshold
        2. Impact of UMI/RMT correction
        3. Impact of disambiguation & error correction
        4. GC bias across cell barcodes
        5. filter correlation heatmap

        args:
        -----
        save_plots (str): folder location where plots should be saved. If False, not plots
          are generated

        returns:
        dictionary of summary statistics. If save_plots is not False, writes the result
          to save_plots/summary.json
        """
        # get some ids from the gtf file
        reader = seqc.gtf.Reader(gtf)
        phix_id = reader.get_phix_id()
        mt_ids = reader.get_mitochondrial_ids()

        tx_aligned = sum(1 for _ in self.features if _) / self.data.shape[0]
        phix_aligned = sum(1 for f in self.features if phix_id in f) / self.data.shape[0]
        mt_aligned = sum(1 for features in self.features if
                         any(f in mt_ids for f in features)) / self.data.shape[0]
        is_cell = np.sum(self.data['valid_cell'])

        # get number contributing to all filters


        # get error number todo | need to wait for Rami's error correction
        # is_error = np.sum(self.data['is_error']) / self.data.shape[0]


class UniqueReadArray:
    """
    Comparable, but faster datastructure to ReadArray that stores only unique alignments,
    in sorted order

    ideally we would append the fields but this copies data, is too slow. :-(
    """

    def __init__(self, data, features, positions):
        """
        args:
        -----
        data: numpy structured array
        features: np.array
        positions: np.array

        """
        self._data = data
        self._features = features
        self._positions = positions

        # pointer to Experiment holding read and molecule counts
        self._experiment = None

        # hold inds for molecule sort. Guarantees correct read order.
        self._sorted = None

    def __getitem__(self, item):
        if isinstance(item, slice):  # return array slice
            return UniqueReadArray(
                self.data[item], self.positions[item], self.features[item])
        else:  # return column
            return self.data[item], self.features[item], self.positions[item]

    def __repr__(self):
        return '<UniqueReadArray object with %d items>' % len(self)

    def __len__(self):
        return np.ravel(self._data.shape[0])[0]

    @property
    def shape(self):
        return self._data.shape

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
        return self.data.nbytes

    def _sort(self):
        """
        Lexicographical argsort (indirect) of self.

        args:
        -----
        molecules: if True, sorts on [cell, rmt, features]. In this case, stores both
         self._sort_molecules and self._sort_reads because molecule sort is also a read
         sort.

        actions:
        --------
        store result of lexsort in either self._sort_reads (molecules=False) or
         self._sort_reads and self._sort_molecules (molecules=True)
        """
        self._sorted = np.lexsort((self.data['rmt'], self.features, self.data['cell']))

    def to_experiment(self, required_support=0):
        """Generate an Experiment containing read and molecule SparseCount objects

        args:
        -----
        required_support (default: 2): required number of observations of a molecule for
         the molecule to be included in the counts matrix

        returns:
        seqc.analyze.Experiment object

        """

        # don't reprocess
        if self._experiment is not None:
            return self._experiment

        # get sorted index
        if self._sorted is None:
            self._sort()

        sort_ord = self._sorted

        # todo could be some 1-off bugs: Mostly working.
        # todo I believe the first molecule ends up missing 1 read
        # todo the final molecule may be restricted to 1 read.
        # I see this as minor, will fix later in the process of refactoring this stuff.

        # find boundaries between cell, rmt, and feature
        all_diff = np.zeros(len(self), dtype=np.bool)
        all_diff[0] = True  # keep the first one; it's different from void (preceding)
        all_diff[1:] |= np.diff(self.data['cell'][sort_ord]).astype(np.bool)  # diff cell
        all_diff[1:] |= np.diff(self.data['rmt'][sort_ord]).astype(np.bool)  # diff rmt
        all_diff[1:] |= np.diff(self.features[sort_ord]).astype(np.bool)  # diff feature

        # get key_index (into original ReadArray) and reads per molecule counts
        # adding required support here will cut any molecules with < required_support from
        # downstream steps that rely upon ra_molecule_idx (such as reads per molecule
        # calculations
        i = np.concatenate((np.where(all_diff)[0], [len(self)]))
        rpm_count = np.diff(i)

        # filter which molecules we want to keep by thresholding ra_molecule_idx
        ra_molecule_idx = sort_ord[i[np.concatenate((rpm_count > required_support, [False]))]]

        # make sure at least one molecule passes, otherwise raise.
        if len(ra_molecule_idx) == 0:
            raise ValueError('No molecules passed support filter')

        # to get reads per cell, I discard the notion of molecular correction by diffing
        # on only cell and feature from the original sort
        # diff is working as intended
        rpc_diff = np.zeros(len(self), dtype=np.bool)
        rpc_diff[0] = True
        rpc_diff[1:] |= np.diff(self.data['cell'][sort_ord]).astype(np.bool)
        rpc_diff[1:] |= np.diff(self.features[sort_ord]).astype(np.bool)
        i = np.concatenate((np.where(rpc_diff)[0], [len(self)]))
        ra_rpc_idx = sort_ord[i[:-1]]
        rpc_count = np.ravel(np.diff(i))

        # because reads per molecule gives me a list of all molecules, molecules per
        # cell can be calculated by re-diffing on the reads per molecule without
        # considering the rmt. This has the effect of counting unique RMTs per molecule
        # and per cell.
        mpc_diff = np.zeros(len(ra_molecule_idx), dtype=np.bool)
        mpc_diff[0] = True
        mpc_diff[1:] |= np.diff(self.data['cell'][ra_molecule_idx]).astype(np.bool)
        mpc_diff[1:] |= np.diff(self.features[ra_molecule_idx]).astype(np.bool)
        i = np.concatenate((np.where(mpc_diff)[0], [len(ra_molecule_idx)]))
        ra_rpm_idx = ra_molecule_idx[i[:-1]]
        mpc_count = np.ravel(np.diff(i))

        def map_to_unique_index(vector):
            """function to map a vector of gene or feature ids to an index"""
            ids = np.unique(np.ravel(vector))
            # key gives the new values you wish original ids to be mapped to.
            key = np.arange(ids.shape[0])
            index = np.digitize(vector.ravel(), ids, right=True)
            return key[index].reshape(vector.shape), ids

        # generate sparse matrices

        # reads per cell
        row, cells = map_to_unique_index(self.data['cell'][ra_rpc_idx])
        col, genes = map_to_unique_index(self.features[ra_rpc_idx])
        shape = (len(cells), len(genes))
        rpc = coo_matrix((rpc_count, (row, col)), shape=shape, dtype=np.int32)
        reads_per_cell = seqc.analyze.SparseCounts(rpc, cells, genes)

        # molecules per cell
        row, cells = map_to_unique_index(self.data['cell'][ra_rpm_idx])
        col, genes = map_to_unique_index(self.features[ra_rpm_idx])
        shape = (len(cells), len(genes))
        mpc = coo_matrix((mpc_count, (row, col)), shape=shape, dtype=np.int32)
        molecules_per_cell = seqc.analyze.SparseCounts(mpc, cells, genes)

        self._experiment = seqc.analyze.Experiment(reads_per_cell, molecules_per_cell)
        return self._experiment

    @staticmethod
    def from_read_array(ra, n_poly_t_required):
        """
        Construct a unique ReadArray. Note that this will copy data. Requires O(2n) memory

        args:
        -----
        ra: ReadArray we are copying from
        filter_records: If true, will eliminate any records that fail filters before
         constructing the UniqueReadArray

        Can think about ways to get counts via searchsorted
        """
        return ra.to_unique(n_poly_t_required)

    def save_h5(self, archive_name):
        """
        efficiently save the UniqueReadArray to a compressed h5 archive; note that this
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
        f.create_table(f.root, 'data', self.data)
        store_carray(f, self.features, f.root, 'features')
        store_carray(f, self.positions, f.root, 'positions')

        f.close()

    @classmethod
    def from_h5(cls, archive_name):
        f = tb.open_file(archive_name, mode='r')

        data = f.root.data.read()
        features = f.root.features.read()
        positions = f.root.positions.read()

        f.close()

        return cls(data, features, positions)


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
