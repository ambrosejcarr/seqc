import numpy as np

# todo required api
# indexing of individual values
# setting of individual values
# length of returned values
# comparison to returned values

# JaggedArray has a slightly awkward API with respect to uniqueness due to numpy slicing
# basics: a slice that would return a single integer returns a scalar, whereas a slice
# that would return a tensor returns a numpy.ndarray

class JaggedArray:

    def __init__(self, data, index):
        diff = index[:, 0] - index[:, 1]
        if np.any(np.equal(diff, 0)):
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

    # todo note that under the current numpy release this is _very_ slow. Need to download
    # at least 1.11.0.dev0-f428bce
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

    # todo add shrink() call to resolve_alignments()
    # todo write tests
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
        # # build a new index; identify any empty space
        # index = np.zeros(self.index.shape, dtype=self.index.dtype)
        # empty_size = 0
        # current_position = 0
        # for i in self._index:
        #     if i[0] > current_position:
        #         empty_size += i[0] - current_position
        #     index[i] = [current_position, len(self._data[i])]
        #
        # # check that there is no excess in the final position
        # last_index = self._index[-1][1]
        # if last_index > current_position:
        #     empty_size += last_index - current_position
        #
        # # if no empty space, return self
        # if not empty_size:
        #     return self
        #
        # # otherwise, build new data
        # new_size = self.data.shape[0] - empty_size
        # data = np.zeros(new_size, dtype=self.data.dtype)
        # for i, v in enumerate(self):
        #     data[index[i][0]:index[i][1]] = v
        #
        # return JaggedArray(data, index)


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

    @staticmethod
    def build_index():
        pass


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



                # class ArraySet:
#
#     def __init__(self, iterable_of_arrays):
#         """identifies a set of unique arrays"""
#         self._seen = set()
#         self._arrays = []
#         for array in iterable_of_arrays:
#             hashed = hash(array.tobytes())
#             if not hashed in self._seen:
#                 self._seen.add(hashed)
#                 self._arrays.append(array)
#
#     def __iter__(self):
#         return iter(self._arrays)
#
#     def __len__(self):
#         return len(self._arrays)
#
#     def __repr__(self):
#         if len(self) > 10:
#             return '<ArraySet {%s ...}>' % ', '.join(repr(a) for a in self._arrays)
#         else:
#             return '<ArraySet {%s}>' % ', '.join(repr(a) for a in self._arrays)
#             # repr(np.array(self._arrays, dtype=object))
#
#     def __contains__(self, array):
#         return True if hash(array.tobytes()) in self._seen else False
#
#     def add(self, array):
#         try:
#             hashed = hash(array.tobytes())
#         except AttributeError:
#             raise ValueError('array must be a np.ndarray object, not %s' % type(array))
#         if array.tobytes() not in self._seen:
#             self._arrays.append(array)
#             self._seen.add(hashed)
#
#     def union(self, other):
#         union_ = copy(self)
#         for v in other:
#             union_.add(v)
#         return union_
#
#     def tolist(self):
#         return self._arrays


# class ArrayCounter:
#
#     def __init__(self, iterable_of_arrays):
#         self.counts = {}
#         self._arrays = []
#         self._keys = set()
#         for a in iterable_of_arrays:
#             hashed = hash(a.tobytes())
#             if hashed not in self._keys:
#                 self._keys.add(hashed)
#                 self._arrays.append(a)
#                 self.counts[hashed] = 1
#             else:
#                 self.counts[hashed] += 1
#
#     def __getitem__(self, item):
#         try:
#             return self.counts[hash(item.tobytes())]
#         except AttributeError:
#             if not isinstance(item, np.ndarray):
#                 raise ValueError('item must be a np.array object')
#             else:
#                 raise
#
#     def __delitem__(self, key):
#         # if key not in self, do nothing.
#         if key in self:
#             del self.counts[key]
#             self._keys.remove(key)
#             for i, v in enumerate(self._arrays):
#                 if v == key:
#                     del self._arrays[i]
#
#     def __contains__(self, key):
#         return key in self._keys
#
#     def __iter__(self):
#         return iter(self._arrays)
#
#     def __repr__(self):
#         sorted_data = sorted(self.items(), key=itemgetter(1), reverse=True)
#         top10 = ', '.join('%s: %d' % (repr(i[0]), i[1]) for i in sorted_data[:10])
#         if len(sorted_data) > 10:
#             return '<ArrayCounter {%s ...}>' % top10
#         else:
#             return '<ArrayCounter {%s}>' % top10
#
#     def __len__(self):
#         return len(self._arrays)
#
#     def keys(self):
#         for k in self._arrays:
#             yield k
#
#     def values(self):
#         for k in self._keys:
#             yield self[k]
#
#     def items(self):
#         for k in self._arrays:
#             yield k, self[k]