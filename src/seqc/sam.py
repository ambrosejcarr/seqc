__author__ = 'ambrose'

from seqc.convert_features import construct_feature_table
from seqc import three_bit
from seqc.qc import multinomial_loglikelihood, likelihood_ratio_test
from array import array
from threading import Thread
from queue import Queue, Full, Empty
from time import sleep
from collections import defaultdict, Counter
from seqc.qc import UnionFind, merge_observations_and_expectations
import numpy as np
import collections
from itertools import islice, tee, chain
import pickle
from subprocess import Popen, PIPE, check_output
import numpy.lib.recfunctions as rfn


def multi_delete(sorted_deque, *lists):
    while sorted_deque:
        i = sorted_deque.pop()
        for l in lists:
            del l[i]
    return lists


def mask_failing_cells(df_or_array, n_poly_t_required):
    return ((df_or_array['cell'] != 0) &
            (df_or_array['rmt'] != 0) &
            (df_or_array['n_poly_t'] >= n_poly_t_required) &
            (df_or_array['is_aligned'])
            )


"""
the sole purpose of the two classes below is to defeat the confounding way that numpy
treats tuples with len == 1 when constructing arrays by obscuring the existence of tuple
data by 1 level. Pass a tuple as a feature or it will break
"""


class ObfuscatedTuple:

    def __init__(self, tuple_):
        if not isinstance(tuple_, tuple):
            tuple_ = tuple(tuple_)
        self._tuple = tuple_

    def __repr__(self):
        return '<ObfuscatedTuple: %s>' % self._tuple.__repr__()

    def __lt__(self, other):
        return self._tuple.__lt__(other)

    def __gt__(self, other):
        return self._tuple.__gt__(other)

    def __eq__(self, other):
        return self._tuple.__eq__(other)

    def __ne__(self, other):
        return self._tuple.__ne__(other)

    def __le__(self, other):
        return self._tuple.__le__(other)

    def __ge__(self, other):
        return self._tuple.__ge__(other)

    def __len__(self):
        return len(self._tuple)

    def __contains__(self, item):
        return self._tuple.__contains__(item)

    def __iter__(self):
        return iter(self._tuple)

    def __hash__(self):
        return hash(self._tuple)

    def to_tuple(self):
        return self._tuple


class Peekable(collections.Iterator):

    def __init__(self, it):
        self.it, self.nextit = tee(iter(it))
        self._advance()

    def _advance(self):
        self.peek = next(self.nextit, None)

    def __next__(self):
        self._advance()
        return next(self.it)


def average_quality(quality_string):
    """calculate the average quality of a sequencing read from and ASCII quality string"""
    n_bases = len(quality_string)
    return (sum(ord(q) for q in quality_string) - n_bases * 33) // n_bases


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


def group_multialignments(alignments):
    iterator = Peekable(alignments)

    # get first non-header alignment
    multialignment = next(iterator)
    while multialignment.startswith('@'):
        multialignment = next(iterator)
    multialignment = [multialignment]

    while True:
        next_alignment = iterator.peek
        if not next_alignment:
            yield multialignment
            break
        if multialignment[0].split('\t')[0] == next_alignment.split('\t')[0]:
            multialignment.append(next(iterator))
        else:
            yield multialignment
            multialignment = [next(iterator)]


def process_multialignment(alignments, feature_positions, feature_table):
    """translate a sam record into a recarray row"""

    # all fields are identical except feature; get from first alignment
    first = alignments[0].strip().split('\t')

    rev_quality = average_quality(first[10])
    alignment_score = int(first[13].split(':')[-1])

    # parse data from name field, previously extracted from forward read
    forward_metadata = (int(f) for f in first[0].strip().split(';')[0].split(':'))
    cell, rmt, n_poly_t, valid_cell, trimmed_bases, fwd_quality = forward_metadata

    # get all features and positions
    if len(alignments) == 1:
        true_position = int(first[3])
        flag = int(first[1])
        strand = '-' if (flag & 16) else '+'
        reference_name = first[2]
        features = [translate_feature(reference_name, strand, true_position,
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
            features.append(translate_feature(reference_name, strand, true_position,
                                              feature_table, feature_positions))
            positions.append(true_position)

    delete = collections.deque()
    for i in range(len(features)):
        if features[i] == 0:
            delete.append(i)

    features, positions = multi_delete(delete, features, positions)

    positions = ObfuscatedTuple(positions)
    features = ObfuscatedTuple(features)
    is_aligned = True if features else False

    rec = (cell, rmt, n_poly_t, valid_cell, trimmed_bases, rev_quality, fwd_quality,
           features, positions, is_aligned, alignment_score)

    return rec


class JaggedArray:

    def __init__(self, data, index):
        self._data = data
        self._index = index

    @property
    def data(self):
        return self._data

    @property
    def index(self):
        return self._index

    def __len__(self):
        return len(self._index)

    def __getitem__(self, item):
        """"""
        if isinstance(item, slice):
            min_index = self.index[slice.start, 0]
            max_index = self.index[slice.stop, 1]
            new_index = self.index[slice, :]
            adj_index = new_index - min_index  # zero the index
            new_data = self.data[min_index:max_index]
            return self.__init__(new_data, adj_index)
        else:
            i = self._index[item, :]
            return self._data[i[0]:i[1]]

    def __setitem__(self, i, value):
        """
        allows a feature set to be re-set by adjusting the index and values. Note that
        increasing the number of features is not supported, and will raise a ValueError.
        Features can only be changed or decreased in number. Note that setting features
        in this way generates unused space in the array.
        """
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

    def __repr__(self):
        if len(self.data) < 10:
            return '<JaggedArray:\n[' + '\n '.join(str(i) for i in self) + ']>'
        else:
            return ('<JaggedArray:\n[' + '\n '.join(str(self[i]) for i in range(5)) +
                    '\n ...\n ' +
                    '\n '.join(str(self[i]) for i in range(len(self) - 5, len(self))) +
                    ']>'
                    )

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

    @classmethod
    def from_samfile(cls, samfile, feature_positions, feature_table):

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
                if not multi_alignment:
                    break
                recd, feat, posn = cls._process_multialignment(
                    multi_alignment, feature_positions, feature_table)
                read_data[i] = recd
                features[current_index_position:current_index_position + s] = feat
                positions[current_index_position:current_index_position + s] = posn
                index[i, :] = (current_index_position, current_index_position + s)
                current_index_position += s

        # eliminate any extra size from the features & positions arrays
        features = features[:current_index_position]
        positions = positions[:current_index_position]

        jagged_features = JaggedArray(features, index)
        jagged_positions = JaggedArray(positions, index)

        return cls(read_data, jagged_features, jagged_positions)

    @staticmethod
    def _process_multialignment(alignments, feature_positions, feature_table):
        """translate a sam record into a recarray row"""

        # all fields are identical except feature; get from first alignment
        first = alignments[0].strip().split('\t')

        # todo remove; debugging
        try:
            rev_quality = average_quality(first[10])
        except IndexError:
            print(first)
            raise
        alignment_score = int(first[13].split(':')[-1])

        # parse data from name field, previously extracted from forward read
        forward_metadata = (int(f) for f in first[0].strip().split(';')[0].split(':'))
        cell, rmt, n_poly_t, valid_cell, trimmed_bases, fwd_quality = forward_metadata

        # get all features and positions
        if len(alignments) == 1:
            true_position = int(first[3])
            flag = int(first[1])
            strand = '-' if (flag & 16) else '+'
            reference_name = first[2]
            features = [translate_feature(reference_name, strand, true_position,
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
                features.append(translate_feature(reference_name, strand, true_position,
                                                  feature_table, feature_positions))
                positions.append(true_position)

        delete = collections.deque()
        for i in range(len(features)):
            if features[i] == 0:
                delete.append(i)

        features, positions = multi_delete(delete, features, positions)

        if not features:
            features, positions = 0, 0

        is_aligned = True if features else False

        rec = (cell, rmt, n_poly_t, valid_cell, trimmed_bases, rev_quality, fwd_quality,
               is_aligned, alignment_score)

        return rec, features, positions

    def group_for_disambiguation(self, required_poly_t):
        """
        generate an heirarchical dictionary that groups the indices of reads according to
        their gene and sequence information

        any indices associated with a failing (masked) read is not included in the
        index dictionary, and so is implicity given a zero value in the results vector
        """
        indices = np.arange(self.data.shape[0])

        mask = mask_failing_cells(self.data, n_poly_t_required=required_poly_t)

        # filter, then merge cell & rmt; this creates a copy of cell + rmt
        # memory overhead = n * 16b
        seq_data = np.vstack([self.data['cell'][mask],
                              self.data['rmt'][mask].astype(np.int64)]).T
        seq = np.apply_along_axis(three_bit.ThreeBit.ints2int, axis=1,
                                  arr=seq_data)
        indices = indices[mask]

        # get indices of reads associated with each putative molecule (rmt/cell pair)
        molecules = defaultdict(list)
        for i, s in zip(indices, seq):
            molecules[s].append(i)
        for k, v in molecules.items():
            molecules[k] = np.array(v)  # covert all lists to arrays for indexing

        return molecules, seq

    def group_for_error_correction(self, required_poly_t):
        """better if disambiguation is done by now!"""
        indices = np.arange(self.data.shape[0])

        mask = mask_failing_cells(self.data, n_poly_t_required=required_poly_t)

        # filter and merge cell/rmt
        seq_data = np.vstack([self.data['cell'][mask],
                              self.data['rmt'][mask].astype(np.int64)]).T
        seq = np.apply_along_axis(three_bit.ThreeBit.ints2int, axis=1,
                                  arr=seq_data)

        indices = indices[mask]

        # get indices of reads associated with each putative molecule (gene/rmt/cell)
        # cannot hash numpy arrays; instead hash the bytes representation of the array
        # this means that the feature (f) is giberish, but that's not important for
        # the experiment.

        molecules = defaultdict(dict)
        for i, s, f in zip(indices, seq, self.features):
            try:
                molecules[f.tobytes()][s].append(i)
            except AttributeError:
                molecules[np.array(f).tobytes()][s].append(i)
            except KeyError:
                try:
                    molecules[f.tobytes()][s] = [i]
                except AttributeError:
                    molecules[np.array(f).tobytes()][s].append(i)
        for f in molecules.keys():
            for s, v in molecules[f].items():
                molecules[f][s] = np.array(v)

        return dict(molecules), seq

    def correct_errors(self):
        raise NotImplementedError

    def disambiguate(self, expectations, required_poly_t=0, alpha=0.1):
        molecules, seqs = self.group_for_disambiguation(required_poly_t)
        results = np.zeros(self.data.shape[0], dtype=np.int8)

        # load the expectations
        if isinstance(expectations, str):
            with open(expectations, 'rb') as f:
                expectations = pickle.load(f)

        # loop over potential molecules
        for read_indices in molecules.values():

            # get a view of the structarray
            putative_molecule = self.data[read_indices]
            putative_features = self.features[read_indices]

            # check if all reads have a single, identical feature
            check_trivial = Counter(putative_features)
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
                disjoint_features = putative_features[set_membership == s]
                disjoint_group_idx = read_indices[set_membership == s]  # track index

                # get observed counts
                obs_counter = Counter(f.to_tuple() for f in putative_features)
                # obs_features = np.array(list(obs_counter.keys()))
                # obs_counts = np.array(list(obs_counter.values()))

                # check that creating disjoint sets haven't made this trivial
                if len(obs_counter) == 1 and len(next(iter(obs_counter))) == 1:
                    results[disjoint_group_idx] = 2
                    continue  # no need to calculate probabilities for single model

                # get all potential molecule identities, create container for lh, df
                possible_models = np.array(
                    list(set(chain(*(f.to_tuple() for f in disjoint_features))))
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
                self.features[disjoint_group_idx] = passing_models

        return results, self.data


def create_structured_array(n):
    """pre-allocate a recarray of size n for sam processing

    note that the complete array will often not be utilized, since multimapping fastq
    records will be compressed to a single recarray row"""

    dtype = [
        ('cell', np.int64),
        ('rmt', np.int32),
        ('n_poly_t', np.uint8),
        ('valid_cell', np.bool),
        ('trimmed_bases', np.uint8),
        ('rev_quality', np.uint8),
        ('fwd_quality', np.uint8),
        ('features', np.object),
        ('positions', np.object),
        ('is_aligned', np.bool),
        ('alignment_score', np.uint8)
    ]

    # trade-off: 50-95% of features will be unique. Therefore, it makes more sense to
    # store them as ints. lists have an overhead of ~ 64b; dictionaries have overhead of
    # 288 bytes (?) so it looks like we should have a VL array.
    # each record in the sam array corresponds to an entry in the VL array; we could know
    # the fixed size. To use a VL array would require two pointer lists, each with int8
    # size, and then the data itself.

    return np.zeros((n,), dtype=dtype)


def process_alignments(samfile, n_threads, gtf, fragment_len):

    """create a record array containing alignment, barcode, and filtering information"""

    def read(in_queue):
        s = open(samfile)
        iterator = Peekable(s)
        try:
            while True:
                try:
                    # get data
                    data = list(islice(iterator, 10000))
                    # make sure we haven't bisected a multialignment
                    final_record = data[-1]
                    while True:
                        next_record = iterator.peek
                        if not next_record:
                            raise StopIteration
                        if final_record.split('\t')[0] == next_record.split('\t')[0]:
                            data.append(next(iterator))
                        else:
                            break
                except StopIteration:
                    if data:
                        while True:  # put the final data chunk on the queue
                            try:
                                in_queue.put(data)
                                break
                            except Full:
                                sleep(0.1)
                                continue
                    break

                # put the read data on the queue
                while True:
                    try:
                        in_queue.put(data)
                        break
                    except Full:
                        sleep(0.1)
                        continue
        finally:
            s.close()

    def process(in_queue, out_queue, feature_positions, feature_table):
        while True:
            # get data and process it into records
            try:
                data = in_queue.get_nowait()

                # process the alignment group
                arr = create_structured_array(10000)
                for i, ma in enumerate(group_multialignments(data)):
                    row = process_multialignment(ma, feature_positions, feature_table)
                    arr[i] = row
                arr = arr[0:i + 1]  # cut any unfilled rows
            except Empty:  # wait a bit for the next set of data
                if not read_thread.is_alive():
                    break
                else:
                    sleep(0.1)
                    continue

            # put the data on the concatenation queue
            while True:
                try:
                    out_queue.put(arr)
                    break
                except Full:
                    sleep(0.1)
                    continue

    # create a feature table todo | load from index
    feature_table_, feature_positions_ = construct_feature_table(gtf, fragment_len)

    # read the files
    alignment_groups = Queue()
    read_thread = Thread(target=read, args=[alignment_groups])
    read_thread.start()

    # process the data
    processed_alignments = Queue()
    process_threads = [Thread(target=process,
                              args=[alignment_groups, processed_alignments,
                                    feature_positions_, feature_table_])
                       for _ in range(n_threads - 3)]

    for t in process_threads:
        t.start()

    # wait for each process to finish
    read_thread.join()
    for t in process_threads:
        t.join()

    # once all processing threads are dead, concatenate all of the arrays

    # get the first array
    while True:
        try:
            arr = processed_alignments.get_nowait()
            break
        except Empty:
            if not any(t.is_alive() for t in process_threads):
                raise ValueError('no processed input received')
            else:
                sleep(0.1)
                continue

    # concatenate all the other arrays
    while True:
        try:
            next_arr = processed_alignments.get_nowait()
            arr = np.hstack([arr, next_arr])
        except Empty:
            if not any(t.is_alive() for t in process_threads):
                break
            else:
                sleep(0.1)
                continue
    return arr


def process_fluidigm_data():
    """best we can do with this data is reads per gene

    process each gene with htseq-count
    """


def assess_outliers(samfile):
    """
    parse samfile, generate reads per cell, & molecules per cell lists and total
    alignments.

    Normalize two ways:
     1: total input reads
     2: aligned reads
    """
    raise NotImplementedError