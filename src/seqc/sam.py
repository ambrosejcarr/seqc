__author__ = 'ambrose'

import numpy as np
import os
import seqc
import tables as tb


class ReadArrayH5Writer:

    _filters = tb.Filters(complevel=5, complib='blosc')

    class _DataTable(tb.IsDescription):
        cell = tb.UInt64Col(pos=0)
        rmt = tb.UInt32Col(pos=1)
        n_poly_t = tb.UInt8Col(pos=2)
        valid_cell = tb.BoolCol(pos=3)
        dust_score = tb.UInt8Col(pos=4)
        rev_quality = tb.UInt8Col(pos=5)
        fwd_quality = tb.UInt8Col(pos=6)
        is_aligned = tb.BoolCol(pos=7)
        alignment_score = tb.UInt8Col(pos=8)

    def __init__(self, h5file):
        self._filename = h5file
        self._fobj = None
        self._is_open = False

    @property
    def file(self):
        return self._filename

    def open(self, mode='w'):
        if self._is_open:
            return self._fobj
        else:
            self._fobj = tb.open_file(self._filename, mode=mode, filters=self._filters,
                                      title='ReadArray data')
            self._is_open = True
            return self._fobj

    def close(self):
        if self._is_open:
            self._fobj.close()
        else:
            return

    def create(self, expectedrows):
        if not self._is_open:
            self.open()
        a = tb.UInt32Atom()
        # create the tables and arrays needed to store data
        self._fobj.create_table(self._fobj.root, 'data', self._DataTable, 'ReadArray.data',
                                filters=self._filters, expectedrows=expectedrows)
        self._fobj.create_earray(self._fobj.root, 'index', a, (0, 2), filters=self._filters,
                                 expectedrows=expectedrows)
        self._fobj.create_earray(self._fobj.root, 'features', a, (0,),
                                 filters=self._filters, expectedrows=expectedrows)
        self._fobj.create_earray(self._fobj.root, 'positions', a, (0,),
                                 filters=self._filters, expectedrows=expectedrows)

    def write(self, data):
        data, index, features, positions = data
        # dtable = self._fobj.root.data
        # dtable.append(data)
        # dtable.flush()
        self._fobj.root.data.append(data)
        self._fobj.root.data.flush()
        self._fobj.root.index.append(index)
        self._fobj.root.features.append(features)
        self._fobj.root.positions.append(positions)


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


def _iterate_chunk(samfile, n=int(1e9)):
    """open a sam file, yield chunks of size n bytes; default ~ 1GB"""
    with open(samfile, 'rb') as f:
        while True:
            d = f.read(n)
            if d:
                yield d
            else:
                break

def _average_quality(quality_string):
    """calculate the average quality of a sequencing read from ASCII quality string"""
    n_bases = len(quality_string)
    return (sum(ord(q) for q in quality_string) - n_bases * 33) // n_bases


def _process_multialignment(alignments, feature_converter):
    """translate a sam multi-alignment into a ReadArray row"""

    # all fields are identical except feature; get from first alignment
    first = alignments[0]

    rev_quality = _average_quality(first[10])
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
            # alignment = alignment.strip().split('\t')
            pos = int(alignment[3])
            flag = int(alignment[1])
            strand = '-' if (flag & 16) else '+'
            reference_name = alignment[2]
            feature = feature_converter.translate(strand, reference_name, pos)
            if feature:
                features += feature
                positions.append(pos)

    is_aligned = True if features else False

    rec = (cell, rmt, n_poly_t, valid_cell, trimmed_bases, rev_quality,
           fwd_quality, is_aligned, alignment_score)

    return rec, features, positions


def _process_chunk(chunk, feature_converter):
    """process the samfile chunk"""
    # decode the data, discarding the first and last record which may be incomplete
    records = [r.split('\t') for r in chunk.decode().strip('\n').split('\n')[1:-1]]

    # discard any headers
    while records[0][0].startswith('@'):
        records = records[1:]

    # over-allocate a numpy structured array, we'll trim n records == number of
    # reads that multi-aligned
    n = len(records)
    data = np.zeros((n,), dtype=_dtype)
    index = np.zeros((n, 2), dtype=np.uint32)
    features = np.zeros(n, dtype=np.uint32)
    positions = np.zeros(n, dtype=np.uint32)

    # process multi-alignments
    i = 0  # index for: data array, features & positions JaggedArray index
    j = 0  # index for: feature & position JaggedArrays
    s = 0  # index for: start of multialignment
    e = 1  # index for: end of multialignment
    while e < len(records):
        # get first unread record name
        name = records[s][0]

        # get record names until the name is different, this is a multialignment
        try:
            while name == records[e][0]:
                e += 1
        except IndexError:
            if e == len(records):
                pass
            else:
                print(e)
                print(len(records))
                print(records[e])

        multialignment = records[s:e]
        s = e
        e = s + 1

        # process the multialignment
        rec, feat, pos = _process_multialignment(multialignment, feature_converter)

        # fill data array
        data[i] = rec

        # fill the JaggedArrays
        if feat:
            n_features = len(feat)
            features[j:j + n_features] = feat
            positions[j:j + n_features] = pos
            index[i, :] = (j, j + n_features)
            j += n_features
        else:
            index[i, :] = (j, j)

        i += 1

    # trim all of the arrays
    data = data[:i]
    index = index[:i]
    features = features[:j]
    positions = positions[:j]

    return data, index, features, positions


def _write_chunk(data, h5_file):
    """write a data chunk to the h5 file"""
    data, index, features, positions = data
    dtable = h5_file.root.data
    dtable.append(data)
    dtable.flush()
    h5_file.root.index.append(np.ravel(index))
    h5_file.root.features.append(features)
    h5_file.root.positions.append(positions)


def to_h5(samfile, h5_name, n_processes, chunk_size, gtf, fragment_length=1000):
    """Process a samfile in parallel, dump results into an h5 database.

    Note that this method uses several shortcuts that have minor adverse consequences for
    the resulting data. We deem them worthwhile, but it is important to be aware.

    (1) it reads in large chunks. Thus, reads on the edge may be: (i) discarded, or (ii)
        improperly broken into multiple alignments despite being a part of a single
        multialignment

    """
    # get file size
    filesize = os.stat(samfile).st_size

    # get a bunch of records to check average size of a record
    with open(samfile) as f:
        records = []
        for line in f:
            if line.startswith('@'):
                continue
            records.append(line)
            if len(records) > 1000:
                break
    average_record_size = np.mean([len(r) for r in records])

    # estimate expected rows for the h5 database
    expectedrows = filesize / average_record_size

    # create a feature_converter object to convert genomic -> transcriptomic features
    fc = seqc.convert_features.ConvertFeatureCoordinates.from_gtf(gtf, fragment_length)

    read_kwargs = dict(samfile=samfile, n=chunk_size)
    process_kwargs = dict(feature_converter=fc)
    write_kwargs = dict(expectedrows=expectedrows)
    seqc.parallel.process_parallel(
        n_processes, h5_name, _iterate_chunk, _process_chunk, ReadArrayH5Writer,
        read_kwargs=read_kwargs, process_kwargs=process_kwargs, write_kwargs=write_kwargs)

    return h5_name


