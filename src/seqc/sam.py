__author__ = 'ambrose'

import numpy as np
import pickle
import os
import io
import seqc
import tables as tb
import gzip
from collections import defaultdict, namedtuple
from scipy.sparse import coo_matrix



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

    # over-allocate a numpy structured array, we'll trim n == number of
    # reads that multi-aligned records at the end
    n = len(records)
    data = np.zeros((n,), dtype=_dtype)
    index = np.zeros((n, 2), dtype=np.uint32)
    features = np.zeros(n, dtype=np.int64)  # changed feature storage for hashed genes
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


class EmptyAligmentFile(Exception):
    pass


def to_h5(samfile, h5_name, n_processes, chunk_size, gtf, fragment_length=1000):
    """Process a samfile in parallel, dump results into an h5 database.

    Note that this method uses several shortcuts that have minor adverse consequences for
    the resulting data. We deem them worthwhile, but it is important to be aware.

    (1) it reads in large chunks. Thus, reads on the edge may be: (i) discarded, or (ii)
        improperly broken into multiple alignments despite being a part of a single
        multialignment

    """
    # check that samefile is non-empty:
    with open(samfile, 'r') as f:
        empty = True
        for line in f:
            if not line.startswith('@'):
                empty = False
                break
    if empty:
        raise EmptyAligmentFile('Alignment may have failed. Sam file has no alignments.')

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
    # todo this is where the conversion method is defined; swapped genes in here.
    # fc = seqc.convert_features.ConvertFeatureCoordinates.from_gtf(gtf, fragment_length)
    fc = seqc.convert_features.ConvertGeneCoordinates.from_gtf(gtf, fragment_length)
    fc_name = h5_name.replace('.h5', '_gene_id_map.p')
    fc.pickle(fc_name)  # save the id map

    read_kwargs = dict(samfile=samfile, n=chunk_size)
    process_kwargs = dict(feature_converter=fc)
    write_kwargs = dict(expectedrows=expectedrows)
    seqc.parallel.process_parallel(
        n_processes, h5_name, _iterate_chunk, _process_chunk, ReadArrayH5Writer,
        read_kwargs=read_kwargs, process_kwargs=process_kwargs, write_kwargs=write_kwargs)

    return h5_name


class GenerateSam:

    @staticmethod
    def in_drop(n, prefix, fasta, gtf, index, barcodes, tag_type='gene_id', replicates=3,
                n_threads=7, *args, **kwargs):
        """generate an in-drop .sam file"""

        with open(barcodes, 'rb') as f:
            cb = pickle.load(f)

        output_dir = '/'.join(prefix.split('/')[:-1]) + '/'
        forward, reverse = seqc.fastq.GenerateFastq.in_drop(
            n, prefix, fasta, gtf, barcodes, tag_type=tag_type, replicates=replicates)

        # merge the generated fastq file
        merged = seqc.fastq.merge_fastq(forward, reverse, 'in-drop', output_dir, cb,
                                        n_threads)

        # generate alignments from merged fastq
        sam = seqc.align.STAR.align(merged, index, n_threads, output_dir)
        return sam


    @staticmethod
    def drop_seq(n, prefix, fasta, gtf, index, tag_type='gene_id', replicates=3,
                 n_threads=7, *args, **kwargs):
        """generate a drop-seq .sam file"""

        output_dir = '/'.join(prefix.split('/')[:-1]) + '/'
        forward, reverse = seqc.fastq.GenerateFastq.drop_seq(
            n, prefix, fasta, gtf, tag_type=tag_type, replicates=replicates)

        # merge the generated fastq file
        merged = seqc.fastq.merge_fastq(forward, reverse, 'drop-seq', output_dir,
                                        n_threads)


        # generate alignments from merged fastq
        sam = seqc.align.STAR.align(merged, index, n_threads, output_dir)
        return sam


def to_count_single_file(sam_file, gtf_file):
    """cannot separate file due to operating system limitations. Instead, implement
    a mimic of htseq-count that uses the default 'union' approach to counting, given the
    same gtf file; for fluidigm/SMART data"""

    # get conversion table, all possible genes for the count matrix
    gt = seqc.convert_features.GeneTable(gtf_file)
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

    # estimate the average read length
    with open(sam_file, 'r') as f:
        sequences = []
        line = f.readline()
        while line.startswith('@'):
            line = f.readline()
        while len(sequences) < 100:
            sequences.append(f.readline().strip().split('\t')[9])
        read_length = round(np.mean([len(s) for s in sequences]))

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

                try:
                    genes = gt.coordinates_to_gene_ids(chromosome, start, end, strand)
                except KeyError:
                    continue  # todo include these non-chromosome scaffolds instead of
                    # discarding
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
    gene_index = np.array(sorted(all_genes) + ['ambiguous', 'no_feature', 'not_aligned'],
                          dtype=object)
    cell_index = np.array(['no_cell'] + list(range(1, cell_number)), dtype=object)

    return coo, gene_index, cell_index


def to_count_multiple_files(sam_files, gtf_file):
    """count genes in each cell; fluidigm data"""
    gt = seqc.convert_features.GeneTable(gtf_file)
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

        # estimate the average read length
        with open(sam_file, 'r') as f:
            sequences = []
            line = f.readline()
            while line.startswith('@'):
                line = f.readline()
            while len(sequences) < 100:
                sequences.append(f.readline().strip().split('\t')[9])
            read_length = round(np.mean([len(s) for s in sequences]))

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

                    try:
                        genes = gt.coordinates_to_gene_ids(chromosome, start, end, strand)
                    except KeyError:
                        continue  # todo count these weird non-chromosome scaffolds
                        # right now, we just throw them out...
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
    gene_index = np.array(sorted(all_genes) + ['ambiguous', 'no_feature', 'not_aligned'],
                          dtype=object)
    cell_index = np.array(['no_cell'] + list(range(1, cell_number)), dtype=object)

    return coo, gene_index, cell_index


class SamRecord:
    """simple record object to use when iterating over sam files"""

    __slots__  = ['qname', 'flag', 'rname', 'pos', 'mapq', 'cigar', 'rnext', 'pnext',
                  'tlen', 'seq', 'qual', 'optional_fields']

    def __init__(self, qname, flag, rname, pos, mapq, cigar, rnext, pnext, tlen, seq,
                 qual, *optional_fields):
        self.qname = qname
        self.flag = flag
        self.rname = rname
        self.pos = pos
        self.mapq = mapq
        self.cigar = cigar
        self.rnext = rnext
        self.pnext = pnext
        self.tlen = tlen
        self.seq = seq
        self.qual = qual
        self.optional_fields = optional_fields

    def __repr__(self):
        return ('<SamRecord:' + ' %s' * 12 + '>') % \
               (self.qname, self.flag, self.rname, self.pos, self.mapq, self.cigar,
                self.rnext, self.pnext, self.tlen, self.seq, self.qual,
                ' '.join(self.optional_fields))

class Reader:
    """simple sam reader, optimized for utility rather than speed"""

    def __init__(self, samfile):

        seqc.util.check_type(samfile, str, 'samfile')
        seqc.util.check_file(samfile, 'samfile')

        self._samfile = samfile
        try:
            samfile_iterator = iter(self)
            next(samfile_iterator)
        except:
            raise ValueError('%s is an invalid samfile. Please check file formatting.' %
                             samfile)

    @property
    def samfile(self):
        return self._samfile

    def _open(self) -> io.TextIOBase:
        """
        seamlessly open self._samfile, whether gzipped or uncompressed

        returns:
        --------
        fobj: open file object
        """
        if self._samfile.endswith('.gz'):
            fobj = gzip.open(self._samfile, 'rt')
        else:
            fobj = open(self._samfile)
        return fobj

    def __len__(self):
        return sum(1 for _ in self)

    def __iter__(self):
        """return an iterator over all non-header records in samfile"""
        fobj = self._open()
        try:
            for line in fobj:
                if line.startswith('@'):
                    continue
                yield SamRecord(*line.strip().split('\t'))
        finally:
            fobj.close()

    def iter_multialignments(self):
        """yields tuples of all alignments for each fastq record"""
        sam_iter = iter(self)
        fq = [next(sam_iter)]
        for record in sam_iter:
            if record.qname == fq[0].qname:
                fq.append(record)
            else:
                yield tuple(fq)
                fq = [record]
        yield record
