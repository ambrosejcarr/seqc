from multiprocessing import Pool
from contextlib import closing
from functools import partial
import numpy as np
from seqc import log
from seqc.alignment import sam
from seqc.sequence.gtf import GeneIntervals
from seqc.sequence.encodings import DNA3Bit
from seqc.sparse_frame import SparseFrame
from scipy.sparse import coo_matrix, csr_matrix, vstack
import seqc.multialignment as mm
import seqc.sequence.barcodes
from itertools import permutations
import time
import tables as tb
import pysam

# Some constants
PASS_FILTERS = 0
NO_ALIGNMENT = 1
GENE_0 = 2
NO_CELL = 3
NO_RMT = 4
POLY_T = 5
SCAFFOLD = 6
NO_GENES = 7
DUST_SCORE = 8

BC_OK = 0
BC_FIXED = 1
BC_BAD = 2

DEFAULT_BASE_CONVERTION_RATE = 0.02


class ReadArray:

    _dtype = [
        ('status', np.uint16),
        ('cell', np.int64),
        ('rmt', np.int32),
        ('n_poly_t', np.uint8)]

    def __init__(self, data, genes, positions, non_unique_alignments=-1):
        """
        Enhanced np.ndarray (structured array) with several companion functions to hold
        filter, sort, and access compressed fastq read information.

        :param np.recarray data: np.ndarray with dtype of self._dtype
        :param scipy.sparse.coo_matrix genes:
        :param scipy.sparse.coo_matrix positions:

        """

        # set genes and positions to csr matrices.
        if not isinstance(genes, csr_matrix):
            try:
                genes = genes.tocsr()
            except AttributeError:
                raise TypeError(
                    'genes must be a scipy sparse matrix that is convertible to '
                    'csr_matrix, not %s' % repr(type(genes)))
        self._genes = genes
        if not isinstance(positions, csr_matrix):
            try:
                positions = positions.tocsr()
            except AttributeError:
                raise TypeError(
                    'positions must be a scipy sparse matrix that is convertible to '
                    'csr_matrix, not %s' % repr(type(positions)))
        self._positions = positions

        self._data = data

        # statistical values that are set by filtering functions
        self._total_alignments = 0
        self._total_reads = 0
        self._non_unique_align = 0
        self._filtered = {
            NO_ALIGNMENT: 0,
            GENE_0: 0,
            NO_CELL: 0,
            NO_RMT: 0,
            POLY_T: 0,
            SCAFFOLD: 0,
            NO_GENES: 0,
            DUST_SCORE: 0}
        self._bad_barcodes = 0
        self._fix_barcodes = 0
        self._ok_barcodes = 0
        self._multialignment = {
            mm.NO_DISAMBIGUATION: 0,
            mm.RESOLVED_GENE: 0,
            mm.NO_GENE_RESOLVED: 0,
            mm.MULTIPLE_MODELS: 0}
        self._total_errors = 0
        
        if non_unique_alignments != -1:
            self._non_unique_align = non_unique_alignments

    def __len__(self):
        return len(self.data)

    def __iter__(self):
        for i in np.arange(self.data.shape[0]):
            return self.data[i], self.genes[i], self.positions[i]

    @property
    def data(self):
        return self._data

    @property
    def genes(self):
        return self._genes

    @property
    def positions(self):
        return self._positions

    # @staticmethod
    # def set_status_active(r):
    #     """Toggles active flag on (byte 15)"""
    #     r['status'] |= 0b1  # todo bug? this isn't byte 15, it's bit 0 ??

    @staticmethod
    def set_status_inactive(r):
        r['status'] &= 0b1111111111111110

    # @staticmethod
    # def set_filt_fail(r, filt):
    #     if filt > 0b1111:
    #         raise ValueError('Filter value must be 4 bits. value given: {}'.format(filt))
    #     r['status'] &= 0b1111111111100001
    #     r['status'] |= filt << 1
        
    @staticmethod
    def set_barcodes(r, bc_status):
        if bc_status > 0b11:
            raise ValueError('Barcodes status must be 0,1 or 2. value given: {}'.format(
                bc_status))
        r['status'] &= 0b1111111110011111
        r['status'] |= bc_status << 5
    
    # @staticmethod
    # def set_mm_status(r, mm_stat):
    #     if mm_stat > 0b11:
    #         raise ValueError(
    #             'Multimapping status must be 0,1,2 or 3. value given: {}'.format(mm_stat))
    #     r['status'] &= 0b1111110011111111
    #     r['status'] |= mm_stat << 8
        
    @staticmethod
    def set_error_status_on(r):
        r['status'] |= 0b0000000010000000
    
    # @staticmethod
    # def set_error_status_off(r):
    #     r['status'] &= 0b1111111101111111
    
    @staticmethod
    def is_active_read(r):
        return (r['status'] & 0b1) == 0b1
    
    # @staticmethod
    # def is_error(r):
    #     return (r['status'] & 0b0000000010000000) == 0b0000000010000000
        
    # def count_active_reads(self):
    #     res = 0
    #     for r in self.data:
    #         if ReadArray.is_active_read(r):
    #             res += 1
    #     return res
    
    # def count_errors(self):
    #     res = 0
    #     for r in self.data:
    #         if ReadArray.is_error(r):
    #             res += 1
    #     return res

    # @staticmethod
    # def _iterate_pysam(reference, open_alignment_file):
    #     """group and return multialignments for each read associated with reference in
    #     bam_file
    #
    #     :param str reference: the reference within the file to parse
    #     :param pysam.AlignmentFile open_alignment_file: an open file descriptor pointing
    #       to the bam file that should be parsed
    #     :yields list: list of all alignments for each read
    #     """
    #     iter_bam = open_alignment_file.fetch(reference=reference, multiple_iterators=True)
    #     ma = [next(iter_bam)]  # get first alignment
    #     for alignment in iter_bam:
    #         if alignment.qname == ma[0].qname:
    #             ma.append(alignment)
    #         else:
    #             yield ma
    #             ma = [alignment]
    #     yield ma  # yield final alignment
    #
    # @classmethod
    # def _from_pysam(cls, reference_number, bam_file, translator):
    #     """
    #
    #     :param reference: the reference within the file to parse
    #     :param bam_file: the alignment file to iterate over
    #     :param translator: a data structure that maps alignments in chromosome space to
    #       genes
    #     :return np.recarray, csr_matrix, csr_matrix: data, genes, and positions input for
    #       ReadArray construction
    #     """
    #
    #     reader = pysam.AlignmentFile(bam_file, 'rb')
    #     reference = reader.references[reference_number]
    #     length = reader.lengths[reference_number]
    #
    #     # pre-allocate arrays
    #     data = np.recarray((length,), cls._dtype)
    #     row = np.zeros(length, dtype=np.int32)
    #     col = np.zeros(length, dtype=np.int8)
    #     position = np.zeros(length, dtype=np.int32)
    #     gene = np.zeros(length, dtype=np.int32)
    #
    #     row_idx = 0  # identifies the read index
    #     arr_idx = 0  # identifies the alignment index across all reads
    #     max_ma = 0  # maximum number of alignments observed for a read
    #     for multialignment in cls._iterate_pysam(reference, reader):
    #
    #         col_idx = 0  # identifies the alignment index within a read
    #         for a in multialignment:
    #             strand = '-' if a.is_reverse else '+'
    #             genes = translator.translate(a.reference_name, strand, a.pos)
    #
    #             # if gene passes filter, store passing genes and increment appropriate
    #             # arr and col indices
    #             if genes is not None and len(genes) == 1:
    #                 row[arr_idx] = row_idx
    #                 col[arr_idx] = col_idx
    #                 gene[arr_idx] = genes[0]
    #                 position[arr_idx] = a.pos
    #                 arr_idx += 1
    #                 col_idx += 1
    #         max_ma = max(max_ma, col_idx)
    #
    #         _, cell, rmt, poly_t = a.qname[:a.qname.index(';')].split(':')
    #         cell = seqc.sequence.encodings.DNA3Bit.encode(cell)
    #         rmt = seqc.sequence.encodings.DNA3Bit.encode(rmt)
    #         n_poly_t = poly_t.count('T') + poly_t.count('N')
    #
    #         data[row_idx] = (1, cell, rmt, n_poly_t)
    #         row_idx += 1
    #
    #     # some reads will not have aligned, throw away excess allocated space before
    #     # creating the ReadArray
    #     row, col, position, gene = (
    #         row[:arr_idx], col[:arr_idx], position[:arr_idx], gene[:arr_idx])
    #
    #     gene = coo_matrix((gene, (row, col)), shape=(row_idx, max_ma), dtype=np.int32).tocsr()
    #     position = coo_matrix((position, (row, col)), shape=(row_idx, max_ma),
    #                           dtype=np.int32).tocsr()
    #
    #     return data, gene, position
    #
    # @classmethod
    # def from_pysam(cls, bam_file, translator):
    #
    #     reader = pysam.AlignmentFile(bam_file, 'rb')
    #
    #     if not reader.has_index():  # todo instead of raising, sort & index bam file
    #         raise RuntimeError('Bam file must be coordinate sorted and indexed')
    #
    #     reference_numbers = np.arange(len(reader.references))
    #
    #     parser = partial(cls._from_pysam, bam_file=bam_file, translator=translator)
    #
    #     with closing(Pool()) as pool:
    #         data, gene, position = zip(*pool.map(parser, reference_numbers))
    #
    #     gene = vstack(gene)
    #     position = vstack(position)
    #     data = np.concatenate(data, axis=0)
    #     cls(data, gene, position)

    @classmethod
    def from_samfile(cls, sam_file: str, gtf_file: str):
        """
        construct a ReadArray object from a samfile containing only uniquely aligned
        records

        :param gtf_file: str, filename of annotations.gtf file
        :param sam_file: str, filename of alignment file.
        :return:
        """

        # todo replace most of this function with the pysam API, and read the coordinate
        # sorted bam file with index by chromosome
        # farm out individual chromosomes to a processor pool
        # if scaffolds are not supported, filter out scaffolds.
        # ignore dust score, it takes 66% of time in normal processing
        # this should speed up RA creation to ~ 1h or less, from 17h on 400M read file.

        reader = sam.Reader(sam_file)

        # todo if this takes a long time, we can read this from the alignment summary!
        # go over the file to find how many unique reads there are
        num_reads = 0
        num_unique = 0
        prev_alignment_name = ''
        for alignment in reader:
            num_reads += 1
            if alignment.qname != prev_alignment_name:
                num_unique += 1
                prev_alignment_name = alignment.qname       
        
        translator = GeneIntervals(gtf_file)

        # pre-allocate arrays
        data = np.recarray((num_unique,), cls._dtype)
        row = np.zeros(num_reads, dtype=np.int32)
        col = np.zeros(num_reads, dtype=np.int8)
        position = np.zeros(num_reads, dtype=np.int32)
        gene = np.zeros(num_reads, dtype=np.int32)

        # loop over multialignments
        row_idx = 0  # identifies the read index
        arr_idx = 0  # identifies the alignment index across all reads
        max_ma = 0  # maximum number of alignments observed for a read
        for ma in reader.iter_multialignments():
            col_idx = 0  # identifies the alignment index within a read
            for a in ma:
                genes = translator.translate(a.rname, a.strand, a.pos)

                # if gene passes filter, store passing genes and increment appropriate
                # arr and col indices
                if genes is not None and len(genes) == 1:
                    row[arr_idx] = row_idx
                    col[arr_idx] = col_idx
                    gene[arr_idx] = genes[0]
                    position[arr_idx] = a.pos
                    arr_idx += 1
                    col_idx += 1
            max_ma = max(max_ma, col_idx)

            cell = seqc.sequence.encodings.DNA3Bit.encode(alignment.cell)
            rmt = seqc.sequence.encodings.DNA3Bit.encode(alignment.rmt)
            n_poly_t = alignment.poly_t.count(b'T') + alignment.poly_t.count(b'N')
            # dust_score = alignment.dust_low_complexity_score
            data[row_idx] = (0, cell, rmt, n_poly_t)
            row_idx += 1

        # some reads will not have aligned, throw away excess allocated space before
        # creating the ReadArray
        row, col, position, gene = (
            row[:arr_idx], col[:arr_idx], position[:arr_idx], gene[:arr_idx])

        gene = coo_matrix((gene, (row, col)), shape=(row_idx, max_ma), dtype=np.int32)
        position = coo_matrix((position, (row, col)), shape=(row_idx, max_ma),
                              dtype=np.int32)

        log.info('ReadArray constructed')
        return cls(data, gene, position)
        # res_ra._total_alignments = len(reader)  # todo I removed these.
        # res_ra._total_reads = num_unique

    def apply_filters(self, required_poly_t=1, max_dust_score=10):
        """Apply different filters to the read array. If a read fails a filter, it is not
        passed to the others and so the counts for each filter is the number of reads that
        failed that one but passed all previous ones Filters are not ordered in any
        particular way.

        :param required_poly_t: minimum number of T nucleotides in the primer tail of a
          valid read
        :param max_dust_score: maximum complexity score of a valid read, where higher
          implies lower complexity
        :return None: method sets the status vector according to the result of filtering.
        """

        failing = np.zeros(self.data.shape[0])
        failing[self.genes.sum(axis=1) == 0] &= 0b1  # no gene
        failing[self.data['rmt'] == 0] &= 0b10  # no rmt
        failing[self.data['cell'] == 0] &= 0b100 # no cell
        failing[self.data['n_poly_t'] < required_poly_t] &= 0b1000
        # failing[self.data['dust_score'] > max_dust_score] &= 0b10000

        self.data['status'] &= failing

    def _disambiguate(self, inds):
        """side-effecting thread-safe function that sets the reads within inds to the
        gene that is identified as the unambiguous source of each read.

        Two-stage function that:
        (1) separates reads into disjoint sets
        (2) takes the intersection of all gene lists in each disjoint set and, if the
            intersection is a single gene, sets all reads equal to that gene alone.

        :param inds: indices localizing all the reads of a molecule, designated as a
          unique (cell, rmt) pair.
        :return None:
        """

        # expects d[(cell, rmt)][(gene(s))] = [inds]

        reads = self._genes[slice(*inds)]


        total_reads = 0
        g_ra = {}  # I added this

        for g_ra_rec in g_ra.values():  # loop over molecules; this is in the outer func
            obs = {}
            for genes in g_ra_rec:  # loop over genes
                if genes == ():
                    continue
                if genes not in obs:
                    obs[genes] = len(g_ra_rec[genes])
                    total_reads += len(g_ra_rec[genes])
                    # TODO: do a sanity check here? I think genes should not be in obs

            if () in g_ra_rec:
                if len(g_ra_rec) == 1:
                    continue
            if () in obs:
                del obs[()]

            obs_s = mm.strip_model(obs)
            ind_s = mm.strip_model(g_ra_rec)

            subsets = mm.split_to_disjoint(obs_s)
            for obs_subset in subsets:
                resolved = False
                models, r = mm.best_fit_model(obs_subset, mat)
                self._multialignment[r] += sum(obs_subset.values())
                if r == mm.NO_DISAMBIGUATION or r == mm.RESOLVED_GENE:
                    model = models[0]
                    g = mm.model_to_gene(model)
                    resolved = True

                indices = mm.get_indices(ind_s, obs_subset)
                for ind in indices:
                    ReadArray.set_mm_status(self.data[ind], r)
                    if resolved:
                        self.data['gene'][ind] = g  # api (indexing gene field)
                        self.data['position'][ind] = (  # api (indexing position field)
                            self.data['ma_pos'][ind][
                                self.data['ma_genes'][ind].index(g)])  # api (indexing)
                    else:
                        ReadArray.set_status_inactive(self.data[ind])

        log.info('Total reads in samfile: {}'.format(len(self)))
        log.info('Reads that passed the filter: {}'.format(total_reads))
        tot_time = time.process_time() - start
        log.info('Multimapping resolution completed in {} seconds.'.format(tot_time))

    def disambiguation(self):
        """

        :return:
        """

        start_time = time.time()

        # todo could prefilter, but would need to copy all data in self. ~ 6g / 200M reads
        # sort data according to cell, rmt
        idx = np.argsort(self.data, order=('cell', 'rmt'))
        self._data = self.data[idx]
        self._genes = self.genes[idx]
        self._positions = self.positions[idx]

        # identify points separating molecules
        molecule_breaks = np.concatenate([
            np.array([0]),  # start point
            np.where(np.diff(self.data['rmt']))[0] + 1,  # molecule break points
            self.data.shape[0]  # end point
        ])

        # zip together two views of the created list for passing to multiprocessing
        molecule_inds = zip(molecule_breaks[:-1], molecule_breaks[1:])

        with closing(Pool()) as pool:
            pool.map(self._disambiguate, molecule_inds)

    def save(self, archive_name: str) -> None:
        """save a ReadArray object as an hdf5 archive

        :param archive_name: filestem for the new archive
        :return: None
        """

        def store_carray(archive, array, name):
            atom = tb.Atom.from_dtype(array.dtype)
            store = archive.create_carray(archive.root, name, atom, array.shape)
            store[:] = array
            store.flush()

        if not archive_name.endswith('.h5'):
            archive_name += '.h5'

        # construct container
        blosc5 = tb.Filters(complevel=5, complib='blosc')
        f = tb.open_file(archive_name, mode='w', title='Data for seqc.ReadArray',
                         filters=blosc5)

        f.create_table(f.root, 'data', self.data)

        # data, indices, indptr
        store_carray(f, self.genes.indices, 'indices')
        store_carray(f, self.genes.indptr, 'indptr')
        store_carray(f, self.genes.data, 'gene_data')
        store_carray(f, self.positions.data, 'positions_data')

        f.close()

    @classmethod
    def load(cls, archive_name: str):
        """load a ReadArray from an hdf5 archive, note that ma_pos and ma_genes are
        discarded.

        :param archive_name: name of a .h5 archive containing a saved ReadArray object
        :return: ReadArray
        """

        f = tb.open_file(archive_name, mode='r')
        data = f.root.data.read()
        indptr = f.root.indptr.read()
        indices = f.root.indices.read()
        genes = f.root.gene_data.read()
        positions = f.root.positions_data.read()
        f.close()

        # todo may need to save and load shape if not properly inferred.
        genes_csr = csr_matrix((genes, indices, indptr))
        positions_csr = csr_matrix((positions, indices, indptr))

        return cls(data, genes_csr, positions_csr)

    def apply_barcode_correction(self, platform, barcode_files, max_ed=2):
        """
        Correct reads with incorrect barcodes according to the correct barcodes files.
        Reads with barcodes that have too many errors are filtered out.
        """

        # todo once the sorted list of barcodes is created, we should be able to
        # parallelize because the error correction requires only reading, not writing

        # Read the barcodes into lists
        correct_barcodes = []
        for barcode_file in barcode_files:
            with open(barcode_file, 'r') as f:
                correct_barcodes.append(
                    set(DNA3Bit.encode(line.strip()) for line in f.readlines()))

        num_barcodes = platform.num_barcodes
        dynamic_codes_table = [{}] * num_barcodes
        correct = [None] * num_barcodes
        edit_dist = [None] * num_barcodes
        err_list = [None] * num_barcodes
        sequencer_errors_cnt = [None] * num_barcodes

        errors = list(DNA3Bit.ints2int([p[0], p[1]]) for p in
                      permutations(DNA3Bit.bin2strdict.keys(), r=2))
        error_table = dict(zip(errors, [0] * len(errors)))
        cor_instance_table = {DNA3Bit.encode(b'A'): 0,
                              DNA3Bit.encode(b'C'): 0,
                              DNA3Bit.encode(b'G'): 0,
                              DNA3Bit.encode(b'T'): 0,
                              DNA3Bit.encode(b'N'): 0}

        for r in self.data:
            if not ReadArray.is_active_read(r):
                continue
            barcodes = platform.extract_barcodes(int(r['cell']))
            for i in range(num_barcodes):
                try:
                    correct[i], edit_dist[i], err_list[i], sequencer_errors_cnt[i] = (
                        dynamic_codes_table[i][barcodes[i]])
                except KeyError:
                    correct[i], edit_dist[i] = (
                        seqc.sequence.barcodes.find_correct_barcode(
                            barcodes[i], correct_barcodes[i]))
                    sequencer_errors_cnt[i] = DNA3Bit.count(
                        barcodes[i], DNA3Bit.encode(b'N'))

                    # No point calculating errors on codes with more errors than allowed
                    if edit_dist[i] > max_ed:
                        err_list[i] = None
                    else:
                        err_list[i] = seqc.sequence.barcodes.list_errors(
                            correct[i], barcodes[i])
                    dynamic_codes_table[i][barcodes[i]] = (
                        correct[i], edit_dist[i], err_list[i], sequencer_errors_cnt[i])

            # Filter reads with too many errors in their barcodes
            if max(edit_dist) > max_ed:
                self._bad_barcodes += 1
                ReadArray.set_barcodes(r, BC_BAD)
                ReadArray.set_status_inactive(r)
                continue

            # Count the different error types
            # For error rate estimation, only reads with at most a single error in both
            # barcodes are counted. note that only base swicth errors are counted towards
            # this threshold, sequencer errors (that cause an 'N') are ignored for this.
            if sum(edit_dist) - sum(sequencer_errors_cnt) <= 1:
                for err_l in err_list:
                    # todo bug here, if err_l is None (see above) this will throw a
                    # TypeError because NoneType is not iterable. Since this does not
                    # happen in practice, I think the above section is not necessary.
                    for err in err_l:
                        try:
                            error_table[err] += 1
                        except KeyError:  # happens for sequencer errors
                            continue

            # count non error bases
            tmp_bc = DNA3Bit.ints2int(barcodes)
            tmp_cor = DNA3Bit.ints2int(correct)
            while tmp_bc > 0:
                if tmp_bc & 0b111 == tmp_cor & 0b111:
                    cor_instance_table[tmp_bc & 0b111] += 1
                tmp_bc >>= 3
                tmp_cor >>= 3

            # correct cell and update read status
            if tmp_bc != tmp_cor:
                r['cell'] = tmp_cor
                ReadArray.set_barcodes(r, BC_FIXED)
                self._fix_barcodes += 1
            else:
                ReadArray.set_barcodes(r, BC_OK)
                self._ok_barcodes += 1

        # Create error rate table
        default_error_rate = DEFAULT_BASE_CONVERTION_RATE
        if sum(error_table.values()) == 0:
            log.info('No errors were detected, using %f uniform error chance.'
                     % default_error_rate)
            # todo bug here <-- this value is never used
            err_rate = dict(zip(errors, [default_error_rate] * len(errors)))
        err_rate = dict(zip(errors, [0.0] * len(errors)))
        for k, v in error_table.items():
            if DNA3Bit.decode(k)[0] in b'Nn':
                continue
            try:
                err_rate[k] = v / (sum(n for err_type, n in error_table.items()
                                       if err_type & 0b111000 == k & 0b111000) +
                                   cor_instance_table[(k & 0b111000) >> 3])
            except ZeroDivisionError:
                log.info('Warning: too few reads to estimate error rate for %s, setting '
                         'default rate of %f'
                         % (str(DNA3Bit.decode(k)), default_error_rate))

                err_rate[k] = default_error_rate

        log.info('Barcodes correction results: Good barcodes: {}, fixed barcodes: {}, '
                 'bad barcodes: {}'.format(
            self._ok_barcodes, self._fix_barcodes, self._bad_barcodes))
        return err_rate

    def to_count_matrix(self, csv_path=None, sparse_frame=False, genes_to_symbols=False):
        """Convert the ReadArray into a count matrix.

        Since the matrix is sparse we represent it with 3 columns: row (cell),
         col (gene), value if not 0 (reads/molecules)

        :param str csv_path: file prefix for .csv output
        :param bool sparse_frame: if True, return sparse frame objects
        :param genes_to_symbols: if SparseFrame is True, integer gene IDs are converted
          to symbols in the returned SparseFrame objects
        :return str, str, dict, dict:
          filename (read counts),
          filename (mol counts),
          dict of (cell, rmt) -> read counts
          dict of (cell, rmt) -> molecule counts
        """
        reads_mat = {}
        mols_mat = {}
        for r in self.data:
            if not ReadArray.is_active_read(r):
                continue
            cell = r['cell']
            gene = r['gene']
            rmt = r['rmt']
            if gene == 0:
                # todo bug here (ma_genes is a list, these will never be true)
                if r['ma_genes'] == -1 or r['ma_genes'] == 0:  # api (comparison to row value of ma_genes)
                    continue
                elif len(r['ma_genes']) == 1:  # api (length)
                    gene = r['ma_genes'][0]  # api (indexing)
                else:
                    continue
            try:
                reads_mat[cell, gene] += 1
            except KeyError:
                reads_mat[cell, gene] = 1
                                
            try:
                if rmt not in mols_mat[cell, gene]:
                    mols_mat[cell, gene].append(rmt)
            except KeyError:
                mols_mat[cell, gene] = [rmt]

        if sparse_frame:
            return (SparseFrame.from_dict(reads_mat, genes_to_symbols=genes_to_symbols),
                    SparseFrame.from_dict(
                        {k: len(v) for k, v in mols_mat.items()},
                        genes_to_symbols=genes_to_symbols))

        if csv_path is None:
            return reads_mat, mols_mat

        f = open(csv_path+'reads_count.csv', 'w')
        for cell, gene in reads_mat:
            f.write('{},{},{}\n'.format(cell, gene, reads_mat[cell, gene]))
        f.close()
        f = open(csv_path+'mols_count.csv', 'w')
        for cell, gene in mols_mat:
            f.write('{},{},{}\n'.format(cell, gene, len(mols_mat[cell, gene])))
        f.close()

        return csv_path + 'reads_count.csv', csv_path + 'mols_count.csv',
