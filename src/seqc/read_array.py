import numpy as np
from seqc import log
from seqc.alignment import sam
from seqc.sequence.gtf import GeneIntervals
from seqc.sequence.encodings import DNA3Bit
from seqc.sparse_frame import SparseFrame
import seqc.multialignment as mm
import seqc.sequence.barcodes
import pickle
from itertools import permutations
import time
# import tables as tb

# Some constants
PASS_FILTERS    = 0
NO_ALIGNMENT    = 1
GENE_0          = 2
NO_CELL         = 3
NO_RMT          = 4
POLY_T          = 5
SCAFFOLD        = 6
NO_GENES        = 7
DUST_SCORE      = 8

BC_OK           = 0
BC_FIXED        = 1
BC_BAD          = 2

DEFAULT_BASE_CONVERTION_RATE = 0.02
    
class ReadArray:

    _dtype = [
        ('status', np.uint16),
        ('cell', np.int64),
        ('rmt', np.int32),
        ('n_poly_t', np.uint8),
        ('dust_score', np.uint8),
        ('gene', np.uint32),
        ('position', np.uint64),
        ('ma_genes', np.dtype(object)),
        ('ma_pos', np.dtype(object))]

    def __init__(self, data, non_unique_alignments=-1):
        """
        Enhanced np.ndarray (structured array) with several companion functions to hold
        filter, sort, and access compressed fastq read information.

        :param data: np.ndarray with dtype of self._dtype

        :property data: stored np.ndarray
        :method save: saves the ReadArray in compressed .h5 format
        :method load: load a saved compressed .h5 representation of a ReadArray
        :method from_samfile: constructs a ReadArray object from a samfile (uniquely
          aligned records only)
        :method reads_passing_filters: Return a ReadArray containing only reads
          that pass all filters (n_poly_t, dust_complexity_score, presence of cell and
          rmt)
        """
        self._data = data
        
        #some stats
        self._total_alignments  = 0
        self._total_reads       = 0
        self._non_unique_align  = 0
        self._filtered          = {NO_ALIGNMENT: 0, GENE_0: 0, NO_CELL: 0, NO_RMT: 0, POLY_T: 0, SCAFFOLD: 0, NO_GENES: 0, DUST_SCORE: 0}
        self._bad_barcodes      = 0
        self._fix_barcodes      = 0
        self._ok_barcodes       = 0
        self._multialignment    = {mm.NO_DISAMBIGUATION: 0, mm.RESOLVED_GENE: 0, mm.NO_GENE_RESOLVED: 0, mm.MULTIPLE_MODELS: 0}
        self._total_errors      = 0
        
        if non_unique_alignments != -1:
            self._non_unique_align = non_unique_alignments

    def __repr__(self):
        s = ""
        s += "Read array of size: {:,}\n".format(len(self._data))
        s += "Non unique alignments (mapped to more than one gene): {:,}\n".format(self._non_unique_align)
        s += "Total reads: {:,}\n".format(self._total_reads)
        s += "Filtered: No alignment:{:,}, gene is 0:{:,}, no cell: {:,}, no rmt: {:,}, poly_t:{:,}, scaffold: {:,}, no genes at all: {:,}\n".\
             format(self._filtered[NO_ALIGNMENT], self._filtered[GENE_0],self._filtered[NO_CELL], self._filtered[NO_RMT], self._filtered[POLY_T], self._filtered[SCAFFOLD], self._filtered[NO_GENES])
        s += "Barcodes: good:{:,}, fixed: {:,}, bad:{:,}\n".format(self._ok_barcodes, self._fix_barcodes, self._bad_barcodes)
        s += "Multialignment results: no disambiguation: {:,}, resolved gene: {:,}, no gene resolved: {:,}, competing models (unhandeled): {:,}\n".\
             format(self._multialignment[mm.NO_DISAMBIGUATION], self._multialignment[mm.RESOLVED_GENE], self._multialignment[mm.NO_GENE_RESOLVED], self._multialignment[mm.MULTIPLE_MODELS])
        s += "Total errors: {:,}\n".format(self._total_errors)
        s += "Total active reads: {:,}\n".format(self.count_active_reads())
        return s
    
    @property
    def data(self):
        return self._data

    @staticmethod
    def active_status():
        return 0
    # Toggle active flag on (byte 15)
    @staticmethod
    def set_status_active(r):
        r['status'] |= 0b1
    @staticmethod
    def set_status_inactive(r):
        r['status'] &= 0b1111111111111110

    @staticmethod
    def set_filt_fail(r, filt):
        if filt > 0b1111:
            raise ValueError('Filter value must be 4 bits. value given: {}'.format(filt))
        r['status'] &= 0b1111111111100001
        r['status'] |= filt<<1
        
    @staticmethod
    def set_barcodes(r, bc_status):
        if bc_status > 0b11:
            raise ValueError('Barcodes status must be 0,1 or 2. value given: {}'.format(bc_status))
        r['status'] &= 0b1111111110011111
        r['status'] |= bc_status<<5
    
    @staticmethod
    def set_mm_status(r, mm_stat):
        if mm_stat > 0b11:
            raise ValueError('Multimapping status must be 0,1,2 or 3. value given: {}'.format(mm_stat))
        r['status'] &= 0b1111110011111111
        r['status'] |= mm_stat<<8
        
    @staticmethod
    def set_error_status_on(r):
        r['status'] |= 0b0000000010000000
    
    @staticmethod
    def set_error_status_off(r):
        r['status'] &= 0b1111111101111111
    
    @staticmethod
    def is_active_read(r):
        return (r['status']&0b1)==0b1
    
    @staticmethod
    def is_error(r):
        return (r['status']&0b0000000010000000)==0b0000000010000000
        
    def count_active_reads(self):
        res = 0
        for r in self.data:
            if ReadArray.is_active_read(r):
                res+=1
        return res
    
    def count_errors(self):
        res = 0
        for r in self.data:
            if ReadArray.is_error(r):
                res += 1
        return res
        
    @classmethod
    def from_samfile(cls, samfile: str, gtf_f: str):
        """
        construct a ReadArray object from a samfile containing only uniquely aligned
        records

        :param gtf_f: str, filename of annotations.gtf file
        :param samfile: str, filename of alignment file.
        :return:
        """
        non_unique_align = 0
        reader = sam.Reader(samfile)
        #go over the file to find how many unique reads there are
        num_reads=0
        prev_alignment_name=''
        for alignment in reader:
            if alignment.qname != prev_alignment_name:
                num_reads +=1
                prev_alignment_name = alignment.qname       
        
        translator = GeneIntervals(gtf_f)
        
        data = np.recarray((num_reads,), cls._dtype)
        prev_alignment_name=''
        read_idx = 0
        for alignment in reader:    
            #new read
            if alignment.qname != prev_alignment_name:
                if read_idx>0:
                    data[read_idx-1] = (ReadArray.active_status(), cell, rmt, n_poly_t, dust_score, 0, 0, genes, ma_pos)
                    cls.set_status_active(data[read_idx-1])
                prev_alignment_name = alignment.qname
                read_idx+=1
                genes = translator.translate(alignment.rname, alignment.strand, alignment.pos)
                ma_pos = []
                if genes != -1:
                    if len(genes) > 1:
                        non_unique_align+=1
                    for i in range(len(genes)):
                        ma_pos.append(alignment.pos)
                cell = seqc.sequence.encodings.DNA3Bit.encode(alignment.cell)
                rmt = seqc.sequence.encodings.DNA3Bit.encode(alignment.rmt)
                n_poly_t = alignment.poly_t.count(b'T') + alignment.poly_t.count(b'N')
                dust_score = alignment.dust_low_complexity_score
                
            #the same read
            else:
                genes_l = translator.translate(alignment.rname, alignment.strand, alignment.pos)
                if genes_l != -1:
                    if len(genes_l) > 1:
                        non_unique_align+=1
                    if genes == -1:
                        genes = genes_l
                    else:
                        genes+=genes_l
                    for i in range(len(genes_l)):
                        ma_pos.append(alignment.pos)
            
        log.info('non unique alignments count: {}'.format(non_unique_align))
        res_ra = cls(data, non_unique_align)
        res_ra._total_alignments = len(reader)
        res_ra._total_reads = num_reads

        return res_ra, res_ra._total_alignments #TODO: It might be nicer to return just the RA and take the total alignments from it

    def apply_filters(self, required_poly_t = 1, max_dust_score=10):
        """
        Apply different filters to the read array.
        If a read fails a filter, it is not passed to the others and so the counts
        for each filter is the number of reads that failed that one but passed all previous ones
        Filters are not ordered in any particular way.
        """
        for r in self.data:
            if r['ma_genes'] == None:
                ReadArray.set_filt_fail(r, NO_GENES)
                self._filtered[NO_GENES] += 1
                ReadArray.set_status_inactive(r)
                continue
            if r['ma_genes'] == -1:
                ReadArray.set_filt_fail(r, SCAFFOLD)
                self._filtered[SCAFFOLD] += 1
                ReadArray.set_status_inactive(r)
                continue
            if len(r['ma_genes']) == 0:
                ReadArray.set_filt_fail(r, NO_ALIGNMENT)
                self._filtered[NO_ALIGNMENT] += 1
                ReadArray.set_status_inactive(r)
                continue
            if r['cell'] == 0:
                ReadArray.set_filt_fail(r, NO_CELL)
                self._filtered[NO_CELL] += 1
                ReadArray.set_status_inactive(r)
                continue
            if r['rmt'] == 0:
                ReadArray.set_filt_fail(r, NO_RMT)
                self._filtered[NO_RMT] += 1
                ReadArray.set_status_inactive(r)
                continue
            if r['n_poly_t'] < required_poly_t:
                ReadArray.set_filt_fail(r, POLY_T)
                self._filtered[POLY_T] += 1
                ReadArray.set_status_inactive(r)
                continue
            if 0 in r['ma_genes']:
                ReadArray.set_filt_fail(r, GENE_0)
                self._filtered[GENE_0] += 1
                ReadArray.set_status_inactive(r)
                continue
            if r['dust_score'] > max_dust_score:
                ReadArray.set_filt_fail(r, DUST_SCORE)
                self._filtered[DUST_SCORE] += 1
                ReadArray.set_status_inactive(r)
                continue
                
        log.info('Filters results: No alignment: {}, scaffold: {}, no cell {}, no rmt: {}, low poly_t: {}, mapped to gene 0: {}, genes is None: {}'.
        format(self._filtered[NO_ALIGNMENT], self._filtered[SCAFFOLD], self._filtered[NO_CELL], self._filtered[NO_RMT], self._filtered[POLY_T], self._filtered[GENE_0], self._filtered[NO_GENES]))
        
            
    def apply_barcode_correction(self, platform, barcode_files, max_ed=2):
        """
        Correct reads with incorrect barcodes according to the correct barcodes files.
        Reads with barcodes that have too many errors are filtered out.
        """
        
        # Read the barcodes into lists
        correct_barcodes = []
        for barcode_file in barcode_files:
            with open(barcode_file, 'r') as f:
                correct_barcodes.append(set(DNA3Bit.encode(line.strip()) for line in f.readlines()))
        
        num_barcodes = platform.num_barcodes
        dynamic_codes_table     = [{}]*num_barcodes
        correct                 = [None]*num_barcodes
        edit_dist               = [None]*num_barcodes
        err_list                = [None]*num_barcodes
        sequencer_errors_cnt    = [None]*num_barcodes
        
        errors = list(DNA3Bit.ints2int([p[0], p[1]]) for p in permutations(DNA3Bit._bin2strdict.keys(), r=2))
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
                    correct[i], edit_dist[i], err_list[i], sequencer_errors_cnt[i] = dynamic_codes_table[i][barcodes[i]]
                except KeyError:
                    correct[i], edit_dist[i] = seqc.sequence.barcodes.find_correct_barcode(barcodes[i], correct_barcodes[i])
                    sequencer_errors_cnt[i] = DNA3Bit.count(barcodes[i], DNA3Bit.encode(b'N'))
                    
                    # No point calculating errors on codes with more errors than allowed
                    if edit_dist[i] > max_ed:
                        err_list[i] = None
                    else:
                        err_list[i] = seqc.sequence.barcodes.list_errors(correct[i], barcodes[i])
                    dynamic_codes_table[i][barcodes[i]] = (correct[i], edit_dist[i], err_list[i], sequencer_errors_cnt[i])
                    

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
                    for err in err_l:
                        try:
                            error_table[err] += 1
                        except KeyError:    # happens for sequencer errors
                            continue
            
            # count non error bases
            tmp_bc = DNA3Bit.ints2int(barcodes)
            tmp_cor = DNA3Bit.ints2int(correct)
            while tmp_bc > 0:
                if tmp_bc & 0b111 == tmp_cor & 0b111:
                    cor_instance_table[tmp_bc & 0b111] += 1
                tmp_bc >>= 3
                tmp_cor >>= 3
            
            #correct cell and update read status
            if tmp_bc!=tmp_cor:
                r['cell'] = tmp_cor
                ReadArray.set_barcodes(r, BC_FIXED)
                self._fix_barcodes += 1
            else:
                ReadArray.set_barcodes(r, BC_OK)
                self._ok_barcodes += 1
                
        # Create error rate table
        default_error_rate = DEFAULT_BASE_CONVERTION_RATE
        if sum(error_table.values()) == 0:
            log.info('No errors were detected, using %f uniform error chance.' % (default_error_rate))
            err_rate = dict(zip(errors, [default_error_rate] * len(errors)))
        err_rate = dict(zip(errors, [0.0] * len(errors)))
        for k, v in error_table.items():
            if DNA3Bit.decode(k)[0] in b'Nn':
                continue
            try:
                err_rate[k] = v / (sum(n for err_type, n in error_table.items() if err_type&0b111000 == k&0b111000) + cor_instance_table[(k & 0b111000) >> 3])
            except ZeroDivisionError:
                log.info('Warning: too few reads to estimate error rate for %s, setting default rate of %f' % (str(DNA3Bit.decode(k)), default_error_rate))
                err_rate[k] = default_error_rate
            
        log.info('Barcodes correction results: Good barcodes: {}, fixed barcodes: {}, bad barcodes: {}'.format(self._ok_barcodes, self._fix_barcodes, self._bad_barcodes))
        return err_rate

    def group_for_disambiguation(self):
        """
        Prepare the RA for multi alignment disambiguation. Apply filters, correct the barcodes and group by cell/RMT pairs.
        Retain the indices corresponding to each read.
        The resulting dictionary has the following structure:
            res[correct cell barcode, rmt] = {(mapped genes): [read1, read2, ...]}
        Note the mapped genes list should be sorted before converted to a tuple
        """
        start = time.process_time()

        res = {}

        for i, r in enumerate(self.data):
            if not ReadArray.is_active_read(r):
                continue

            genes = tuple(set(sorted(r['ma_genes'])))
            # group according to the correct barcodes and gene
            cell = r['cell']
            rmt = r['rmt']
            try:
                res[cell, rmt][genes].append(i)
            except KeyError:
                try:
                    res[cell, rmt][genes] = [i]
                except KeyError:
                    res[cell, rmt] = {genes: [i]}

        tot_time=time.process_time()-start
        log.info('Grouping for multimapping completed in {} seconds.'.format(tot_time))
        return res
    
    def save(self, archive_name: str) -> None:
        """save a ReadArray in .h5 format

        :param archive_name: filename of a new .h5 archive in which to save the ReadArray
        :return: None
        """

        # create table
#        blosc5 = tb.Filters(complevel=5, complib='blosc')
#        f = tb.open_file(archive_name, mode='w', title='Data for seqc.ReadArray',
#                         filters=blosc5)

        # store data
#        f.create_table(f.root, 'data', self._data)
#        f.close()
        f = open(archive_name, 'wb')
        pickle.dump(self._data, f, protocol=4)
        f.close()

    @classmethod
    def load(cls, archive_name: str):
        """load a ReadArray from a .h5 archive

        :param archive_name: name of a .h5 archive containing a saved ReadArray object
        :return: ReadArray
        """

        #f = tb.open_file(archive_name, mode='r')
        #data = f.root.data.read()
        #f.close()
        #return cls(data)
        
        f = open(archive_name, 'rb')
        data = pickle.load(f)
        f.close()
        return cls(data)
        
    def resolve_alignments(self, index):
        """
        Resolve ambiguously aligned molecules and edit the ReadArray data structures
        in-place to reflect the more specific gene assignments.
        args:
        index: path were p_coalignment_array.p is found - . 
        
        After loading the co alignment matrix we group the reads of the ra by cell/rmt. In each group we look at the different disjoint subsetes of genes reads are aligned to.
        For each subset there are 4 options:
        1. It includes one gene only (no disambiguation).
        2. There is a single gene that explains all (resolved gene)
        3. There is more than one gene that can explain - can't be resolved
        4. There are a few models - for now no single model is being picked.
        """
              
        start = time.process_time()
        total_reads = 0

        log.info('grouping for disambiguation')
        #Group according to cell/RMT
        g_ra = self.group_for_disambiguation()
        
        #TODO - for now we're not using the expectations so I commented this part. -> Rami
        # After uncommenting, we need a new p_coalignment_array.p file to be created with the same mosule hirerachy
        # load the expectations
        #expectations = index + 'p_coalignment_array.p'
        #if os.path.isfile(expectations):
        #    with open(expectations, 'rb') as f:
        #        expectations_mat = pickle.load(f)
        #elif isinstance(expectations, str):
        #    raise FileNotFoundError('could not locate serialized expectations object: %s'
        #                            % expectations)
        #elif isinstance(expectations, dict):
        #    pass  # expectations already loaded
        #else:
        #    raise TypeError('invalid expectation object type, must be a dict expectation'
        #                    ' object or a string filepath')
        #mat = mm.reduce_coalignment_array(expectations_mat)
        mat=[]
        
        log.info('disambiguating')
        for g_ra_rec in g_ra.values():
            obs = {}
            for genes in g_ra_rec:
                if genes == ():
                    continue
                if genes not in obs:
                    obs[genes] = len(g_ra_rec[genes])
                    total_reads += len(g_ra_rec[genes])
                #TODO: do a sanity check here? I think genes should not be in obs
                    
            
            res = []
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
                if r==mm.NO_DISAMBIGUATION or r==mm.RESOLVED_GENE:
                        model=models[0]
                        g = mm.model_to_gene(model)
                        resolved = True
                        
                indices = mm.get_indices(ind_s, obs_subset)
                for ind in indices:
                    ReadArray.set_mm_status(self.data[ind], r)
                    if resolved:
                        self.data['gene'][ind] = g
                        self.data['position'][ind] = self.data['ma_pos'][ind][self.data['ma_genes'][ind].index(g)]
                    else:
                        ReadArray.set_status_inactive(self.data[ind])        
        
        log.info('Total reads in samfile: {}'.format(len(self)))
        log.info('Reads that passed the filter: {}'.format(total_reads))
        tot_time=time.process_time()-start
        log.info('Multimapping resolution completed in {} seconds.'.format(tot_time))

        return

    def __len__(self):
        return len(self.data)
        
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
                if r['ma_genes'] == -1 or r['ma_genes'] == 0:
                    continue
                elif len(r['ma_genes']) == 1:
                    gene = r['ma_genes'][0]
                else:
                    continue
            try:
                reads_mat[cell, gene] = reads_mat[cell, gene]+1
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

