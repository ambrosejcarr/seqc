import numpy as np
import tables as tb
from seqc.alignment import sam
from seqc.sequence.gtf import GeneIntervals
from seqc.sequence.encodings import DNA3Bit
import os
import seqc.correct_errors
import seqc.multialignment as mm
import pickle
    

class ReadArray:

    _dtype = [
        ('status', np.uint16),
        ('cell', np.int64),
        ('rmt', np.int32),
        ('n_poly_t', np.uint8),
        ('dust_score', np.uint8),
        ('gene', np.uint32),
        ('position', np.uint64),
        ('ma_genes', np.dtype(object))]
    
    # Some constants
    PASS_FILTERS    = 0
    NO_ALIGNMENT    = 1
    GENE_0          = 2
    CELL_0          = 3
    RMT_0           = 4
    POLY_T_0        = 5
    SCAFFOLD        = 6

    def __init__(self, data):
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
        self._no_alignment      = 0
        self._non_unique_align  = 0
        self._filtered          = {NO_ALIGNMENT: 0, GENE_0: 0, NO_CELL: 0, NO_RMT: 0, POLY_T: 0, SCAFFOLD: 0}

    
    @property
    def data(self):
        return self._data

    # Toggle active flag on (byte 15)
    @staticmethod
    def set_status_active(r):
        r['status'] |= 0b1
    
    def set_status_inactive(r):
        r['status'] &= 0b1111111111111110

    @staticmethod
    def set_filt_fail(r, filt):
        if filt > 0b1111:
            raise ValueError('Filter value must be 4 bits. value given: {}'.format(filt))
        r['status'] &= 0b1111111111100001
        r['status'] |= filt<<1
        
    @classmethod
    def from_samfile(cls, samfile: str, gtf_f: str):
        """
        construct a ReadArray object from a samfile containing only uniquely aligned
        records

        :param gtf: str, filename of annotations.gtf file
        :param samfile: str, filename of alignment file.
        :return:
        """
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
                    data[read_idx-1] = (cell, rmt, n_poly_t, dust_score, 0, pos, genes, 0)
                    cls.set_status_active(data[read_idx-1])
                prev_alignment_name = alignment.qname
                read_idx+=1
                genes = translator.translate(alignment.rname, alignment.strand, alignment.pos)
                if genes != -1:
                    if len(genes) > 1:
                        self._non_unique_align+=1
                cell = seqc.sequence.encodings.DNA3Bit.encode(alignment.cell)
                rmt = seqc.sequence.encodings.DNA3Bit.encode(alignment.rmt)
                n_poly_t = alignment.poly_t.count(b'T') + alignment.poly_t.count(b'N')
                dust_score = alignment.dust_low_complexity_score
                pos = alignment.pos
            #the same read
            else:
                genes_l = translator.translate(alignment.rname, alignment.strand, alignment.pos)
                if genes_l != -1:
                    if len(genes_l) > 1:
                        self._non_unique_align+=1
                    if genes == -1:
                        genes = genes_l
                    else:
                        genes+=genes_l
            
        seqc.log.info('non unique alignments count: {}'.format(self._non_unique_align))
        #seqc.log.info('focus gene {}. Total alignmens: {}, total reads: {}, total uniquely aligned reads: {}'.format(focus_gene, fgene_tot_reads, fgene_tot_reads, fgene_tot_uniq))
        
        res_ra = cls(data)
        res_ra._total_alignments = len(reader)
        res_ra._total_reads = num_reads

        return res_ra, res_ra._total_alignments #TODO: It might be nicer to return just the RA and take the total alignments from it

    def apply_filters(self, required_poly_t = 1):
        """
        Apply different filters to the read array.
        If a read fails a filter, it is not passed to the others and so the counts
        for each filter is the number of reads that failed that one but passed all previous ones
        Filters are not ordered in any particular way.
        """
        for r in self.data:
            if len(r['ma_genes']) == 0:
                self.set_filt_fail(r, self.NO_ALIGNMENT)
                self._filtered[NO_ALIGNMENT] += 1
                self.set_status_inactive(r)
                continue
            if r['ma_genes'] == -1:
                self.set_filt_fail(r, self.SCAFFOLD)
                self._filtered[SCAFFOLD] += 1
                self.set_status_inactive(r)
                continue
            if r['cell'] == 0:
                self.set_filt_fail(r, self.NO_CELL)
                self._filtered[NO_CELL] += 1
                self.set_status_inactive(r)
                continue
            if r['rmt'] == 0:
                self.set_filt_fail(r, self.NO_RMT)
                self._filtered[NO_RMT] += 1
                self.set_status_inactive(r)
                continue
            if r['n_poly_t'] < required_poly_t:
                self.set_filt_fail(r, self.POLY_T)
                self._filtered[POLY_T] += 1
                self.set_status_inactive(r)
                continue
            if 0 in r['ma_genes']:
                self.set_filt_fail(r, self.GENE_0)
                self._filtered[GENE_0] += 1
                self.set_status_inactive(r)
                continue
                
        seqc.log.info('Filters results: No alignment: {}, scaffold: {}, no cell {}, no rmt: {}, low poly_t: {}, mapped to gene 0: {}'.
        format(self._filtered[NO_ALIGNMENT], self._filtered[SCAFFOLD], self._filtered[NO_CELL], self._filtered[NO_RMT], self._filtered[POLY_T], self._filtered[GENE_0]))
        
            
                
        
    def reads_passing_filters(self, min_poly_t: int, max_dust_score: int):
        """
        Subset the ReadArray returning a new ReadArray containing eads that passed all
        filters

        :param min_poly_t: int, minimum number of T nucleotides that defined a valid
          capture primer
        :param max_dust_score: int, higher scores indicate increasingly degenerate
          sequences. Typically sequences with dust_score > 10 may be discarded.
        :return: ReadArray
        """
        phix_genes = np.array(range(1, 7)) * 111111111
        data = self.data[((self.data['n_poly_t'] >= min_poly_t) &
                          (self.data['dust_score'] <= max_dust_score) &
                          (self.data['gene'] != 0) &
                          (self.data['cell'] != 0) &
                          (self.data['rmt'] != 0))]

        # filter out phiX genes
        not_phix = ~np.in1d(data['gene'], phix_genes)
        data = data[not_phix]

        # filter out N's in rmt
        res = np.zeros(len(data), dtype=np.bool)
        rmt = data['rmt'].copy()
        while np.any(rmt):
            n_filter = rmt & 0b111 == 0b111
            res[n_filter] = True
            rmt >>= 3
        data = data[~res]
        return ReadArray(data)

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
        
    def resolve_alignments(self, index, barcode_files, required_poly_t=1, reverse_complement = True, max_ed=2):
        """
        Resolve ambiguously aligned molecules and edit the ReadArray data structures
        in-place to reflect the more specific gene assignments.
        args:
        index: path were p_coalignment_array.p is found - . 
        required_poly_t (default: 0): the number of poly_t required for a read to be
          considered a valid input for resolution. Reads lacking a poly_t tail are likely
          to randomly prime into non-mRNA sequences, so requiring the presence of poly_t
          in the 5' end of the forward read can reduce unwanted contamination
        
        After loading the co alignment matrix we group the reads of the ra by cell/rmt. In each group we look at the different disjoint subsetes of genes reads are aligned to.
        For each subset there are 4 options:
        1. It includes one gene only (no disambiguation).
        2. There is a single gene that explains all (resolved gene)
        3. There is more than one gene that can explain - can't be resolved
        4. There are a few models - for now no single model is being picked.
        """
        
        #for debugging
        focus_gene = 3
        fg_gra_cnt = 0
        #
                
        reads_res = {mm.NO_DISAMBIGUATION:0,mm.RESOLVED_GENE:0,mm.NO_GENE_RESOLVED:0,mm.MULTIPLE_MODELS:0}
        total_reads = 0

        seqc.log.info('grouping for disambiguation')
        #Group according to cell/RMT
        g_ra = seqc.correct_errors.group_for_disambiguation(self, barcode_files, required_poly_t, reverse_complement, max_ed=2)
        
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
        
        seqc.log.info('disambiguating')
        for g_ra_rec in g_ra.values():
            obs = {}
            for genes in g_ra_rec:
                if genes == ():
                    continue
                if focus_gene in genes:
                    fg_gra_cnt += len(g_ra_rec[genes])
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
                models, r = mm.best_fit_model(obs_subset, mat)
                reads_res[r] += sum(obs_subset.values())
                if r==mm.NO_DISAMBIGUATION or r==mm.RESOLVED_GENE:
                    model=models[0]
                    indices = mm.get_indices(ind_s, obs_subset)
                    g = mm.model_to_gene(model)
                    for ind in indices:
                        self.data['gene'][ind] = g
                        
        
        seqc.log.info('Total reads in samfile: {}'.format(len(self)))
        seqc.log.info('Reads that passed the filter: {}'.format(total_reads))
        seqc.log.info('Aligned reads: {}, resolved reads: {}, can\'t be resolved: {}, competing models (unhandeled): {}'.format(reads_res[mm.NO_DISAMBIGUATION], 
                reads_res[mm.RESOLVED_GENE], reads_res[mm.NO_GENE_RESOLVED], reads_res[mm.MULTIPLE_MODELS]))
        
        
        #for debugging
        fg_mm_resolve_cnt = 0
        for v in self.data:
            if v['gene'] == focus_gene:
                fg_mm_resolve_cnt += 1
        seqc.log.info('Focus gene: {}. total reads in grouped ra: {}, total reads after mm resolving: {}'.format(focus_gene, fg_gra_cnt, fg_mm_resolve_cnt))
        #
        return

    def __len__(self):
        return len(self.data)
        
    def to_count_matrix(self):
        #TODO - write this
