import numpy as np
from seqc import log
from seqc.alignment import sam
from seqc.sequence.gtf import GeneIntervals
from seqc.sequence.encodings import DNA3Bit
from seqc.sparse_frame import SparseFrame
import seqc.multialignment as mm
import seqc.sequence.barcodes
from itertools import permutations
import tables as tb
import pandas as pd
from statsmodels.sandbox.stats.multicomp import multipletests as mt

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
LOW_COVERAGE    = 9

BC_OK           = 0
BC_FIXED        = 1
BC_BAD          = 2

DEFAULT_BASE_CONVERTION_RATE = 0.02


class ReadArray:

    _read_dtype = [
        ('status', np.uint16),
        ('cell', np.int64),
        ('rmt', np.int32),
        ('n_poly_t', np.uint8),
        ('dust_score', np.uint8),
        ('gene', np.uint32),
        ('position', np.uint64)
    ]

    _alignment_dtype = [
        ('read_idx', np.uint32),
        ('pos', np.uint64),
        ('gene', np.uint32)
    ]

    def __init__(self, data, alignments, non_unique_alignments=-1):
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
        self._alignments = alignments
        
        #some stats
        self._total_alignments  = 0
        self._total_reads       = 0
        self._total_active_reads = 0
        self._non_unique_align  = 0
        self._filtered          = {NO_ALIGNMENT: 0, GENE_0: 0, NO_CELL: 0, NO_RMT: 0, POLY_T: 0, SCAFFOLD: 0, NO_GENES: 0, DUST_SCORE: 0, LOW_COVERAGE: 0}
        self._bad_barcodes      = 0
        self._fix_barcodes      = 0
        self._ok_barcodes       = 0
        self._multialignment    = {mm.NO_DISAMBIGUATION: 0, mm.RESOLVED_GENE: 0, mm.NO_GENE_RESOLVED: 0, mm.MULTIPLE_MODELS: 0}
        self._total_errors      = 0

        self._resolved_alignments = False
        
        if non_unique_alignments != -1:
            self._non_unique_align = non_unique_alignments

    def __repr__(self):
        s = ""
        s += "Read array of size: {:,}\n".format(len(self._data))
        s += "Non unique alignments (mapped to more than one gene): {:,}\n".format(self._non_unique_align)
        s += "Total reads: {:,}\n".format(self._total_reads)
        s += "Filtered: No alignment:{:,}, gene is 0:{:,}, no cell: {:,}, no rmt: {:,}, poly_t:{:,}, scaffold: {:,}," \
             " no genes at all: {:,}, low coverage(marseq only): {:,}\n".\
             format(self._filtered[NO_ALIGNMENT], self._filtered[GENE_0],self._filtered[NO_CELL],
                    self._filtered[NO_RMT], self._filtered[POLY_T], self._filtered[SCAFFOLD], self._filtered[NO_GENES],
                    self._filtered[LOW_COVERAGE])
        s += "Barcodes: good:{:,}, fixed: {:,}, bad:{:,}\n".\
            format(self._ok_barcodes, self._fix_barcodes, self._bad_barcodes)
        s += "Multialignment results: no disambiguation: {:,}, resolved gene: {:,}, no gene resolved: {:,}, " \
             "competing models (unhandeled): {:,}\n".\
             format(self._multialignment[mm.NO_DISAMBIGUATION], self._multialignment[mm.RESOLVED_GENE],
                    self._multialignment[mm.NO_GENE_RESOLVED], self._multialignment[mm.MULTIPLE_MODELS])
        s += "Total errors: {:,}\n".format(self._total_errors)
        s += "Total active reads: {:,}\n".format(self.count_active_reads())
        return s
    
    @property
    def data(self):
        return self._data

    @property
    def alignments(self):
        return self._alignments

    @property
    def resolved_alignments(self):
        return self._resolved_alignments

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
    def from_samfile(cls, samfile, gtf_f):
        """
        construct a ReadArray object from a samfile containing only uniquely aligned
        records

        :param gtf_f: str, filename of annotations.gtf file
        :param samfile: str, filename of alignment file.
        :return:
        """
        #non_unique_align = 0
        no_alignment_count = 0
        reader = sam.Reader(samfile)
        #go over the file to find how many unique reads there are
        num_reads=0
        prev_alignment_name=''
        for alignment in reader:
            if alignment.qname != prev_alignment_name:
                num_reads +=1
                prev_alignment_name = alignment.qname

        translator = GeneIntervals(gtf_f)

        data = np.recarray((num_reads,), ReadArray._read_dtype)
        alignments = np.recarray((len(reader),), ReadArray._alignment_dtype)
        prev_alignment_name=''
        read_idx = 0
        alignment_ctr = 0
        genes = []

        max_alignments = 100
        for alignment in reader:
            ########## debug
            #if alignment_ctr > max_alignments:
            #    break
            ##############
            # new read
            if alignment.qname != prev_alignment_name:
                if read_idx>0:
                    data[read_idx-1] = (ReadArray.active_status(), cell, rmt, n_poly_t, dust_score, 0, 0)
                    cls.set_status_active(data[read_idx-1])
                    # todo can count the gene=0 and no genes seperately
                    if len(genes)>0 and 0 not in genes:    # Add all valid alignments to alignment table
                        for gene in genes:
                            alignments[alignment_ctr] = (read_idx-1, alignment.pos, gene)
                            alignment_ctr += 1
                    else:
                        ReadArray.set_filt_fail(data[read_idx-1], NO_ALIGNMENT)
                        no_alignment_count += 1
                        ReadArray.set_status_inactive(data[read_idx-1])

                prev_alignment_name = alignment.qname
                read_idx+=1
                genes = []

                cell = seqc.sequence.encodings.DNA3Bit.encode(alignment.cell)
                rmt = seqc.sequence.encodings.DNA3Bit.encode(alignment.rmt)
                n_poly_t = alignment.poly_t.count(b'T') + alignment.poly_t.count(b'N')
                dust_score = alignment.dust_low_complexity_score
                genes_mapped = translator.translate(alignment.rname, alignment.strand, alignment.pos)
                if genes_mapped != -1:
                    if len(genes_mapped) == 1:
                        genes.append(genes_mapped[0])


            # the same read
            else:
                genes_mapped = translator.translate(alignment.rname, alignment.strand, alignment.pos)
                if genes_mapped != -1:
                    if len(genes_mapped) == 1:
                        genes.append(genes_mapped[0])

        #log.info('non unique alignments count: {}'.format(non_unique_align))
        res_ra = cls(data, np.delete(alignments, range(alignment_ctr,len(alignments))))
        res_ra._total_alignments = alignment_ctr
        res_ra._total_reads = num_reads
        res_ra._filtered[NO_ALIGNMENT] += no_alignment_count    #todo remove these and add them to the ctr

        return res_ra

    def apply_filters(self, required_poly_t = 1, max_dust_score=10):
        """
        Apply different filters to the read array.
        If a read fails a filter, it is not passed to the others and so the counts
        for each filter is the number of reads that failed that one but passed all previous ones
        Filters are not ordered in any particular way.0
        """
        for r in self.data:
            if not ReadArray.is_active_read(r):
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
        The result is 3 dictionaries: reads, genes and positions. Each grouped by cell,rmt.
            r_dic[cell, rmt] = list of indices of reads with that cell,rmt
            g_dic[cell, rmt] = list of lists of genes corresponding to the reads in r of the genes the read is aligned to
            p_dic[cell, rmt] = the same with positions
        Note the mapped genes list should be sorted before converted to a tuple
        """

        r_dic = {}
        g_dic = {}
        p_dic = {}

        for i, alignment in enumerate(self.alignments):
            read_idx = alignment['read_idx']
            if not ReadArray.is_active_read(self.data[read_idx]):
                continue
            cell = self.data[read_idx]['cell']
            rmt = self.data[read_idx]['rmt']
            gene = alignment['gene']
            pos = alignment['pos']

            try:
                if read_idx in r_dic[cell,rmt]:  #This is an alignment of a read we've already seen
                    inner_idx = r_dic[cell,rmt].index(read_idx)
                    g_dic[cell,rmt][inner_idx].append(gene)
                    p_dic[cell, rmt][inner_idx].append(pos)
                else:
                    r_dic[cell, rmt].append(read_idx)
                    g_dic[cell, rmt].append([gene])
                    p_dic[cell, rmt].append([pos])
            except KeyError:
                r_dic[cell, rmt] = [read_idx]
                g_dic[cell, rmt] = [[gene]]
                p_dic[cell, rmt] = [[pos]]

        return r_dic, g_dic, p_dic

    def save(self, archive_name):
        """save a ReadArray object as a .h5 archive: note that ma_genes and ma_pos are
        discarded

        :param archive_name: filestem for the new archive
        :return: None
        """

        if not archive_name.endswith('.h5'):
            archive_name += '.h5'

        # construct container
        blosc5 = tb.Filters(complevel=5, complib='blosc')
        f = tb.open_file(archive_name, mode='w', title='Data for seqc.ReadArray',
                         filters=blosc5)

        # select the subset of the readarray that is storable and put that in h5
        data = self._data[
            ['status', 'cell', 'rmt', 'n_poly_t', 'dust_score', 'gene', 'position']]
        f.create_table(f.root, 'data', data)
        f.close()

    @classmethod
    def load(cls, archive_name):
        """load a ReadArray from an hdf5 archive

        :param archive_name: name of a .h5 archive containing a saved ReadArray object
        :return: ReadArray
        """

        f = tb.open_file(archive_name, mode='r')
        data = f.root.data.read()
        f.close()

        return cls(data)

    def resolve_alignments(self):
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


        log.info('grouping for disambiguation')
        #Group according to cell/RMT
        r_dic, g_dic, p_dic = self.group_for_disambiguation()
        
        log.info('disambiguating')

        for cell,rmt in g_dic:
            gene_lists = g_dic[cell, rmt]
            reads_list = r_dic[cell,rmt]
            uf = mm.UnionFind()
            uf.union_all(gene_lists)
            set_membership, sets = uf.find_all(gene_lists)

            for s in sets:
                relevant_reads = list(np.array(reads_list)[set_membership == s])
                res_gene, mm_status = mm.resolve_gene(list(np.array(gene_lists)[set_membership == s]))
                self._multialignment[mm_status] += len(relevant_reads)

                for read_idx in relevant_reads:
                    ReadArray.set_mm_status(self.data[read_idx], mm_status)

                    if mm_status == mm.NO_DISAMBIGUATION or mm_status == mm.RESOLVED_GENE:
                        i = reads_list.index(read_idx)
                        res_pos = p_dic[cell, rmt][i][gene_lists[i].index(res_gene)]
                        self.data['gene'][read_idx] = res_gene
                        self.data['position'][read_idx] = res_pos
                    else:
                        ReadArray.set_status_inactive(self.data[read_idx])


        self.resolve_alignments = True

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
        if not self.resolved_alignments:
            log.info('Read array converted to count matrix before alignments were resolved. {} will be empty.'.format(csv_path))

        reads_mat = {}
        mols_mat = {}
        for r in self.data:
            if not ReadArray.is_active_read(r):
                continue
            cell = r['cell']
            gene = r['gene']
            rmt = r['rmt']
            if gene == 0:
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


    def filter_low_coverage(self, alpha_value = 0.25):
        cell = self.data['cell']
        position = self.data['position']
        rmt = self.data['rmt']
        genes = self.data['gene']

        # A triplet is a (cell, position, rmt) triplet in each gene
        df = pd.DataFrame({'gene': genes, 'cell': cell, 'position': position,
                           'rmt': rmt})
        grouped = df.groupby(['gene', 'position'])
        # This gives the gene followed by the number of triplets at each position
        # Summing across each gene will give the number of total triplets in gene
        num_per_position = (grouped['position'].agg({
            'Num Triplets at Pos': np.count_nonzero})).reset_index()

        # Total triplets in each gene
        trips_in_gene = (num_per_position.groupby(['gene'])
                         )['Num Triplets at Pos'].agg({'Num Triplets at Gene': np.sum})

        trips_in_gene = trips_in_gene.reset_index()

        num_per_position = num_per_position.merge(trips_in_gene, how='left')

        # for each (c,rmt) in df check in grouped2 if it is lonely
        # determine number of lonely triplets at each position
        grouped2 = df.groupby(['gene', 'cell', 'rmt'])
        # lonely_triplets = grouped2["position"].apply(lambda x: len(x.unique()))
        # This is a list of each gene, cell, rmt combo and the positions with that criteria
        lonely_triplets = grouped2['position'].apply(np.unique)
        lonely_triplets = pd.DataFrame(lonely_triplets)

        # if the length is one, this is a lonely triplet
        lonely_triplets_u = lonely_triplets['position'].apply(len)
        lonely_triplets_u = pd.DataFrame(lonely_triplets_u)

        lonely_triplets_u = lonely_triplets_u.reset_index()
        lonely_triplets = lonely_triplets.reset_index()

        # Rename the columns
        lonely_triplets = lonely_triplets.rename(columns=lambda x: x.replace(
            'position', 'lonely position'))
        lonely_triplets_u = lonely_triplets_u.rename(columns=lambda x: x.replace(
            'position', 'num'))

        # merge the column that is the length of the positions array
        # take the ones with length 1
        lonely_triplets = lonely_triplets.merge(lonely_triplets_u, how='left')
        lonely_triplets = lonely_triplets.loc[lonely_triplets.loc[:, 'num'] == 1, :]

        # This is the gene, cell, rmt combo and the position that is lonely
        # We need to convert the array to a scalar
        scalar = lonely_triplets["lonely position"].apply(np.asscalar)
        lonely_triplets["lonely position"] = scalar
        # Now if we group as such, we can determine how many (c, rmt) paris exist at each position
        # This would be the number of lonely pairs at a position
        grouped3 = lonely_triplets.groupby(["gene", "lonely position"])
        l_num_at_position = (grouped3["cell"].agg(['count'])).reset_index()
        l_num_at_position = l_num_at_position.rename(columns=lambda x: x.replace(
            'count', 'lonely triplets at pos'))
        l_num_at_position = l_num_at_position.rename(columns=lambda x: x.replace(
            'lonely position', 'position'))
        # lonely pairs in each gene
        l_num_at_gene = (lonely_triplets.groupby(["gene"]))['lonely position'].agg(
            ['count'])
        l_num_at_gene = l_num_at_gene.reset_index()
        l_num_at_gene = l_num_at_gene.rename(columns=lambda x: x.replace(
            'count', 'lonely triplets at gen'))

        # aggregate
        total = l_num_at_position.merge(l_num_at_gene, how='left')
        total = total.merge(num_per_position, how='left')

        # scipy hypergeom
        p = total.apply(ReadArray.hypergeom_wrapper, axis=1)
        p = 1 - p

        adj_p = mt(p, alpha=alpha_value, method='fdr_bh')

        keep = pd.DataFrame(adj_p[0])
        total['keep'] = keep

        remove = total[total['keep'] == True]
        d = df.merge(remove, how="left")
        final = d[d["keep"] == True]

        for idx in final.index:
            ReadArray.set_filt_fail(self.data[idx], LOW_COVERAGE)
            self._filtered[LOW_COVERAGE]+=1
            ReadArray.set_status_inactive(self.data[idx])

        return

    @staticmethod
    def hypergeom_wrapper(x):
        from scipy.stats import hypergeom
        p = hypergeom.cdf(x['lonely triplets at pos'], x['Num Triplets at Gene'],
                          x['lonely triplets at gen'], x['Num Triplets at Pos'])
        return p
