from scipy.special import gammainc
from itertools import permutations
from sys import maxsize
import time
from scipy.sparse import coo_matrix
from seqc.sequence.encodings import ThreeBit as BinRep
import numpy as np
import seqc
import random

high_value = maxsize  # Used for sorting, needs to be longer than any sequence

# todo: increase number of error correction methods by 2
NUM_OF_ERROR_CORRECTION_METHODS = 3
#ERROR_CORRECTION_AJC = 0

ERROR_CORRECTION_BC_FILTERS = 3
ERROR_CORRECTION_LIKELIHOOD_NOT_ERROR = 0
#ERROR_CORRECTION_jaitin = 1
ERROR_CORRECTION_JAIT_LIKELIHOOD = 2

DEFAULT_BASE_CONVERTION_RATE = 0.02


def prepare_for_ec(ra, barcode_files, required_poly_t=1, reverse_complement=True,
                   max_ed=2, err_correction_mat=''):
    """
    Prepare the RA for error correction. Apply filters, estimate error correction and
    correct the barcodes
    :param err_correction_mat:
    :param max_ed:
    :param reverse_complement:
    :param required_poly_t:
    :param barcode_files:
    :param ra:
    """
    res = {}
    tot = 0
    bc_filter = 0
    filtered = {'gene_0': 0, 'phi_x': 0, 'cell_0': 0, 'rmt_0': 0, 'poly_t': 0, 'rmt_N': 0}

    dynamic_codes_table_c1 = {}
    dynamic_codes_table_c2 = {}

    errors = list(BinRep.ints2int([p[0], p[1]]) for p in permutations(BinRep.bin_nums,
                                                                      r=2))
    error_table = dict(zip(errors, [0] * len(errors)))
    cor_instance_table = {BinRep._str2bindict['A']: 0,
                          BinRep._str2bindict['C']: 0,
                          BinRep._str2bindict['G']: 0,
                          BinRep._str2bindict['T']: 0}

    # Read the barcodes into lists
    correct_barcodes = []
    if reverse_complement:
        for barcode_file in barcode_files:
            with open(barcode_file) as f:
                correct_barcodes.append(set(BinRep.str2bin(rev_comp(line.strip()))
                                            for line in f.readlines()))
    else:
        for barcode_file in barcode_files:
            with open(barcode_file) as f:
                correct_barcodes.append(set(BinRep.str2bin(line.strip())
                                            for line in f.readlines()))

    for i, v in enumerate(ra.data):
        tot += 1
        if not pass_filter(v, filtered, required_poly_t):
            continue

        gene = v['gene']

        # correct and filter barcodes
        c1 = BinRep.c1_from_codes(int(v['cell']))
        try:
            cor_c1, ed_1, err_l_1, seq_err_1 = dynamic_codes_table_c1[c1]
        except KeyError:
            cor_c1, ed_1 = find_correct_barcode(c1, correct_barcodes[0])
            seq_err_1 = count_seqeuncer_errors(c1)
            # No point calculating errors on codes with more errors than allowed
            if ed_1 > max_ed:
                err_l_1 = None
            else:
                err_l_1 = list_errors(cor_c1, c1)
            dynamic_codes_table_c1[c1] = (cor_c1, ed_1, err_l_1, seq_err_1)

        c2 = BinRep.c2_from_codes(int(v['cell']))
        try:
            cor_c2, ed_2, err_l_2, seq_err_2 = dynamic_codes_table_c2[c2]
        except KeyError:
            cor_c2, ed_2 = find_correct_barcode(c2, correct_barcodes[1])
            seq_err_2 = count_seqeuncer_errors(c2)
            # No point calculating errors on codes with more errors than allowed
            if ed_2 > max_ed:
                err_l_2 = None
            else:
                err_l_2 = list_errors(cor_c2, c2)
            dynamic_codes_table_c2[c2] = (cor_c2, ed_2, err_l_2, seq_err_2)

        # Filter reads with too many errors in their barcodes
        if ed_1 > max_ed or ed_2 > max_ed:
            bc_filter += 1
            # Uncomment for debugging
#            if err_correction_mat != '':
#                err_correction_mat[i,ERROR_CORRECTION_BC_FILTERS] = 1
            continue
        
        # For error rate estimation, only reads with at most a single error in both barcodes are counted.
        # note that only base swicth errors are counted towards this threshold, sequencer errors (that cause an 'N') are ignored for this.
        if ed_1 + ed_2 - seq_err_1 - seq_err_2 <= 1:
            for er in err_l_1 + err_l_2:
                try:
                    error_table[er] += 1
                except KeyError:    #happens for sequencer errors
                    continue

            # count non error bases
            tmp_c = BinRep.ints2int([c1, c2])
            tmp_cor = BinRep.ints2int([cor_c1, cor_c2])
            while tmp_c > 0:
                if tmp_c & 0b111 == tmp_cor & 0b111:
                    cor_instance_table[tmp_c & 0b111] += 1
                tmp_c >>= 3
                tmp_cor >>= 3

        # group according to the correct barcodes and gene
        cell = BinRep.ints2int([cor_c1, cor_c2])
        rmt = v['rmt']
        try:
            res[gene, cell][rmt].append(i)
        except KeyError:
            try:
                res[gene, cell][rmt] = [i]
            except KeyError:
                res[gene, cell] = {}
                res[gene, cell][rmt] = [i]

    # convert to error rates
    default_error_rate = DEFAULT_BASE_CONVERTION_RATE
    err_rate = dict(zip(errors, [0.0] * len(errors)))
    if sum(error_table.values()) == 0:
        seqc.log.info('No errors were detected, using %f uniform error chance.' % (
            default_error_rate))
        err_rate = dict(zip(errors, [default_error_rate] * len(errors)))
    for k, v in error_table.items():
        try:
            err_rate[k] = v / (sum(n for err_type, n in error_table.items()
                               if err_type & 0b111000 == k & 0b111000) +
                               cor_instance_table[(k & 0b111000) >> 3])
        except ZeroDivisionError:
            seqc.log.info('Warning: too few reads to estimate error rate for %r '
                          'setting default rate of %f' % (k, default_error_rate))
            err_rate[k] = default_error_rate

    seqc.log.info('Error correction filtering results: total reads: {}; '
                  'did not pass preliminary filters: {}; cell barcodes are wrong: '
                  '{}'.format(tot, sum(filtered.values()), bc_filter))
    # print('error_table: ', error_table, ' cor_instance_table: ', cor_instance_table)

    return res, err_rate, filtered


# Return the hamming distance between two numbers representing a sequence (3 bits / base)
def hamming_dist_bin(c1, c2):
    """
    :param c1:
    :param c2:
    :return:
    """
    if BinRep.seq_len(c1) != BinRep.seq_len(c2):
        return high_value
    d = 0
    while c1 > 0:
        if c1 & 0b111 != c2 & 0b111:
            d += 1
        c1 >>= 3
        c2 >>= 3
    return d

def count_seqeuncer_errors(c):
    err=0
    while c>0:
        if BinRep._bin2strdict[c&0b111] == 'N':
            err+=1
        c>>=3
    return err

def generate_close_seq(seq):
    """ Return a list of all sequences that are up to 2 hamm distance from seq
    :param seq:
    """
    res = []
    l = BinRep.seq_len(seq)

    # generate all sequences that are dist 1
    for i in range(l):
        mask = 0b111 << (i * 3)
        cur_chr = (seq&mask) >> (i * 3)
        res += [seq & (~mask) | (new_chr << (i * 3))
                for new_chr in BinRep.bin_nums if new_chr != cur_chr]
    # generate all sequences that are dist 2
    for i in range(l):
        mask_i = 0b111 << (i * 3)
        chr_i = (seq & mask_i) >> (i * 3)
        for j in range(i + 1, l):
            mask_j = 0b111 << (j * 3)
            chr_j = (seq & mask_j) >> (j * 3)
            mask = mask_i | mask_j
            res += [seq & (~mask) | (new_chr_i << (i * 3)) | (new_chr_j << (j * 3)) for
                    new_chr_i in BinRep.bin_nums if new_chr_i != chr_i for
                    new_chr_j in BinRep.bin_nums if new_chr_j != chr_j]

    return res


def prob_d_to_r_bin(d_seq, r_seq, err_rate):
    """
    Return the probability of d_seq turning into r_seq based on the err_rate table
    (all binary)

    :param err_rate:
    :param r_seq:
    :param d_seq:
    """

    if BinRep.seq_len(d_seq) != BinRep.seq_len(r_seq):
        return 1

    p = 1.0
    while d_seq > 0:
        if d_seq & 0b111 != r_seq & 0b111:
            p *= err_rate[BinRep.ints2int([d_seq & 0b111, r_seq & 0b111])]
        d_seq >>= 3
        r_seq >>= 3
    return p


_complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}  # for reverse complement


def rev_comp(s):
    """Return the reverse complement of an ACGT string s

    :param s:
    :returns:
    """
    ret = []
    for i in range(len(s)):
        ret.append(_complement[(s[-1 - i])])
    return ''.join(ret)


def list_errors(s1, s2):
    """
    Return the list of nucleotide transformations that turn s1 to s2.
    An error is a six bit int representing a two chr string of type "AG","CT", etc.

    :param s2:
    :param s1:

    :returns:
    """

    # return the actual error
    err_list = []
    while s1 > 0:
        if s1 & 0b111 != s2 & 0b111:
            err_list.append(BinRep.ints2int([s1 & 0b111, s2 & 0b111]))
        s1 >>= 3
        s2 >>= 3
    return err_list


def find_correct_barcode(code, barcodes_list):
    """
    For a given barcode find the closest correct barcode to it from the list (limited to
    one ED), a string representing the error and the edit distance
    NOTE: for now this function looks for a barcode with ED==1 and does not bother
    looking for the minimum

    :param barcodes_list:
    :param code:
    :returns:
    """
    if code in barcodes_list:
        return code, 0

    min_ed = high_value
    cor_code = 0
    for bc in barcodes_list:
        hamm_d = hamming_dist_bin(code, bc)
        if hamm_d == 1:
            min_ed = 1
            cor_code = bc
            break
        if hamm_d < min_ed:
            min_ed = hamm_d
            cor_code = bc

    return cor_code, min_ed
        

def pass_filter(read, filters_counter, required_poly_t):
    """
    return true if a read pass the filters.
    """
    phix_genes = np.array(range(1, 7)) * 111111111
    N = BinRep._str2bindict['N']

    if read['gene'] == 0:
        filters_counter['gene_0'] += 1
        return False
    if read['gene'] in phix_genes:
        filters_counter['phi_x'] += 1
        return False
    # todo: collapse cell_0 and rmt_0 into one dictionary key
    if read['cell'] == 0:
        filters_counter['cell_0'] += 1
        return False
    if read['rmt'] == 0:
        filters_counter['rmt_0'] += 1
        return False
    if read['n_poly_t'] < required_poly_t:
        filters_counter['poly_t'] += 1
        return False
    if BinRep.contains(int(read['rmt']), N):
        filters_counter['rmt_N'] += 1
        return False
    return True


def group_for_dropseq(ra, required_poly_t=4):
    """
    Prepare the RA for the dropSeq error correction. Apply filters and group by gene/cell
    :param err_correction_mat:
    :param max_ed:
    :param reverse_complement:
    :param required_poly_t:
    :param barcode_files:
    :param ra:
    """
    res = {}
    tot = 0
    filtered = {'gene_0': 0, 'phi_x': 0, 'cell_0': 0, 'rmt_0': 0, 'poly_t': 0, 'rmt_N': 0}

    for i, v in enumerate(ra.data):
        tot += 1
        if not pass_filter(v, filtered, required_poly_t):
            continue
        
        # Build a dictionary of {cell11b:{rmt1:{gene1,cell1:#reads,
        # gene2,cell2:#reads},rmt2:{...}},cell211b:{...}}
        cell = v['cell']
        cell_header = cell >> 3   # we only look at the first 11 bases
        rmt = v['rmt']
        gene = v['gene']
        try:
            res[cell_header][rmt][gene,cell] += 1
        except KeyError:
            try:
                res[cell_header][rmt][gene,cell] = 1
            except KeyError:
                try:
                    res[cell_header][rmt]={(gene,cell):1}
                except KeyError:
                    res[cell_header] = {rmt:{(gene,cell):1}}

    seqc.log.info('Error correction filtering results: total reads: {}; did not pass '
                  'preliminary filters: {}'.format(tot, sum(filtered.values())))
    return res, filtered


def drop_seq(alignments_ra, *args, **kwargs):
    """pass-through function that groups the read_array for count matrix construction but
    does not perform any error correction. To be replaced by Ashley's function.

    :param alignments_ra: seqc.core.ReadArray object
    :param args:
    :param kwargs:
    :return:
    """
    
    start = time.process_time()
    pos_filter_threshold = 0.8
    barcode_base_shift_threshold = 0.9
    umi_len = 8
    min_umi_cutoff = 10
    # TODO:
    # 4. collapse UMI's that are 1 ED from one another (This requires more thought as there are plenty of edge cases)
    
    res_dic = {}
    
    if 'grouped_ra' in kwargs:
        grouped_ra = kwargs['grouped_ra']
        summary = None
    else:
        grouped_ra, summary = group_for_dropseq(alignments_ra)
    if 'min_umi_cutoff' in kwargs:
        min_umi_cutoff = kwargs['min_umi_cutoff']
        
    base_shift_count = 0
    pos_bias_count = 0
    small_cell_groups = 0
    
    #close_pairs = 0
    
    for cell in grouped_ra:
        retain = True
        base_shift = False
        correct_cell = cell
        size = len(grouped_ra[cell])
        
        #close_pairs += count_close_umi(grouped_ra[cell])
        
        #Check minimum size of cell group
        if size < min_umi_cutoff:   #if the number of UMI's for this cell don't meet the threshold, we have to retain them
            retain = True
            small_cell_groups += 1
            
        else:
            # Check for bias in any single base of the UMI (1-7)
            base_mat = base_count(grouped_ra[cell])
            for pos in range(umi_len-1):
                for base in base_mat:
                    if base_mat[base][pos]>pos_filter_threshold:
                        retain = False
                        pos_bias_count += 1
                        continue
            
            # Check for base shift
            if base_mat['T'][-1] > barcode_base_shift_threshold:
                retain = True
                base_shift = True
                # Retain, but insert an 'N'
                base_shift_count += 1            
                
        if retain:
            for rmt in grouped_ra[cell]:
                for gene, full_cell in grouped_ra[cell][rmt]:
                    if base_shift:
                        correct_cell = BinRep.ints2int([cell, BinRep._str2bindict['N']])    #replace the last base of cell with 'N'
                        # todo next line unnecessary now that rmts are not being stored
                        correct_rmt = BinRep.ints2int([full_cell&0b111, rmt]) >> 3         # correct the rmt with the last base of the cell, and dump the last 'T'
                    else:
                        correct_cell = full_cell
                        # todo next line unnecessary now that rmts are not being stored
                        correct_rmt = rmt
                    try:
                        res_dic[gene, correct_cell][0] += 1
                        res_dic[gene, correct_cell][1] += (
                            grouped_ra[cell][rmt][gene, full_cell])

                    except KeyError:
                        res_dic[gene, correct_cell] = [
                            1, grouped_ra[cell][rmt][gene, full_cell]]

    seqc.log.info('base shift: {}, pos_bias: {}, small cell groups: {}'.format(base_shift_count, pos_bias_count, small_cell_groups))
    #print('base shift: {}, pos_bias: {}, small cell groups: {}, close pairs: {}'.format(base_shift_count, pos_bias_count, small_cell_groups, close_pairs))
    tot_time=time.process_time()-start
    #print('tot time: {}'.format(tot_time))
    return res_dic, grouped_ra, summary

def base_count(seq_dic, umi_len=8):
    count_mat = {'A':np.zeros(umi_len), 'C':np.zeros(umi_len), 'G':np.zeros(umi_len), 'T':np.zeros(umi_len)}
    for seq in seq_dic:
        for i, base in enumerate(BinRep.bin2str(seq)):
            if base=='N':
              continue
            count_mat[base][i] += 1
    tot = len(seq_dic)
    for base in count_mat:
        count_mat[base] = count_mat[base]/tot
    return count_mat

#Used for research only
def count_close_umi(seq_dic):
    count=0
    for seq1 in seq_dic:
        for seq2 in seq_dic:
            if hamming_dist_bin(seq1, seq2) <= 1:
                count+=1
    return (count-len(seq_dic))/2


def mars1_seq(*args, **kwargs):
    """very simple pass-through wrapper for mars1_seq error correction, needs to be
    replaced with the actual error correction method instead of just drop_seq()

    :param args:
    :param kwargs:
    :return:
    """
    return in_drop(*args, **kwargs)


def mars2_seq(*args, **kwargs):
    """very simple pass-through wrapper for mars1_seq error correction, designed so that
    getattr() on seqc.correct_errors will find the correct error function for mars2_seq

    :param args:
    :param kwargs:
    :return:
    """
    return in_drop(*args, **kwargs)


def in_drop_v2(*args, **kwargs):
    """very simple pass-through wrapper for in_drop error correction, designed so that
    getattr() on seqc.correct_errors will find the correct error function for in_drop_v2

    :param args:
    :param kwargs:
    :return:
    """
    return in_drop(*args, **kwargs)


# TODO: check this. clean other ec methods, comments and prob_d_to_r. push.
def in_drop(alignments_ra, barcode_files=list(), apply_likelihood=True,
            reverse_complement=True, donor_cutoff=1, alpha=0.05,
            required_poly_t=1, max_ed=2, singleton_weight=1):
    """
    Recieve an RA and return a bool matrix of identified errors according to each
    method

    :param alignments_ra:
    :param barcode_files:
    :param apply_likelihood:
    :param reverse_complement:
    :param donor_cutoff:
    :param alpha:
    :param required_poly_t:
    :param max_ed:
    :return:
    """

    res_time_cnt = {}
    err_correction_res = ''#np.zeros((len(alignments_ra), NUM_OF_ERROR_CORRECTION_METHODS))
    ra_grouped, error_rate, summary = prepare_for_ec(
            alignments_ra, barcode_files, required_poly_t, reverse_complement, max_ed,
            err_correction_mat='')
    grouped_res_dic, error_count = correct_errors(
            alignments_ra, ra_grouped, error_rate, err_correction_res='', p_value=alpha,
            singleton_weight=singleton_weight)

    return grouped_res_dic, err_correction_res, summary


def correct_errors(ra, ra_grouped, err_rate, err_correction_res='', donor_cutoff=1,
                   p_value=0.05, apply_likelihood=True, singleton_weight=1, fname=''):
    """
    Calculate and correct errors in barcode sets

    :param ra:
    :param ra_grouped:
    :param err_rate:
    :param err_correction_res:
    :param donor_cutoff:
    :param p_value:
    :param apply_likelihood:
    :param fname:
    :return:
    """
    start = time.process_time()
    d = ra_grouped
    grouped_res_dic = {}
    error_count = 0

    tot_feats = len(ra_grouped)
    cur_f = 0
    N = BinRep._str2bindict['N']
    for_removal = []
    for feature, cell in d.keys():
        # sys.stdout.write('\r' + str(cur_f) + '/' + str(tot_feats) +
        #                  ' groups processed. (' + str((100 * cur_f) / tot_feats) + '%)')
        cur_f += 1
        if feature == 0:
            continue

        retained_rmts = 0
        retained_reads = 0

        for r_seq in d[feature, cell].keys():
            if BinRep.contains(r_seq, N):  # todo This throws out more than we want?
                continue

            gene = feature
            r_num_occurences = len(d[gene , cell][r_seq])
            r_pos_list = np.hstack(ra.data['position'][d[feature, cell][r_seq]])

            expected_errors = 0
            # actual_donors = 0
            p_val_not_err = 0
            # P(x<=r_num_occurences) Note: these are different distributions
            # p_val_not_correct = poisson.cdf(r_num_occurences, avg_rmt_size)
            jait = False
            for d_rmt in generate_close_seq(r_seq):
                try:
                    d_num_occurences = len(d[gene, cell][d_rmt])
                    # actual_donors += d_num_occurences
                except KeyError:
                    continue
                # This can cut running time but will have a less accurate likelihood model. It is currently not used.
                # if d_num_occurences<=donor_cutoff:
                #    continue

                p_dtr = prob_d_to_r_bin(d_rmt, r_seq, err_rate)
                expected_errors += d_num_occurences * p_dtr

                # do Jaitin
                if not jait:
                    d_pos_list = np.hstack(ra.data['position'][d[feature, cell][d_rmt]])
                    if set(r_pos_list).issubset(set(d_pos_list)):
                        #err_correction_res[ra_grouped[gene, cell][r_seq],
                        #                   ERROR_CORRECTION_jaitin] = 1
                        jait = True

            # P(x>=r_num_occurences)
            p_val_not_err = gammainc(r_num_occurences, expected_errors)
            # Used for error correction methods comparison
            #if p_val_not_err <= p_value:
            #    err_correction_res[ra_grouped[gene, cell][r_seq],
            #                       ERROR_CORRECTION_LIKELIHOOD_NOT_ERROR] = 1
            #elif jait:
            #    err_correction_res[ra_grouped[gene, cell][r_seq],
            #                       ERROR_CORRECTION_JAIT_LIKELIHOOD] = 1

            if apply_likelihood:
                if not jait or p_val_not_err <= p_value:
                    if r_num_occurences == 1:
                        retained_rmts += singleton_weight
                    else:
                        retained_rmts += 1
                    retained_reads += r_num_occurences
            elif not jait:
                if r_num_occurences == 1:
                    retained_rmts += singleton_weight
                else:
                    retained_rmts += 1
                retained_reads += r_num_occurences

        grouped_res_dic[feature, cell] = retained_rmts, retained_reads

    # print('\nLikelihood model error_count: ', error_count)
    tot_time=time.process_time()-start
    # seqc.log.info('total error correction runtime: {}'.format(tot_time))
    # f.close()
    return grouped_res_dic, error_count


def convert_to_matrix(counts_dictionary):
    """

    :param counts_dictionary:
    :return:
    """
    # Set up entries for sparse matrix
    cols = [k[0] for k in counts_dictionary.keys()]
    rows = [k[1] for k in counts_dictionary.keys()]

    molecules = np.array(list(v[0] for v in counts_dictionary.values()))
    reads = np.array(list(v[1] for v in counts_dictionary.values()))

    # Map row and cell to integer values for indexing
    unq_row = np.unique(rows)  # these are the ids for the new rows / cols of the array
    unq_col = np.unique(cols)
    row_map = dict(zip(unq_row, np.arange(unq_row.shape[0])))
    col_map = dict(zip(unq_col, np.arange(unq_col.shape[0])))
    row_ind = np.array([row_map[i] for i in rows])
    col_ind = np.array([col_map[i] for i in cols])

    # change dtype, set shape
    molecules = molecules.astype(np.uint32)
    reads = reads.astype(np.uint32)
    shape = (unq_row.shape[0], unq_col.shape[0])

    # return a sparse array
    read_coo = coo_matrix((reads, (row_ind, col_ind)), shape=shape, dtype=np.uint32)
    mol_coo = coo_matrix((molecules, (row_ind, col_ind)), shape=shape, dtype=np.uint32)
    return {'molecules': {'matrix': mol_coo, 'row_ids': unq_row, 'col_ids': unq_col},
            'reads': {'matrix': read_coo, 'row_ids': unq_row, 'col_ids': unq_col}}


# For research use only.
def plot_ed_dist(ra, iter):
    dist_dic = {'umi':{},'cell':{}}
    for i in range(iter):
        read1 = ra.data[random.randint(0,len(ra)-1)]
        read2 = ra.data[random.randint(0,len(ra)-1)]
        ed = hamming_dist_bin(read1['rmt'],read2['rmt'])
        try:
            dist_dic['umi'][ed] += 1
        except KeyError:
            dist_dic['umi'][ed] = 1
            
        ed = hamming_dist_bin(read1['cell'],read2['cell'])
        try:
            dist_dic['cell'][ed] += 1
        except KeyError:
            dist_dic['cell'][ed] = 1
    return dist_dic
