from scipy.special import gammainc
from itertools import permutations
from sys import maxsize
from scipy.sparse import coo_matrix
from seqc.sequence.encodings import DNA3Bit
from seqc.sequence.barcodes import revcomp_bytes
import numpy as np
from seqc import log
import random
from seqc.read_array import ReadArray
import time


DEFAULT_BASE_CONVERTION_RATE = 0.02

def prepare_for_ec(ra):
    """
    Prepare the RA for error correction by grouping it according to cell,gene and rmt

    :param ra:
    """
    res = {}

    for i, v in enumerate(ra.data):
        if not ReadArray.is_active_read(v):
            continue

        gene = v['gene']
        cell = v['cell']
        rmt = v['rmt']
        try:
            res[gene, cell][rmt].append(i)
        except KeyError:
            try:
                res[gene, cell][rmt] = [i]
            except KeyError:
                res[gene, cell] = {}
                res[gene, cell][rmt] = [i]

    return res


def generate_close_seq(seq):
    """ Return a list of all sequences that are up to 2 hamm distance from seq
    :param seq:
    """
    res = []
    l = DNA3Bit.seq_len(seq)

    # generate all sequences that are dist 1
    for i in range(l):
        mask = 0b111 << (i * 3)
        cur_chr = (seq & mask) >> (i * 3)
        res += [seq & (~mask) | (new_chr << (i * 3))
                for new_chr in DNA3Bit._bin2strdict.keys() if new_chr != cur_chr]
    # generate all sequences that are dist 2
    for i in range(l):
        mask_i = 0b111 << (i * 3)
        chr_i = (seq & mask_i) >> (i * 3)
        for j in range(i + 1, l):
            mask_j = 0b111 << (j * 3)
            chr_j = (seq & mask_j) >> (j * 3)
            mask = mask_i | mask_j
            res += [seq & (~mask) | (new_chr_i << (i * 3)) | (new_chr_j << (j * 3)) for
                    new_chr_i in DNA3Bit._bin2strdict.keys() if new_chr_i != chr_i for
                    new_chr_j in DNA3Bit._bin2strdict.keys() if new_chr_j != chr_j]

    return res


def prob_d_to_r_bin(d_seq, r_seq, err_rate):
    """
    Return the probability of d_seq turning into r_seq based on the err_rate table
    (all binary)

    :param err_rate:
    :param r_seq:
    :param d_seq:
    """

    if DNA3Bit.seq_len(d_seq) != DNA3Bit.seq_len(r_seq):
        return 1

    p = 1.0
    while d_seq > 0:
        if d_seq & 0b111 != r_seq & 0b111:
            p *= err_rate[DNA3Bit.ints2int([d_seq & 0b111, r_seq & 0b111])]
        d_seq >>= 3
        r_seq >>= 3
    return p

def group_for_dropseq(ra):
    """
    Prepare the RA for the dropSeq error correction. Apply filters and group by gene/cell
    :param required_poly_t:
    :param max_dust:
    :param ra:
    """
    res = {}

    for i, v in enumerate(ra.data):
        if not ReadArray.is_active_read(v):
            continue
        
        # Build a dictionary of {cell11b:{rmt1:{gene1,cell1:#reads,
        # gene2,cell2:#reads},rmt2:{...}},cell211b:{...}}
        cell = v['cell']
        cell_header = cell >> 3   # we only look at the first 11 bases
        rmt = v['rmt']
        gene = v['gene']
        try:
            res[cell_header][rmt][gene, cell].append(i)
        except KeyError:
            try:
                res[cell_header][rmt][gene, cell] = [i]
            except KeyError:
                try:
                    res[cell_header][rmt] = {(gene, cell): [i]}
                except KeyError:
                    res[cell_header] = {rmt: {(gene, cell): [i]}}
    return res


def drop_seq(alignments_ra, *args, **kwargs):
    """pass-through function that groups the read_array for count matrix construction but
    does not perform any error correction. To be replaced by Ashley's function.

    :param alignments_ra: seqc.read_array.ReadArray object
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
    # 4. collapse UMI's that are 1 ED from one another (This requires more thought as
    # there are plenty of edge cases)
    
    res_dic = {}
    
    grouped_ra = group_for_dropseq(alignments_ra)
    if 'min_umi_cutoff' in kwargs:
        min_umi_cutoff = kwargs['min_umi_cutoff']
        
    base_shift_count = 0
    pos_bias_count = 0
    small_cell_groups = 0
    
    # close_pairs = 0
    
    for cell in grouped_ra:
        retain = True
        base_shift = False
        # correct_cell = cell
        size = len(grouped_ra[cell])
        # close_pairs += count_close_umi(grouped_ra[cell])
        
        # Check minimum size of cell group
        # if the number of UMI's for this cell don't meet the threshold, we have to
        # retain them
        if size < min_umi_cutoff:
            retain = True
            small_cell_groups += 1
            
        else:
            # Check for bias in any single base of the UMI (1-7)
            base_mat = base_count(grouped_ra[cell])
            for pos in range(umi_len-1):
                for base in base_mat:
                    if base_mat[base][pos] > pos_filter_threshold:
                        retain = False
                        pos_bias_count += 1
                        continue
            
            # Check for base shift
            if base_mat['T'][-1] > barcode_base_shift_threshold:
                # Retain, but insert an 'N'
                retain = True
                base_shift = True
                base_shift_count += 1            
                
        if retain:
            for rmt in grouped_ra[cell]:
                for gene, full_cell in grouped_ra[cell][rmt]:
                    if base_shift:
                        # replace the last base of cell with 'N'
                        correct_cell = DNA3Bit.ints2int([cell, DNA3Bit.encode(b'N')])
                        # todo next line unnecessary now that rmts are not being stored
                        # correct the rmt with the last base of cb, remove the last 'T'
                        # correct_rmt = BinRep.ints2int([full_cell & 0b111, rmt]) >> 3
                    else:
                        correct_cell = full_cell
                    for read_idx in grouped_ra[cell][rmt][gene, full_cell]:
                        alignments_ra.data[read_idx]['cell'] = correct_cell
        
        else:
            for rmt in grouped_ra[cell]:
                for gene, full_cell in grouped_ra[cell][rmt]:
                    for read_idx in grouped_ra[cell][rmt][gene, full_cell]:
                        ReadArray.set_error_status_on(alignments_ra.data[read_idx])     
                    
    log.info('base shift: {}, pos_bias: {}, small cell groups: {}'.format(base_shift_count, pos_bias_count, small_cell_groups))
    log.info('Error correction finished in {}'.format(time.process_time() - start))
    return

def base_count(seq_dic, umi_len=8):
    count_mat = {'A': np.zeros(umi_len), 'C': np.zeros(umi_len), 'G': np.zeros(umi_len),
                 'T': np.zeros(umi_len)}
    for seq in seq_dic:
        for i, base in enumerate(DNA3Bit.decode(seq).decode()):
            if base == 'N':
                continue
            count_mat[base][i] += 1
    tot = len(seq_dic)
    for base in count_mat:
        count_mat[base] /= tot
    return count_mat

def in_drop(alignments_ra, error_rate, apply_likelihood=True, alpha=0.05, singleton_weight=1):
    """
    Recieve an RA and return a bool matrix of identified errors according to each
    method

    :param alignments_ra:
    :param apply_likelihood:
    :param alpha:
    :return:
    """

    ra_grouped = prepare_for_ec(alignments_ra)
    grouped_res_dic = correct_errors(alignments_ra, ra_grouped, error_rate, p_value=alpha, apply_likelihood = apply_likelihood, singleton_weight=singleton_weight)

    return grouped_res_dic


def correct_errors(ra, ra_grouped, err_rate, p_value=0.05, apply_likelihood=True, singleton_weight=1):
    """
    Calculate and correct errors in barcode sets

    :param ra:
    :param ra_grouped:
    :param err_rate: A table of the estimate base switch error rate, produced ReadArray.applu_barcode_correction().
    :param p_value: The p_value used for the likelihood method
    :param apply_likelihood: apply the likelihood method in addition to the Jaitin et al method (Allons)
    :return:
    """
    start = time.process_time()
    d = ra_grouped
    grouped_res_dic = {}
    rmt_N = 0
    N = DNA3Bit.encode(b'N')

    for feature, cell in d.keys():
        retained_rmts = 0.0
        retained_reads = 0

        for r_seq in d[feature, cell].keys():
            if DNA3Bit.contains(r_seq, N):  # todo This throws out more than we want?
                rmt_N += 1
                continue

            gene = feature
            r_num_occurences = len(d[gene , cell][r_seq])
            r_pos_list = np.hstack(ra.data['position'][d[feature, cell][r_seq]])

            expected_errors = 0
            jait = False
            for d_rmt in generate_close_seq(r_seq):
                try:
                    d_num_occurences = len(d[gene, cell][d_rmt])
                except KeyError:
                    continue
                # This can cut running time but will have a less accurate likelihood
                # model. It is currently not used.
                # if d_num_occurences<=donor_cutoff:
                #    continue

                p_dtr = prob_d_to_r_bin(d_rmt, r_seq, err_rate)
                expected_errors += d_num_occurences * p_dtr

                # do Jaitin
                if not jait:
                    d_pos_list = np.hstack(ra.data['position'][d[feature, cell][d_rmt]])
                    if set(r_pos_list).issubset(set(d_pos_list)):
                        # err_correction_res[ra_grouped[gene, cell][r_seq],
                        #                    ERROR_CORRECTION_jaitin] = 1
                        jait = True

            # P(x>=r_num_occurences)
            p_val_not_err = gammainc(r_num_occurences, expected_errors)

            if apply_likelihood:
                if not jait or p_val_not_err <= p_value:
                    #not error
                    if r_num_occurences == 1:
                        retained_rmts += singleton_weight
                    else:
                        retained_rmts += 1
                    retained_reads += r_num_occurences
                else:
                    #error
                    for ind in d[gene, cell][r_seq]:
                        ReadArray.set_error_status_on(ra.data[i])
            elif not jait:
                #not error
                if r_num_occurences == 1:
                    retained_rmts += singleton_weight
                else:
                    retained_rmts += 1
                retained_reads += r_num_occurences
            else:
                #error
                for ind in d[gene, cell][r_seq]:
                    ReadArray.set_error_status_on(ra.data[i])

        grouped_res_dic[feature, cell] = retained_rmts, retained_reads
    
    log.info('Molecules with N that were ignored: {}'.format(rmt_N))
    tot_time=time.process_time() - start
    log.info('total error correction runtime: {}'.format(tot_time))    
    
    return grouped_res_dic
