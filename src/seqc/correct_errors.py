from scipy.special import gammainc
from seqc.sequence.encodings import DNA3Bit
import numpy as np
from seqc import log
from seqc.read_array import ReadArray
import time

# todo PROPOSAL: RENAME MODULE TO rmt_correction


# todo document me
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
                for new_chr in DNA3Bit.bin2strdict.keys() if new_chr != cur_chr]
    # generate all sequences that are dist 2
    for i in range(l):
        mask_i = 0b111 << (i * 3)
        chr_i = (seq & mask_i) >> (i * 3)
        for j in range(i + 1, l):
            mask_j = 0b111 << (j * 3)
            chr_j = (seq & mask_j) >> (j * 3)
            mask = mask_i | mask_j
            res += [seq & (~mask) | (new_chr_i << (i * 3)) | (new_chr_j << (j * 3)) for
                    new_chr_i in DNA3Bit.bin2strdict.keys() if new_chr_i != chr_i for
                    new_chr_j in DNA3Bit.bin2strdict.keys() if new_chr_j != chr_j]

    return res


# todo document me
def probability_for_convert_d_to_r(d_seq, r_seq, err_rate):
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
            p *= err_rate[(d_seq & 0b111, r_seq & 0b111)]
        d_seq >>= 3
        r_seq >>= 3
    return p

# def group_for_dropseq(ra):
#     """
#
#     This function is likely deprecated and will be replaced by proper grouping method
#     based on array sorting.
#
#     Prepare the RA for the dropSeq error correction. Apply filters and group by gene/cell
#     :param seqc.read_array.ReadArray ra:  ReadArray object to be grouped
#     """
#     res = {}
#
#     for i, v in enumerate(ra.data):
#         if not ReadArray.is_active_read(v):
#             continue
#
#         # Build a dictionary of {cell11b:{rmt1:{gene1,cell1:#reads,
#         # gene2,cell2:#reads},rmt2:{...}},cell211b:{...}}
#         cell = v['cell']
#         cell_header = cell >> 3   # we only look at the first 11 bases
#         rmt = v['rmt']
#         gene = v['gene']
#         try:
#             res[cell_header][rmt][gene, cell].append(i)
#         except KeyError:
#             try:
#                 res[cell_header][rmt][gene, cell] = [i]
#             except KeyError:
#                 try:
#                     res[cell_header][rmt] = {(gene, cell): [i]}
#                 except KeyError:
#                     res[cell_header] = {rmt: {(gene, cell): [i]}}
#     return res


# todo document me
def drop_seq(alignments_ra, *args, **kwargs):
    raise NotImplementedError  # todo port to use gene/position csr

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
    # todo:
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

# todo document me
def base_count(seq_dic, umi_len=8):
    """

    :param seq_dic:
    :param umi_len:
    :return:
    """
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


# todo document me
def in_drop(ra, error_rate, alpha=0.05):
    """ Tag any RMT errors

    :param ra: Read array
    :param error_rate: Sequencing error rate determined during barcode correction
    :param alpha: Tolerance for errors
    """

    # Group read array by cells and submit indices of each cell to a different processor
    indices_grouped_by_cells = ra.group_indices_by_cell()
    _correct_errors(indices_grouped_by_cells, ra, error_rate, alpha)


# todo document me
def _correct_errors(indices_grouped_by_cells, ra, err_rate, p_value=0.05):
    """Calculate and correct errors in RMTs

    :param ra:
    :param err_rate: A table of the estimate base switch error rate, produced
      ReadArray.apply_barcode_correction().
    :param p_value: The p_value used for the likelihood method
    :return None:
    """

    for cell_group in indices_grouped_by_cells:

        # Breaks for each gene
        gene_inds = cell_group[np.argsort(ra.genes[cell_group])]
        breaks = np.where(np.diff(ra.genes[gene_inds]))[0] + 1
        splits = np.split(gene_inds, breaks)

        for inds in splits:
            # RMT groups
            rmt_groups = {}
            for ind in inds:
                rmt = ra.data['rmt'][ind]
                try:
                    rmt_groups[rmt].append(ind)
                except KeyError:
                    rmt_groups[rmt] = [ind]

            if len(rmt_groups) == 1:
                continue

            # This logic retains RMTs with N if no donor is found and contributes to the
            # molecule count
            for rmt in rmt_groups.keys():

                # Enumerate all possible RMTs with hamming distances 1 and/or 2
                # to build a probablitiy that this particular RMT was not an error
                # Simulatenously, check if Jaitin error correction can be applied
                jaitin_corrected = False
                expected_errors = 0
                for donor_rmt in generate_close_seq(rmt):

                    # Check if donor is detected
                    try:
                        donor_count = len(rmt_groups[donor_rmt])
                    except KeyError:
                        continue

                    # Build likelihood
                    # Probability of converting donor to target
                    p_dtr = probability_for_convert_d_to_r(donor_rmt, rmt, err_rate)
                    # Number of occurances
                    expected_errors += donor_count * p_dtr

                    # Check if jaitin correction is feasible
                    if not jaitin_corrected: 
                        ref_positions = ra.positions[rmt_groups[rmt]]
                        donor_positions = ra.positions[rmt_groups[donor_rmt]]

                        # Is reference a subset of the donor ? 
                        if (set(ref_positions)).issubset(donor_positions):
                            jaitin_corrected = True
                            jaitin_donor = donor_rmt

                # Probability that the RMT is an error
                p_val_err = gammainc(len(rmt_groups[rmt]), expected_errors)

                # Remove Jaitin corrected reads if probability of RMT == error is high
                if p_val_err > p_value and jaitin_corrected:
                    # Mark as error
                    ra.data['status'][rmt_groups[rmt]] |= ra.filter_codes['rmt_error']
                    ra.data['rmt'][rmt_groups[rmt]] = jaitin_donor
