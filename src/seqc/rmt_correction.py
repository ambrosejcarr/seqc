from scipy.special import gammainc
from seqc.sequence.encodings import DNA3Bit
import numpy as np
from seqc import log
from seqc.read_array import ReadArray
import time
import pandas as pd
import multiprocessing_on_dill as multi
from multiprocessing_on_dill import Manager

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
            if isinstance(err_rate,float):
                p *= err_rate
            else:
                p *= err_rate[(d_seq & 0b111, r_seq & 0b111)]
        d_seq >>= 3
        r_seq >>= 3
    return p


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


def _correct_errors(indices_grouped_by_cells, ra, err_rate, p_value=0.05):
    """Calculate and correct errors in RMTs

    :param ra:
    :param err_rate: A table of the estimate base switch error rate, produced
      ReadArray.apply_barcode_correction().
    :param p_value: The p_value used for the likelihood method
    :return None:
    """
    
    # A shared array among the parallel processes
    # Each entry represent a molecule RMT barcode with
    # -1: can't be corrected or already valid, >=0: represent the donor RMT 
    manager = Manager()
    shared_rmt_array = manager.Array('i',[-1 for i in range(len(ra.data['rmt']))])

        
    # a method called by each process to correct RMT for each cell
    def _correct_errors_by_cell_group(cell_group):
        
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
                    # Number of occurrences
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
                    # Save the RMT donor
                    shared_rmt_array[rmt_groups[rmt]]=jaitin_donor
                    #ra.data['status'][rmt_groups[rmt]] |= ra.filter_codes['rmt_error']
                    #ra.data['rmt'][rmt_groups[rmt]] = jaitin_donor
    
    # create a pool of workers and let each work on each single cell
    p=multi.Pool(processes=multi.cpu_count())
    #p.map(partial(_correct_errors_by_cell_group, shared_bc_correction_array=shared_rmt_array), indices_grouped_by_cells)
    p.map(_correct_errors_by_cell_group,indices_grouped_by_cells)
    #for cell_group in indices_grouped_by_cells:
        #p.apply_async(_correct_errors_by_cell_group, cell_group)
    p.close()
    p.join()
        
    # iterate through each RMT and do correction 
    # if it is indicated by the shared_rmt_array
    for i in range(len(shared_rmt_array)):
        if shared_rmt_array[i]>=0:
            ra.data['rmt'][i]=shared_rmt_array[i]
            ra.data['status'][i]|= ra.filter_codes['rmt_error']
