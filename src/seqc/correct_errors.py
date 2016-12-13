from scipy.special import gammainc
from seqc.sequence.encodings import DNA3Bit
import numpy as np
from seqc import log
from seqc.read_array import ReadArray
import time
import pandas as pd

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


# todo document me
def drop_seq(ra, *args, min_rmt_cutoff=10, rmt_error_frequency=0.8, barcode_base_shift_threshold=0.9):

    """Drop-seq barcode correction suggested by Ashley
    1. Barcodes can be truncated to 11 bases because of synthesis error. Therefore a single
       barcode can be potentially be split to 4 barcodes
       Solution: Fix barcode: At the 8th position of RMT, if the fraction of T > 80%,
                 replace the 12th position of the cell barcode with N
                 Fix RMT: Remove the T in the last position of the RMT and prepend the 
                 first base from the uncorrected cell barcode
    2. If a particular base dominates any of the positions of the RMT,
       remove that cell barcode
    3. TODO: Primer match

    :param alignments_ra: seqc.read_array.ReadArray object
    :param min_rmt_cutoff: Minimum number of RMTs to apply barcode correction
    :param rmt_error_frequency: If a base appears with this frequency across the RMTs associated with the barcode
           in any position, the barcode is removed
    :param barcode_base_shift_threshold: Thresholds for detecting barcode shift 
    :return:
    """
    
    # Cell header [First 11 bases only - this should be parametrized]
    cell_header = ra.data['cell'] >> 3
    idx = np.argsort( cell_header )
    # Active reads
    passing = ra.data['status'][idx] == 0
    idx = idx[passing]
    breaks = np.where(np.diff(cell_header[idx]))[0] + 1
    indices_grouped_by_cell_headers = np.split(idx, breaks)

    # RMT length
    rmt_length = DNA3Bit.seq_len( ra.data['rmt'][idx[0]] )

    # 1. Barcode synthesis errors
    for header_group in indices_grouped_by_cell_headers:

        # RMT set
        # todo this could potentially be used in RMT correction / barcode correction in indrop
        all_rmts = list(set(ra.data['rmt'][header_group]))
        if len(all_rmts) < min_rmt_cutoff:
            continue

        # Count Ts in the last RMT position
        nuc_counts = dict(zip(DNA3Bit.bin2strdict.keys(), np.zeros(len(DNA3Bit.bin2strdict))))
        for rmt in all_rmts:
            nuc_counts[rmt & 0b0111] += 1

        # Correct cell barcode if necessary
        if nuc_counts[DNA3Bit.str2bindict['T']] > barcode_base_shift_threshold * len(all_rmts):

            # Correct the RMTs [This needs to done for each cell/RMT combination]
            idx = header_group[np.argsort(ra.data['cell'][header_group])]
            breaks = np.where(np.diff(ra.data['cell'][idx]))[0] + 1
            cell_groups = np.split(idx, breaks)

            for cell_group in cell_groups:
                last_base = ra.data['cell'][cell_group[0]] & 0b111

                # Correct the RMTs
                idx = cell_group[np.argsort(ra.data['rmt'][cell_group])]
                breaks = np.where(np.diff(ra.data['rmt'][cell_group]))[0] + 1
                rmt_groups = np.split(idx, breaks)

                for rmt_group in rmt_groups:
                    # Skip the last base
                    new_rmt = ra.data['rmt'][rmt_group[0]] >> 3
                    # Get the last base from the cell barcode
                    new_rmt = DNA3Bit.ints2int([last_base, new_rmt ])
                    ra.data['rmt'][rmt_group] = new_rmt

            # Append N to the cell header
            correct_barcode = DNA3Bit.ints2int([cell_header[header_group[0]], DNA3Bit.str2bindict['N']])
            ra.data['cell'][header_group] = correct_barcode

    # 2. Single UMI error
    indices_grouped_by_cells = ra.group_indices_by_cell()
    for cell_group in indices_grouped_by_cells:

        # RMT set
        all_rmts = list(set(ra.data['rmt'][cell_group]))
        if len(all_rmts) < min_rmt_cutoff:
            continue

        # RMT nucleotide frequency per position
        base_frequencies = dict()
        for i in DNA3Bit.bin2strdict.keys():
            base_frequencies[i] = np.zeros(rmt_length)
        for i in range(len(all_rmts)):
            rmt = all_rmts[i]
            position = rmt_length-1
            while rmt > 0:
                base_frequencies[rmt & 0b111][position] += 1
                rmt >>= 3
                position -= 1

        # Chuck N
        base_frequencies = pd.DataFrame(base_frequencies).T
        base_frequencies.ix[DNA3Bit.str2bindict['N']] = 0

        # Identify incorrect UMIs
        if any( base_frequencies.iloc[:, 0:(rmt_length-1)].max() > rmt_error_frequency * len(all_rmts)):
            ra.data['status'][cell_group] |= ra.filter_codes['cell_error']


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
