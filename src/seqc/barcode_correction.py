import numpy as np
from seqc.sequence.encodings import DNA3Bit
import pandas as pd
from itertools import permutations
import seqc.sequence.barcodes
from seqc import log

# todo document me
def in_drop(ra, platform, barcode_files, max_ed=2,
            default_error_rate=0.02):
    """
    Correct reads with incorrect barcodes according to the correct barcodes files.
    Reads with barcodes that have too many errors are filtered out.
    """

    # Read the barcodes into lists
    valid_barcodes = []
    for barcode_file in barcode_files:
        with open(barcode_file, 'r') as f:
            valid_barcodes.append(set(DNA3Bit.encode(line.strip()) for line in
                                      f.readlines()))
    
    # Containers         
    num_barcodes = platform.num_barcodes
    correct = [None] * num_barcodes
    edit_dist = [None] * num_barcodes

    # Error table container
    errors = [p for p in permutations(DNA3Bit.bin2strdict.keys(), r=2)]
    error_table = dict(zip(errors, np.zeros(len(errors))))
    cor_instance_table = dict(zip(DNA3Bit.bin2strdict.keys(),
                                  np.zeros(len(DNA3Bit.bin2strdict))))
    
    # Check if the barcode has to be an exact match
    exact_match = False
    if max_ed == 0:
        exact_match = True

    # Group reads by cells
    indices_grouped_by_cells = ra.group_indices_by_cell(multimapping=True)

    for inds in indices_grouped_by_cells:

        # Extract barcodes for one of the reads
        barcodes = platform.extract_barcodes(ra.data['cell'][inds[0]])

        # Identify correct barcode
        for i in range(num_barcodes):
            correct[i], edit_dist[i] = seqc.sequence.barcodes.find_correct_barcode(
                barcodes[i], valid_barcodes[i], exact_match)

        # 1. If all edit distances are 0, barcodes are correct,
        #    update the correct instance table
        # 2. Correct any barcodes within permissible edit distance,
        #    update the correct instance table for non-errored bases,
        #    update error table for the errored bases
        # 3. Mark the uncorrectable barcodes as cell errors

        if all(np.array(edit_dist) == 0):
            # Temp container to increment the correct instance counter
            tmp_bc = DNA3Bit.ints2int(barcodes)
            while tmp_bc > 0:
                cor_instance_table[tmp_bc & 0b111] += 1
                tmp_bc >>= 3

        elif max(edit_dist) > max_ed:
            ra.data['status'][inds] |= ra.filter_codes['cell_error']
            continue

        else:
            # These barcodes can be corrected, Count the number of correct bases
            # Update the error table if there was only one error across the barcodes                
            tmp_bc = DNA3Bit.ints2int(barcodes)
            tmp_cor = DNA3Bit.ints2int(correct)

            # Update the read array with the correct barcode
            ra.data['cell'][inds] = tmp_cor

            # Iterating through the sequences
            while tmp_bc > 0:
                if tmp_bc & 0b111 == tmp_cor & 0b111:
                    cor_instance_table[tmp_bc & 0b111] += 1
                elif sum(edit_dist) == 1:
                    error_table[(tmp_cor & 0b111, tmp_bc & 0b111)] += 1
                tmp_bc >>= 3
                tmp_cor >>= 3

    # Create error rate table
    if sum(error_table.values()) == 0:
        log.info('No errors were detected or barcodes do not support error '
                 'correction, using %f uniform error chance.' % default_error_rate)
        err_rate = dict(zip(errors, [default_error_rate] * len(errors)))
    # todo @Manu bug here, we're always setting the error rate even if there are
    # no detected errors. should the following line be in an "else" clause?
    err_rate = dict(zip(errors, [0.0] * len(errors)))
    for k, v in error_table.items():
        if DNA3Bit.decode(k[0]) in b'Nn':
            continue
        try:
            err_rate[k] = v / (sum(n for err_type, n in error_table.items()
                                   if err_type[0] == k[0]) + cor_instance_table[k[0]])
        except ZeroDivisionError:
            log.info('Warning: too few reads to estimate error rate for %s, setting '
                     'default rate of %f' %
                     (str(DNA3Bit.decode(k)), default_error_rate))
            err_rate[k] = default_error_rate

    return err_rate


def drop_seq(ra, min_rmt_cutoff=10, rmt_error_frequency=0.8, barcode_base_shift_threshold=0.9):

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

    :param ra: seqc.read_array.ReadArray object
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


