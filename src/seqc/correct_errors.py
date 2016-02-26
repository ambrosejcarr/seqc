import random
from scipy.special import gammaincinv
from itertools import permutations
from sys import maxsize
import time
from three_bit import ThreeBit as bin_rep       #Using Ambrose coder for unity
import numpy as np
import seqc
import sys

high_value = maxsize  # Used for sorting, needs to be longer than any sequence  

NUM_OF_ERROR_CORRECTION_METHODS = 4
ERROR_CORRECTION_AJC = 0
ERROR_CORRECTION_STEN = 1
ERROR_CORRECTION_FILTERS = 1        # recycling the column in the matrix to save room
ERROR_CORRECTION_jaitin = 2
ERROR_CORRECTION_ALLON = 3

def group_for_ec(ra, required_poly_t=1):
    res = {}
    for i, v in enumerate(ra.data):
        if len(ra.features[i]) != 1:
            continue
            
        seq = bin_rep.ints2int([int(v['cell']),int(v['rmt'])])
        gene = ra.features[i][0]
        try:
            res[gene][seq] += 1
        except KeyError:
            try:
                res[gene][seq] = 1
            except KeyError:
                res[gene]={}
                res[gene][seq] = 1
    return res
    
# Return the hamming distance bewteen two numbers representing a sequence (3 bits per base)
def hamming_dist_bin(c1, c2):
    if (bin_rep.seq_len(c1) != bin_rep.seq_len(c2)):
        return high_value
    d=0
    while c1>0:
        if c1&0b111 != c2&0b111:
            d+=1
        c1>>=3
        c2>>=3
    return d
    
def generate_close_seq(seq):
    ''' Return a list of all sequences that are up to 2 hamm distance from seq'''
    res = []
    l = bin_rep.seq_len(seq)
    
    #genereta all sequences that are dist 1
    for i in range(l):
        mask = 0b111<<(i*3)
        cur_chr = (seq&mask)>>(i*3)
        res += [seq&(~mask)|(new_chr<<(i*3)) for new_chr in bin_rep.bin_nums if new_chr!=cur_chr]
    #genereta all sequences that are dist 2
    for i in range(l):
        mask_i = 0b111<<(i*3)
        chr_i = (seq&mask_i)>>(i*3)
        for j in range(i+1,l):
            mask_j = 0b111<<(j*3)
            chr_j = (seq&mask_j)>>(j*3)
            mask = mask_i|mask_j
            res += [seq&(~mask)|(new_chr_i<<(i*3))|(new_chr_j<<(j*3)) for new_chr_i in bin_rep.bin_nums if new_chr_i!=chr_i for new_chr_j in bin_rep.bin_nums if new_chr_j!=chr_j]
    
    return res 
    
def prob_d_to_r_bin(d_seq, r_seq, err_rate):
    '''Return the probability of d_seq turning into r_seq based on the err_rate table (all binary)'''
    
    if (bin_rep.seq_len(d_seq) != bin_rep.seq_len(r_seq)):
        return 1
    
    p=1.0
    while d_seq>0:
        if d_seq&0b111 != r_seq&0b111:
            p*= err_rate[bin_rep.ints2int([d_seq&0b111, r_seq&0b111])]
        d_seq>>=3
        r_seq>>=3
    return p
    
complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}  # for reverse complement
def rev_comp(s):
    """Return the reverse complement of an ACGT string s"""
    ret = []
    for i in range(len(s)):
        ret.append(complement[(s[-1-i])])
    return ''.join(ret)
            

def estimate_error_rate(barcode_files, grouped_ra, reverse_complement=True):
    '''
    Estimate the error rate based on the barcodes in the data and the correct barcodes in the barcode file.
    Return an error_rate table.
    ''' 
    correct_barcodes = []
    if reverse_complement:
        for barcode_file in barcode_files:
            with open(barcode_file) as f:
                correct_barcodes.append(set(bin_rep.str2bin(rev_comp(line.strip()))
                                            for line in f.readlines()))
    else:
        for barcode_file in barcode_files:
            with open(barcode_file) as f:
                correct_barcodes.append(set(bin_rep.str2bin(line.strip())
                                            for line in f.readlines()))

    # go over the sequences in the file to calculate the error rate
    correct_instances = 0
    errors = list(bin_rep.ints2int([p[0],p[1]]) for p in permutations(bin_rep.bin_nums, r=2))
    error_table = dict(zip(errors, [0] * len(errors)))
    
    N = bin_rep._str2bindict['N']     
    
    
    dynamic_codes_table = {}
    #tot = len(ra.data)
    tot = 0
    print('\n')
    ignored=0
    repeated = 0
    new = 0
    #for i, read in enumerate(ra.data):
    #for read in ra.data:
    for gene in grouped_ra.keys():
        tot+=len(grouped_ra[gene])
    #    if i%100000==0:
    #        sys.stdout.write('\rDoing read: '+str(i)+'/'+str(tot)+' ('+str(i/tot)+'%)')
        for seq in grouped_ra[gene].keys():
                    
            #if bin_rep.contains(int(read['cell']), N):
            if bin_rep.contains(seq, N):
                    ignored+=1
                    continue
                    
            #c1 = bin_rep.c1_from_codes(int(read['cell']))
            #c2 = bin_rep.c2_from_codes(int(read['cell']))
            c1 = bin_rep.c1_from_int(seq)
            c2 = bin_rep.c2_from_int(seq)
            #print(str(c1), str(c2))
            
            try:
                #print('a.')
                cor_c1, err_c1, ed_1 = dynamic_codes_table[c1]
                repeated+=1
            except KeyError:
                #print('b')
                new+=1
                cor_c1, err_c1, ed_1 = find_correct_barcode(c1, correct_barcodes[0])
                dynamic_codes_table[c1] = cor_c1, err_c1, ed_1
            try:
                #print('c')
                cor_c2, err_c2, ed_2 = dynamic_codes_table[c2]
                repeated+=1
            except KeyError:
                #print ('d')
                #return correct_barcodes
                cor_c2, err_c2, ed_2 = find_correct_barcode(c2, correct_barcodes[1])
                #print('d1')
                dynamic_codes_table[c2] = cor_c2, err_c2, ed_2
                new+=1
            
            if ed_1+ed_2 == 0:
                #print('e')
                correct_instances += ((bin_rep.seq_len(c1) + bin_rep.seq_len(c2)) / 4)*len(grouped_ra[gene][seq])
            elif ed_1+ed_2 > 1:
                #print('f')
                continue    # Ignore codes that are too noisy in the error correction calculation
            else:
                try:
                    if ed_1 == 1:
                        #print('g')
                        error_table[err_c1] += len(grouped_ra[gene][seq])
                    elif ed_2 == 1:
                        #print('h')
                        error_table[err_c2] += len(grouped_ra[gene][seq])
                except TypeError:   # TODO: The 'N' in the codes create a key error. A good idea might be to just ignore them or even filter them completely.
                    #print('i')
                    pass
                except KeyError:
                    pass
        

    
    # convert to error rates    
    default_error_rate = 0.02
    err_rate = dict(zip(errors, [0.0] * len(errors)))
    if sum(error_table.values()) == 0:
        print('No errors were detected, using %f uniform error chance.' % (
            default_error_rate))
        err_rate = dict(zip(errors, [default_error_rate] * len(errors)))
    for k, v in error_table.items():
        try:
            err_rate[k] = v / (sum(n for err_type, n in error_table.items()
                               if err_type&0b111000 == k&0b111000) + correct_instances)
        except ZeroDivisionError:
            print('Warning: too few reads to estimate error rate for %r '
                  'setting default rate of %f' % (k, default_error_rate))
            err_rate[k] = default_error_rate
    print('total reads: ',str(tot),', ignored: ',str(ignored), ', new entries: ',str(new),', repeated entries: ',str(repeated))
    return err_rate


def list_errors(code, correct_barcodes):
    """
    For a given code and a list of correct barcodes - find the correct barcode
    that is closest to code and return the list of errors that turned it into
    code. An error is a six bit int representing a two chr string of type "AG","CT", etc.
    """
    # find the closest correct barcode
    min_dist = high_value
    donor = 0
    for cor_code in correct_barcodes:
        hamm_d = hamming_dist_bin(code, cor_code)
        if hamm_d < min_dist:
            min_dist = hamm_d
            donor = cor_code
            if hamm_d == 1:
                break
    
    if donor==0:
        print('Error: no donor code was found to be closest. code = ', bin_rep.bin2str(code))
    # return the actual error
    err_list = []
    while code > 0:
        if code&0b111 != donor&0b111:
            err_list.append(bin_rep.ints2int([donor&0b111, code&0b111]))
        code>>=3
    return err_list
    
def find_correct_barcode(code, barcodes_list):
    """
    For a given barcode find the closest correct barcode to it from the list (limited to one ED), a string representing the error and the edit distance
    NOTE: for now this function looks for a barcode with ED==1 and does not bother looking for the minimum
    """
    if code in barcodes_list:
        return code, 0, 0
        
    found = False
    for cor_code in barcodes_list:
        hamm_d = hamming_dist_bin(code, cor_code)
        if hamm_d == 1:
            found = True
            break
        
    if not found:
        return '', 0, 10        # code is more than one ED away

    ret_code = cor_code
    # return the actual error
    while code > 0:
        if code&0b111 != cor_code&0b111:
            err = bin_rep.ints2int([cor_code&0b111, code&0b111])
        code>>=3
        cor_code>>=3
    return ret_code, err, 1
    
def correct_errors(alignments_ra, barcode_files = [], reverse_complement = True, donor_cutoff=1, P_VALUE=0.05, compare_methods = False):
    '''Recieve an RA and return a bool matrix of identified errors according to each method'''
    res_time_cnt = {}
    err_correction_res = np.zeros((alignments_ra.shape[0], NUM_OF_ERROR_CORRECTION_METHODS))
    ra_grouped = alignments_ra.group_for_error_correction(required_poly_t = 1)
    print ('doing error correction version 2.2 - rmt only, ignore gene 0 for ec only, er based on group dynamic')
    res_time_cnt[ERROR_CORRECTION_AJC] = correct_errors_AJC(alignments_ra, ra_grouped, err_correction_res, barcode_files, reverse_complement, donor_cutoff, P_VALUE)
    if compare_methods:
        #print ('doing Sten')
        #res_time_cnt[ERROR_CORRECTION_STEN] = correct_errors_sten(ra_grouped, err_correction_res)
        print ('doing jaitin')
        res_time_cnt[ERROR_CORRECTION_jaitin] = correct_errors_jaitin(alignments_ra, ra_grouped, err_correction_res)
        #print ('doing Allon')
        #res_time_cnt[ERROR_CORRECTION_ALLON] = correct_errors_allon(ra_grouped, err_correction_res, barcode_files, reverse_complement)
        
    return err_correction_res, res_time_cnt
    
        
def correct_errors_AJC(ra, ra_grouped, err_correction_res, barcode_files, reverse_complement=True, donor_cutoff=1, P_VALUE=0.1):
    """calculate and correct errors in barcode sets"""
    start = time.process_time()
    d = ra_grouped

    error_count = 0
    err_rate = estimate_error_rate(barcode_files, ra_grouped, reverse_complement)
    
    tot_feats = len(ra_grouped)
    cur_f = 0
    
    N = bin_rep._str2bindict['N']
    for_removal = []
    for feature in d.keys():
        sys.stdout.write('\r' + str(cur_f) + '/' + str(tot_feats) + ' features processed. ('+str((100*cur_f)/tot_feats)+'%)')
        cur_f += 1
        if feature==0:  
            continue
        for r_seq in d[feature].keys():
            if bin_rep.contains(r_seq, N):
                continue
            
            gene = feature
            r_c1 = bin_rep.c1_from_int(r_seq)
            r_c2 = bin_rep.c2_from_int(r_seq)
            #TODO: get the correct codes according to barcodes
            r_rmt = bin_rep.rmt_from_int(r_seq)
            r_num_occurences = d[gene][r_seq].shape[0]
            
            """
            We want to calculate the error rate (L) required such that we expect
            to see k observations of r_seq or more 5% of the time.

            This is:
               Pois(K >= k) = 0.1; P(r_seq == error)
             = Pois(K <= k) = 0.9

             Pois(K <= k) is the poisson cdf function, which is numerically
             approximated by the regularized incomplete gamma function:

            Pois(K <= k | L) = gammainc(k + 1, L), x >= 0
                             = 0, otherwise

            thus L = gammaincinv(k + 1, 0.9)
            
            Rami: after checking it should be L = gammaincinv(k + 1, 0.1)
            summary - threshold is the mean of the dist for which the chance of getting r_num_occurences is P_VALUE
            if the expected errors is bigger than that, than the chance of getting r_num_occurences by chance of error is higher than P_VALUE
            and so we need to regect r_seq as an error. lower P_VALUE will increase the number of errors (we require higher level of certainty)         
            """

            threshold = gammaincinv(r_num_occurences, P_VALUE)
            
            expected_errors = 0
            #for d_seq in generate_close_seq(r_seq):
            for d_rmt in generate_close_seq(r_rmt):
                d_seq = bin_rep.ints2int([r_c1,r_c2,d_rmt])
                try:
                    d_num_occurences = d[gene][d_seq].shape[0]
                except KeyError:
                    continue
                if d_num_occurences<=donor_cutoff:
                    continue
                #d_c1 = bin_rep.c1_from_int(d_seq)
                #d_c2 = bin_rep.c2_from_int(d_seq)
                #d_rmt = bin_rep.rmt_from_int(d_seq)

                #p_dtr = prob_d_to_r_bin(d_seq, r_seq, err_rate)
                p_dtr = prob_d_to_r_bin(d_rmt, r_rmt, err_rate)
                
                expected_errors += d_num_occurences * p_dtr
                if expected_errors > threshold:
                    for_removal.append((gene, r_seq))
                    error_count+=1
                    break
    for (gene, r_seq) in for_removal:
        err_correction_res[ra_grouped[gene][r_seq],[ERROR_CORRECTION_AJC]] = 1
    
    print ('\nAJC error_count: ', error_count)
    tot_time=time.process_time()-start
    print('total AJC error_correction runtime: ',tot_time)
    return error_count, tot_time
    
def correct_errors_sten(ra_grouped, err_correction_res):
    """Correct errors using the method in Sten's paper.
       Remove any molecule supported by only a single read"""
    
    # for python 3
    start = time.process_time() 

    # comment out for testing with from_thin_air
    d = ra_grouped
    
    #   print('version 6.0 - using readArray')

    error_count = 0
    N = bin_rep._str2bindict['N']
    for_removal = []

    for feature in d.keys():
        for r_seq in d[feature].keys():
            if bin_rep.contains(r_seq, N):
                continue

            gene = feature
            r_num_occurences = d[gene][r_seq].shape[0]

            if(r_num_occurences <= 1):
                for_removal.append((gene, r_seq))
                error_count+=1

    for (gene, r_seq) in for_removal:
        err_correction_res[ra_grouped[gene][r_seq],[ERROR_CORRECTION_STEN]] = 1        #TODO: check that this is the correct way to address the array
    print ('Sten error_count: ', error_count)
    tot_time=time.process_time()-start
    print('total Sten error_correction runtime: ',tot_time)
    return error_count, tot_time
    
def correct_errors_jaitin(alignment_ra, ra_grouped, err_correction_res):
    """Correct errors according to jaitin method """
    
    #sort reads by number of positions
    #go from least to most:
    #   go over all d_seq:
    #       if d_seq is 1 dist from r_seq and covers the same positions:
    #           remove r_seq.
    #
    # Note: the sorting isn't really important, it just saves time.
    
    """calculate and correct errors in barcode sets"""
    start = time.process_time()
    d = ra_grouped
#   print('version 6.0 - using readArray')

    error_count = 0
    
    tot_feats = len(ra_grouped)
    cur_f = 0
    
    N = bin_rep._str2bindict['N']
    for_removal = []
    tot=0
    for feature in d.keys():
        sys.stdout.write('\r' + str(cur_f) + '/' + str(tot_feats) + ' features processed. ('+str((100*cur_f)/tot_feats)+'%)')
        cur_f += 1
        if feature==0:
            continue
        cur_seq=0
        tot_seq=len(d[feature].keys())
        sorted_seq_l = sorted([(seq, len(set(np.hstack(alignment_ra.positions[d[feature][seq]])))) for seq in d[feature].keys()], key=lambda x:x[1])
        for idx, r_seqs in enumerate(sorted_seq_l):
            r_seq = r_seqs[0]
            cur_seq+=1
            sys.stdout.write('\rfeature: '+str(cur_f) + '/' + str(tot_feats) + ', seq: ' + str(cur_seq) + '/' + str(tot_seq))
            if bin_rep.contains(r_seq, N):
                continue
            
            gene = feature
            r_rmt = bin_rep.rmt_from_int(r_seq)
            r_pos_list = np.hstack(alignment_ra.positions[d[feature][r_seq]])

            for d_idx in range(idx-1, -1, -1):
                d_seq = sorted_seq_l[d_idx][0]
                d_rmt = bin_rep.rmt_from_int(d_seq)
                d_pos_list = np.hstack(alignment_ra.positions[d[feature][r_seq]])

                if hamming_dist_bin(r_rmt, d_rmt) == 1 and set(r_pos_list).issubset(set(d_pos_list)):
                    for_removal.append((gene, r_seq))
                    error_count+=1
                    break
    for (gene, r_seq) in for_removal:
        err_correction_res[ra_grouped[gene][r_seq],[ERROR_CORRECTION_jaitin]] = 1
    
    print ('\njaitin error_count: ', error_count)
    tot_time=time.process_time()-start
    print('total jaitin error_correction runtime: ',tot_time)
    return error_count, tot_time

def correct_errors_allon(ra_grouped, err_correction_res, barcode_files, reverse_complement=True):
    """Correct errors using the method in Allon paper.
       Compare barcodes to list and discard any reads that have more
       than two mismatches from all barcodes on the list"""

    start = time.process_time()
    d = ra_grouped
    
    error_count = 0
    N = bin_rep._str2bindict['N']
    for_removal = []
    
    tot_feats = len(ra_grouped)
    cur_f = 0

    #create barcode list from barcode_files
    correct_barcodes = []
    if reverse_complement:
        for barcode_file in barcode_files:
            with open(barcode_file) as f:
                correct_barcodes.append(set(bin_rep.str2bin(rev_comp(line.strip()))
                                            for line in f.readlines()))
    else:
        for barcode_file in barcode_files:
            with open(barcode_file) as f:
                correct_barcodes.append(set(bin_rep.str2bin(line.strip())
                                            for line in f.readlines()))
    for feature in d.keys():
        sys.stdout.write('\r' + str(cur_f) + '/' + str(tot_feats) + ' features processed. ('+str((100*cur_f)/tot_feats)+'%)')
        cur_f += 1
        for r_seq in d[feature].keys():
            if bin_rep.contains(r_seq, N):
                continue

            gene = feature
            r_c1 = bin_rep.c1_from_int(r_seq)
            r_c2 = bin_rep.c2_from_int(r_seq)
            r_err_cnt = 0
            if r_c1 not in correct_barcodes[0]:
                r_err_cnt += len(list_errors(r_c1, correct_barcodes[0]))
            if r_c2 not in correct_barcodes[1]:
                r_err_cnt += len(list_errors(r_c2, correct_barcodes[1]))
            if(r_err_cnt > 2):
                for_removal.append((gene, r_seq))
                error_count+=1

    for (gene, r_seq) in for_removal:
        err_correction_res[ra_grouped[gene][r_seq],[ERROR_CORRECTION_ALLON]] = 1

    print ('\nAllon error count: ', error_count)
    tot_time=time.process_time()-start
    print('total Allon error_correction runtime: ',tot_time)
    return error_count, tot_time



    