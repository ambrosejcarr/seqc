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

NUM_OF_ERROR_CORRECTION_METHODS = 3
ERROR_CORRECTION_AJC = 0
#ERROR_CORRECTION_STEN = 1
ERROR_CORRECTION_BC_FILTERS = 1        # recycling the column in the matrix to save room
ERROR_CORRECTION_jaitin = 2
#ERROR_CORRECTION_ALLON = 3

    
#This is used just for running the Jaitin method for comparison.
def group_for_ec_pos(ra, scid_to_gene_map, required_poly_t=1):
    res = {}
    no_single_gene = 0
    cell_0 = 0
    rmt_0 = 0
    small_poly_t = 0
    tot = 0
    
    for i, v in enumerate(ra.data):
        tot += 1
        #Apply various perliminary filters on the data
        if ra.features[i][0] == 0:
            no_single_gene += 1
            continue
        if v['cell']==0:
            cell_0 += 1
            continue
        if v['rmt']==0:
            rmt_0 += 1
            continue
        if v['n_poly_t']<=required_poly_t:
            small_poly_t += 1
            continue
            
        
        #map the feature to a gene
        try:
            gene = scid_to_gene_map[ra.features[i][0]]
        except KeyError:
            print ('No gene mapped to feature ', ra.features[i][0], ' ignoring.')
            continue
        
        seq = bin_rep.ints2int([int(v['cell']),int(v['rmt'])])
        try:
            res[gene][seq].append(i)
        except KeyError:
            try:
                res[gene][seq] = [i]
            except KeyError:
                res[gene]={}
                res[gene][seq] = [i]
    print('total reads: ',tot,' no single gene: ',no_single_gene,' cell is zero: ',cell_0, ' rmt is zero: ',rmt_0,' small poly_t: ',small_poly_t)
    return res
    

def prepare_for_ec(ra, scid_to_gene_map, barcode_files, required_poly_t=1, reverse_complement=True, max_ed=2, err_correction_mat=''):
    ''' Prepare the RA for error correction. Apply filters, estimate error correction and correct the barcodes'''
    res = {}
    tot = 0
    filtered = 0
    bc_filter = 0
    
    N = bin_rep._str2bindict['N']
    
    dynamic_codes_table_c1 = {}
    dynamic_codes_table_c2 = {}
    
    errors = list(bin_rep.ints2int([p[0],p[1]]) for p in permutations(bin_rep.bin_nums, r=2))
    error_table = dict(zip(errors, [0] * len(errors)))
    cor_instance_table = {bin_rep._str2bindict['A']:0,
                        bin_rep._str2bindict['C']:0,
                        bin_rep._str2bindict['G']:0,
                        bin_rep._str2bindict['T']:0}
    
    #Read the barcodes into lists
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
    
    for i, v in enumerate(ra.data):
        tot += 1
        #Apply various perliminary filters on the data
        if ra.features[i][0] == 0:
            filtered += 1
            continue
        if v['cell']==0:
            filtered += 1
            continue
        if v['rmt']==0:
            filtered += 1
            continue
        if v['n_poly_t']<=required_poly_t:
            filtered += 1
            continue
        if bin_rep.contains(int(v['cell']), N):
            filtered += 1
            continue
        if bin_rep.contains(int(v['rmt']), N):
            filtered += 1
            continue
            
        
        #map the feature to a gene
        try:
            gene = scid_to_gene_map[ra.features[i][0]]
        except KeyError:
            print ('No gene mapped to feature ', ra.features[i][0], ' ignoring.')
            filtered += 1
            continue
        
        #correct and filter barcodes
        c1 = bin_rep.c1_from_codes(int(v['cell']))
        try:
            cor_c1, ed_1, err_l_1 = dynamic_codes_table_c1[c1]
        except KeyError:
            cor_c1, ed_1 = find_correct_barcode(c1, correct_barcodes[0])
            if ed_1 > max_ed:   #No point calculating errors on codes with more errors than allowed
                err_l_1 = None
            else:
                err_l_1 = list_errors(cor_c1, c1)
            dynamic_codes_table_c1[c1] = (cor_c1, ed_1, err_l_1)
            
        c2 = bin_rep.c2_from_codes(int(v['cell']))
        try:
            cor_c2, ed_2, err_l_2 = dynamic_codes_table_c2[c2]
        except KeyError:
            cor_c2, ed_2 = find_correct_barcode(c2, correct_barcodes[1])
            if ed_2 > max_ed:   #No point calculating errors on codes with more errors than allowed
                err_l_2 = None
            else:
                err_l_2 = list_errors(cor_c2, c2)
            dynamic_codes_table_c2[c2] = (cor_c2, ed_2, err_l_2)
        
        # Filter reads with too many errors in their barcodes
        if ed_1+ed_2 > max_ed:
            bc_filter += 1
            if err_correction_mat!='':
                err_correction_mat[i,ERROR_CORRECTION_BC_FILTERS] = 1
            continue
        
        for er in err_l_1 + err_l_2:
            error_table[er] += 1        
        
        #count non error bases
        tmp_c = bin_rep.ints2int([c1,c2])
        tmp_cor = bin_rep.ints2int([cor_c1,cor_c2])
        while tmp_c>0:
            if tmp_c&0b111 == tmp_cor&0b111:
                cor_instance_table[tmp_c&0b111]+=1
            tmp_c >>= 3
            tmp_cor >>= 3
            
        #group according to the correct barcodes and gene
        seq = bin_rep.ints2int([cor_c1,cor_c2,int(v['rmt'])])
        try:
            res[gene][seq].append(i)
        except KeyError:
            try:
                res[gene][seq] = [i]
            except KeyError:
                res[gene]={}
                res[gene][seq] = [i]
                
    # convert to error rates    #TODO: This is ambrose code, low priority to test it and /or replace it with mine.
    default_error_rate = 0.02
    err_rate = dict(zip(errors, [0.0] * len(errors)))
    if sum(error_table.values()) == 0:
        print('No errors were detected, using %f uniform error chance.' % (
            default_error_rate))
        err_rate = dict(zip(errors, [default_error_rate] * len(errors)))
    for k, v in error_table.items():
        try:
            err_rate[k] = v / (sum(n for err_type, n in error_table.items()
                               if err_type&0b111000 == k&0b111000) + cor_instance_table[(k&0b111000)>>3])
        except ZeroDivisionError:
            print('Warning: too few reads to estimate error rate for %r '
                  'setting default rate of %f' % (k, default_error_rate))
            err_rate[k] = default_error_rate            
    #print('total reads: ',tot,' no single gene: ',no_single_gene,' cell is zero: ',cell_0, ' rmt is zero: ',rmt_0,' small poly_t: ',small_poly_t)
    print('total reads: ',tot,' filtered: ',filtered, ' bc_filter: ', bc_filter)
    #print('error_table: ',error_table,' cor_instance_table: ',cor_instance_table)
    return res, err_rate

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
    
    
    dynamic_codes_table_c1 = {}
    dynamic_codes_table_c2 = {}
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
                cor_c1, err_c1, ed_1 = dynamic_codes_table_c1[c1]
                repeated+=1
            except KeyError:
                #print('b')
                new+=1
                cor_c1, err_c1, ed_1 = find_correct_barcode(c1, correct_barcodes[0])
                dynamic_codes_table_c1[c1] = cor_c1, err_c1, ed_1
            try:
                #print('c')
                cor_c2, err_c2, ed_2 = dynamic_codes_table_c2[c2]
                repeated+=1
            except KeyError:
                #print ('d')
                #return correct_barcodes
                cor_c2, err_c2, ed_2 = find_correct_barcode(c2, correct_barcodes[1])
                #print('d1')
                dynamic_codes_table_c2[c2] = cor_c2, err_c2, ed_2
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

#deprecated
#def list_errors(code, correct_barcodes):
    """
    For a given code and a list of correct barcodes - find the correct barcode
    that is closest to code and return the list of errors that turned it into
    code. An error is a six bit int representing a two chr string of type "AG","CT", etc.
    """
    # find the closest correct barcode
#    min_dist = high_value
#    donor = 0
#    for cor_code in correct_barcodes:
#        hamm_d = hamming_dist_bin(code, cor_code)
#        if hamm_d < min_dist:
#            min_dist = hamm_d
#            donor = cor_code
#            if hamm_d <= max_ed:
#                break
    
#    if donor==0:
#        print('Error: no donor code was found to be closest. code = ', bin_rep.bin2str(code))
    # return the actual error
#    err_list = []
#    while code > 0:
#        if code&0b111 != donor&0b111:
#            err_list.append(bin_rep.ints2int([donor&0b111, code&0b111]))
#        code>>=3
#    return err_list

def list_errors(s1, s2):
    """
    Return the list of nucleotide transformations that turn s1 to s2.
    An error is a six bit int representing a two chr string of type "AG","CT", etc.
    """
    
    # return the actual error
    err_list = []
    while s1 > 0:
        if s1&0b111 != s2&0b111:
            err_list.append(bin_rep.ints2int([s1&0b111, s2&0b111]))
        s1>>=3
        s2>>=3
    return err_list
    
def find_correct_barcode(code, barcodes_list):
    """
    For a given barcode find the closest correct barcode to it from the list (limited to one ED), a string representing the error and the edit distance
    NOTE: for now this function looks for a barcode with ED==1 and does not bother looking for the minimum
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

    
def correct_errors(alignments_ra, scid_to_gene_map, barcode_files = [], reverse_complement=True, donor_cutoff=1, alpha=0.05, AJC=True, jaitin=False, required_poly_t=1, max_ed=2):
    '''Recieve an RA and return a bool matrix of identified errors according to each method'''
    res_time_cnt = {}
    err_correction_res = np.zeros((len(alignments_ra), NUM_OF_ERROR_CORRECTION_METHODS))
    if AJC:
        print ('Starting likelihood model error detection')
        print ('Applying filterring and barcode correction...')
        ra_grouped, error_rate = prepare_for_ec(alignments_ra, scid_to_gene_map, barcode_files, required_poly_t, reverse_complement, max_ed, err_correction_res)
        print ('Applying likelihood model...')
        res_time_cnt[ERROR_CORRECTION_AJC] = correct_errors_AJC(ra_grouped, error_rate, err_correction_res, donor_cutoff, alpha)
    if jaitin:
        print ('Starting jaitin error detection')
        print ('Grouping by genes...')
        ra_grouped = group_for_ec_pos(alignments_ra, scid_to_gene_map, required_poly_t)
        print ('Applying error detection...')
        res_time_cnt[ERROR_CORRECTION_jaitin] = correct_errors_jaitin(alignments_ra, ra_grouped, err_correction_res)
        
    return err_correction_res, res_time_cnt
    
        
def correct_errors_AJC(ra_grouped, err_rate, err_correction_res, donor_cutoff=1, alpha=0.025):
    """calculate and correct errors in barcode sets"""
    start = time.process_time()
    d = ra_grouped

    error_count = 0
    
    tot_feats = len(ra_grouped)
    cur_f = 0
    
    N = bin_rep._str2bindict['N']
    for_removal = []
    for feature in d.keys():
        #sys.stdout.write('\r' + str(cur_f) + '/' + str(tot_feats) + ' features processed. ('+str((100*cur_f)/tot_feats)+'%)')
        cur_f += 1
        if feature==0:  
            continue
        
        for r_seq in d[feature].keys():
            if bin_rep.contains(r_seq, N):
                continue
                
        for r_seq in d[feature].keys():
            if bin_rep.contains(r_seq, N):
                continue
            
            gene = feature
            r_c1 = bin_rep.c1_from_int(r_seq)
            r_c2 = bin_rep.c2_from_int(r_seq)
            r_rmt = bin_rep.rmt_from_int(r_seq)
            r_num_occurences = len(d[gene][r_seq])


            threshold = gammaincinv(r_num_occurences, alpha)
            
            expected_errors = 0
            tot_donors = 0
            for d_rmt in generate_close_seq(r_rmt):
                d_seq = bin_rep.ints2int([r_c1,r_c2,d_rmt])
                try:
                    d_num_occurences = len(d[gene][d_seq])
                    tot_donors += d_num_occurences
                except KeyError:
                    continue
                if d_num_occurences<=donor_cutoff:
                    continue
                d_rmt = bin_rep.rmt_from_int(d_seq)

                p_dtr = prob_d_to_r_bin(d_rmt, r_rmt, err_rate)
                
                expected_errors += d_num_occurences * p_dtr
                if expected_errors > threshold:
                    for_removal.append((gene, r_seq))
                    error_count+=r_num_occurences
                    break
                    
    for (gene, r_seq) in for_removal:
        err_correction_res[ra_grouped[gene][r_seq],[ERROR_CORRECTION_AJC]] = 1
    
    print ('\nLikelihood model error_count: ', error_count)
    tot_time=time.process_time()-start
    print('total error_correction runtime: ',tot_time)
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

    error_count = 0
    
    tot_feats = len(ra_grouped)
    cur_f = 0
    
    N = bin_rep._str2bindict['N']
    for_removal = []
    tot=0
    for feature in d.keys():
        #sys.stdout.write('\r' + str(cur_f) + '/' + str(tot_feats) + ' features processed. ('+str((100*cur_f)/tot_feats)+'%)')
        #cur_f += 1
        if feature==0:
            continue
        #cur_seq=0
        #tot_seq=len(d[feature].keys())
        sorted_seq_l = sorted([(seq, len(set(np.hstack(alignment_ra.positions[d[feature][seq]])))) for seq in d[feature].keys()], key=lambda x:x[1])
        for idx, r_seqs in enumerate(sorted_seq_l):
            r_seq = r_seqs[0]
            #cur_seq+=1
            #sys.stdout.write('\rfeature: '+str(cur_f) + '/' + str(tot_feats) + ', seq: ' + str(cur_seq) + '/' + str(tot_seq))
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
                    error_count+=len(d[feature][r_seq])
                    break
    for (gene, r_seq) in for_removal:
        err_correction_res[ra_grouped[gene][r_seq],[err.ERROR_CORRECTION_jaitin]] = 1
    
    print ('\nJaitin error_count: ', error_count)
    tot_time=time.process_time()-start
    print('total Jaitin error_correction runtime: ',tot_time)
    return error_count, tot_time

# This is old and probably needs to be updated before use
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

# This is old and probably needs to be updated before use
def correct_errors_sten(ra_grouped, err_correction_res):
    """Correct errors using the method in Sten's paper.
       Remove any molecule supported by only a single read"""
    
    # for python 3
    start = time.process_time() 

    d = ra_grouped
    

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

    