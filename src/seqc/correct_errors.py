import random
#import bin_rep
from scipy.special import gammaincinv
from itertools import permutations
from sys import maxsize
import time
from seqc.three_bit import ThreeBit as bin_rep		#Using Ambrose coder for unity

high_value = maxsize  # Used for sorting, needs to be longer than any sequence  


def from_thinair(barcodes_files, num_reads = 1000000, err_rate = 0.02, num_features = 1, reverse_complement=True):
    '''
    USED FOR DEVELOPMENT ONLY
    Simulate aligment reads. This is used for testing error_correction, so positions are for now ignored.
    '''
    bases = ['A','C','G','T']
    alignment_d = {}
    
    def apply_error(str, err_rate):
        res = list(str)
        for i in range(len(str)):
            if random.uniform(0,1)<err_rate:
                res[i] = random.choice(bases)   #There's a 1/4 chance for the same base to be elected so the err_rate is a bit smaller
        return ''.join(res)
    
    def generate_correct_rmt(rmt_length = 6):
        #return 'AAAAAA'
        res = ['A','A','A','A'] # To make sure erroneous rmt's are not as common as correct ones, we limit the range of correct rmt's to be only those starting with 'AA...'
        for i in range(rmt_length-len(res)):
            res.append(random.choice(bases))
        return ''.join(res)
    
    #TODO: this can be moved outside and merged with the snippet inside estimate error rate, also change to binary
    def get_codes(barcode_files, reverse_complement):
        correct_barcodes = []
        if reverse_complement:
            for barcode_file in barcode_files:
                with open(barcode_file) as f:
                    correct_barcodes.append(list(set(rev_comp(line.strip())
                                                for line in f.readlines())))
        else:
            for barcode_file in barcode_files:
                with open(barcode_file) as f:
                    correct_barcodes.append(list(set(line.strip()
                                                for line in f.readlines())))
        return correct_barcodes
        
    feature_list = range(num_features)
    codes = get_codes(barcode_files=barcodes_files, reverse_complement=reverse_complement)
    
    for i in range(num_reads):
        c1 = bin_rep.str2bin(apply_error(random.choice(codes[0]), err_rate))
        c2 = bin_rep.str2bin(apply_error(random.choice(codes[1]), err_rate))
        rmt = bin_rep.str2bin(apply_error(generate_correct_rmt(), err_rate))
        flist = tuple([random.choice(feature_list)])
        seq = bin_rep.ints2int([c1, c2, rmt])
        if flist in alignment_d.keys():
            try:
                alignment_d[flist][seq]=alignment_d[flist][seq]+1
            except KeyError:
                alignment_d[flist][seq] = 1
        else:
            alignment_d[flist] = {}
            alignment_d[flist][seq] = 1
    return alignment_d
    
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
            

def estimate_error_rate(barcode_files, alignment_dic, reverse_complement=True):
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

    for gene in alignment_dic.keys():
        for seq in alignment_dic[gene].keys():
            c1 = bin_rep.c1_from_int(seq)
            c2 = bin_rep.c2_from_int(seq)
            code_errors = []
            
            if c1 not in correct_barcodes[0]:
                code_errors += list_errors(c1, correct_barcodes[0])
            if c2 not in correct_barcodes[1]:
                code_errors += list_errors(c2, correct_barcodes[1])
                    
            if code_errors==[]:
                correct_instances += (bin_rep.seq_len(c1) + bin_rep.seq_len(c2)) / 4 
            else:
                for err_type in code_errors:
                    try:
                        error_table[err_type] += 1
                    except TypeError:
                        pass
                    except KeyError:        # TODO: The 'N' in the codes create a key error. A good idea might be to just ignore them or even filter them completely.
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
            err_list.append(bin_rep.ints2int([code&0b111, donor&0b111]))
        code>>=3
    
    return err_list
        
def correct_errors(alignment_dic, barcode_files, reverse_complement=True, donor_cutoff=1, P_VALUE=0.1):
    """calculate and correct errors in barcode sets"""
    start = time.process_time()
    d = prepare_dic(alignment_dic)
#   print('version 6.0 - using readArray')

    error_count = 0
    err_rate = estimate_error_rate(barcode_files, d, reverse_complement)
    
    N = bin_rep._str2bindict['N']
    for_removal = []
    tot=0
    for feature in d.keys():
        for r_seq in d[feature].keys():
            
            if bin_rep.contains(r_seq, N):
                continue
            
            gene = feature
            r_c1 = bin_rep.c1_from_int(r_seq)
            r_c2 = bin_rep.c2_from_int(r_seq)
            r_rmt = bin_rep.rmt_from_int(r_seq)
            r_num_occurences = d[gene][r_seq]


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
a           and so we need to regect r_seq as an error. lower P_VALUE will increase the number of errors (we require higher level of certainty)         
            """

            threshold = gammaincinv(r_num_occurences, P_VALUE)
            
            expected_errors = 0
            for d_seq in generate_close_seq(r_seq):
                try:
                    d_num_occurences = d[gene][d_seq]
                except KeyError:
                    continue
                if d_num_occurences<=donor_cutoff:
                    continue
                tot+=1
                d_c1 = bin_rep.c1_from_int(d_seq)
                d_c2 = bin_rep.c2_from_int(d_seq)
                d_rmt = bin_rep.rmt_from_int(d_seq)

                p_dtr = prob_d_to_r_bin(d_seq, r_seq, err_rate)
                
                expected_errors += d_num_occurences * p_dtr
                if expected_errors > threshold:
                    for_removal.append((gene, r_seq))
                    error_count+=1
                    break
    for (gene, r_seq) in for_removal:
        del d[gene][r_seq]
    
    print (error_count, ' errors were found')
    tot_time=time.process_time()-start
    print('total error_correction runtime: ',tot_time)
    return
    
#An interface function to convert the read array into a sequence dic used by error_correction.
def prepare_dic(read_array):
    tmp = read_array.group_for_error_correction(required_poly_t = 0)
    d={}
    for k_lv1 in tmp.keys():
        d[k_lv1] = tmp[k_lv1]
        for k_lv2 in tmp[k_lv1].keys():
            d[k_lv1][k_lv2] = tmp[k_lv1][k_lv2].shape[0]
    return d
    
