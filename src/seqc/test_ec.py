import correct_errors as err
import seqc
import pickle
import time
from three_bit import ThreeBit as bin_rep       #Using Ambrose coder for unity
import numpy as np
import sam_reader


def from_thinair(barcode_files = ['/mnt/gel_barcode1_list_10.txt','/mnt/gel_barcode2_list_10.txt'], num_reads = 10000, err_rate = 0.01, num_features = 1, reverse_complement=True):
    '''
    USED FOR DEVELOPMENT ONLY
    Return a simulated grouped ReadArray for testing
    '''
    bases = ['A','C','G','T']
    alignment_d = {}
    
    def apply_error(str, err_rate):
        res = list(str)
        error = False
        for i in range(len(str)):
            if random.uniform(0,1)<err_rate:
                res[i] = random.choice(bases)   #There's a 1/4 chance for the same base to be elected so the actual err_rate is a bit smaller
                error = True
        return ''.join(res), error
    
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
    codes = get_codes(barcode_files=barcode_files, reverse_complement=reverse_complement)
    validation_arr = np.zeros(num_reads)
    
    for i in range(num_reads):
        error = False
        c1, err = apply_error(random.choice(codes[0]), err_rate)
        c1 = bin_rep.str2bin(c1)
        error |= err
        c2, err = apply_error(random.choice(codes[1]), err_rate)
        c2 = bin_rep.str2bin(c2)
        error|= err
        rmt, err = apply_error(generate_correct_rmt(), err_rate)
        rmt = bin_rep.str2bin(rmt)
        error |= err
        feature = random.choice(feature_list)
        seq = bin_rep.ints2int([c1, c2, rmt])
        if feature in alignment_d.keys():
            try:    #TODO: try and use a list instead of an ndarray to save memory
                alignment_d[feature][seq] = np.append(alignment_d[feature][seq],i)    # Since we don't have a real RA to link back to, we use the index i itself
            except KeyError:
                alignment_d[feature][seq] = np.array([i])
        else:
            alignment_d[feature] = {}
            alignment_d[feature][seq] = np.array([i])
        if error:
            validation_arr[i] = 1
    return alignment_d, validation_arr


def compare_methods(err_res_mat):
    """Return a matrix containig the Jaccard score between each two methods"""
    jac_mat = np.zeros((err.NUM_OF_ERROR_CORRECTION_METHODS,err.NUM_OF_ERROR_CORRECTION_METHODS))
    for i in range(jac_mat.shape[0]):
        for j in range(jac_mat.shape[1]):
            jac_mat[i,j] = sum(err_res_mat[:,i]*err_res_mat[:,j]) / (sum(err_res_mat)[i]+sum(err_res_mat)[j]-sum(err_res_mat[:,i]*err_res_mat[:,j]))
    return jac_mat

def compare_methods_venn(err_res_mat):
    ''' Return the size of each area in the 4 venn diagram'''
    r = {'none': 0, 
          'AJC': 0,
          'sten': 0,
          'allon': 0,
          'jaitin': 0,
          'AJC-sten': 0,
          'AJC-allon': 0,
          'AJC-jaitin':0,
          'sten-allon': 0,
          'sten-jaitin': 0,
          'allon-jaitin': 0,
          'AJC-sten-allon': 0,
          'AJC-allon-jaitin': 0,
          'sten-allon-jaitin': 0,
          'all': 0}
          
    #Doing the calculations beforehand saves time by storing the result and avoiding duplicate calculations
    sums = sum(err_res_mat)         # area of all 4 methods
    aj_s = sum(err_res_mat[:,err.ERROR_CORRECTION_AJC]*err_res_mat[:,err.ERROR_CORRECTION_STEN])    #|AJC /\ sten|
    aj_al = sum(err_res_mat[:,err.ERROR_CORRECTION_AJC]*err_res_mat[:,err.ERROR_CORRECTION_ALLON])  #|AJC /\ allon|
    aj_j = sum(err_res_mat[:,err.ERROR_CORRECTION_AJC]*err_res_mat[:,err.ERROR_CORRECTION_jaitin]) #|AJC /\ jaitin|
    s_al = sum(err_res_mat[:,err.ERROR_CORRECTION_ALLON]*err_res_mat[:,err.ERROR_CORRECTION_STEN])  #|sten /\ allon|
    s_j = sum(err_res_mat[:,err.ERROR_CORRECTION_jaitin]*err_res_mat[:,err.ERROR_CORRECTION_STEN]) #|sten /\ jaitin|
    al_j = sum(err_res_mat[:,err.ERROR_CORRECTION_ALLON]*err_res_mat[:,err.ERROR_CORRECTION_jaitin])   #|allon /\ jaitin|
    aj_s_a = sum(err_res_mat[:,err.ERROR_CORRECTION_ALLON]*err_res_mat[:,err.ERROR_CORRECTION_STEN]*err_res_mat[:,err.ERROR_CORRECTION_AJC])    #|allon /\ sten /\ AJC|
    aj_s_j = sum(err_res_mat[:,err.ERROR_CORRECTION_AJC]*err_res_mat[:,err.ERROR_CORRECTION_STEN]*err_res_mat[:,err.ERROR_CORRECTION_jaitin])  #|AJC /\ sten /\ jaitin|
    s_a_j = sum(err_res_mat[:,err.ERROR_CORRECTION_ALLON]*err_res_mat[:,err.ERROR_CORRECTION_STEN]*err_res_mat[:,err.ERROR_CORRECTION_jaitin]) #|allon /\ sten /\ jaitin|
    aj_a_j = sum(err_res_mat[:,err.ERROR_CORRECTION_ALLON]*err_res_mat[:,err.ERROR_CORRECTION_jaitin]*err_res_mat[:,err.ERROR_CORRECTION_AJC]) #|allon /\ jaitin /\ AJC|

    r['all'] = sum(err_res_mat[:,err.ERROR_CORRECTION_AJC]*err_res_mat[:,err.ERROR_CORRECTION_STEN]*err_res_mat[:,err.ERROR_CORRECTION_jaitin]*err_res_mat[:,err.ERROR_CORRECTION_ALLON]) #|AJC /\ sten /\ allon /\ jaitin|

    r['AJC-sten-allon'] = aj_s_a - r['all']
    r['AJC-sten-jaitin'] = aj_s_j - r['all']
    r['AJC-allon-jaitin'] = aj_a_j - r['all']
    r['sten-allon-jaitin'] = s_a_j - r['all']

    r['AJC-allon'] = aj_al - r['AJC-sten-allon'] - r['AJC-allon-jaitin'] - r['all']
    r['AJC-sten'] = aj_s - r['AJC-sten-allon'] - r['AJC-sten-jaitin'] - r['all']
    r['AJC-jaitin'] = aj_j - r['AJC-sten-jaitin'] - r['AJC-allon-jaitin'] - r['all']
    r['sten-allon'] = s_al - r['AJC-sten-allon'] - r['sten-allon-jaitin'] - r['all']
    r['sten-jaitin'] = s_j - r['AJC-sten-jaitin'] - r['sten-allon-jaitin'] - r['all']
    r['allon-jaitin'] = al_j - r['AJC-allon-jaitin'] - r['sten-allon-jaitin'] - r['all']

    r['AJC'] = sums[err.ERROR_CORRECTION_AJC] - r['AJC-sten-allon'] - r['AJC-sten-jaitin'] - r['AJC-allon-jaitin'] - r['AJC-allon'] - r['AJC-sten'] - r['AJC-jaitin'] - r['all']
    r['sten'] = sums[err.ERROR_CORRECTION_STEN] - r['AJC-sten-allon'] - r['AJC-sten-jaitin'] -  r['sten-allon-jaitin'] - r['AJC-sten'] - r['sten-allon'] - r['sten-jaitin'] - r['all']
    r['allon'] = sums[err.ERROR_CORRECTION_ALLON] - r['AJC-sten-allon'] - r['AJC-allon-jaitin'] - r['sten-allon-jaitin'] - r['AJC-allon'] - r['sten-allon'] - r['allon-jaitin'] - r['all']
    r['jaitin'] = sums[err.ERROR_CORRECTION_jaitin] - r['AJC-sten-jaitin'] - r['AJC-allon-jaitin'] - r['sten-allon-jaitin'] - r['AJC-jaitin'] - r['sten-jaitin'] - r['allon-jaitin'] - r['all']

    r['none']  = err_res_mat.shape[0] - sum(r.values())
    
    return r

def test_err_correction(fname, read_array = None, barcodes = ['/mnt/gel_barcode1_list.txt','/mnt/gel_barcode2_list.txt'], annotations = '/mnt/annotations.gtf', fl = 1000, mat_file=''):
    '''run the error correction on all 4 methods and analyse the results'''
    
    #fout = open(fname+'_filter_out.txt','w+')
    ra = read_array
    if ra == None:
        print('converting features')
        cf = seqc.convert_features.ConvertFeatureCoordinates.from_gtf(annotations, fl)
        print ('reading from samfile')
        ra = seqc.arrays.ReadArray.from_samfile(fname,cf)
    #fout.write('len ra: '+str(len(ra))+'\n\n')
    print('correcting errors')  
    res_full, res_sum = err.correct_errors(ra, barcodes, compare_methods=True, P_VALUE=0.05)
    #fout.write(str(res_sum)+'\n\n')
    print('Applying filters')
    apply_filters(ra, res_full)
    #print('analyzing results - jaccard score')
    #r1 = compare_methods(res_full)
    #fout.write(str(r1)+'\n\n')
    #print('analyzing results - venn')
    #r2 = compare_methods_venn(res_full)
    #fout.write(str(r2)+'\n\n')
    if mat_file != '':
        print('Dumping result matrix to file')
        matf = open(mat_file, 'wb')
        pickle.dump(res_full, matf, protocol=4)
        matf.close()
    #fout.close()
    return res_full
    
def calc_read_dist(fname, barcodes = ['/mnt/gel_barcode1_list.txt','/mnt/gel_barcode2_list.txt'], annotations = '/mnt/annotations.gtf', fl = 1000):
    '''return a histogram of numbers of reads per gene/seq'''
    
    fout = open(fname+'_dist.txt','w+')
    print('converting features')
    cf = seqc.convert_features.ConvertFeatureCoordinates.from_gtf(annotations, fl)
    print ('reading from samfile')
    ra = seqc.arrays.ReadArray.from_samfile(fname,cf)
    
    print('grouping ReadArray')  
    ra_grouped = ra.group_for_error_correction(required_poly_t = 1)
    
    print('Counting read appearence frequency')
    hist = {}
    for feature in ra_grouped.keys():
        for seq in ra_grouped[feature].keys():
            q = ra_grouped[feature][seq].shape[0]
            if q in hist.keys():
                hist[q] += 1
            else:
                hist[q] = 1
    for k,v in hist.items():
        fout.write(str(k)+'\t'+str(v)+'\n')
    
    fout.close()

def validate_correct_errors(num_reads = 10000, err_rate = 0.01, num_features = 1, barcode_files = ['/mnt/gel_barcode1_list_10.txt','/mnt/gel_barcode2_list_10.txt'], reverse_complement = True, donor_cutoff=1, P_VALUE=0.1):
    ''' Use a simulation model to validate error_correction. The jaitin method is not validated and instead its column in the matrix is used for the
        true errors.'''
    
    fname = 'r'+str(num_reads)+'_er'+str(err_rate)+'_f'+str(num_features)
    fout = open(fname+'_simul.txt','w+')
    
    res_time_cnt = {}
    
    err_correction_res = np.zeros((num_reads, NUM_OF_ERROR_CORRECTION_METHODS))
    print('simulationg data')
    ra_grouped, validation = from_thinair(barcode_files, num_reads, err_rate, num_features)
    # The ra is not needed for these 3 methods and should also be removed as a parameter...
    print ('doing AJC')
    res_time_cnt[ERROR_CORRECTION_AJC] = correct_errors_AJC(ra_grouped, err_correction_res, barcode_files, reverse_complement, donor_cutoff, P_VALUE)
    print ('doing Sten')
    res_time_cnt[ERROR_CORRECTION_STEN] = correct_errors_sten(ra_grouped, err_correction_res)
    print ('doing Allon')
    res_time_cnt[ERROR_CORRECTION_ALLON] = correct_errors_allon(ra_grouped, err_correction_res, barcode_files, reverse_complement)
    err_correction_res[:,ERROR_CORRECTION_jaitin] = validation
    
    fout.write(str(res_time_cnt)+'\n\n')
    print('analyzing results - jaccard score')
    r1 = compare_methods(err_correction_res)
    fout.write(str(r1)+'\n\n')
    print('analyzing results - venn')
    r2 = compare_methods_venn(err_correction_res)
    fout.write(str(r2)+'\n\n')
        
    fout.close()
    
def apply_filters(alignments_ra, err_correction_res):
    """Correct errors using the method in Sten's paper.
       Remove any molecule supported by only a single read"""
    
    # for python 3
    start = time.process_time() 

    error_count = 0
    N = bin_rep._str2bindict['N']

    for idx, read in enumerate(alignments_ra.data):
        if bin_rep.contains(read['rmt'], N):
                continue
        if read['dust_score'] >= 10:
            err_correction_res[idx][err.ERROR_CORRECTION_FILTERS] = 1
            error_count+=1
    
    print ('filter count: ', error_count)
    tot_time=time.process_time()-start
    print('total filters runtime: ',tot_time)
    return error_count, tot_time
    
def gc_content(res_mat, samfile, col1, col2):
    ''' return the average gc content of reads from a sam_file that correspond to col in res_mat'''
    
    tot_len = {'all':0, 'col1': 0, 'col2': 0, 'both': 0}
    tot_gc = {'all':0, 'col1': 0, 'col2': 0, 'both': 0}
    sr = sam_reader.Reader(samfile)
    for i, r in enumerate(sr.iter_multialignments()):
        l = len(r[0].seq)
        gc = gc_count(r[0].seq)
        if i < res_mat.shape[0]:    # For a weird reason the ra is one record short on my test files
            tot_len['all'] += l
            tot_gc['all'] += gc
            if res_mat[i][col1] == 1 and res_mat[i][col2] == 1:
                tot_len['both'] += l
                tot_gc['both'] += gc
            elif res_mat[i][col1] == 1 and res_mat[i][col2] == 0:
                tot_len['col1'] += l
                tot_gc['col1'] += gc
            elif res_mat[i][col1] == 0 and res_mat[i][col2] == 1:
                tot_len['col2'] += l
                tot_gc['col2'] += gc
                
    res = {}
    for k in tot_len.keys():
        res[k] = tot_gc[k]/tot_len[k]
    return res
    
def triple_entropy(res_mat, samfile, col1, col2, fname=''):
    ''' return the average gc content of reads from a sam_file that correspond to col in res_mat'''
    if fname != '':
        fout = open(fname, 'w+')
    dist_all = {}
    dist_col1 = {}
    dist_col2 = {}
    dist_both = {}
    sr = sam_reader.Reader(samfile)
    l = len(sr)
    for i, r in enumerate(sr.iter_multialignments()):
        if 'N' in r[0].seq:
            continue
        if i%10000 == 0:
            print('\r'+str(i)+'/'+str(l)+' done')
        if i < res_mat.shape[0]:    # For a weird reason the ra is one record short on my test files
            dist_all = triple_entropy_count(r[0].seq, dist_all)
            if res_mat[i][col1] == 1 and res_mat[i][col2] == 1:
                dist_both = triple_entropy_count(r[0].seq, dist_both)
            elif res_mat[i][col1] == 1 and res_mat[i][col2] == 0:
                dist_col1 = triple_entropy_count(r[0].seq, dist_col1)
            elif res_mat[i][col1] == 0 and res_mat[i][col2] == 1:
                dist_col2 = triple_entropy_count(r[0].seq, dist_col2)                
    
    if fname != '':
        fout.write('triplet\tall\tcol1\tcol2\tboth\n')
        for k in sorted(dist_all.keys()):
            fout.write(k + '\t')
            try:
                fout.write(str(dist_all[k])+'\t')
            except KeyError:
                fout.write('0\t')
            try:
                fout.write(str(dist_col1[k])+'\t')
            except KeyError:
                fout.write('0\t')
            try:
                fout.write(str(dist_col2[k])+'\t')
            except KeyError:
                fout.write('0\t')
            try:
                fout.write(str(dist_both[k])+'\n')
            except KeyError:
                fout.write('0\n')
        fout.close()
    
def gc_count(seq):
    count=0
    for base in seq:
        if base=='G' or base=='C':
            count+=1
    return count

def triple_entropy_count(seq, dist_d):
    for i in range(3,len(seq)+1):
        try:
            dist_d[seq[i-3:i]] += 1
        except KeyError:
            dist_d[seq[i-3:i]] = 1
    return dist_d
    
def cell_dist(res_mat, ra, col1, col2, fname=''):
    '''return the average cell distribution for reads that correspond to col1, col2, col1/\col2 and all in res_mat'''
    if fname != '':
        fout = open(fname, 'w+')
    dist_all = {}
    dist_col1 = {}
    dist_col2 = {}
    dist_both = {}
    
    for i, r in enumerate(ra):
        if bin_rep.contains(r[0]['rmt'], 'N'):
            continue
        cell = r[0]['cell']
        try:
            dist_all[cell] += 1
        except KeyError:
            dist_all[cell] = 1
        if res_mat[i][col1] == 1 and res_mat[i][col2] == 1:
            try:
                dist_both[cell] += 1
            except KeyError:
                dist_both[cell] = 1
        elif res_mat[i][col1] == 1 and res_mat[i][col2] == 0:
            try:
                dist_col1[cell] += 1
            except KeyError:
                dist_col1[cell] = 1
        elif res_mat[i][col1] == 0 and res_mat[i][col2] == 1:
            try:
                dist_col2[cell] += 1
            except KeyError:
                dist_col2[cell] = 1

    if fname != '':
        fout.write('cell\tall\tcol1\tcol2\tboth\n')
        for k in dist_all.keys():
            fout.write(str(k) + '\t')
            try:
                fout.write(str(dist_all[k])+'\t')
            except KeyError:
                fout.write('0\t')
            try:
                fout.write(str(dist_col1[k])+'\t')
            except KeyError:
                fout.write('0\t')
            try:
                fout.write(str(dist_col2[k])+'\t')
            except KeyError:
                fout.write('0\t')
            try:
                fout.write(str(dist_both[k])+'\n')
            except KeyError:
                fout.write('0\n')
        fout.close()
        
        
def num_pos_dist(ra, fname = ''):
    dist_dic = {}
    res_dic = {}
    for i, r in enumerate(ra):
        if ra.features[i].shape[0] != 1:
            continue
        pos = ra.positions[i]
        if len(pos) < 1:
            continue
        gene = ra.features[i][0]
        try:
            l = dist_dic[(r[0]['cell'], r[0]['rmt'], gene)]
            dist_dic[(r[0]['cell'], r[0]['rmt'], gene)] = np.concatenate((l,pos))
        except KeyError:
            dist_dic[(r[0]['cell'], r[0]['rmt'], gene)] = pos
    
    for k in dist_dic.keys():
        try:
            res_dic[len(set(dist_dic[k]))]+=1
        except KeyError:
            res_dic[len(set(dist_dic[k]))] = 1
    if fname!='':
        fout = open(fname,'w')
        fout.write('num of positions\tnum reads\n')
        for num_pos in res_dic.keys():
            fout.write(str(num_pos)+'\t'+str(res_dic[num_pos])+'\n')    
        fout.close()
    return res_dic

def pos_instances(ra, fname = ''):
    pos_dic = {}
    for i, r in enumerate(ra):
        if ra.features[i].shape[0] != 1:
            continue
        for pos in ra.positions[i]:
            try:
                pos_dic[pos]+=1
            except KeyError:
                pos_dic[pos] = 1
    if fname!='':
        fout = open(fname, 'w')
        fout.write('position\tnum reads\n')
        for pos in pos_dic.keys():
            fout.write(str(pos)+'\t'+str(pos_dic[pos])+'\n')
        fout.close()
        

def pos_gene_instances(ra, fname = ''):
    pos_dic = {}
    for i, r in enumerate(ra):
        if ra.features[i].shape[0] != 1:
            continue
        gene = ra.features[i][0]
        for pos in ra.positions[i]:
            try:
                pos_dic[gene][pos]+=1
            except KeyError:
                try:
                    pos_dic[gene][pos] = 1
                except KeyError:
                    pos_dic[gene] = {}
                    pos_dic[gene][pos] = 1
    if fname != '':
        fout = open(fname, 'w')
        fout.write('gene\tposition\tnum reads\n')
        for gene in pos_dic.keys():
            for pos in pos_dic[gene].keys():
                fout.write(str(gene)+'\t'+str(pos)+'\t'+str(pos_dic[gene][pos])+'\n')
        fout.close()

    
def group_hamm_dist(ra, reads):
    '''Gets an ra and a list on indices representing reads and returns:
        1. Average hamming distance between all pairs
        2. fraction of pairs with hamming distance of 1 or less'''

    cnt = 0
    sum = 0
    low_cnt = 0
    for r1 in reads:
        for r2 in reads:
            d_codes = err.hamming_dist_bin(int(ra[r1][0]['cell']), int(ra[r2][0]['cell']))
            d_rmt = err.hamming_dist_bin(int(ra[r1][0]['rmt']), int(ra[r2][0]['rmt']))
            if d_codes!=err.high_value:
                cnt+=1
                sum+=d_codes+d_rmt
                if d_codes+d_rmt<= 1:
                    low_cnt+=1
    return sum/cnt, low_cnt/cnt
    
def test_ham_dist(ra, pos, gene=None):
    
    reads = []
    for i, r in enumerate(ra):
        if pos in ra.positions[i]:
            if gene==None:
                reads.append(i)
            elif gene in ra.features[i]:
                reads.append(i)
    avg, frac = group_hamm_dist(ra, reads)
    return avg, frac

def n_random_reads(ra, fname='', n=100):
    '''Returns indexes of n random reads from the readArray'''
    iter = 1000
    res = []
    for i in range(iter):
        res.append(group_hamm_dist(ra, np.random.randint(len(ra), size=n)))
    if fname!='':
        f = open(fname,'w')
        f.write('average distance\tfraction of close reads\n')
        for avg, frac in res:
            f.write(str(avg) + '\t' + str(frac) + '\n')
        f.close()