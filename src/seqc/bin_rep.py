#coding of bases into binary. Note the first bit is on for a valid base and off for "NA". the representation 0 is not used and means "empty".
str2bindict = {'A':0b100, 'C':0b101, 'G':0b110, 'T':0b111, 'N':0b001,
'a':0b100, 'c':0b101, 'g':0b110, 't':0b111, 'n':0b001}
bin2strdict = {0b100:'A', 0b101:'C', 0b110:'G', 0b111:'T', 0b001:'N'}
bin_nums = [0b100, 0b101, 0b110, 0b111]

def str2bin(s):
	'''turn a binary representation into a string'''
    res = 0
    for c in s:
        res<<=3
        res+=str2bindict[c]
    return res

def bin2str(bin_):
	'''Return the binary representation of a sequence'''
    r=''
    while bin_>0:
        r=bin2strdict[bin_ & 0b111]+r
        bin_>>=3
    return r
    
def ints2int(l):
	''' Turn a list of ints into a long int comprising all of them in order of l[0] on the msb side'''
    res = 0
    for num in l:
        tmp = num
        while tmp>0:    # Get the length of the next number to concatenate (with enough room for leading zeros)
            res<<=3
            tmp>>=3
        res+=num
    return res

def rmt_from_int(seq):
	'''Extract the rmt (6 chrs) from an encoded sequecne'''
    return seq&0o777777

def c2_from_int(seq):
	'''Extract barcode2 from a sequence'''
    return (seq&0o77777777000000)>>(6*3)

def c1_from_int(seq):
	'''Extract barcode1 from a sequence'''
    return seq>>((8+6)*3)


def contains(seq, chr):
	'''return true if the chr (bin representation) is contained in seq (binary representation)	'''
	while seq>0:
		if chr==(seq&0b111):
			return True
		seq>>=3
	return False

def seq_len(seq):
	'''Return the length of a sequence based on its binary representation'''
	l=0
	while seq>0:
		l+=1
		seq>>=3
	return l
