__author__ = 'Ambrose J. Carr'


class ThreeBit:
    """
    method to extract cell barcodes and rmts given variable lengths of barcodes
    """

    _str2bindict = {'A': 0b100, 'C': 0b110, 'G': 0b101, 'T': 0b011, 'N': 0b111,
                    'a': 0b100, 'c': 0b110, 'g': 0b101, 't': 0b011, 'n': 0b111}
    _bin2strdict = {0b100: 'A', 0b110: 'C', 0b101: 'G', 0b011: 'T', 0b111: 'N'}
    bin_nums = [0b100, 0b110, 0b101, 0b011]

    def __init__(self, cell_len, rmt_len):
        """
        Class to allow interconversion between python string and three-bit bitstring
        encoded python integer objects.

        The standard form for this structure is a bitfield with subfields of lengths
         defined by [cell (can be variable) | rmt ]

        args:
        -----
        cell: length of the cell barcode from the start of the forward read
        rmt: length of the rmt from the end of the cell barcode
        """
        self._rmt_len = rmt_len
        # self._rmt_mask = int('0b' + '111' * rmt_len, base=2)
        self._rmt_mask = (1 << (rmt_len * 3)) - 1
        self._cell_len = cell_len
        self._cell_mask = (1 << (cell_len * 3)) - 1
        # self._cell_mask = int('0b' + '111' * cell_len, base=2)

    @classmethod
    def str2bin(cls, s):
        """Convert string nucleotide sequence into binary, note: string is reversed so
        that the first nucleotide is in the LSB position"""
        res = 0
        for c in s:
            res <<= 3
            res += cls._str2bindict[c]
        return res

    @classmethod
    def bin2str(cls, i):
        """Convert binary nucleotide sequence into string"""
        if i < 0:
            raise ValueError('i must be an unsigned (positive) integer, not %d' % i)
        r = ''
        while i > 0:
            r = cls._bin2strdict[i & 0b111] + r
            i >>= 3
        return r

    @staticmethod
    def ints2int(ints):
        """convert [i1, i2, i3] -> 0bi1i2i3"""
        if not all(i > 0 for i in ints):
            raise ValueError('all inputs must be unsigned (positive) integers')
        res = 0
        for num in ints:
            tmp = num
            # Get length of next number to concatenate (with enough room for leading 0's)
            while tmp > 0:
                res <<= 3
                tmp >>= 3
            res += num
        return res
    
    @staticmethod
    def rmt_from_int(seq):
<<<<<<< HEAD
        '''Extract the rmt (6 chrs) from an encoded sequecne'''
        return seq&0o777777

    @staticmethod
    def c2_from_int(seq):
        '''Extract barcode2 from a sequence'''
        return (seq&0o77777777000000)>>(6*3)

    @staticmethod
    def c1_from_int(seq):
        '''Extract barcode1 from a sequence'''
        return seq>>((8+6)*3)
=======
        """Extract the rmt (6 chrs) from an encoded sequecne"""
        return seq & 0o777777

    @staticmethod
    def c2_from_int(seq):
        """Extract barcode2 from a sequence"""
        return (seq & 0o77777777000000) >> (6 * 3)

    @staticmethod
    def c1_from_int(seq):
        """Extract barcode1 from a sequence"""
        return seq >> ((8 + 6) * 3)
>>>>>>> origin/master

    @staticmethod
    def gc_content(s):
        """calculate percentage of s that is G or C"""
        gc = 0
        length = 0
        while s > 0:
            length += 1
            masked = s & 111
            if masked == 0b100 or masked == 0b100:
                gc += 1
            s >>= 3
        return gc / length

    @staticmethod
    def _num_poly_t(s):
        n = 0
        while s:
            if (s & 0b111) == 0b011:
                n += 1
            s >>= 3
        return n

    @staticmethod
    def _num_poly_a(s):
        n = 0
        while s:
            if (s & 0b111) == 0b100:
                n += 1
            s >>= 3
        return n

    @classmethod
    def default_processors(cls, key):
        if key == 'mars-seq':
            raise NotImplementedError
        proc = {
            'in-drop': (8, 6),
            'drop-seq': (12, 8),
            'cel-seq': (6, 6),
            'avo-seq': (8, 4),
        }
        try:
            p = proc[key]
            if key == 'in-drop':
                return ThreeBitInDrop(*p)
            else:
                return cls(*p)
        except KeyError:
            raise ValueError('unknown processor')

    @staticmethod
    def split_cell(s):
        """return all of the cell barcodes present in s"""
        return tuple(s)

    @staticmethod
    def seq_len(s):
        """Return the length of a sequence based on its binary representation"""
        l = 0
        while s > 0:
            l += 1
            s >>= 3
        return l

    @staticmethod
    def contains(s, char):
        """
        return true if the char (bin representation) is contained in seq (binary
        representation)
        """
        while s > 0:
            if char == (s & 0b111):
                return True
            s >>= 3
        return False

    @staticmethod
    def gc_content(s):
        t = 0
        gc = 0
        while s > 0:
            if s & 0b110 == 0b110 or s & 0b101 == 0b101:
                gc += 1
            t += 1
            s >>= 3
        return gc / t

    def process_forward_sequence(self, record):
        s = self.str2bin(record)

        # get length of sequence
        bitlen = s.bit_length()

        # correct for leading T-nucleotide (011) whose leading 0 gets trimmed
        if bitlen % 3:
            bitlen += 1

        # get expected primer length
        primer_bitsize = 3 * (self._cell_len + self._rmt_len)
        if 3 * self._cell_len <= bitlen < primer_bitsize:  # truncated read
            cell = s >> (bitlen - (3 * self._cell_len))
            return cell, 0, 0
        elif bitlen < 3 * self._cell_len:
            return 0, 0, 0
        else:
            # count T in poly-T tail
            poly_t_mask = (1 << (bitlen - primer_bitsize)) - 1
            n_poly_t = self._num_poly_t(s & poly_t_mask)

            # get RMT
            s >>= (bitlen - primer_bitsize)
            rmt = s & self._rmt_mask

            # get cell
            cell = s >> (3 * self._rmt_len)

            return cell, rmt, n_poly_t

    @staticmethod
    def bitlength(s):
        bitlen = s.bit_length()
        # correct for leading T-nucleotide (011) whose leading 0 gets trimmed
        if bitlen % 3:
            bitlen += 1
        return bitlen


class ThreeBitInDrop(ThreeBit):
    """
    method to extract cell barcodes and rmts given variable lengths of barcodes
    """

    _w1 = 'GAGTGATTGCTTGTGACGCCTT'

    def process_forward_sequence(self, s):
        """the 1st cell barcode is 8-11bp, and the spacer is found directly afterwards.
        find the cell barcode in

        for our data, position 24 should have a different value for each of the c1
         lengths. We check this position first, then check p23 and p24 to ensure that
         the spacer is correctly localized.
        """

        s = self.str2bin(s)

        # get bit length of sequence
        bitlen = self.bitlength(s)

        # use to check for spacer positioning
        tetramer = (s >> (bitlen - (3 * 28))) & 0o7777

        # check that the bases match expected spacer locations, then extract c1 position
        if tetramer == 3446:  # CGCC
            cb1_length = 8 * 3
        elif tetramer == 2478:  # ACGC
            cb1_length = 9 * 3
        elif tetramer == 2869:  # GACG
            cb1_length = 10 * 3
        elif tetramer == 1894:  # TGAC
            cb1_length = 11 * 3
        else:
            return 0, 0, 0

        # set some variables
        cb2_length = 8 * 3
        spacer_length = 22 * 3
        rmt_length = 6 * 3

        # set expected primer length
        primer_bitsize = cb1_length + spacer_length + cb2_length + rmt_length

        # some reads are truncated, find the cell barcode then return
        truncated_length = cb1_length + spacer_length + cb2_length
        if truncated_length <= bitlen < primer_bitsize:

            # get CB2
            s >>= (bitlen - truncated_length)
            cb2 = s & self._cell_mask

            # get CB1
            cb1 = s >> cb2_length + spacer_length

            return self.ints2int([cb1, cb2]), 0, 0
        elif truncated_length > bitlen:
            return 0, 0, 0  # can't resolve barcode
        else:
            # count T in poly-T tail
            poly_t_mask = (1 << (bitlen - primer_bitsize)) - 1
            num_poly_t = self._num_poly_t(s & poly_t_mask)

            # get RMT
            s >>= (bitlen - primer_bitsize)
            rmt = s & self._rmt_mask

            # get CB2
            s >>= rmt_length
            cb2 = s & self._cell_mask

            # get CB1
            cb1 = s >> cb2_length + spacer_length

            return self.ints2int([cb1, cb2]), rmt, num_poly_t

    @staticmethod
    def split_cell(s):
        """return all of the cell barcodes present in s"""
        c2 = s & 0o77777777
        c1 = s >> (3 * 8)
        return c1, c2
