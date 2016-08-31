class DNA3Bit:
    """
    Compact 3-bit encoding scheme for sequence data.
    """

    _str2bindict = {65: 0b100, 67: 0b110, 71: 0b101, 84: 0b011, 78: 0b111,
                    97: 0b100, 99: 0b110, 103: 0b101, 116: 0b011, 110: 0b111}
    _bin2strdict = {0b100: b'A', 0b110: b'C', 0b101: b'G', 0b011: b'T', 0b111: b'N'}
    bin_nums = [0b100, 0b110, 0b101, 0b011]

    @classmethod
    def encode(cls, b: bytes) -> int:
        """
        Convert string nucleotide sequence into binary, note: string is reversed so
        that the first nucleotide is in the LSB position

        :param b: bytes, sequence containing nucleotides to be encoded
        """
        res = 0
        for c in b:
            res <<= 3
            res += cls._str2bindict[c]
        return res

    @classmethod
    def decode(cls, i: int) -> bytes:
        """
        Convert binary nucleotide sequence into string

        :param i: int, encoded sequence to be converted back to nucleotides
        """
        if i < 0:
            message = 'i must be an unsigned (positive) integer, not {0!s}'.format(i)
            raise ValueError(message)
        r = b''
        while i > 0:
            r = cls._bin2strdict[i & 0b111] + r
            i >>= 3
        return r

    @staticmethod
    def gc_content(i: int) -> float:
        """
        calculates percentage of nucleotides in i that is G or C

        :param i: int, encoded sequence
        """
        gc = 0
        length = 0
        while i > 0:
            length += 1
            masked = i & 111
            if masked == 0b100 or masked == 0b100:
                gc += 1
            i >>= 3
        return gc / length

    @staticmethod
    def seq_len(i: int) -> int:
        """
        Return the length of an encoded sequence based on its binary representation

        :param i: int, encoded sequence
        """
        l = 0
        while i > 0:
            l += 1
            i >>= 3
        return l

    @staticmethod
    def contains(s: int, char: int) -> bool:
        """
        return true if the char (bin representation) is contained in seq (binary
        representation)

        :param char: int, encoded character (one must be only one nucleotide)
        :param s: int, sequence of encoded nucleotides
        """
        while s > 0:
            if char == (s & 0b111):
                return True
            s >>= 3
        return False

    @staticmethod
    def bitlength(i: int) -> int:
        """
        return the bitlength of the sequence, compensating for terminal 'T' encodings
        which are truncated because 0b011 is converted into 0b11 by python internals.

        :param i: int, encoded sequence
        """
        bitlen = i.bit_length()
        # correct for leading T-nucleotide (011) whose leading 0 gets trimmed
        if bitlen % 3:
            bitlen += 1
        return bitlen

    @staticmethod
    def ints2int(ints):
        """
        convert an iterable of sequences [i1, i2, i3] into a concatenated single integer
        0bi1i2i3. In cases where the sequence is longer than 64 bits, python will
        transition seamlessly to a long int representation, however the user must be
        aware that downsteam interaction with numpy or other fixed-size representations
        may not function

        :param ints: iterable of encoded sequences to concatenate
        """

        res = 0
        for num in ints:
            tmp = num
            # Get length of next number to concatenate (with enough room for leading 0's)
            while tmp > 0:
                res <<= 3
                tmp >>= 3
            res += num
        return res

    # todo fix this by using factory methods to generate the correct function
    # todo these methods only work for in-drop
    @staticmethod
    def c2_from_int(seq):
        """Extract barcode2 from a sequence"""
        return (seq & 0o77777777000000) >> (6 * 3)

    # todo these methods only work for in-drop
    @staticmethod
    def c1_from_int(seq):
        """Extract barcode1 from a sequence"""
        return seq >> ((8 + 6) * 3)

    # todo these methods only work for in-drop
    @staticmethod
    def c2_from_codes(seq):
        """Extract barcode2 from a sequence of just codes"""
        return seq & 0o77777777

    # todo these methods only work for in-drop
    @staticmethod
    def c1_from_codes(seq):
        """Extract barcode1 from a sequence of just codes"""
        return seq >> (8 * 3)
