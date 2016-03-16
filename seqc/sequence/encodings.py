class DNA3Bit:
    """
    Compact encoding scheme for sequence data.
    """

    _str2bindict = {65: 0b100, 67: 0b110, 71: 0b101, 84: 0b011, 78: 0b111,
                    97: 0b100, 99: 0b110, 103: 0b101, 116: 0b011, 110: 0b111}
    _bin2strdict = {0b100: b'A', 0b110: b'C', 0b101: b'G', 0b011: b'T', 0b111: b'N'}
    bin_nums = [0b100, 0b110, 0b101, 0b011]

    @classmethod
    def encode(cls, b: bytes) -> int:
        """Convert string nucleotide sequence into binary, note: string is reversed so
        that the first nucleotide is in the LSB position
        :param b: """
        res = 0
        for c in b:
            res <<= 3
            res += cls._str2bindict[c]
        return res

    @classmethod
    def decode(cls, i: int) -> bytes:
        """Convert binary nucleotide sequence into string
        :param i:
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
        """calculate percentage of i that is G or C
        :param i:
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
        """Return the length of a sequence based on its binary representation
        :param i:
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
        :param char:
        :param s:
        """
        while s > 0:
            if char == (s & 0b111):
                return True
            s >>= 3
        return False

    @staticmethod
    def bitlength(i: int) -> int:
        """return the bitlength of the sequence
        :param i:
        """
        bitlen = i.bit_length()
        # correct for leading T-nucleotide (011) whose leading 0 gets trimmed
        if bitlen % 3:
            bitlen += 1
        return bitlen

    @staticmethod
    def ints2int(ints):
        """convert [i1, i2, i3] -> 0bi1i2i3
        :param ints:
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
