import random
import nose2
import unittest
import seqc



seqc_dir = '/'.join(seqc.__file__.split('/')[:-3]) + '/'


# todo this function is defunct
class TestMergeFunctionsRegex(unittest.TestCase):

    def test_merge_functions_specific_spacer(self):
        name = b'@something\n'
        sequence = b'GTATGTTGACGGCCATAAGGCTGCTTCTGACGTTCGTGATGAGTTTGTAT\n'
        name2 = b'+\n'
        quality = b'GTATGTTGACGGCCATAAGGCTGCTTCTGACGGTCGTGATGAGTTTGTAT\n'
        record = seqc.sequence.fastq.FastqRecord([name, sequence, name2, quality])
        merged = seqc.sequence.merge_functions.in_drop(record, record)
        print(merged)
        print(record.sequence[28:32])

    def test_merge_functions_regex_spacer(self):
        name = b'@something\n'
        sequence = b'GTATGTTGACGGCCATAAGGCTGCTTCTGACNTTCGTGATGAGTTTGTAT\n'
        name2 = b'+\n'
        quality = b'GTATGTTGACGGCCATAAGGCTGCTTCTGACNTTCGTGATGAGTTTGTAT\n'
        record = seqc.sequence.fastq.FastqRecord([name, sequence, name2, quality])
        merged = seqc.sequence.merge_functions.in_drop(record, record)
        print(merged)
        print(record.sequence[28:32])
        b'AAAAAAAAGAGTGATTGCTTGTGACGCCAACCCCCCCCNNNNNNTTTTT'
        b'AAAAAAAAAGAGTGATTGCTTGTGACGCCAACCCCCCCCNNNNNNTTTTT'
        b'AAAAAAAAAAGAGTGATTGCTTGTGACGCCAACCCCCCCCNNNNNNTTTTT'
        b'AAAAAAAAAAAGAGTGATTGCTTGTGACGCCAACCCCCCCCNNNNNNTTTTT'


class TestThreeBitEquivalence(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cb1_file = seqc_dir + 'data_test/in_drop/barcodes/cb1.txt'
        cb2_file = seqc_dir + 'data_test/in_drop/barcodes/cb2.txt'

        with open(cb1_file) as f1, open(cb2_file) as f2:
            cls.cb1_codes = [c.strip() for c in f1.readlines()]
            cls.cb2_codes = [c.strip() for c in f2.readlines()]

    def test_three_bit_version_1(self):
        """
        :return:
        """

        # test conversion of single cb1
        for _ in range(100):
            case = random.choice(self.cb1_codes).encode()
            coded = seqc.sequence.encodings.DNA3Bit.encode(case)
            decoded = seqc.sequence.encodings.DNA3Bit.decode(coded)
            self.assertEqual(case, decoded)

        # test conversion of single cb2
        for _ in range(100):
            case = random.choice(self.cb2_codes).encode()
            coded = seqc.sequence.encodings.DNA3Bit.encode(case)
            decoded = seqc.sequence.encodings.DNA3Bit.decode(coded)
            self.assertEqual(case, decoded)

        # test conversion of merged barcodes
        for _ in range(100):
            case1 = random.choice(self.cb1_codes).encode()
            case2 = random.choice(self.cb2_codes).encode()
            coded = seqc.sequence.encodings.DNA3Bit.encode(case1 + case2)
            decoded = seqc.sequence.encodings.DNA3Bit.decode(coded)
            self.assertEqual(case1 + case2, decoded)

    def test_three_bit_version_2(self):
        """
        :return:
        """

        # test conversion of single cb1
        for _ in range(100):
            case = random.choice(self.cb1_codes)
            coded = seqc.sequence.encodings.ThreeBit.str2bin(case)
            decoded = seqc.sequence.encodings.ThreeBit.bin2str(coded)
            self.assertEqual(case, decoded)

        # test conversion of single cb2
        for _ in range(100):
            case = random.choice(self.cb2_codes)
            coded = seqc.sequence.encodings.ThreeBit.str2bin(case)
            decoded = seqc.sequence.encodings.ThreeBit.bin2str(coded)
            self.assertEqual(case, decoded)

        # test conversion of merged barcodes
        for _ in range(100):
            case1 = random.choice(self.cb1_codes)
            case2 = random.choice(self.cb2_codes)
            coded = seqc.sequence.encodings.ThreeBit.str2bin(case1 + case2)
            decoded = seqc.sequence.encodings.ThreeBit.bin2str(coded)
            self.assertEqual(case1 + case2, decoded)

    def test_mixed_encoding(self):

        # test conversion of single cb1
        for _ in range(100):
            case = random.choice(self.cb1_codes).encode()
            coded = seqc.sequence.encodings.DNA3Bit.encode(case)
            decoded = seqc.sequence.encodings.ThreeBit.bin2str(coded).encode()
            self.assertEqual(case, decoded)

        # test conversion of single cb2
        for _ in range(100):
            case = random.choice(self.cb2_codes).encode()
            coded = seqc.sequence.encodings.DNA3Bit.encode(case)
            decoded = seqc.sequence.encodings.ThreeBit.bin2str(coded).encode()
            self.assertEqual(case, decoded)

        # test conversion of merged barcodes
        for _ in range(100):
            case1 = random.choice(self.cb1_codes).encode()
            case2 = random.choice(self.cb2_codes).encode()
            coded = seqc.sequence.encodings.DNA3Bit.encode(case1 + case2)
            decoded = seqc.sequence.encodings.ThreeBit.bin2str(coded).encode()
            self.assertEqual(case1 + case2, decoded)

        # test conversion and extraction of cb1, cb2
        for _ in range(100):  # todo fails
            case1 = random.choice(self.cb1_codes)
            case2 = random.choice(self.cb2_codes)
            coded = seqc.sequence.encodings.DNA3Bit.encode(
                    case1.encode() + case2.encode())
            case1_int = seqc.sequence.encodings.ThreeBit.c1_from_codes(coded)
            case2_int = seqc.sequence.encodings.ThreeBit.c2_from_codes(coded)
            case1_decoded = seqc.sequence.encodings.ThreeBit.bin2str(case1_int)
            case2_decoded = seqc.sequence.encodings.ThreeBit.bin2str(case2_int)
            self.assertEqual(case1, case1_decoded)
            self.assertEqual(case2, case2_decoded)


if __name__ == "__main__":
    nose2.main()
