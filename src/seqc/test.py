import nose2
import unittest
import seqc

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


if __name__ == "__main__":
    nose2.main()
