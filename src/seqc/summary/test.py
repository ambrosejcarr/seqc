import nose2
import unittest
from seqc.summary import summary
from collections import OrderedDict


class TestSummary(unittest.TestCase):

    def test_render_section(self):
        s1 = summary.Section.from_alignment_summary(
            '/var/folders/y3/ysxvl2w921d881nfpvx5ypvh0000gn/T/seqc/test_no_aws_in_drop_v2'
            '/alignment_summary.txt')
        s1.render('./src/seqc/summary/test_summary.html')

if __name__ == "__main__":
    nose2.main()
