import nose2
import unittest
from seqc.summary import summary
from collections import OrderedDict


class TestSummary(unittest.TestCase):

    def test_render_section(self):
        s1 = summary.KVSection('test1', OrderedDict(zip(range(5, 10), 'fghij')))
        s2 = summary.KVSection('test2', OrderedDict(zip(range(5), 'abcde')))
        summary.create_summary_section((s1, s2))


if __name__ == "__main__":
    nose2.main()
