#!/usr/local/bin/python3

import pickle
import sys
import numpy as np


def main(old, new):
    with open(old, 'rb') as pickle_object:
        old_counts = pickle.load(pickle_object)
    with open(new, 'rb') as pickle_object:
        new_counts = pickle.load(pickle_object)

    assert np.array_equal(old_counts['molecules']['matrix'].data,
                          new_counts['molecules']['matrix'].data)
    assert np.array_equal(old_counts['molecules']['matrix'].row,
                          new_counts['molecules']['matrix'].row)
    assert np.array_equal(old_counts['molecules']['matrix'].col,
                          new_counts['molecules']['matrix'].col)
    assert np.array_equal(old_counts['molecules']['col_ids'],
                          new_counts['molecules']['col_ids'])
    assert np.array_equal(old_counts['molecules']['row_ids'],
                          new_counts['molecules']['row_ids'])

    assert np.array_equal(old_counts['reads']['matrix'].data,
                          new_counts['reads']['matrix'].data)
    assert np.array_equal(old_counts['reads']['matrix'].row,
                          new_counts['reads']['matrix'].row)
    assert np.array_equal(old_counts['reads']['matrix'].col,
                          new_counts['reads']['matrix'].col)
    assert np.array_equal(old_counts['reads']['col_ids'],
                          new_counts['reads']['col_ids'])
    assert np.array_equal(old_counts['reads']['row_ids'],
                          new_counts['reads']['row_ids'])

if __name__ == "__main__":
    if not len(sys.argv) == 3:
        print('usage: test_result_equal old_counts, new_counts')
    else:
        main(*sys.argv[1:])