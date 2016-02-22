#!/usr/local/bin/python

import pickle
import hashlib
import sys
import numpy as np

# import numpy object

def main(_, old, new):
    with open(old, 'rb') as pickle_object:
        old_counts = pickle.load(pickle_object)
    with open(new, 'rb') as pickle_object:
        new_counts = pickle.load(pickle_object)

    # expect a coo matrix

    old_digest = hashlib.sha256(old_counts).hexdigest()
    new_digest = hashlib.sha256(new_counts).hexdigest()


if __name__ == "__main__":
    main(*sys.argv)