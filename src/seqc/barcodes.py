__author__ = 'ambrose'

import seqc
from itertools import product
from collections import defaultdict, Mapping
import pickle
import io


class CellBarcodes:
    """"""

    def __init__(self, perfect_codes: list, error_codes: list):
        self._perfect_codes = perfect_codes
        self._error_codes = error_codes

    @property
    def perfect_codes(self):
        return self._perfect_codes

    @property
    def error_codes(self):
        return self._error_codes

    @classmethod
    def from_files(cls, *barcode_files, reverse_complement: bool=False):
        """Create a CellBarcodes object from text files containing barcodes"""
        tbp = seqc.three_bit.ThreeBit

        perfect_codes = cls.load_barcodes(barcode_files, reverse_complement)
        current_codes = cls.add_single_error(perfect_codes)

        # todo
        # perfectly possible to store all 2-error codes in int form, but need to build
        # it in int form as well; building in string form takes way too much memory.
        # going to put this off for a later time. for now, using only 1-base errors

        # error_codes = dict(deepcopy(current_codes))
        # for code, error_data in current_codes.items():
        #     try:
        #         error, pos = error_data[0]
        #     except IndexError:
        #         continue  # this is an original code, all errors already entered
        #     for i, original in enumerate(code):
        #         for e in alphabet:
        #             temp = list(code)  # mutable string
        #             if e == error[-1] and i == pos:  # don't revert original error
        #                 continue
        #             temp[i] = e
        #             try:
        #                 error_codes[(''.join(temp))].append((original + e, i))
        #             except KeyError:
        #                 error_codes[(''.join(temp))] = [(original + e, i)]


        perfect_codes = cls.codes2bin(perfect_codes, tbp)
        error_codes = cls.codes2bin(current_codes, tbp)

        return cls(perfect_codes, error_codes)

    @staticmethod
    def from_dtype(data_type, *barcode_files, reverse_complement: bool=False):
        if data_type == 'drop_seq':
            return DropSeqCellBarcodes()
        else:
            return CellBarcodes.from_files(
                *barcode_files, reverse_complement=reverse_complement)

    def pickle(self, fname: str):
        """save a serialized version of this object to file"""
        with open(fname, 'wb') as f:
            data = {'perfect': self.perfect_codes, 'current': self.current_codes}
            pickle.dump(data, f)

    @staticmethod
    def from_pickle(fname: str):
        """recreate a barcodes file from serialized data saved using self.pickle()"""
        with open(fname, 'rb') as f:
            data = pickle.load(f)
        if not all(k in data for k in ['perfect', 'current']):
            raise ValueError('file "%s" does not contain CellBarcodes data. Must contain '
                             'both perfect and current barcode lists.' % fname)
        if not data['perfect'] and not data['current']:  # drop-seq data; no barcodes
            return DropSeqCellBarcodes()
        else:
            return CellBarcodes(data['perfect'], data['current'])

    def perfect_match(self, barcode):
        return 1 if barcode in self.perfect_codes else 0

    def close_match(self, barcode):
        return 1 if barcode in self.error_codes else 0

    def map_errors(self, barcode):
        """
        Note: return of None is somewhat counter-intuitive when the barcode has
        TOO MANY errors, but that's what it does.
        """
        try:
            return self.error_codes[barcode]
        except KeyError:
            return None

    _revcomp = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}

    @classmethod
    def revcomp(cls, s):
        return ''.join(cls._revcomp[n] for n in s[::-1])

    @classmethod
    def load_barcodes(cls, barcode_files, reverse_complement=False):
        lists = []
        for blist in barcode_files:
            if isinstance(blist, io.TextIOBase):
                lists.append([bc.strip() for bc in blist.readlines()])
            elif isinstance(blist, list):
                lists.append([bc.strip() for bc in blist])
            elif isinstance(blist, str):
                with open(blist) as f:
                    lists.append([bc.strip() for bc in f.readlines()])

        if reverse_complement:
            rc_lists = []
            for l in lists:
                rc_lists.append([cls.revcomp(bc) for bc in l])
            lists = rc_lists

        # if only one list in lists, make a set
        if len(lists) == 1:
            perfect_codes = set(lists[  0])
        else:
            perfect_codes = set(''.join(r) for r in product(*lists))

        return perfect_codes

    @staticmethod
    def add_single_error(perfect_codes):
        alphabet = 'ACGTN'
        current_codes = defaultdict(list)
        # add existing barcodes
        for c in perfect_codes:
            current_codes[c] = []

        for code in perfect_codes:
            for i, original in enumerate(code):
                for e in alphabet:
                    if original == e:
                        continue
                    temp = list(code)  # mutable
                    temp[i] = e
                    current_codes[(''.join(temp))].append((original + e, i))

        return current_codes

    @staticmethod
    def codes2bin(codes, tbp):
        if isinstance(codes, set):
            bin_codes = []
            for code in codes:
                bin_codes.append(tbp.str2bin(code))
            bin_codes = set(bin_codes)
        elif isinstance(codes, Mapping):
            bin_codes = {}
            for k, v in codes.items():  # current_codes --> error_codes when doing 2e
                bin_codes[tbp.str2bin(k)] = tuple(e for (e, pos) in v)
        else:
            raise TypeError('codes must be a dict or set of barcodes')
        return bin_codes


class DropSeqCellBarcodes(CellBarcodes):

    def __init__(self):
        super().__init__([], [])

    def close_match(self, barcode):
        return 1

    def perfect_match(self, barcode):
        return 1

    def pickle(self, fname: str):
        """save a serialized version of this object to file"""
        with open(fname, 'wb') as f:
            data = {'perfect': [], 'current': []}
            pickle.dump(data, f)
