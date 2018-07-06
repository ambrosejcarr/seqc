from abc import ABCMeta, abstractmethod
from seqc import rmt_correction
from seqc import barcode_correction
import regex as re
from seqc.sequence.encodings import DNA3Bit
from seqc import log
import itertools
from seqc.sequence.encodings import DNA3Bit

# todo REMOVE POOL FROM ALL THESE CLASSES


class AbstractPlatform:

    __metaclass__ = ABCMeta

    def __init__(self, barcodes_len, filter_lonely_triplets=False, filter_low_count=True):
        """
        Ctor for the abstract class. barcodes_len is a list of barcodes lengths,
        check_barcodes is a flag signalling whether or not the barcodes are known apriori
        and so can be filtered and corrected. correct_errors_func points to the correct
        function for error correction. barcode_files is a list of valid barcodes files

        """
        self._barcodes_lengths = barcodes_len
        self._filter_lonely_triplets = filter_lonely_triplets
        self._filter_low_count = filter_low_count

    def factory(type):
        if type == "in_drop":
            return in_drop()
        if type == "in_drop_v2":
            return in_drop_v2()
        if type == "in_drop_v3":
            return in_drop_v3()
        if type == "in_drop_v4":
            return in_drop_v4()
        if type == "in_drop_v5":
            return in_drop_v5()
        if type == "drop_seq":
            return drop_seq()
        if type == "mars1_seq":
            return mars1_seq()
        if type == "mars2_seq":
            return mars2_seq()
        if type == "mars_germany":
            return mars_germany()
        if type == "ten_x":
            return ten_x()
        if type == "ten_x_v2":
            return ten_x_v2()

    @property
    def num_barcodes(self):
        """
        return the number of barcodes used by this platform
        """
        return len(self._barcodes_lengths)

    @property
    def filter_lonely_triplets(self):
        return self._filter_lonely_triplets

    @property
    def filter_low_count(self):
        return self._filter_low_count

    @property
    def resolve_alignments(self):
        """Gets the resolve_alignments function sharing the same name as the platform.
        If self.name matches no resolve alignments function, it is assumed that there is
        no function implemented and calls pass_through(), which does not resolve any
        alignments

        :return function:
        """
        # todo: pass-through method for resolve_alignments not implemented yet
        raise NotImplementedError

    @abstractmethod
    def merge_function(self, g, b):
        """Defines an abstract method for the merge function, which is specific to each
        platform. The merge function must be implemented separately for each instance
        of the AbstractPlatform class.
        """
        pass

    @abstractmethod
    def primer_length(self):
        """Defines an abstract method for the primer length, which is specific to each
        platform. This value must be defined for each instance of the AbstractPlatform
        class to be used for estimating the min_poly_t value at each SEQC run.
        """
        pass

    @abstractmethod
    def extract_barcodes(self, seq):
        """
        Return a list of barcodes from the sequence. The number of barcodes and the way
        they are stored in 'cell' is platform dependant.
        """
        res = []
        for bc_len in reversed(self._barcodes_lengths):
            if bc_len == -1:  # Return the rest of the MSb's
                res.insert(0, seq)
                return res
            res.insert(0, seq & ((1 << bc_len * DNA3Bit.bits_per_base()) - 1))
            seq >>= bc_len * DNA3Bit.bits_per_base()
        return res

    @abstractmethod
    def apply_barcode_correction(self, ra, barcode_files):
        """
        Apply platform specific barcode correction. The read array should be modified in place
        with status of reads with incorrect cell barcodes updated to 'cell_error'
        """
        pass

    @abstractmethod
    def apply_rmt_correction(self, ra, error_rate):
        """
        Apply platform specific RMT correction. The read array should be modified in place
        with status of reads with incorrect cell barcodes updated to 'rmt_error'
        """
        pass

class in_drop(AbstractPlatform):

    def __init__(self):
        AbstractPlatform.__init__(self, [-1, 8])

    @classmethod
    def check_spacer(cls, sequence):
        """a fast, in-drop-v1 specific command to find a spacer sequence and cb length

        :param sequence: fastq sequence data
        :returns: (cell_barcode, rmt, poly_t)
        """
        assert ~sequence.endswith(b'\n')
        identifier = sequence[24:28]
        if identifier == b'CGCC':
            cb1 = sequence[:8]
            cb2 = sequence[30:38]
            rmt = sequence[38:44]
            poly_t = sequence[44:]
        elif identifier == b'ACGC':
            cb1 = sequence[:9]
            cb2 = sequence[31:39]
            rmt = sequence[39:45]
            poly_t = sequence[45:]
        elif identifier == b'GACG':
            cb1 = sequence[:10]
            cb2 = sequence[32:40]
            rmt = sequence[40:46]
            poly_t = sequence[46:]
        elif identifier == b'TGAC':
            cb1 = sequence[:11]
            cb2 = sequence[33:41]
            rmt = sequence[41:47]
            poly_t = sequence[47:]
        else:
            return b'', b'', b''
        cell = cb1 + cb2
        return cell, rmt, poly_t

    def primer_length(self):
        """The appropriate value is used to approximate the min_poly_t for each platform.
        :return: appropriate primer length for in_drop_v1
        """
        return 47

    def merge_function(self, g, b):
        """
        merge forward and reverse in-drop v1 reads, annotating the reverse read
        (containing genomic information) with the cell barcode, rmt, and number of
        poly_t. Pool is left empty.

        :param g: genomic fastq sequence data
        :param b: barcode fastq sequence data
        :return: annotated genomic sequence.
        """
        pattern = re.compile(b'(.{8,11}?)(GAGTGATTGCTTGTGACGCCTT){s<=2}(.{8})(.{6})(.*?)')
        cell, rmt, poly_t = self.check_spacer(b.sequence[:-1])
        if not cell:
            try:
                cell1, spacer, cell2, rmt, poly_t = re.match(
                    pattern, b.sequence[:-1]).groups()
                cell = cell1 + cell2
            except AttributeError:
                cell, rmt, poly_t = b'', b'', b''
        g.add_annotation((b'', cell, rmt, poly_t))
        return g

    def apply_barcode_correction(self, ra, barcode_files):
        """
        Apply barcode correction and return error rate

        :param ra: Read array
        :param barcode_files: Valid barcodes files
        :returns: Error rate table

        """
        error_rate = barcode_correction.in_drop(ra, self, barcode_files, max_ed=1)
        return error_rate

    def apply_rmt_correction(self, ra, error_rate):
        """
        Apply RMT correction

        :param ra: Read array
        :param error_rate: Error rate table from apply_barcode_correction

        """
        rmt_correction.in_drop(ra, error_rate)


class in_drop_v2(AbstractPlatform):

    def __init__(self):
        AbstractPlatform.__init__(self, [-1, 8])

    @classmethod
    def check_spacer(cls, sequence):
        """a fast, in-drop-v2 specific command to find a spacer sequence and cb length

        :param sequence: fastq sequence data
        :returns: (cell_barcode, rmt, poly_t)
        """
        assert ~sequence.endswith(b'\n')
        identifier = sequence[24:28]
        if identifier == b'CGCC':
            cb1 = sequence[:8]
            cb2 = sequence[30:38]
            rmt = sequence[38:46]
            poly_t = sequence[46:]
        elif identifier == b'ACGC':
            cb1 = sequence[:9]
            cb2 = sequence[31:39]
            rmt = sequence[39:47]
            poly_t = sequence[47:]
        elif identifier == b'GACG':
            cb1 = sequence[:10]
            cb2 = sequence[32:40]
            rmt = sequence[40:48]
            poly_t = sequence[48:]
        elif identifier == b'TGAC':
            cb1 = sequence[:11]
            cb2 = sequence[33:41]
            rmt = sequence[41:49]
            poly_t = sequence[49:]
        else:
            return b'', b'', b''
        cell = cb1 + cb2
        return cell, rmt, poly_t

    def primer_length(self):
        """The appropriate value is used to approximate the min_poly_t for each platform.
        :return: appropriate primer length for in_drop_v2
        """
        return 49

    def merge_function(self, g, b):
        """
        merge forward and reverse in-drop v2 reads, annotating the reverse read
        (containing genomic information) with the cell barcode, rmt, and number of poly_t.
        Pool is left empty.

        :param g: genomic fastq sequence data
        :param b: barcode fastq sequence data
        :return: annotated genomic sequence.
        """
        pattern = re.compile(b'(.{8,11}?)(GAGTGATTGCTTGTGACGCCAA){s<=2}(.{8})(.{8})(.*?)')
        cell, rmt, poly_t = self.check_spacer(b.sequence[:-1])
        if not cell:
            try:
                cell1, spacer, cell2, rmt, poly_t = re.match(
                    pattern, b.sequence[:-1]).groups()
                cell = cell1 + cell2
            except AttributeError:
                cell, rmt, poly_t = b'', b'', b''
        g.add_annotation((b'', cell, rmt, poly_t))
        return g

    def apply_barcode_correction(self, ra, barcode_files):
        """
        Apply barcode correction and return error rate

        :param ra: Read array
        :param barcode_files: Valid barcodes files
        :returns: Error rate table

        """
        error_rate = barcode_correction.in_drop(ra, self, barcode_files, max_ed=2)
        return error_rate

    def apply_rmt_correction(self, ra, error_rate):
        """
        Apply RMT correction

        :param ra: Read array
        :param error_rate: Error rate table from apply_barcode_correction

        """
        rmt_correction.in_drop(ra, error_rate)


class in_drop_v3(AbstractPlatform):

    def __init__(self):
        AbstractPlatform.__init__(self, [-1, 8])

    def primer_length(self):
        """The appropriate value is used to approximate the min_poly_t for each platform.
        :return: appropriate primer length for in_drop_v3
        """
        return 16

    def merge_function(self, g, b):
        """
        merge forward and reverse in-drop v3 reads, annotating the reverse read
        (containing genomic information) with the rmt and number of poly_t from the
        forward read. Pool is left empty, and the cell barcode is reconstructed from the
        second index and the second barcode.

        Please note that R1 is genomic, and R2 is barcode, unlike previous iterations

        :param g: genomic fastq sequence data
        :param b: barcode fastq sequence data
        :return: annotated genomic sequence.
        """
        seq = b.sequence.strip()
        cell2 = seq[:8]
        rmt = seq[8:16]
        poly_t = seq[16:]
        # bc is in a fixed position in the name; assumes 8bp indices.
        cell1 = g.name.strip()[-17:-9]
        g.add_annotation((b'', cell1 + cell2, rmt, poly_t))
        return g

    def apply_barcode_correction(self, ra, barcode_files):
       raise NotImplementedError

    def apply_rmt_correction(self, ra, error_rate):
       raise NotImplementedError


class in_drop_v4(AbstractPlatform):

    def __init__(self):
        AbstractPlatform.__init__(self, [-1, 8])

    def primer_length(self):
        """The appropriate value is used to approximate the min_poly_t for each platform.
        :return: appropriate primer length for in_drop_v3
        """
        return 28

    def merge_function(self, g, b):
        """
        merge forward and reverse in-drop v3 reads, annotating the reverse read
        (containing genomic information) with the rmt and number of poly_t from the
        forward read. Pool is left empty, and the cell barcode is reconstructed from the
        second index and the second barcode.

        Please note that R1 is genomic, and R2 is barcode, unlike previous iterations

        :param g: genomic fastq sequence data
        :param b: barcode fastq sequence data
        :return: annotated genomic sequence.
        """
        seq = b.sequence.strip()
        cell1 = seq[:8]
        cell2 = seq[12:20]
        rmt = seq[20:28]
        poly_t = seq[28:]
        g.add_annotation((b'', cell1 + cell2, rmt, poly_t))
        return g

    def apply_barcode_correction(self, ra, barcode_files):
        """
        Apply barcode correction and return error rate

        :param ra: Read array
        :param barcode_files: Valid barcodes files
        :returns: Error rate table
        """
        error_rate = barcode_correction.in_drop(ra, self, barcode_files, max_ed=2)
        return error_rate

    def apply_rmt_correction(self, ra, error_rate):
        """
        Apply RMT correction

        :param ra: Read array
        :param error_rate: Error rate table from apply_barcode_correction
        """
        rmt_correction.in_drop(ra, error_rate)


class in_drop_v5(AbstractPlatform):

    def __init__(self, potential_barcodes=None):
        AbstractPlatform.__init__(self, [-1, 8])
        self.potential_barcodes = potential_barcodes
        if self.potential_barcodes is not None:
            self.potential_encoded_bcs = set(DNA3Bit.encode(pb) for pb in self.potential_barcodes)

    @classmethod
    def check_spacer(cls, sequence):
        """a fast, in-drop-v5 specific command to find a spacer sequence and cb1 length

        :param sequence: fastq sequence data
        :returns: (cb1, rest), where rest includes cb2, rmt, poly_t
        """
        assert ~sequence.endswith(b'\n')
        identifier = sequence[24:28]
        if identifier == b'CGCC':
            cb1 = sequence[:8]
            rest = sequence[30:]
        elif identifier == b'ACGC':
            cb1 = sequence[:9]
            rest = sequence[31:]
        elif identifier == b'GACG':
            cb1 = sequence[:10]
            rest = sequence[32:]
        elif identifier == b'TGAC':
            cb1 = sequence[:11]
            rest = sequence[33:]
        else:
            return b'', b''

        return cb1, rest

    def check_cb2(self, rest):
        """a fast, in-drop-v5 specific command to find determine len of cb2

        :param rest: sequence remaining after spacer returned by check_spacer
        :returns: (cb2, rmt, poly_t)
        """

        # Check for cb2 length
        if rest[:8] in self.potential_barcodes:
            cb2 = rest[:8]
            rmt = rest[8:16]
            poly_t = rest[16:]
        elif rest[:9] in self.potential_barcodes:
            cb2 = rest[:9]
            rmt = rest[9:17]
            poly_t = rest[17:]
        else:
            return b'', b'', b''

        return cb2, rmt, poly_t

    @classmethod
    def build_cb2_barcodes(cls, barcode_files, max_ed=1):
        """
        build a set of valid and invalid barcodes used to determine
        length of cb2 in self.check_cb2

        :param barcode_files: Valid barcodes files
        :param max_ed: number of allowable mismatches
        :returns: new class with potential barcodes set
        """

        # Build set of all potential correct and incorrect cb2
        potential_barcodes = set()
        cb2_file = barcode_files[1]
        with open(cb2_file, 'r') as f:
            valid_barcodes = set([line.strip() for line in f.readlines()])
        # This will work for any number of allowable mismatches
        for bc in valid_barcodes:
            potential_barcodes.add(bc)
            for inds in itertools.combinations(range(len(bc)), max_ed):
                invalid_bc = [[nt] for nt in bc]
                for ind in inds:
                    valid_nt = bc[ind]
                    invalid_bc[ind] = [nt for nt in ['A', 'C', 'G', 'T', 'N'] if nt != valid_nt]
                for mut in itertools.product(*invalid_bc):
                    potential_barcodes.add(''.join(mut))
        potential_barcodes = set([pb.encode() for pb in potential_barcodes])

        return cls(potential_barcodes=potential_barcodes)

    def primer_length(self):
        """The appropriate value is used to approximate the min_poly_t for each platform.
        :return: appropriate primer length for in_drop_v5
        """
        return 28

    def merge_function(self, g, b):
        """
        merge forward and reverse in-drop v2 reads, annotating the reverse read
        (containing genomic information) with the cell barcode, rmt, and number of poly_t.
        Pool is left empty.

        :param g: genomic fastq sequence data
        :param b: barcode fastq sequence data
        :return: annotated genomic sequence.
        """
        cb1, rest = self.check_spacer(b.sequence[:-1])
        if not cb1:
             cell, rmt, poly_t = b'', b'', b''
        else:
            cb2, rmt, poly_t = self.check_cb2(rest)
            if not cb2:
                cell = b''
            else:
                cell = cb1 + cb2

        g.add_annotation((b'', cell, rmt, poly_t))
        return g

    def extract_barcodes(self, seq):
        """
        Return a list of barcodes from the sequence. A bit hacky right now.
        Specific to v5 platform.
        """
        res = []
        for bc_len in reversed(self._barcodes_lengths):
            if bc_len == -1:  # Return the rest of the MSb's
                res.insert(0, seq)
                return res
            potent_cb2 = seq & ((1 << bc_len * DNA3Bit.bits_per_base()) - 1)
            # First assume it is length 8 through self._barcodes_lengths, then
            # if it isn't in potentials, assume 9.
            if potent_cb2 not in self.potential_encoded_bcs:
                    potent_cb2 = seq & ((1 << 9 * DNA3Bit.bits_per_base()) - 1)
                    bc_len = 9
            res.insert(0, seq & ((1 << bc_len * DNA3Bit.bits_per_base()) - 1))
            seq >>= bc_len * DNA3Bit.bits_per_base()

        return res

    def apply_barcode_correction(self, ra, barcode_files):
        """
        Apply barcode correction and return error rate

        :param ra: Read array
        :param barcode_files: Valid barcodes files
        :returns: Error rate table

        """
        error_rate = barcode_correction.in_drop(ra, self, barcode_files, max_ed=1)
        return error_rate

    def apply_rmt_correction(self, ra, error_rate):
        """
        Apply RMT correction

        :param ra: Read array
        :param error_rate: Error rate table from apply_barcode_correction

        """
        rmt_correction.in_drop(ra, error_rate)


class drop_seq(AbstractPlatform):

    def __init__(self):
        AbstractPlatform.__init__(self, [12])

    def primer_length(self):
        """The appropriate value is used to approximate the min_poly_t for each platform.
        :return: appropriate primer length for drop_seq.
        """
        return 20

    def merge_function(self, g, b):
        """
        merge forward and reverse drop-seq reads, annotating the reverse read (containing
        genomic information) with the cell barcode and rmt. Number of poly_t and pool
        fields are left empty.

        :param g: genomic fastq sequence data
        :param b: barcode fastq sequence data
        :return: annotated genomic sequence.
        """
        cell = b.sequence[:12]
        rmt = b.sequence[12:20]
        poly_t = b.sequence[20:-1]
        g.add_annotation((b'', cell, rmt, poly_t))
        return g

    def apply_barcode_correction(self, ra, barcode_files):
        """
        Apply barcode correction and return error rate

        :param ra: Read array
        :param barcode_files: Valid barcodes files
        :returns: Error rate table

        """
        barcode_correction.drop_seq(ra)
        return None

    def apply_rmt_correction(self, ra, error_rate):
        """
        Apply RMT correction

        :param ra: Read array
        :param error_rate: Error rate table from apply_barcode_correction

        """
        log.info('Drop-seq barcodes do not support RMT correction')


class mars1_seq(AbstractPlatform):

    def __init__(self):
        AbstractPlatform.__init__(self, [4, 8], True, False)

    def primer_length(self):
        """The appropriate value is used to approximate the min_poly_t for each platform.
        :return: appropriate primer length for mars1_seq.
        """
        # todo: we don't have pre-demultiplexed data for mars1-seq
        # raise NotImplementedError
        return None

    def merge_function(self, g, b):
        """
        re-annotate reverse mars-seq reads in a format consistent with other SEQC platforms,
        annotating the reverse read (containing genomic information) with the pool, cell
        barcode, rmt, number of poly_t.

        :param g: genomic fastq sequence data
        :param b: barcode fastq sequence data
        :return: annotated genomic sequence.
        """

        *name_fields, pool, cell, rmt = g.name[1:-1].split(b':')
        g.name = (b'@' + b':'.join((pool, cell, rmt, b'')) + b';' +
                  b':'.join(name_fields) + b'\n')
        return g

    def apply_barcode_correction(self, ra, barcode_files):
        """
        Apply barcode correction and return error rate

        :param ra: Read array
        :param barcode_files: Valid barcodes files
        :returns: Error rate table

        """
        # todo: verify max edit distance
        error_rate = barcode_correction.in_drop(ra, self, barcode_files, max_ed=0)
        return error_rate

    def apply_rmt_correction(self, ra, error_rate):
        """
        Apply RMT correction

        :param ra: Read array
        :param error_rate: Error rate table from apply_barcode_correction

        """
        log.info('Mars-seq barcodes do not support RMT correction')


class mars2_seq(AbstractPlatform):

    def __init__(self):
        AbstractPlatform.__init__(self, [4, 8], True, False)

    def primer_length(self):
        """The appropriate value is used to approximate the min_poly_t for each platform.
        :return: appropriate primer length for mars2_seq.
        """
        return 15

    def merge_function(self, g, b):
        """
        re-annotate reverse mars-seq v2 reads in a format consistent with other SEQC
        platforms, annotating the reverse read (containing genomic information) with
        the pool, cell barcode, rmt, number of poly_t.

        :param g: genomic fastq sequence data
        :param b: barcode fastq sequence data
        :return: annotated genomic sequence.
        """
        pool = g.sequence.strip()[5:9]
        seq = b.sequence.strip()
        cell = seq[:7]
        rmt = seq[7:15]
        poly_t = seq[15:]
        g.add_annotation((b'', pool + cell, rmt, poly_t))
        return g

    def apply_barcode_correction(self, ra, barcode_files):
        """
        Apply barcode correction and return error rate

        :param ra: Read array
        :param barcode_files: Valid barcodes files
        :returns: Error rate table

        """
        # todo: verify max edit distance
        error_rate = barcode_correction.in_drop(ra, self, barcode_files, max_ed=0)
        return error_rate

    def apply_rmt_correction(self, ra, error_rate):
        """
        Apply RMT correction

        :param ra: Read array
        :param error_rate: Error rate table from apply_barcode_correction

        """
        log.info('Mars-seq barcodes do not support RMT correction')


class mars_germany(AbstractPlatform):

    def __init__(self):
        AbstractPlatform.__init__(self, [10], True, False)

    def primer_length(self):
        return 15

    def merge_function(self, g, b):
        pool = g.sequence.strip()[3:7]  # 4 bp
        g.sequence = g.sequence.strip()[7:] + b'\n'  # strip() is necessary in case there is a truncated read. \n=good, \n\n=bad
        # Need to skip over the quality as well
        g.quality = g.quality.strip()[7:] + b'\n'
        seq = b.sequence.strip()
        cell = seq[:6]  # 6 bp
        rmt = seq[6:12]  # 6 bp
        g.add_annotation((b'', pool + cell, rmt, b''))
        return g

    def apply_barcode_correction(self, ra, barcode_files):
        """
        Apply barcode correction and return error rate

        :param ra: Read array
        :param barcode_files: Valid barcodes files
        :returns: Error rate table

        """
        # todo: verify max edit distance
        error_rate = barcode_correction.in_drop(ra, self, barcode_files, max_ed=0)
        return error_rate

    def apply_rmt_correction(self, ra, error_rate):
        """
        Apply RMT correction

        :param ra: Read array
        :param error_rate: Error rate table from apply_barcode_correction

        """
        log.info('Mars-seq barcodes do not support RMT correction')


class ten_x(AbstractPlatform):
    # 10X version 1 chemistry

    def __init__(self):
        AbstractPlatform.__init__(self, [14])

    def primer_length(self):
        """The appropriate value is used to approximate the min_poly_t for each platform.
        :return: appropriate primer length for 10X
        """
        return 10

    def merge_function(self, g, b):
        """
        merge forward and reverse 10x reads, annotating the reverse read
        (containing genomic information) with the rmt from the forward read.
        Pool is left empty, and the cell barcode is obtained from the
        name field of the forward read.

        Please note that R1 is genomic, and R2 is RMT, unlike previous iterations

        :param g: genomic fastq sequence data
        :param b: barcode fastq sequence data
        :return: annotated genomic sequence.
        """
        combined = b.sequence.strip()
        rmt = combined[:10]
        # bc is in a fixed position in the name; assumes 10bp indices.
        cell = g.name.strip()[-23:-9]
        poly_t = combined[10:]
        g.add_annotation((b'', cell, rmt, poly_t))
        return g

    def apply_barcode_correction(self, ra, barcode_files):
        error_rate = barcode_correction.ten_x_barcode_correction(ra, self, barcode_files, max_ed=0)
        return error_rate

    def apply_rmt_correction(self, ra, error_rate):
        rmt_correction.in_drop(ra, error_rate=0.02)


class ten_x_v2(AbstractPlatform):
    # 10X version 2 chemistry

    def __init__(self):
        AbstractPlatform.__init__(self, [16])

    def primer_length(self):
        """The appropriate value is used to approximate the min_poly_t for each platform.
        :return: appropriate primer length for 10X
        """
        return 26

    def merge_function(self, g, b):
        """
        merge forward and reverse 10x reads, annotating the reverse read
        (containing genomic information) with the rmt from the forward read.
        Pool is left empty, and the cell barcode is obtained from the
        forward read.

        :param g: genomic fastq sequence data
        :param b: barcode fastq sequence data
        :return: annotated genomic sequence.
        """
        combined = b.sequence.strip()
        cell = combined[0:16]  # v2 chemistry has 16bp barcodes
        rmt = combined[16:26]  # 10 baselength RMT
        poly_t = combined[26:]
        g.add_annotation((b'', cell, rmt, poly_t))
        return g

    def apply_barcode_correction(self, ra, barcode_files):
        """
        Apply barcode correction and return error rate

        :param ra: Read array
        :param barcode_files: Valid barcodes files
        :returns: Error rate table

        """
        # todo: verify max edit distance
        error_rate = barcode_correction.ten_x_barcode_correction(ra, self, barcode_files, max_ed=0)
        return error_rate

    def apply_rmt_correction(self, ra, error_rate):
        """
        Apply RMT correction

        :param ra: Read array
        :param error_rate: Error rate table from apply_barcode_correction

        """
        rmt_correction.in_drop(ra, error_rate=0.02)
