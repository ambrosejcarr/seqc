from abc import ABCMeta, abstractmethod
from seqc import correct_errors
import regex as re
from seqc.sequence.encodings import DNA3Bit


class AbstractPlatform:

    __metaclass__ = ABCMeta

    def __init__(self, barcodes_len, check_barcodes, correct_errors_func):
        """
        Ctor for the abstract class. barcodes_len is a list of barcodes lengths, 
        check_barcodes is a flag signalling wether or not the barcodes are known apriori and so can be filtered and corrected
        correct_errors_func points to the correct function for error correction
        """
        self._correct_errors = correct_errors_func
        self._barcodes_lengths = barcodes_len
        self._check_barcodes = check_barcodes

    def factory(type):
        if type == "in_drop":
            return in_drop()
        if type == "in_drop_v2":
            return in_drop_v2()
        if type == "in_drop_v3":
            return in_drop_v3()
        if type == "drop_seq":
            return drop_seq()
        if type == "mars1_seq":
            return mars1_seq()
        if type == "mars2_seq":
            return mars2_seq()
        if type == "ten_x":
            return ten_x()
    
    @property
    def num_barcodes(self):
        """
        return the number of barcodes used by this platform
        """
        return len(self._barcodes_lengths)
        
    @property
    def check_barcodes(self):
        return self._check_barcodes
    
    @property
    def correct_errors(self):
        """Gets the error correction function sharing the same name as the platform.
        If self.name matches no error correction function, it is assumed that there is
        no function implemented and calls pass_through(), which does not correct any
        errors

        :return function:
        """
        return self._correct_errors

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
        Return a list of barcodes from the sequence. The number of barcodes and the way they are stored in 'cell' is platform dependant.
        """
        res = []
        for bc_len in reversed(self._barcodes_lengths):
            if bc_len == -1:    #Return the rest of the MSb's
                res.insert(0,seq)
                return res
            res.insert(0,seq&((1<<bc_len*DNA3Bit.bits_per_base())-1))
            seq>>=bc_len*DNA3Bit.bits_per_base()
        return res


class in_drop(AbstractPlatform):

    def __init__(self):
        AbstractPlatform.__init__(self, [-1,8], True, correct_errors.in_drop)
        
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


class in_drop_v2(AbstractPlatform):

    def __init__(self):
        AbstractPlatform.__init__(self, [-1,8], True, correct_errors.in_drop)
        
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
        merge forward and reverse in-drop v2 reads, annotating the reverse read (containing
        genomic information) with the cell barcode, rmt, and number of poly_t. Pool is left
        empty.

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


class in_drop_v3(AbstractPlatform):
    
    def __init__(self):
        AbstractPlatform.__init__(self, [-1,8], True, correct_errors.in_drop)

    def primer_length(self):
        """The appropriate value is used to approximate the min_poly_t for each platform.
        :return: appropriate primer length for in_drop_v3
        """
        return 16

    def merge_function(self, g, b):
        """
        merge forward and reverse in-drop v3 reads, annotating the reverse read (containing
        genomic information) with the rmt and number of poly_t from the
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


class drop_seq(AbstractPlatform):

    def __init__(self):
        AbstractPlatform.__init__(self, [12], False, correct_errors.drop_seq)

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


class mars1_seq(AbstractPlatform):

    def __init__(self):
        AbstractPlatform.__init__(self, [4,8], True, correct_errors.in_drop)
        
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


class mars2_seq(AbstractPlatform):

    def __init__(self):
        AbstractPlatform.__init__(self, [4,8], True, correct_errors.in_drop)
        
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


class ten_x(AbstractPlatform):

    def __init__(self):
        AbstractPlatform.__init__(self, [14], True, correct_errors.in_drop)
        
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
        rmt = b.sequence.strip()
        # bc is in a fixed position in the name; assumes 10bp indices.
        cell = g.name.strip()[-23:-9]
        g.add_annotation((b'', cell, rmt, b''))
        return g
