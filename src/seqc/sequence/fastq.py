from io import StringIO
import os
import gzip
import random
from collections.abc import Iterable, Mapping
from collections import Counter
import numpy as np
import bz2
from functools import lru_cache
import seqc


class FastqRecord:

    __slots__ = ['_data']

    def __init__(self, record: [bytes, bytes, bytes, bytes]):
        self._data = list(record)

    @property
    def name(self) -> bytes:
        return self._data[0]

    @name.setter
    def name(self, value: bytes):
        self._data[0] = value

    @property
    def sequence(self) -> bytes:
        return self._data[1]

    @sequence.setter
    def sequence(self, value: bytes):
        self._data[1] = value

    @property
    def name2(self) -> bytes:
        return self._data[2]

    @name2.setter
    def name2(self, value: bytes):
        self._data[2] = value

    @property
    def quality(self) -> bytes:
        return self._data[3]

    @quality.setter
    def quality(self, value: bytes):
        self._data[3] = value

    def __bytes__(self) -> bytes:
        return b''.join(self._data)

    def __str__(self) -> str:
        return bytes(self).decode()

    def __len__(self) -> int:
        return len(self.sequence)

    @property
    def annotations(self) -> list:
        """
        returns:
        --------
        list of annotations present in the fastq header
        """
        try:
            end = self.name.index(b';')
            return self.name[:end].split(b':')
        except ValueError:
            return []

    @property
    def metadata(self) -> dict:
        """
        returns:
        --------
        dictionary of annotations and fields, if any are present"""
        try:
            start = self.name.rindex(b'|')
        except ValueError:
            return {}
        fields = {}
        for field in self.name[start + 1:].split(b':'):
            k, v = field.split(b'=')
            fields[k] = v
        return fields

    def add_annotation(self, values: Iterable) -> None:
        """prepends a list of annotations to the name field of self.name"""
        self._data[0] = b'@' + b':'.join(values) + b';' + self.name[1:]

    def add_metadata(self, values: Mapping) -> None:
        """appends a list of metadata fields to the name field of self.name"""
        self.name += b'|' + b':'.join(k + '=' + v for k, v in values.items())

    def average_quality(self) -> int:
        """"""
        return np.mean(np.frombuffer(self.quality, dtype=np.int8, count=len(self)))\
            .astype(int) - 33


class Reader(seqc.reader.Reader):
    """simple Reader Class, designed for inheritance across data types"""

    @staticmethod
    def record_grouper(iterable):
        args = [iter(iterable)] * 4
        return zip(*args)

    def __iter__(self):
        for f in self._files:
            if f.endswith('.gz'):
                file_input = gzip.open(f, 'rb')
            elif f.endswith('.bz2'):
                file_input = bz2.open(f, 'rb')
            else:
                file_input = open(f, 'rb')
            for record in self.record_grouper(file_input):
                yield FastqRecord(record)

    @lru_cache(maxsize=1)
    def __len__(self):
        """
        return the length of the Reader object. This depends on the implementation of
        self.__iter__(); it does not necessarily represent the length of the file in
        lines.
        """
        return sum(1 for _ in self) / 4


def merge_paired(merge_function, fout, genomic, barcode=None):
    """
    annotate genomic fastq with barcode information from reverse read

    args:
    -----
    merge_function: function from merge_functions.py
    fout: merged output file name
    genomic: fastq containing genomic data
    barcode: fastq containing barcode data
    """
    directory, filename = os.path.split(fout)
    if not os.path.isdir(directory):
        os.makedirs(directory, exist_ok=True)
    genomic = Reader(genomic)
    if barcode:
        barcode = Reader(barcode)
        with open(fout, 'wb') as f:
            for g, b in zip(genomic, barcode):
                r = merge_function(g, b)
                f.write(bytes(r))
    else:
        with open(fout, 'wb') as f:
            for g in genomic:
                r = merge_function(g)
                f.write(bytes(r))

    return fout


def truncate_read_length(fastq_files, length):
    for f in fastq_files:
        rd = Reader(f)
        with open(f.replace('.fastq', '_truncated_{}.fastq'.format(length)), 'wb') as o:
            for record in rd:
                record.sequence = record.sequence[:length] + b'\n'
                record.quality = record.quality[:length] + b'\n'
                o.write(bytes(record))


class GenerateSynthetic:

    # define some general constants
    # _alphabet = ['A', 'C', 'G', 'T']

    # _platforms = {
    #     'in_drop': {'cell': None, 'rmt': 6, 'pool': 0, 'blength': 50, 'glength': 100},
    #     'in_drop_v2': {'cell': None, 'rmt': 8, 'pool': 0, 'blength': 55, 'glength': 95}
    # }

    def __init__(self, gtf: str, fasta: str):
        """
        :param gtf:  annotations.gtf filename
        :param fasta:  genome.fa filename
        :return: None
        """
        self._gtf = gtf
        self._fasta = fasta

        self._genome = None
        self._annotation = None

    def add_alignable(self, barcode, genomic, length: int, n: int
                      ) -> (str, str):
        """
        generate n random sequence with length=length from annotation. Please not that
        if any restriction on the location of sequences is desired, the interval object
        must be sliced a priori

        :param gtf: annotations.gtf filename
        :param barcode: barcodes.fastq file
        :param genomic: genomic.fastq file
        :param n: number of sequences to generate
        :return: seq: str, scid: str
        """

        if not self._genome:
            self._genome = seqc.sequence.fasta.Genome.from_file(self._fasta)
        if not self._annotation:
            self._annotation = seqc.gtf.SCID_set(self._gtf)

        # restrict to final 1000 bases
        self._annotation.slice_SCIDs(-1000, None)

        all_genes = list(self._annotation.values())
        n_genes = len(all_genes) - 1
        s = 0
        genes = []
        seqs = []
        while s < n:
            gene = all_genes[random.randint(0, n_genes)]
            ivs = gene.intervals()
            i = random.randint(0, len(ivs) - 1)
            istart, istop, ichr, istrand = ivs[i]

            assert isinstance(istrand, bytes)  # make sure not string (silent failures)

            iend = istop - length
            if istart < iend:
                pos = random.randint(istart, iend - 1)
                sequence = self._genome[ichr].sequence[pos:pos + length]
                if istrand == b'-':
                    sequence = seqc.sequence.revcomp(sequence)
                seqs.append(sequence)
                genes.append(gene.id)
                s += 1

        # create a pair of dummy forward fastq record
        forward_fastq1 = (
            '@{index}:{gene}\n'
            'AAAAAAAAAAGAGTGATTGCTTGTGACGCCAATTTTTTTT{rmt}TTTTT\n'
            '+\n'
            'FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF{rmt_len}FFFFF\n')
        forward_fastq2 = (
            '@{index}:{gene}\n'
            'AAAAACCCCCGAGTGATTGCTTGTGACGCCAATTTTGGGG{rmt}TTTTT\n'
            '+\n'
            'FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF{rmt_len}FFFFF\n')

        # create a dummy reverse fastq record
        reverse_fastq = (
            '@{index}:{gene}\n'
            '{sequence!s}\n'
            '+\n' +
            ('F' * length) + '\n'
        )

        # create directories
        for filepath in (barcode, genomic):
            stem, leaf = os.path.split(filepath)
            if stem:
                os.makedirs(stem, exist_ok=True)

        # write fastq
        threshold = length // 2
        with open(barcode, 'a') as barcode, open(genomic, 'a') as genomic:
            for i, (s, g) in enumerate(zip(seqs, genes)):
                if i < threshold:
                    barcode.write(forward_fastq1.format(index=i, gene=g, rmt='ACGTCAAC',
                                                        rmt_len='FFFFFFFF'))
                else:
                    barcode.write(forward_fastq2.format(index=i, gene=g, rmt='ACGTTTTT',
                                                        rmt_len='FFFFFFFF'))
                genomic.write(reverse_fastq.format(index=i, gene=g, sequence=s.decode()))

    @classmethod
    def genomic(cls, platform: str, n: int, fout: str, gtf: str):
        """
        generates n genomic fastq reads from gtf

        :param platform: the platform from which to generate reads
        :param n: the number of reads to generate
        :param fout: name of the output file to be generated
        :param gtf: name of the annotations.gtf file from which to generate data
        :return: None
        """
        raise NotImplementedError


    @classmethod
    def barcode(cls):
        pass

    @classmethod
    def _forward_in_drop(cls, n, barcodes_):
        barcodes_ = seqc.barcodes.CellBarcodes.from_pickle(barcodes_)
        read_length = 50
        names = range(n)
        name2 = '+'
        quality = 'I' * read_length
        records = []
        umi_len = 6
        codes = list(barcodes_.perfect_codes)
        for name in names:
            # for now, translate barcode back into string code
            cb = random.choice(codes)
            c1, c2 = seqc.encodings.ThreeBitInDrop.split_cell(cb)
            c1, c2 = [seqc.encodings.ThreeBit.bin2str(c) for c in [c1, c2]]
            w1 = 'GAGTGATTGCTTGTGACGCCTT'
            cb = ''.join([c1, w1, c2])
            umi = ''.join(np.random.choice(cls._alphabet, umi_len))
            poly_a = (read_length - len(cb) - len(umi)) * 'T'
            records.append('\n'.join(['@%d' % name, cb + umi + poly_a, name2, quality]))
        forward_fastq = StringIO('\n'.join(records) + '\n')
        return forward_fastq

    @classmethod
    def _forward_drop_seq(cls, n, *args):  # args is for unused barcodes parameters
        names = range(n)
        name2 = '+'
        quality = 'I' * 20
        records = []
        for name in names:
            cb = ''.join(np.random.choice(cls._alphabet, 12))
            umi = ''.join(np.random.choice(cls._alphabet, 8))
            records.append('\n'.join(['@%d' % name, cb + umi, name2, quality]))
        forward_fastq = StringIO('\n'.join(records) + '\n')
        return forward_fastq

    @staticmethod
    def _reverse(n: int, read_length: int, fasta: str, gtf: str, tag_type='gene_name'):

        # read gtf
        reader = seqc.gtf.Reader(gtf)
        intervals = []
        for r in reader.iter_exons():
            end = int(r.end) - read_length
            start = int(r.start)
            if end > start:
                intervals.append((r.attribute[tag_type].upper(), start, end, r.strand))

        # pick intervals
        exon_selections = np.random.randint(0, len(intervals), (n,))

        # fasta information:
        with open(fasta) as f:
            fasta = f.readlines()[1:]
            fasta = ''.join(fasta)

        # generate sequences
        sequences = []
        tags = []
        for i in exon_selections:
            tag, start, end, strand = intervals[i]

            # get position within interval
            start = random.randint(start, end)
            end = start + read_length
            seq = fasta[start:end]
            if strand == '-':
                seq = revcomp(seq)
            sequences.append(seq)
            tags.append(tag)

        prefixes = range(n)
        name2 = '+'
        quality = 'I' * read_length
        records = []
        for name, tag, seq in zip(prefixes, tags, sequences):
            records.append('\n'.join(['@%d:%s' % (name, tag), seq, name2, quality]))
        reverse_fastq = StringIO('\n'.join(records) + '\n')
        return reverse_fastq

    @staticmethod
    def _reverse_three_prime(n: int, read_length: int, fasta: str, gtf: str,
                             tag_type='gene_name', fragment_length=1000):

        # read gtf
        reader = seqc.gtf.Reader(gtf)
        intervals = []
        for r in reader.iter_genes_final_nbases(fragment_length):
            for iv in r.intervals:
                start, end = map(int, iv)
                end -= read_length
            # note that this means I never test small intervals
            if end > start:
                intervals.append((r.attribute[tag_type].upper(), start, end, r.strand))

        # pick intervals
        exon_selections = np.random.randint(0, len(intervals), (n,))

        # fasta information:
        with open(fasta) as f:
            fasta = f.readlines()[1:]
            fasta = ''.join(fasta)

        # todo
        # filtering here with a converter. This guarantees that no successful
        # testing of the converter can be done with any data generated by this function.
        fc = seqc.convert_features.ConvertGeneCoordinates.from_gtf(gtf)

        # generate sequences
        sequences = []
        tags = []
        for i in exon_selections:
            tag, start, end, strand = intervals[i]

            # todo
            # remove this check when removing the converter.
            if not fc.translate(strand, 'chr19', start):
                continue

            # get position within interval
            start = random.randint(start, end)
            end = start + read_length
            seq = fasta[start:end]
            if strand == '-':
                seq = revcomp(seq)
            sequences.append(seq)
            tags.append(tag)

        prefixes = range(n)
        name2 = '+'
        quality = 'I' * read_length
        records = []
        for name, tag, seq in zip(prefixes, tags, sequences):
            records.append('\n'.join(['@%d:%s' % (name, tag), seq, name2, quality]))
        reverse_fastq = StringIO('\n'.join(records) + '\n')
        return reverse_fastq, len(records)

    @staticmethod
    def _reverse_split(localized: int, unlocalized: int, unalignable: int, index: str,
                       fasta: str, gtf: str, fragment_length=1000, sequence_length=100,
                       n_threads=7):

        # containers for generated data
        sequences = []
        donor_genes = []

        # read gtf
        annotation = seqc.gtf_new.Annotation(gtf, fasta)

        # generate sequences localized in expected locations
        seqs, genes = annotation.random_sequences(
            localized, sequence_length, return_genes=True,
            start=-fragment_length, stop=None)
        sequences.extend(seqs)
        donor_genes.extend([g.string_gene_id for g in genes])

        # generate sequences outside of expected locations
        seqs, genes = annotation.random_sequences(
            unlocalized, sequence_length, return_genes=True,
            start=None, stop=-fragment_length)
        sequences.extend(seqs)
        donor_genes.extend([g.string_gene_id for g in genes])

        # generate unalignable sequences
        seqs = seqc.align.STAR.generate_unalignable_sequences(
            unalignable, sequence_length, index, n_threads)
        sequences.extend(seqs)
        donor_genes.extend([b'unalignable'] * len(seqs))

        # generate the fastq file
        assert(len(sequences) == len(donor_genes))

        name1 = [b'@' + str(i).encode() + b':' + gene
                 for i, gene in enumerate(donor_genes)]
        name2 = b'+'
        quality = b'I' * sequence_length
        records = []
        for n, s in zip(name1, sequences):
            records.append(b'\n'.join([n, s, name2, quality]))
        data = (b'\n'.join(records) + b'\n').decode()
        reverse_fastq = StringIO(data)  # todo not efficient to operate on string
        return reverse_fastq, len(records)

    @classmethod
    def in_drop(cls, localized, unlocalized, unalignable, prefix, fasta, gtf, index,
                barcodes, replicates=3, n_threads=7, fragment_length=1000,
                sequence_length=100, *args, **kwargs):

        if not replicates >= 0:
            raise ValueError('Cannot generate negative replicates')

        # create directory if it doesn't exist
        directory = '/'.join(prefix.split('/')[:-1])
        if not os.path.isdir(directory):
            os.makedirs(directory)

        reverse, n = cls._reverse_split(localized, unlocalized, unalignable, index, fasta,
                                        gtf, fragment_length, sequence_length, n_threads)
        reverse = reverse.read()  # consume the StringIO object
        forward = cls._forward_in_drop(n, barcodes)
        forward = forward.read()  # consume the StringIO object
        with open(prefix + '_r1.fastq', 'w') as f:
            f.write(''.join([forward] * (replicates + 1)))
        with open(prefix + '_r2.fastq', 'w') as r:
            r.write(''.join([reverse] * (replicates + 1)))
        return prefix + '_r1.fastq', prefix + '_r2.fastq'

    # @classmethod
    # def in_drop(cls, n, prefix, fasta, gtf, barcodes, tag_type='gene_name', replicates=3,
    #             *args, **kwargs):
    #
    #     if not replicates >= 0:
    #         raise ValueError('Cannot generate negative replicates')
    #
    #     # create directory if it doesn't exist
    #     directory = '/'.join(prefix.split('/')[:-1])
    #     if not os.path.isdir(directory):
    #         os.makedirs(directory)
    #
    #     fwd_len = 50
    #     rev_len = 100
    #
    #     reverse, n = cls._reverse_three_prime(n, rev_len, fasta, gtf, tag_type=tag_type)
    #     reverse = reverse.read()  # consume the StringIO object
    #     forward = cls._forward_in_drop(n, barcodes)
    #     forward = forward.read()  # consume the StringIO object
    #     with open(prefix + '_r1.fastq', 'w') as f:
    #         f.write(''.join([forward] * (replicates + 1)))
    #     with open(prefix + '_r2.fastq', 'w') as r:
    #         r.write(''.join([reverse] * (replicates + 1)))
    #     return prefix + '_r1.fastq', prefix + '_r2.fastq'

    @classmethod
    def drop_seq(cls, localized, unlocalized, unalignable, prefix, fasta, gtf, index,
                 replicates=3, n_threads=7, fragment_length=1000, sequence_length=100,
                 *args, **kwargs):

        if not replicates >= 0:
            raise ValueError('Cannot generate negative replicates')

        # create directory if it doesn't exist
        directory = '/'.join(prefix.split('/')[:-1])
        if not os.path.isdir(directory):
            os.makedirs(directory)

        reverse, n = cls._reverse_split(localized, unlocalized, unalignable, index, fasta,
                                        gtf, fragment_length, sequence_length, n_threads)
        reverse = reverse.read()
        forward = cls._forward_drop_seq(n)
        forward = forward.read()
        with open(prefix + '_r1.fastq', 'w') as f:
            f.write(''.join([forward] * (replicates + 1)))
        with open(prefix + '_r2.fastq', 'w') as r:
            r.write(''.join([reverse] * (replicates + 1)))
        return prefix + '_r1.fastq', prefix + '_r2.fastq'

    # @classmethod
    # def drop_seq(cls, n, prefix, fasta, gtf, tag_type='gene_name', replicates=3, *args,
    #              **kwargs):
    #
    #     if not replicates >= 0:
    #         raise ValueError('Cannot generate negative replicates')
    #
    #     # create directory if it doesn't exist
    #     directory = '/'.join(prefix.split('/')[:-1])
    #     if not os.path.isdir(directory):
    #         os.makedirs(directory)
    #
    #     rev_len = 100
    #     reverse, n = cls._reverse_three_prime(n, rev_len, fasta, gtf, tag_type=tag_type)
    #     reverse = reverse.read()  # consume the StringIO object
    #     forward = cls._forward_drop_seq(n)
    #     forward = forward.read()  # consume the StringIO object
    #     with open(prefix + '_r1.fastq', 'w') as f:
    #         f.write(''.join([forward] * (replicates + 1)))
    #     with open(prefix + '_r2.fastq', 'w') as r:
    #         r.write(''.join([reverse] * (replicates + 1)))
    #     return prefix + '_r1.fastq', prefix + '_r2.fastq'

    @classmethod
    def simple_fastq(cls, n, length):
        sequences = np.random.choice(list('ACGT'), n * length)
        sequences = np.reshape(sequences, (n, length))
        qualities = np.random.choice(np.arange(30, 40), n * length)
        qualities = np.reshape(qualities, (n, length))

        fastq_data = ''
        for i in range(n):
            name = '@simple_fastq:i\n'
            seq = ''.join(sequences[i, :]) + '\n'
            name2 = '+\n'
            qual = ''.join(chr(s) for s in qualities[i, :]) + '\n'
            fastq_data += ''.join([name, seq, name2, qual])

        fastq_data = StringIO(fastq_data)
        fastq_data.seek(0)

        return fastq_data