__author__ = "Ambrose J. Carr"

from collections import defaultdict
from intervaltree import IntervalTree
from more_itertools import first
import seqc
import pickle
import re
import os


class Exon:
    """
    Dummy class to store exon data for ConvertFeatureCoordinates. Replaceable with
    NamedTuple
    """

    __slots__ = ['start', 'end']

    def __init__(self, start, end):
        self.start = start
        self.end = end


class Transcript:
    """
    Dummy class to store transcript data for ConvertFeatureCoordinates. Not replaceable
    with NamedTuple because the .exons property must be mutable.
    """

    __slots__ = ['start', 'end', 'strand', 'feature', 'exons', 'scid']

    def __init__(self, start, end, strand, feature, exons, scid):
        self.start = start
        self.end = end
        self.strand = strand
        self.feature = feature
        self.exons = exons
        self.scid = scid


class ConvertFeatureCoordinates:
    """
    creates a lookup object whose .translate() method rapidly converts chromosome
    coordinates to the associated transcripts if the alignment falls in the final n bases
    of the transcript.

    feature_table = dict() (chr, strand position) -> {scids}
    feature_positions = dict() scid -> (start, end)

    """

    def __init__(self, feature_table, feature_positions):
        # todo type checks for these; this has bugged out multiple times...
        self._table = feature_table
        self._positions = feature_positions

    @classmethod
    def from_gtf(cls, gtf, fragment_length):
        """
        Create a ConvertFeatureCoordinates object from a gtf file. Include only
        features that fall in the final <fragment_length> of each transcript.
        """
        feature_table = defaultdict(set)
        feature_positions = defaultdict(list)
        pattern = re.compile(r'(^.*?scseq_id "SC)(.*?)(";.*?$)')
        # pile up exons of each scid
        transcripts = {}
        with open(gtf, 'r') as f:
            for record in f:
                # get rid of comments
                if record.startswith('#'):
                    continue
                record = record.split('\t')
                if record[2] == 'transcript':
                    scid = int(re.match(pattern, record[-1]).group(2))
                    transcripts[scid] = Transcript(
                        start=int(record[3]), end=int(record[4]), strand=record[6],
                        feature=record[0], exons=[], scid=scid
                    )
                elif record[2] == 'exon':
                    scid = int(re.match(pattern, record[-1]).group(2))
                    exon = Exon(start=int(record[3]), end=int(record[4]))
                    transcripts[scid].exons.append(exon)

        # process the exons, deleting duplicates where they are found
        sort_key = lambda x: (x.start, x.end)
        for scid, transcript in transcripts.items():
            # eliminate duplicate exons, then regenerate sorted order among exons
            if transcript.strand == '+':
                exons = sorted((s for s in set(transcript.exons)), key=sort_key)
                transcript.exons = exons
            else:
                exons = sorted((s for s in set(transcript.exons)), key=sort_key,
                               reverse=True)
                transcript.exons = exons
            # process_transcript modifies feature_table and feature_positions in place
            # each time the function is called
            cls._process_transcript(transcript, feature_table, feature_positions,
                                    fragment_length)

        # get rid of defaultdict
        feature_table = dict(feature_table)
        feature_positions = dict(feature_positions)
        return cls(feature_table, feature_positions)

    def save(self, fout):
        """save a serialized version of self as a pickle file"""
        dout = {'feature_table': self._table,
                'feature_positions': self._positions,
                'fragment_length': self._fragment_length}
        with open(fout, 'wb') as f:
            pickle.dump(dout, f)

    @classmethod
    def load(cls, fin):
        """load a ConvertFeatureCoordinates file from a serialized pickle"""
        with open(fin, 'rb') as f:
            din = pickle.load(f)
            return cls(**din)

    @staticmethod
    def _round_up(x):
        """round x up to the nearest hundred"""
        return x if x % 100 == 0 else x + 100 - x % 100

    @staticmethod
    def _round_down(x):
        """round x down to the nearest hundred"""
        return x // 100 * 100

    @classmethod
    def _process_transcript(cls, transcript, feature_table, feature_positions,
                            fragment_len):
        """updates feature_table and feature_positions for the passed transcript"""
        remaining_fragment_size = fragment_len
        if transcript.strand == '+':
            for exon in transcript.exons[::-1]:  # iterate backwards through exons
                feature_end = cls._round_up(exon.end)
                if exon.start < exon.end - remaining_fragment_size:

                    # break exon into 100 base pieces, round down start, up end
                    feature_start = cls._round_down(exon.start)

                    # add each 100-base position in interval to feature_table
                    for position in range(feature_start, feature_end + 100, 100):
                        key = (transcript.strand, transcript.feature, position)
                        feature_table[key].add(transcript.scid)

                    # add feature_positions to second look-up table
                    feature_positions[transcript.scid].append((
                        exon.start, exon.end))

                    remaining_fragment_size -= exon.end - exon.start
                else:
                    # get unrounded start and add to feature_positions
                    feature_start = exon.end - remaining_fragment_size
                    feature_positions[transcript.scid].append((
                        feature_start, exon.end))
                    # round and add positions to feature_table
                    feature_start = cls._round_down(feature_start)
                    for position in range(feature_start, feature_end + 100, 100):
                        key = (transcript.strand, transcript.feature, position)
                        feature_table[key].add(transcript.scid)
                    break  # move on to next set of exons
        else:
            for exon in transcript.exons[::-1]:  # iterate backwards
                feature_start = cls._round_down(exon.start)
                if exon.start + remaining_fragment_size < exon.end:

                    # break exon into 100 base pieces, round down start, up end
                    feature_end = cls._round_up(exon.end)

                    # add each 100-base position in interval to feature_table
                    for position in range(feature_start, feature_end + 100, 100):
                        key = (transcript.strand, transcript.feature, position)
                        feature_table[key].add(transcript.scid)

                    # add feature_positions to second look-up table
                    feature_positions[transcript.scid].append(
                        (exon.start, exon.end))

                    remaining_fragment_size -= exon.end - exon.start
                else:
                    # get unrounded start and add to feature_positions
                    feature_end = exon.start + remaining_fragment_size
                    feature_positions[transcript.scid].append((
                        exon.start, feature_end))
                    # round and add positions to feature_table
                    feature_end = cls._round_up(exon.start + remaining_fragment_size)
                    for position in range(feature_start, feature_end + 100, 100):
                        key = (transcript.strand, transcript.feature, position)
                        feature_table[key].add(transcript.scid)
                    break  # move on to next set of exons

    def translate(self, strand: str, chromosome: str, position: int) -> list:
        """
        translate a strand, chromosome, and position into all associated SCIDs, which
        correspond to groups of overlapping transcripts.

        If no feature, returns None to reduce memory footprint.
        """
        rounded_position = position // 100 * 100
        try:
            potential_ids = self._table[(strand, chromosome, rounded_position)]
        except KeyError:
            return []

        # purge end cases
        for scid in potential_ids:
            if any(s < position < e for (s, e) in self._positions[scid]):
                return [scid]  # todo write test to ensure there is never more than one feat
            else:
                return []

    def add_mtDNA(self, gtf):
        """
        need to associate all positions of the MT, chromosome 'MT', and either
        strand with an SCID designated as mitochondrial. Probably need a fasta file for
        this, since it is a DNA object, not cDNA.
        """
        raise NotImplementedError


class GeneTable:

    def __init__(self, gtf):
        self.interval_table = defaultdict(dict)
        self._all_genes = set()
        pattern = re.compile(r'(^.*?gene_name ")(.*?)(";.*?$)')
        with open(gtf, 'r') as f:
            for record in f:
                # get rid of comments
                if record.startswith('#'):
                    continue
                record = record.split('\t')
                if record[2] == 'gene':
                    chromosome = record[0]
                    start = int(record[3])
                    end = int(record[4])
                    strand = record[6]
                    gene_id = re.match(pattern, record[-1]).group(2)
                    self._all_genes.add(gene_id)
                    try:
                        self.interval_table[chromosome][strand].addi(
                            start, end, gene_id)
                    except KeyError:  # interval does not exist for this strand
                        self.interval_table[chromosome][strand] = IntervalTree()
                        self.interval_table[chromosome][strand].addi(
                            start, end, gene_id)

        self.interval_table = dict(self.interval_table)  # get rid of defaultdict

    def coordinates_to_gene_ids(self, chromosome, start, end, strand):
        intervals = self.interval_table[chromosome][strand].search(start, end)
        return [i.data for i in intervals]

    def all_genes(self):
        return self._all_genes


def construct_gene_table(gtf):
    """
    construct a feature table for n(3) conversions of genomic coordinates to genes
    features
    """

    # define two helper functions to round to the nearest hundred

    interval_table = defaultdict(dict)
    pattern = re.compile(r'(^.*?gene_name ")(.*?)(";.*?$)')
    with open(gtf, 'r') as f:
        for record in f:
            # get rid of comments
            if record.startswith('#'):
                continue
            record = record.split('\t')
            if record[2] == 'gene':
                chromosome = record[0]
                start = int(record[3])
                end = int(record[4])
                strand = record[6]
                gene_id = int(re.match(pattern, record[-1]).group(2))
                try:
                    interval_table[chromosome][strand].addi(start, end, gene_id)
                except KeyError:  # interval does not exist for this strand
                    interval_table[chromosome][strand] = IntervalTree()
                    interval_table[chromosome][strand].addi(start, end, gene_id)
    return interval_table


class ConvertGeneCoordinates:
    """Converts alignments in chromosome coordinates to gene coordinates"""

    def __init__(self, dict_of_interval_trees, id_map):
        """
        see from_gtf() method for in-depth documentation
        """

        seqc.util.check_type(dict_of_interval_trees, dict, 'dict_of_interval_trees')
        seqc.util.check_type(id_map, dict, 'id_map')

        # check that input dictionary isn't empty
        if not dict_of_interval_trees:
            raise ValueError('Cannot create an empty ConvertGeneCoordinates object. '
                             'Please pass a non-empty dict_of_interval_trees input.')

        # check type of individual trees
        err_msg = 'dict_of_interval_trees must contain only IntervalTree leaves, not %s'
        for tree in dict_of_interval_trees.values():
            if not isinstance(tree, IntervalTree):
                    raise TypeError('all dictionary values must be IntervalTrees not %s'
                                    % type(tree))

        # set self.data; self.id_map
        self._data = dict_of_interval_trees
        self._id_map = id_map

    def translate(self, strand: str, chromosome: str, position: int) -> list:
        """
        translate an alignment in chromosome coordinates to gene coordinates

        Note that there are some cases where genomic coordinates do not map to single
        genes due to double utilization of exons. In this case, the method returns None.

        args:
        -----
        strand (+, -): the strand of the record to be translated
        chromosome: the chromosome of the record to be translated
        position: the position of the record to be translated

        returns:
        --------
        records: 1-item list of gene overlapping the given position. If not unique,
         returns [].

        """
        try:
            ivs = self._data[(chromosome, strand)].search(position)
        except KeyError:  # should only happen with malformed input

            # check for bad input
            if not isinstance(chromosome, str):
                raise TypeError('chromosome must be <class "str">, not %s' %
                                type(chromosome))
            elif not strand in ('-', '+'):
                raise ValueError('strand must be one of "-" or "+"')
            else:  # not sure what the problem is here, raise original exception
                raise

        if len(ivs) >= 1:
            return [first(ivs).data]
        else:
            return []

    def int2str_id(self, identifier: int) -> str:
        """accessory function to convert integer ids back into string gene names"""
        return self._id_map[identifier]

    @staticmethod
    def str2int_id(identifier: str) -> int:
        """accessory function to convert string ids into integers"""
        return hash(identifier)

    def pickle(self, fname: str) -> None:
        """Serialize self and save it to disk as fname"""

        # check that fname is a file:
        if fname.endswith('/'):
            raise ValueError('fname must be the name of a file, not a directory')

        # check that directory exists
        *dir, name = fname.split('/')
        if not os.path.isdir('/'.join(dir)):
            raise FileNotFoundError('the directory fname should be saved in does not '
                                    'exist')
        with open(fname, 'wb') as f:
            pickle.dump({'dict_of_interval_trees': self._data, 'id_map': self._id_map}, f)

    @classmethod
    def from_pickle(cls, fname: str) -> None:
        """load a ConvertGeneCoordinates object from a serialized file"""
        with open(fname, 'rb') as f:
            data = pickle.load(f)
        return cls(**data)

    @classmethod
    def from_gtf(cls, gtf: str, fragment_length: int=1000):
        """
        construct a ConvertGeneCoordinates object from a gtf file. Also creates a map
        of integer ids to genes, which can be saved with using pickle

        # todo improve this
        The map of strings to ints can be done by assigning sequential integers to values
        as they are detected. This means that runs using different gtf files will
        obtain different integer values. Hashing is another option, but the methods I've
        looked up cannot generate integers compatible with uint32, which is preferred
        downstream. A superior method would deterministically map gene ids to
        uint32s.

        Current methods use python hash, require the gene field be int64, and that
        a map be saved.

        args;
        -----
        gtf: str file identifier corresponding to the gtf file to construct the
         ConvertGeneCoordinates object from.

        returns:
        --------
        ConvertGeneCoordinates object built from gtf

        """
        data = {}
        id_map = {}
        gtf_reader = seqc.gtf.Reader(gtf)

        for record in gtf_reader.iter_genes_final_nbases(fragment_length):

            # check if gene is in map
            gene = record.attribute['gene_id']
            int_id = hash(gene)
            if not int_id in id_map:
                id_map[int_id] = gene

            for iv in record.intervals:
                try:
                    data[(record.seqname, record.strand)].addi(
                        iv[0], iv[1], int_id)
                except KeyError:
                    data[record.seqname, record.strand] = IntervalTree()
                    data[record.seqname, record.strand].addi(
                        iv[0], iv[1], int_id)
        return cls(data, id_map)
