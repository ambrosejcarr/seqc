__author__ = "Ambrose J. Carr"

from collections import defaultdict
from intervaltree import IntervalTree
import re
from sys import argv, exit


def _round_up(x):
    return x if x % 100 == 0 else x + 100 - x % 100


def _round_down(x):
    return x // 100 * 100


def _process_transcript(transcript, fragment_len, feature_table, feature_positions):
    """updates feature_table and feature_positions (defined below)"""
    remaining_fragment_size = fragment_len
    if transcript.strand == '+':
        for exon in transcript.exons[::-1]:  # iterate backwards through exons
            feature_end = _round_up(exon.end)
            if exon.start < exon.end - remaining_fragment_size:

                # break exon into 100 base pieces, round down start, up end
                feature_start = _round_down(exon.start)

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
                feature_start = _round_down(feature_start)
                for position in range(feature_start, feature_end + 100, 100):
                    key = (transcript.strand, transcript.feature, position)
                    feature_table[key].add(transcript.scid)
                break  # move on to next set of exons
    else:
        for exon in transcript.exons[::-1]:  # iterate backwards
            feature_start = _round_down(exon.start)
            if exon.start + remaining_fragment_size < exon.end:

                # break exon into 100 base pieces, round down start, up end
                feature_end = _round_up(exon.end)

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
                feature_end = _round_up(exon.start + remaining_fragment_size)
                for position in range(feature_start, feature_end + 100, 100):
                    key = (transcript.strand, transcript.feature, position)
                    feature_table[key].add(transcript.scid)
                break  # move on to next set of exons


# define some simple container classes for Exons and Transcripts
class Exon:
    def __init__(self, start, end):
        self.start = start
        self.end = end


class Transcript:
    def __init__(self, start, end, strand, feature, exons, scid):
        self.start = start
        self.end = end
        self.strand = strand
        self.feature = feature
        self.exons = exons
        self.scid = scid


def construct_feature_table(gtf, fragment_len=1000):
    """
    construct a feature table for n(3) conversions of genomic -> transcriptomic
    coordinates
    """

    # define two helper functions to round to the nearest hundred

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
                exon = Exon(int(record[3]), int(record[4]))
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
        _process_transcript(transcript, fragment_len, feature_table, feature_positions)

    return dict(feature_table), dict(feature_positions)  # get rid of defaultdict


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
                        self.interval_table[chromosome][strand].addi(start, end, gene_id)
                    except KeyError:  # interval does not exist for this strand
                        self.interval_table[chromosome][strand] = IntervalTree()
                        self.interval_table[chromosome][strand].addi(start, end, gene_id)

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



if __name__ == "__main__":
    if not len(argv) == 3:
        print('usage: python convert_features.py gtf_file (string) insert_size (int)')
        exit(1)
    else:
        construct_feature_table(argv[1], argv[2])