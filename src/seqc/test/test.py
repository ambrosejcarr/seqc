import nose2
import unittest
from copy import deepcopy
import os
import seqc
import numpy as np
import re
import pickle
from more_itertools import first
from itertools import islice
from operator import attrgetter
import numpy as np
# from nose2.tools import params


# noinspection PyPep8Naming
class config:
    """
    Dummy class to hold state configuration variables that is importable into other
    scopes for testing
    """

    # this is the universal data dir for these tests
    seqc_dir = '/'.join(seqc.__file__.split('/')[:-3]) + '/'

    # set of dtypes
    data_types = ['in_drop', 'drop_seq']

    # universal file patterns
    samfile_pattern = seqc_dir + 'test_data/%s/sam/alignments.sam'
    forward_pattern = seqc_dir + 'test_data/%s/fastq/test_seqc_r1.fastq'
    reverse_pattern = seqc_dir + 'test_data/%s/fastq/test_seqc_r2.fastq'
    merged_pattern = seqc_dir + 'test_data/%s/fastq/merged.fastq'
    barcode_serial_pattern = seqc_dir + 'test_data/%s/barcodes/barcodes.p'
    barcode_partial_serial_pattern = seqc_dir + 'test_data/%s/barcodes/cb_partial.p'
    barcode_prefix_pattern = seqc_dir + 'test_data/%s/barcodes/'
    barcode_serialized_link_pattern = 's3://dplab-data/barcodes/%s/serial/barcodes.p'
    barcode_partial_serial_link_pattern = ('s3://dplab-data/barcodes/%s/serial/'
                                           'cb_partial.p')
    barcode_files_link_prefix_pattern = 's3://dplab-data/barcodes/%s/flat/'
    h5_name_pattern = seqc_dir + 'test_data/%s/h5/test_seqc.h5'
    gene_map_pattern = seqc_dir + 'test_data/%s/h5/test_seqc_gene_id_map.p'
    experiment_pattern = seqc_dir + 'test_data/%s/results/test_seqc.npz'

    # universal index files
    gtf = seqc_dir + 'test_data/genome/annotations.gtf'
    fasta = seqc_dir + 'test_data/genome/mm38_chr19.fa'
    index = seqc_dir + 'test_data/genome/'
    index_link = 's3://dplab-data/genome/mm38_chr19/'

    # config parameters
    n_threads = 7

    # s3 test file download links
    s3_forward_fastq_pattern = 's3://dplab-data/seqc/test_seqc/%s/fastq/forward/'
    s3_reverse_fastq_pattern = 's3://dplab-data/seqc/test_seqc/%s/fastq/reverse/'
    s3_merged_fastq_pattern = 's3://dplab-data/seqc/test_seqc/%s/fastq/merged/merged.fastq'
    s3_sam_pattern = 's3://dplab-data/seqc/test_seqc/%s/sam/alignments.sam'


def check_index():
    """ensure that there is an index present. If not, download it."""
    gfiles = ['Genome', 'SA', 'SAindex', 'annotations.gtf']
    if not os.path.isdir(config.index) or not all(os.path.isfile(config.index + f) for f
                                                  in gfiles):
        index_bucket = 'dplab-data'
        index_prefix = 'genomes/mm38_chr19/'
        seqc.io.S3.download_files(
            bucket=index_bucket, key_prefix=index_prefix, output_prefix=config.index,
            cut_dirs=2)


def check_barcodes(data_type: str):
    """ensure that all necessary barcodes are present and downloaded"""

    # make directory if not present
    barcode_dir = config.barcode_prefix_pattern % data_type
    if not os.path.isdir(barcode_dir):
        os.makedirs(barcode_dir)

    if data_type == 'drop_seq':
        return  # drop-seq has no barcodes

    # download flat files if not present
    if not any(f.endswith('.txt') and 'partial' not in f
               for f in os.listdir(barcode_dir)):
        bucket, key_prefix = seqc.io.S3.split_link(
            config.barcode_files_link_prefix_pattern % data_type)
        seqc.io.S3.download_files(bucket, key_prefix, output_prefix=barcode_dir,
                                  cut_dirs=3)

    # check for serialized barcode file; download if missing.
    serial_barcodes = config.barcode_serial_pattern % data_type
    if not os.path.isfile(serial_barcodes):
        bucket, key = seqc.io.S3.split_link(
            config.barcode_serialized_link_pattern % data_type)
        seqc.io.S3.download_file(bucket, key, fout=serial_barcodes)

    # get all downloaded flat text files; create partial versions for rapid testing
    n_barcodes = 10
    flat_files = [barcode_dir + f for f in os.listdir(barcode_dir) if f.endswith('.txt')
                  and 'partial' not in f]
    for i, fin in enumerate(flat_files):
        with open(fin) as barcode_file:
            data = islice(barcode_file, n_barcodes)
            with open(barcode_dir + 'cb%d_partial.txt' % (i + 1), 'w') as fout:
                for line in data:
                    fout.write(line)

    # create a serialized partial file
    partial_files = [barcode_dir + f for f in os.listdir(barcode_dir) if
                     f.endswith('_partial.txt')]
    cb = seqc.barcodes.CellBarcodes(*partial_files)
    with open(barcode_dir + 'cb_partial.p', 'wb') as f:
        pickle.dump(cb, f)
    assert (config.barcode_partial_serial_pattern % data_type == barcode_dir +
            'cb_partial.p')


def check_fastq(data_type: str) -> (str, str):
    """check if the required fastq files are present, else generate them
    """

    # replace any dashes with underscores
    data_type = data_type.replace('-', '_')

    check_barcodes(data_type)

    # get forward and reverse file ids
    forward = config.forward_pattern % data_type
    reverse = config.reverse_pattern % data_type
    if data_type == 'drop_seq':
        barcodes = seqc.barcodes.DropSeqCellBarcodes()
    else:
        barcodes = config.barcode_partial_serial_pattern % data_type

    prefix = config.seqc_dir + 'test_data/%s/fastq/test_seqc' % data_type
    # if not os.path.isdir(prefix):
    #     os.makedirs(prefix)

    if not all(os.path.isfile(f) for f in [forward, reverse]):
        gen_func = getattr(seqc.sequence.fastq.GenerateFastq, data_type)
        loc, unloc, unalign = 3000, 3000, 3000
        seqlen, fraglen, rep = 100, 1000, 3
        gen_func(loc, unloc, unalign, prefix, config.fasta, config.gtf, config.index,
                 replicates=rep, n_threads=config.n_threads, fragment_length=fraglen,
                 sequence_length=seqlen, barcodes=barcodes)

    return forward, reverse


def check_merged_fastq(data_type: str) -> str:
    """check if the required merged fastq files are present, else generate them"""
    forward, reverse = check_fastq(data_type)

    if not os.path.isfile(config.merged_pattern % data_type):
        forward, reverse = [forward], [reverse]
        exp_type = data_type.replace('_', '-')
        n_processors = 4
        output_dir = config.seqc_dir + 'test_data/%s/fastq' % data_type
        if not os.path.isdir(output_dir):
            os.makedirs(output_dir)
        if data_type == 'drop_seq':
            cb = seqc.barcodes.DropSeqCellBarcodes()
        else:
            with open(config.barcode_partial_serial_pattern % data_type, 'rb') as f:
                cb = pickle.load(f)

        # create merged file
        merged = seqc.sequence.fastq.merge_fastq(
            forward, reverse, exp_type, output_dir, cb, n_processors)
        msg = ('output_file "%s" does not match the expected pattern: "%s"' %
               (merged, config.merged_pattern % data_type))
        assert merged == config.merged_pattern % data_type, msg

    return config.merged_pattern % data_type


def check_sam(data_type: str) -> str:
    """check if the required sam files are present, else generate them"""

    # replace any dashes with underscores
    data_type = data_type.replace('-', '_')

    merged = check_merged_fastq(data_type)

    # generation params
    samfile = config.samfile_pattern % data_type
    sam_dir = '/' + '/'.join(samfile.strip('/').split('/')[:-1])
    if not os.path.isdir(sam_dir):
        os.makedirs(sam_dir)

    if not os.path.isfile(samfile):
        gen_func = seqc.sam.GenerateSam.from_merged_fastq
        gen_func(samfile, merged, config.index, n_threads=config.n_threads)

    assert os.path.isfile(samfile)

    return samfile


def check_h5(data_type: str) -> str:

    samfile = check_sam(data_type)

    expected_h5 = config.h5_name_pattern % data_type
    h5_dir = '/'.join(expected_h5.split('/')[:-1])
    if not os.path.isdir(h5_dir):
        os.makedirs(h5_dir)

    if not os.path.isfile(config.h5_name_pattern % data_type):
        h5 = seqc.arrays.ReadArray.from_samfile(
                samfile, config.h5_name_pattern % data_type, 4, int(1e7), config.gtf)
        assert os.path.isfile(config.h5_name_pattern % data_type)
        assert h5 == config.h5_name_pattern % data_type
        assert os.path.isfile(config.gene_map_pattern % data_type)
    else:
        h5 = config.h5_name_pattern % data_type

    return h5


def check_experiment(data_type: str) -> str:

    h5_file = check_h5(data_type)

    experiment_file = config.experiment_pattern % data_type
    exp_dir = '/'.join(experiment_file.split('/')[:-1])
    if not os.path.isdir(exp_dir):
        os.makedirs(exp_dir)

    if not os.path.isfile(experiment_file):
        ra = seqc.arrays.ReadArray.from_h5(h5_file)
        if data_type == 'drop_seq':
            exp = seqc.arrays.Experiment.from_read_array(ra, 0, 2)
        else:
            exp = seqc.arrays.Experiment.from_read_array(ra, 3, 2)
        exp.to_npz(experiment_file)

    return experiment_file


class GTFTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        check_index()
        cls.gtf = config.gtf
        cls._gene_id_pattern = b'(.*?gene_id ")(\D+)(\d+)(\.?)(\d*)(.*)'
        cls.fasta = config.fasta

    def test_reader_init(self):
        seqc.gtf.Reader(self.gtf)

    def test_reader_iter(self):
        rd = seqc.gtf.Reader(self.gtf)
        first(rd)  # check iterable

    def test_reader_iter_genes(self):
        rd = seqc.gtf.Reader(self.gtf)
        for i in rd.iter_gene_sets():
            break

        # make sure the gene_id is identical for each value
        gene_ids = set()
        for r in i:
            try:
                gene_ids.add(re.match(self._gene_id_pattern, r[-1]).groups()[2])
            except AttributeError:
                print(r)
                raise
        self.assertEqual(len(gene_ids), 1)

    def test_exon(self):
        rd = seqc.gtf.Reader(self.gtf)

        gene_records = first(rd.iter_gene_sets())

        # get an exon from this gene record
        for record in gene_records:
            if record[2] != b'exon':
                continue
            exon = seqc.gtf.Exon(record)
            break

        # test that all the properties are working
        self.assertIsInstance(exon.seqname, bytes)
        self.assertIsInstance(exon.chromosome, bytes)
        self.assertIsInstance(exon.source, bytes)
        self.assertIsInstance(exon.feature, bytes)
        self.assertIsInstance(exon.start, int)
        self.assertIsInstance(exon.end, int)
        self.assertIsInstance(exon.score, bytes)
        self.assertIn(exon.strand, [b'+', b'-'])
        self.assertIsInstance(exon.frame, bytes)

        # test that attribute parsing works
        try:
            self.assertIsInstance(exon.attribute(b'exon_id'), bytes)
        except:
            print(exon._fields[-1])
            raise

        # test repr & str
        self.assertIsInstance(repr(exon), str)
        self.assertIsInstance(str(exon), str)

        # test slicing
        plus_exon = deepcopy(exon)
        plus_exon._fields[6] = b'+'
        minus_exon = deepcopy(exon)
        minus_exon._fields[6] = b'-'
        size = exon.end - exon.start

        # one of these tests is yielding incorrect output!
        self.assertEqual(plus_exon[:round(size * .25)].end,
                         plus_exon.start + round(size * .25))
        self.assertEqual(plus_exon[-round(size * .25):].start,
                         plus_exon.end - round(size * .25))
        self.assertEqual(plus_exon[10:50].end, plus_exon.start + 50)
        self.assertEqual(plus_exon[10:50].start, plus_exon.start + 10)

        self.assertEqual(minus_exon[:round(size * .25)].end, minus_exon.end)
        self.assertEqual(minus_exon[:round(size * .25)].start,
                         minus_exon.end - round(size * .25))
        self.assertEqual(minus_exon[-round(size * .25):].start, minus_exon.start)
        self.assertEqual(minus_exon[-round(size * .25):].end, minus_exon.start +
                         round(size * .25))
        self.assertEqual(minus_exon[10:50].end, minus_exon.end - 10)
        self.assertEqual(minus_exon[10:50].start, minus_exon.end - 50)

    def test_transcript(self):
        rd = seqc.gtf.Reader(self.gtf)
        gene_records = first(rd.iter_gene_sets())

        iter_gene = iter(gene_records)

        # get the transcript record
        transcript_record = next(iter_gene)
        while transcript_record[2] != b'transcript':
            transcript_record = next(iter_gene)

        # now get all the UTRs
        exons = []
        record = next(iter_gene)  # get the first non-transcript record
        while record[2] not in [b'transcript', b'gene']:
            if record[2] in [b'exon', b'UTR']:
                exons.append(record)
            record = next(iter_gene)

        transcript = seqc.gtf.Transcript.from_records(transcript_record, *exons)

        # create a minus and plus transcript for testing
        if transcript._fields[6] == b'+':
            plus_transcript = transcript

            minus_transcript = deepcopy(transcript)
            minus_transcript._fields[6] = b'-'
            minus_transcript._exons = sorted(minus_transcript._exons,
                                             key=attrgetter('start'), reverse=True)
            for exon in minus_transcript:
                exon._fields[6] = b'-'
        else:
            minus_transcript = transcript
            plus_transcript = deepcopy(transcript)
            plus_transcript._fields[6] = b'+'
            plus_transcript._exons = sorted(plus_transcript._exons,
                                            key=attrgetter('start'))
            for exon in plus_transcript:
                exon._fields[6] = b'+'

        # test that all the properties are working
        self.assertIsInstance(transcript.seqname, bytes)
        self.assertIsInstance(transcript.chromosome, bytes)
        self.assertIsInstance(transcript.source, bytes)
        self.assertIsInstance(transcript.feature, bytes)
        self.assertIsInstance(transcript.start, int)
        self.assertIsInstance(transcript.end, int)
        self.assertIsInstance(transcript.score, bytes)
        self.assertIn(transcript.strand, [b'+', b'-'])
        self.assertIsInstance(transcript.frame, bytes)

        # test transcript slicing
        self.assertEqual(plus_transcript.size, minus_transcript.size)
        subset = round(plus_transcript.size * .25)

        self.assertEqual(plus_transcript[:subset].size, subset)
        self.assertEqual(plus_transcript[-subset:].size, subset)
        self.assertEqual(plus_transcript.start, plus_transcript[:subset].start)
        self.assertEqual(plus_transcript[-subset:].end, plus_transcript.end)

        self.assertEqual(minus_transcript[-subset:].size, subset)
        self.assertEqual(minus_transcript[:subset].size, subset)
        self.assertEqual(minus_transcript.end, minus_transcript[:subset].end)
        self.assertEqual(minus_transcript.start, minus_transcript[-subset:].start)

        self.assertEqual(plus_transcript[:subset].size, minus_transcript[:subset].size)
        self.assertEqual(plus_transcript[:-subset].size, minus_transcript[:-subset].size)
        self.assertEqual(plus_transcript[-subset:].size, minus_transcript[-subset:].size)
        self.assertEqual(plus_transcript[subset:].size, minus_transcript[subset:].size)
        self.assertEqual(plus_transcript[1:5].size, 4)
        self.assertEqual(minus_transcript[1:5].size, 4)
        self.assertEqual(plus_transcript[1:5].start, plus_transcript.start + 1)
        self.assertEqual(plus_transcript[1:5].end, plus_transcript.start + 5)
        self.assertEqual(minus_transcript[1:5].end, minus_transcript.end - 1)
        self.assertEqual(minus_transcript[1:5].start, minus_transcript.end - 5)

    def test_gene(self):
        rd = seqc.gtf.Reader(self.gtf)
        minus_gene = None
        plus_gene = None
        gene_iterator = rd.iter_gene_sets()
        while not all((minus_gene, plus_gene)):
            gene_records = list(next(gene_iterator))
            gene = seqc.gtf.Gene.from_records(*gene_records)
            if len(gene) > 1 and gene.strand == b'+':
                plus_gene = gene
            if len(gene) > 1 and gene.strand == b'-':
                minus_gene = gene

        plus_subset = int(np.mean([tx.size * .75 for tx in plus_gene]))
        minus_subset = int(np.mean([tx.size * .75 for tx in minus_gene]))

        for gene in (plus_gene, minus_gene):
            subset = int(np.mean([tx.size * .75 for tx in gene]))
            sliced_gene = gene[-subset:]
            self.assertTrue(all(tx.size == min(subset, original.size) for
                                tx, original in zip(sliced_gene, gene)))

        # test merge intervals
        res = list(plus_gene._merge_intervals([(1, 5), (4, 9)]))
        self.assertEqual(len(res), 1)
        self.assertEqual(res[0], (1, 9))

        res = list(plus_gene._merge_intervals([(1, 5), (4, 9), (-1, 3)]))
        self.assertEqual(len(res), 1)
        self.assertEqual(res[0], (-1, 9))

        res = list(plus_gene._merge_intervals([(4, 9), (-1, 3)]))
        self.assertEqual(len(res), 2)
        self.assertEqual(res, [(-1, 3), (4, 9)])

        # test intervals
        for gene in (plus_gene, minus_gene):
            subset = int(np.mean([tx.size * .75 for tx in gene]))
            res = gene.intervals(-1000)

    def test_gene_slicing(self):

        anno = seqc.gtf.Annotation(self.gtf, self.fasta)

        # test slicing
        sliced = [g[-1000:] for g in anno.genes]
        self.assertEqual(sum(1 for s in sliced if s.start and s.end), len(anno.genes))

        # test intervals
        ivs = [g.intervals(-1000, None) for g in anno.genes]
        self.assertTrue(all(len(iv) != 0 for iv in ivs))


    def test_annotation(self):

        anno = seqc.gtf.Annotation(self.gtf, self.fasta)

        # test properties
        self.assertIsInstance(len(anno), int)
        self.assertIsInstance(repr(anno), str)
        self.assertIsInstance(anno.chromosomes, dict)

        # test __getitem__
        sample_gene = anno.genes[0]
        chromosome = sample_gene.chromosome
        gene_id = sample_gene.string_gene_id
        self.assertIsInstance(anno[chromosome][gene_id], seqc.gtf.Gene)

        prefix = sample_gene.organism_prefix
        int_id = sample_gene.integer_gene_id
        string_id = sample_gene.string_gene_id.split(b'.')[0]
        self.assertEqual(string_id, seqc.gtf.Record.int2str_gene_id(int_id, prefix))

    def test_random_sequences(self):

        anno = seqc.gtf.Annotation(self.gtf, self.fasta)

        n_seqs = 1000
        seqlen = 98

        # check that we can generate sequences inside expected region
        seqs, genes = anno.random_sequences(n_seqs, seqlen, True, -1000, None)
        self.assertEqual(len(seqs), n_seqs)
        self.assertTrue(all(len(s) == seqlen for s in seqs))
        self.assertEqual(len(genes), n_seqs)

        # check that we can generate sequences outside of expected region
        seqs, genes = anno.random_sequences(n_seqs, seqlen, True, None, -1000)
        self.assertEqual(len(seqs), n_seqs)
        self.assertTrue(all(len(s) == seqlen for s in seqs))
        self.assertEqual(len(genes), n_seqs)

        # finally, want to be able to generate giberish sequences; this is in align


class TestResolveAlignments(unittest.TestCase):

    # unique genes:
    # -------------
    # 0-7 are all 100% unique
    #
    # ambiguous genes:
    # ----------------
    # gene_ids 8, 12, 22, 26, and 49 are good values. 8 has the form:
    # {HashableArray([   8, 1100]): 0.035040466553386274,
    #  HashableArray([  8, 808, 898]): 0.23603050954550936,
    #  HashableArray([8]): 0.69332330791897157,
    #  HashableArray([   8,  808,  898, 1404]): 0.035605715982132781}

    @classmethod
    def setUpClass(cls):
        expectations = config.seqc_dir + 'data_test/p_coalignment_array.p'
        with open(expectations, 'rb') as f:
            cls.expectations = pickle.load(f)

    def construct_dummy_array(self):

        # unique
        features, probabilities = zip(*self.expectations[5].items())
        indices = np.random.multinomial(np.array(list(self.expectations[5])), 50)
        reads1 = [features[i] for i in indices]  # list of array features

        features, probabilities = zip(*self.expectations[5].items())
        indices = np.random.multinomial(np.array(list(self.expectations[5])), 1)
        reads2 = [features[i] for i in indices]  # list of array features

        # ambiguous
        features, probabilities = zip(*self.expectations[8].items())
        indices = np.random.multinomial(np.array(list(self.expectations[8])), 100)
        reads3 = [features[i] for i in indices]  # list of array features

        # 2nd ambiguous
        features, probabilities = zip(*self.expectations[49].items())
        indices = np.random.multinomial(np.array(list(self.expectations[49])), 20)
        reads4 = [features[i] for i in indices]  # list of array features

        # make the rest of the ReadArray fields.

        # same cell
        # different RMT (later, allow errors/conflation of molecules)
        # passes all other filters

    def test_case0(self):
        raise NotImplementedError


if __name__ == "__main__":
    nose2.main()