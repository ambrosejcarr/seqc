__author__ = 'ambrose'

import seqc
import nose2
import os
import random
from io import StringIO
import pickle
import numpy as np
from copy import copy
import json


class GenerateFastq:

    # define some general constants
    _alphabet = ['A', 'C', 'G', 'T']

    def __init__(self):
        pass

    @classmethod
    def _forward_in_drop(cls, n, barcodes_):
        with open(barcodes_, 'rb') as f:
            barcodes_ = pickle.load(f)
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
            c1, c2 = seqc.three_bit.ThreeBitInDrop.split_cell(cb)
            c1, c2 = [seqc.three_bit.ThreeBit.bin2str(c) for c in [c1, c2]]
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
    def _reverse(n: int, read_length: int, fasta: str, gtf: str, tag_type='gene_id'):

        # read gtf
        reader = seqc.gtf.Reader(gtf)
        intervals = []
        for r in reader.iter_exons():
            end = int(r.end) - read_length
            start = int(r.start)
            if end > start:
                intervals.append((r.attribute[tag_type], start, end))

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
            tag, start, end = intervals[i]
            # get position within interval
            start = random.randint(start, end)
            end = start + read_length
            seq = fasta[start:end]
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

    @classmethod
    def in_drop(cls, n, prefix_, fasta, gtf, barcodes, tag_type='gene_id', replicates=3,
                *args, **kwargs):

        if not replicates >= 0:
            raise ValueError('Cannot generate negative replicates')

        fwd_len = 50
        rev_len = 100
        forward = cls._forward_in_drop(n, barcodes)
        forward = forward.read()  # consume the StringIO object
        reverse = cls._reverse(n, rev_len, fasta, gtf, tag_type=tag_type)
        reverse = reverse.read()  # consume the StringIO object
        with open(prefix_ + '_r1.fastq', 'w') as f:
            f.write(''.join([forward] * (replicates + 1)))
        with open(prefix_ + '_r2.fastq', 'w') as r:
            r.write(''.join([reverse] * (replicates + 1)))

    @classmethod
    def drop_seq(cls, n, prefix, fasta, gtf, tag_type='gene_id', replicates=3, *args,
                 **kwargs):
        rev_len = 100
        forward = cls._forward_drop_seq(n)
        forward = forward.read()  # consume the StringIO object
        reverse = cls._reverse(n, rev_len, fasta, gtf, tag_type=tag_type)
        reverse = reverse.read()  # consume the StringIO object
        with open(prefix + '_r1.fastq', 'w') as f:
            f.write(''.join([forward] * (replicates + 1)))
        with open(prefix + '_r2.fastq', 'w') as r:
            r.write(''.join([reverse] * (replicates + 1)))


def test_seqc_raw_fastq_input():
    for data_type in ['drop_seq', 'in_drop']:

        seqc.log.info('SEQC functest for %s' % data_type)

        package_dir = '/'.join(seqc.__file__.split('/')[:-3]) + '/'
        output_prefix = package_dir + 'test/%s/seqc_test' % data_type
        output_dir = package_dir + 'test/%s/' % data_type
        # output_dir = testing_dir + data_type + '/'
        index = package_dir + 'test/genome/'
        n_threads = 7

        # check if output_directory exists
        if not os.path.isdir(output_dir):
            seqc.log.info('SEQC test: creating testing directory')
            os.mkdir(output_dir)

        # check for necessary genome index files
        if not os.path.isdir(index):
            seqc.log.info('SEQC test: downloading mouse chr19 genome index files')
            bucket = 'dplab-data'
            key_prefix = 'genomes/mm38_chr19/'
            output_prefix = index
            seqc.io_lib.S3.download_files(bucket=bucket, key_prefix=key_prefix,
                                          output_prefix=output_prefix, cut_dirs=2)

        forward = output_dir + 'fastq/seqc_test_r1.fastq'
        reverse = output_dir + 'fastq/seqc_test_r2.fastq'
        barcodes = output_dir + 'barcodes/barcodes.p'

        ########################### correct special barcode cases ###########################

        if data_type == 'drop_seq':
            barcodes = None

        # get barcode data
        if barcodes:  # skip for types lacking barcode files
            # get directory
            *path, filename = barcodes.strip('\n').split('/')
            barcode_dir = '/'.join(path) + '/'

            # create directory if missing
            seqc.log.info('Creating %s barcode directory.' % data_type)
            if not os.path.isdir(barcode_dir):
                os.makedirs(barcode_dir)

            # download file if missing
            if not os.path.isfile(barcodes):
                seqc.log.info('Downloading %s barcode files' % data_type)
                bucket = 'dplab-data'
                key = 'barcodes/%s/barcodes.p' % data_type
                seqc.io_lib.S3.download_file(bucket=bucket, key=key, fout=barcodes)

        # define some genome files
        gtf = index + 'annotations.gtf'
        fasta = index + 'mm38_chr19.fa'

        # get fastq data
        fastq_files = [forward, reverse]

        # check if there are any fastq files for the given experiment; if not skip
        # this section
        if any(fastq_files):

            # check if files must be created, or if they already exist in the directory
            if not all([True if os.path.isfile(f) or f is None else False
                        for f in fastq_files]):

                # check that fastq directory exists; get the first filename that exists
                # and use it to check/create necessary path.
                for fastq in fastq_files:
                    if fastq:
                        break
                if not fastq:
                    raise ValueError('no fastq file detected')
                *path, _ = fastq.strip('\n').split('/')
                fastq_dir = '/'.join(path) + '/'
                if not os.path.isdir(fastq_dir):
                    os.makedirs(fastq_dir)

                # get generator for correct fastq data type and create fastq files
                generate_fastq = getattr(GenerateFastq, data_type)
                prefix = fastq_dir + 'seqc_test'
                generate_fastq(10000, prefix, fasta, gtf, barcodes=barcodes,
                               tag_type='gene_id')

        args = [
            data_type.replace('_', '-'),
            '-i', index,
            '-n', str(n_threads),
            '-o', output_prefix,
        ]
        if barcodes:
            args.extend(['-b', barcodes])
        if forward:
            args.extend(['-f', forward,])
        if reverse:
            args.extend(['-r', reverse])

        # yield a test
        yield seqc_raw_fastq_input, (args,)


def seqc_raw_fastq_input(args):
    parser = seqc.core.create_parser()
    kwargs = seqc.core.parse_args(parser, args)
    seqc.log.setup_logger('seqc_test.log')
    try:
        # log command line arguments for debugging
        arg_copy = copy(kwargs)
        del arg_copy['func']  # function is not serializable
        seqc.log.info('SEQC version: %s' % seqc.__version__)
        seqc.log.info('SEQC working directory: %s' % os.getcwd())
        seqc.log.info('Passed command line arguments: %s' %
                      json.dumps(arg_copy, separators=(',', ': '), indent=4))

        func = kwargs['func']
        func(**kwargs)

    except:
        seqc.log.exception()
        raise


if __name__ == "__main__":
    seqc.log.setup_logger()
    try:
        nose2.main(exit=False, module=__name__, verbosity=3)
    except:
        seqc.log.exception()
        raise
