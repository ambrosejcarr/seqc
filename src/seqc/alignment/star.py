import os
from subprocess import Popen, PIPE


def default_alignment_args(fastq_records, n_threads, index, output_dir):
    """default arguments for STAR alignment

    To report unaligned reads, add '--outSAMunmapped': 'Within',

    :param fastq_records:
    :param n_threads:
    :param index:
    :param output_dir:
    :return: dict, default alignment arguments
    """
    default_align_args = {
        '--runMode': 'alignReads',
        '--runThreadN': str(n_threads),
        '--genomeDir': index,
        '--outFilterType': 'BySJout',
        '--outFilterMultimapNmax': '1',  # require unique alignments
        '--limitOutSJcollapsed': '2000000',  # deal with many splice variants
        '--alignSJDBoverhangMin': '8',
        '--outFilterMismatchNoverLmax': '0.04',
        '--alignIntronMin': '20',
        '--alignIntronMax': '1000000',
        '--readFilesIn': fastq_records,
        '--outSAMprimaryFlag': 'AllBestScore',  # all equal-scoring reads are primary
        '--outFileNamePrefix': output_dir,
    }
    return default_align_args


def align(fastq_file, index, n_threads, alignment_dir, reverse_fastq_file=None,
          **kwargs):
    """
    :param fastq_file:
    :param index:
    :param n_threads:
    :param alignment_dir:
    :param reverse_fastq_file:
    :param kwargs:
    :return:
    """

    # check if file exists; if it does, return the filename
    # if os.path.isfile(alignment_dir + 'Aligned.out.sam'):
    #     if os.path.getsize(alignment_dir + 'Aligned.out.sam') > 0:
    #         return alignment_dir + 'Aligned.out.sam'

    runtime_args = default_alignment_args(
        fastq_file, n_threads, index, alignment_dir)

    for k, v in kwargs.items():  # overwrite or add any arguments passed from cmdline
        if not isinstance(k, str):
            try:
                k = str(k)
            except ValueError:
                raise ValueError('arguments passed to STAR must be strings')
        if not isinstance(v, str):
            try:
                v = str(v)
            except ValueError:
                raise ValueError('arguments passed to STAR must be strings')
        runtime_args[k] = v

    # construct command line arguments for STAR
    cmd = ['STAR']
    if reverse_fastq_file:
        for key, value in runtime_args.items():
            if key == '--readFilesIn':
                cmd.extend((key, value))
                cmd.append(reverse_fastq_file)
            else:
                cmd.extend((key, value))
    else:
        for pair in runtime_args.items():
            cmd.extend(pair)

    aln = Popen(cmd, stderr=PIPE, stdout=PIPE)
    out, err = aln.communicate()
    if err:
        raise ChildProcessError(err)

    return alignment_dir + 'Aligned.out.sam'
