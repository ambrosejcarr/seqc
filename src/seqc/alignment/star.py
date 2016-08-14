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
    align a fastq file, or a paired set of fastq files

    :param fastq_file: str, location of a fastq file
    :param index: str, folder containing the STAR index
    :param n_threads: int, number of parallel alignment processes to spawn
    :param alignment_dir: directory for output data
    :param reverse_fastq_file: optional, location of reverse paired-end fastq file
    :param kwargs: additional kwargs for STAR, passed without the leading '--'
    :return: str, .sam file location
    """

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
        runtime_args['--' + k] = v

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
