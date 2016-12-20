from subprocess import Popen, PIPE
from multiprocessing import cpu_count
from os import makedirs
import shlex


def default_alignment_args(
        fastq_records: str, n_threads: int or str, index: str, output_dir: str) -> dict:
    """default arguments for STAR alignment

    To report unaligned reads, add '--outSAMunmapped': 'Within',

    :param fastq_records: str, name of fastq file
    :param n_threads: int or str, number of threads to allocate when calling STAR
    :param index: str, location of the STAR index
    :param output_dir: str, prefix for output files
    :return: dict, default alignment arguments
    """
    default_align_args = {
        '--runMode': 'alignReads',
        '--runThreadN': str(n_threads),
        '--genomeDir': index,
        '--outFilterType': 'BySJout',
        '--outFilterMultimapNmax': '100',  # require unique alignments
        '--limitOutSJcollapsed': '2000000',  # deal with many splice variants
        '--alignSJDBoverhangMin': '8',
        '--outFilterMismatchNoverLmax': '0.04',
        '--alignIntronMin': '20',
        '--alignIntronMax': '1000000',
        '--readFilesIn': fastq_records,
        '--outSAMprimaryFlag': 'AllBestScore',  # all equal-scoring reads are primary
        '--outSAMtype': 'BAM Unsorted',
        '--outFileNamePrefix': output_dir,
    }
    if fastq_records.endswith('.gz'):
        default_align_args['--readFilesCommand'] = 'gunzip -c'
    if fastq_records.endswith('.bz2'):
        default_align_args['--readFilesCommand'] = 'bunzip2 -c'
    return default_align_args


def align(fastq_file: str, index: str, n_threads: int, alignment_dir: str,
          reverse_fastq_file: str or bool=None, **kwargs) -> str:
    """align a fastq file, or a paired set of fastq files

    :param fastq_file: str, location of a fastq file
    :param index: str, folder containing the STAR index
    :param n_threads: int, number of parallel alignment processes to spawn
    :param alignment_dir: str, directory for output data
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

    cmd = shlex.split(' '.join(cmd))
    aln = Popen(cmd, stderr=PIPE, stdout=PIPE)
    out, err = aln.communicate()
    if err:
        raise ChildProcessError(err)

    return alignment_dir + 'Aligned.out.bam'


def create_index(
        fasta: str,
        gtf: str,
        genome_dir: str,
        read_length: int=75, **kwargs) -> None:
    """Create a new STAR index

    :param fasta: complete filepath to fasta file
    :param gtf: complete filepath to gtf file
    :param genome_dir: directory in which new index should be constructed
    :param read_length: length of reads that will be aligned against this index
    :param kwargs: additional keyword arguments to pass to the genome construction call.
      to pass --sjdbFileChrStartEnd filename, pass sjdbFileChrStartEnd=filename (no --)
    :return: None
    """
    ncpu = str(cpu_count())
    makedirs(genome_dir, exist_ok=True)
    overhang = str(read_length - 1)

    cmd = (
        'STAR '
        '--runMode genomeGenerate '
        '--runThreadN {ncpu} '
        '--genomeDir {genome_dir} '
        '--genomeFastaFiles {fasta} '
        '--sjdbGTFfile {gtf} '
        '--sjdbOverhang {overhang} '.format(
            ncpu=ncpu, genome_dir=genome_dir, fasta=fasta, gtf=gtf, overhang=overhang)
    )

    for k, v in kwargs.items():
        cmd += '--{k} {v} '.format(k=k, v=v)

    p = Popen(cmd, stderr=PIPE, stdout=PIPE)
    out, err = p.communicate()
    if err:
        raise ChildProcessError(err)
