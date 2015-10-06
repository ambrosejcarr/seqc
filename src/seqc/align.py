__author__ = 'ambrose'

from subprocess import Popen, PIPE, call
import os
from shutil import rmtree, copyfileobj
from collections import defaultdict
# todo eliminate these dependencies
from seqc.sa_preprocess import prepare_fasta
from seqc.sa_postprocess import ordering_out_stacks, create_gtf_reduced
import seqc
import ftplib
import gzip
import bz2


def star(index, n_threads, temp_dir):
    cmd = [
        'STAR', '--runMode', 'alignReads',
        '--runThreadN', str(n_threads),
        '--genomeDir', index,
        '--outFilterType', 'BySJout',
        '--outFilterMultimapNmax', '50',  # require <= 50 matches
        '--alignSJDBoverhangMin', '8',
        '--outFilterMismatchNoverLmax', '0.04',
        '--alignIntronMin', '20',
        '--alignIntronMax', '1000000',
        '--readFilesIn', '%smerged_temp.fastq' % temp_dir,
        '--outSAMunmapped', 'Within',
        '--outSAMprimaryFlag', 'AllBestScore',  # all equal-scoring reads are primary
        '--outFileNamePrefix', temp_dir,
        # '--outSAMreadID', 'Number'  # test if this is output in order
    ]
    aln = Popen(cmd, stderr=PIPE, stdout=PIPE)
    out, err = aln.communicate()
    if err:
        raise ChildProcessError(err)

    # get the number of aligned reads
    with open(temp_dir + 'Log.final.out', 'r') as f:
        data = f.readlines()
    # * 2 heuristic for multiple alignments, in practice we get about x2 reads
    n_alignments = int(data[5].strip().split('\t')[-1]) * 2
    return n_alignments


# todo | make sure organism is provided as a "choice" in the argparser
# define download locations for mouse and human genome files and annotations
_download_links = dict(hg38={
    'genome': ('ftp://ftp.ensembl.org/pub/release-80/fasta/homo_sapiens/dna/'
               'Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz'),
    'gtf': ('ftp://ftp.ensembl.org/pub/release-80/gtf/homo_sapiens/'
            'Homo_sapiens.GRCh38.80.gtf.gz'),
    'cdna': [('ftp://ftp.ensembl.org/pub/release-80/fasta/homo_sapiens/'
              'ncrna/Homo_sapiens.GRCh38.ncrna.fa.gz'),
             ('ftp://ftp.ensembl.org/pub/release-80/fasta/homo_sapiens/cdna/'
              'Homo_sapiens.GRCh38.cdna.all.fa.gz')]
}, mm38={
    'genome': ('ftp://ftp.ensembl.org/pub/release-76/fasta/mus_musculus/dna/'
               'Mus_musculus.GRCm38.dna.primary_assembly.fa.gz'),
    'gtf': ('ftp://ftp.ensembl.org/pub/release-76/gtf/mus_musculus/'
            'Mus_musculus.GRCm38.76.gtf.gz'),
    'cdna': [('ftp://ftp.ensembl.org/pub/release-76/fasta/mus_musculus/ncrna/'
              'Mus_musculus.GRCm38.ncrna.fa.gz'),
             ('ftp://ftp.ensembl.org/pub/release-76/fasta/mus_musculus/cdna/'
              'Mus_musculus.GRCm38.cdna.all.fa.gz')]
}, ci2={
    'genome': ('ftp://ftp.ensembl.org/pub/release-81/fasta/ciona_intestinalis/dna/'
               'Ciona_intestinalis.KH.dna.toplevel.fa.gz'),
    'gtf': ('ftp://ftp.ensembl.org/pub/release-81/gtf/ciona_intestinalis/'
            'Ciona_intestinalis.KH.81.gtf.gz'),
    'cdna': [('ftp://ftp.ensembl.org/pub/release-81/fasta/ciona_intestinalis/cdna/'
              'Ciona_intestinalis.KH.cdna.all.fa.gz')]
}, cs2={
    'genome': ('ftp://ftp.ensembl.org/pub/release-81/fasta/ciona_savignyi/dna/'
               'Ciona_savignyi.CSAV2.0.dna.toplevel.fa.gz'),
    'gtf': ('ftp://ftp.ensembl.org/pub/release-81/gtf/ciona_savignyi/'
            'Ciona_savignyi.CSAV2.0.81.gtf.gz'),
    'cdna': [('ftp://ftp.ensembl.org/pub/release-81/fasta/ciona_savignyi/cdna/'
              'Ciona_savignyi.CSAV2.0.cdna.all.fa.gz')]
}, phix={  # use genome for both DNA and cDNA
    'genome': ('ftp://ftp.ncbi.nlm.nih.gov/genomes/Viruses/'
               'Enterobacteria_phage_phiX174_sensu_lato_uid14015/NC_001422.fna'),
    'cdna': ('ftp://ftp.ncbi.nlm.nih.gov/genomes/Viruses/'
             'Enterobacteria_phage_phiX174_sensu_lato_uid14015/NC_001422.fna'),
    # note GFF not GTF, this is not implemented yet
    'gff': ('ftp://ftp.ncbi.nlm.nih.gov/genomes/Viruses/'
            'Enterobacteria_phage_phiX174_sensu_lato_uid14015/NC_001422.gff')
    # this is a potential gtf record:
    # create a gtf using stringIO? NC_001422.1	RefSeq	region	1	5386	.	+	.	ID=id0;Dbxref=taxon:374840;Is_circular=true;gbkey=Src;genome=genomic;mol_type=genomic DNA;nat-host=Escherichia coli
})


# define expected file names: efn[organism][filetype] = [names]  (list)
def _define_file_names(download_links):
    file_names = defaultdict(dict)
    for organism, file_types in download_links.items():
        for file_type, link in file_types.items():
            if isinstance(link, list):  # deal with lists of cDNA files
                file_name = [l.split('/')[-1] for l in link]
            else:
                file_name = link.split('/')[-1]
            file_names[organism][file_type] = file_name
    return file_names

_file_names = _define_file_names(_download_links)


def _check_type(obj, type_):
    if not isinstance(obj, type_):
        raise TypeError('%s must be of type %s' % (obj, type_))


class STAR:

    def __init__(self, temp_dir, n_threads, index, *organism):
        """
        classdocs
        """

        # type-check arguments
        _check_type(temp_dir, str)
        _check_type(n_threads, int)
        _check_type(index, str)

        # create temp_dir if it does not exist
        if not os.path.isdir(temp_dir):
            os.makedirs(temp_dir)
        if not temp_dir.endswith('/'):
            temp_dir += '/'
        self.temp_dir = temp_dir

        # create index if it does not exist
        if not os.path.isdir(index):
            os.makedirs(index)
        if not index.endswith('/'):
            index += '/'
        self.index = index

        # check that organism is a valid argument
        valid_organisms = ['mm38', 'hg38', 'mmhg38', 'cs2', 'ci2']
        if organism:
            if not all(org in valid_organisms for org in organism):
                raise ValueError('Invalid organism value. Supported organisms: %r' %
                                 valid_organisms)
        self.organism = organism

        self.n_threads = n_threads

    @classmethod
    def default_alignment_args(cls, fastq_records, n_threads, index, temp_dir):
        default_align_args = {
            '--runMode': 'alignReads',
            '--runThreadN': str(n_threads),
            '--genomeDir': index,
            '--outFilterType': 'BySJout',
            '--outFilterMultimapNmax': '10',  # require <= 10 matches
            '--limitOutSJcollapsed': '2000000',  # deal with many splice variants
            '--alignSJDBoverhangMin': '8',
            '--outFilterMismatchNoverLmax': '0.04',
            '--alignIntronMin': '20',
            '--alignIntronMax': '1000000',
            '--readFilesIn': fastq_records,
            '--outSAMunmapped': 'Within',
            '--outSAMprimaryFlag': 'AllBestScore',  # all equal-scoring reads are primary
            '--outFileNamePrefix': temp_dir,
        }
        return default_align_args

    @staticmethod
    def _download_ftp_file(link, prefix, clobber=False):
        """downloads ftp_file available at 'link' into the 'prefix' directory"""
        # check link validity
        if not link.startswith('ftp://'):
            raise ValueError(
                'link must start with "ftp://". Provided link is not valid: %s' % link)

        # create prefix directory if it does not exist
        if not os.path.isdir(prefix):
            os.makedirs(prefix)

        # make sure prefix has a trailing '/':
        if not prefix.endswith('/'):
            prefix += '/'

        ip, *path, file_name = link.split('/')[2:]  # [2:] -- eliminate leading 'ftp://'
        path = '/'.join(path)

        # check if file already exists
        if os.path.isfile(prefix + file_name):
            if not clobber:
                return  # file is already present, nothing to do.

        ftp = ftplib.FTP(ip)
        try:
            ftp.login()
            ftp.cwd(path)
            with open(prefix + file_name, 'wb') as fout:
                ftp.retrbinary('RETR %s' % file_name, fout.write)
        finally:
            ftp.close()

    @staticmethod
    def _gunzip_file(gzipped_file):
        # todo replace with gzip module
        # python is too damn slow.
        # with gzip.open(gzipped_file, 'rb') as f:
        #     fout = gzipped_file.replace('.gz', '')
        #     data = f.read()
        #     with open(fout, 'wb') as fo:
        #         fo.write(data)
        unzipped = gzipped_file.replace('.gz', '')
        call(['gunzip -c %s > %s' % (gzipped_file, unzipped)], shell=True)

    @staticmethod
    def _gzip_files(*file_names):
        # todo replace with gzip module
        for f in file_names:
            call(['gzip', f])

    @staticmethod
    def _merge_files(fout, *file_names):
        with open(fout, 'wb') as wfd:
            for f in file_names:
                if f.endswith('.gz'):
                    fd = gzip.open(f, 'rb')
                elif f.endswith('.bz2'):
                    fd = bz2.open(f, 'rb')
                else:
                    fd = open(f, 'rb')
                copyfileobj(fd, wfd, 1024**2 * 128)

    def _generate_multiorganism_coalignment(self):

        # keys for coalignment relevant files
        keys = ['cdna', 'gtf']

        # set names for merged files
        names = {
            'cdna': self.index + '%s_cdna_all.fa' % '_'.join(self.organism),
            'gtf': self.index + '%s.gtf' % '_'.join(self.organism)
        }

        # merge gtf files
        files = []
        for organism in self.organism:
            try:
                files.append(self.index + _file_names[organism]['gtf'])
            except KeyError:
                pass  # no gtf file for this organism
        if not len(files) > 0:
            raise ValueError('user must supply at least one organism with a gtf '
                             'annotation file')
        self._merge_files(names['gtf'], *files)

        # merge cdna files
        cdna_files = []
        for org in self.organism:
            for file_ in _file_names[org]['cdna']:
                cdna_files.append(self.index + file_)
        self._merge_files(names['cdna'], *cdna_files)

        # label the fasta file
        labeled_fasta = self.index + '%s_cdna_all_labeled.fa' % '_'.join(self.organism)
        prepare_fasta.standard_run(names['gtf'], names['cdna'], labeled_fasta)

        # locate and execute the julia script
        julia_script = ('/'.join(seqc.__file__.split('/')[:-2]) +
                        '/scripts/complete_sa.jl')

        # process with julia script
        call(['julia', julia_script, labeled_fasta, self.index, '50'])

        # complete processing
        scid_to_tx = self.index + 'scid_to_feature.txt'
        ordering_out_stacks.standard_run(
            self.index + 'jelly_output_stack_50.txt', scid_to_tx, labeled_fasta,
            self.index + 'p_coalignment.pckl')
        create_gtf_reduced.create_gtf(names['gtf'], scid_to_tx, self.index +
                                      'annotations.gtf')

    def _generate_coalignment(self):
        # todo | some files will need to be ungzipped

        # organism was of length 1, convert to string
        organism = self.organism[0]

        # merge cDNA files
        cdna_files = [self.index + f for f in _file_names[organism]['cdna']]
        merged_cdna = self.index + organism + '_cdna_all.fa'
        self._merge_files(merged_cdna, *cdna_files)

        # # merge ncrna and transcriptome files
        # ncrna = self.index + _file_names[organism]['ncrna']
        # tx = self.index + _file_names[organism]['tx']
        # merged_tx = self.index + organism + '_cdna_all.fa'
        # self._merge_files(merged_tx, ncrna, tx)

        # label the fasta file
        gtf = self.index + _file_names[organism]['gtf']
        labeled_fasta = self.index + organism + '_cnda_all_labeled.fa'
        prepare_fasta.standard_run(gtf, merged_cdna, labeled_fasta)

        # find julia script
        julia_script = ('/'.join(seqc.__file__.split('/')[:-2]) +
                        '/scripts/complete_sa.jl')

        # process with julia script
        call(['julia', julia_script, labeled_fasta, self.index, '50'])

        # complete processing
        scid_to_tx = self.index + 'scid_to_feature.txt'
        ordering_out_stacks.standard_run(
            self.index + 'jelly_output_stack_50.txt', scid_to_tx, labeled_fasta,
            self.index + 'p_coalignment.pckl')
        create_gtf_reduced.create_gtf(gtf, scid_to_tx, self.index + 'annotations.gtf')

    def _build_multiorganism_index(self):
        """"""
        # check that genome has been unzipped
        genome_files = [self.index + _file_names[org]['genome'] for org in self.organism]
        unzipped_files = [f.replace('.gz', '') for f in genome_files]
        for i, f in enumerate(unzipped_files):
            if not os.path.isfile(f):
                self._gunzip_file(genome_files[i])

        # merge genome files
        merged_genome = self.index + '%s_genome.fa' % '_'.join(self.organism)
        self._merge_files(merged_genome, *unzipped_files)

        # check that coalignment files have been generated
        gtf = self.index + 'annotations.gtf'
        if not os.path.isfile(gtf):
            self._generate_coalignment()

        # make index
        star_args = [
            'STAR',
            '--genomeDir', self.index,
            '--runMode', 'genomeGenerate',
            '--runThreadN', str(self.n_threads),
            '--genomeFastaFiles', merged_genome,
            '--sjdbGTFfile', self.index + 'annotations.gtf',
            '--sjdbOverhang', '75']
        star = Popen(star_args, stderr=PIPE)
        _, err = star.communicate()
        if err:
            raise ChildProcessError(err)

    def _build_index(self):
        """"""
        # organism was singular, convert to string
        organism = self.organism[0]

        # check that genome has been unzipped
        genome = self.index + _file_names[organism]['genome']
        unzipped = genome.replace('.gz', '')
        if not os.path.isfile(unzipped):
            self._gunzip_file(genome)

        # check that coalignment files have been generated
        gtf = self.index + 'annotations.gtf'
        if not os.path.isfile(gtf):
            self._generate_coalignment()

        # make index
        star_args = [
            'STAR',
            '--genomeDir', self.index,
            '--runMode', 'genomeGenerate',
            '--runThreadN', str(self.n_threads),
            '--genomeFastaFiles', unzipped,
            '--sjdbGTFfile', self.index + 'annotations.gtf',
            '--sjdbOverhang', '75']
        star = Popen(star_args, stderr=PIPE)
        _, err = star.communicate()
        if err:
            raise ChildProcessError(err)

    def verify_index(self):

        # check for expected index files. This list is not comprehensive
        all_files = ['annotations.gtf', 'p_coalignment.pckl', 'SA', 'SAIndex', 'Genome',
                     'scid_to_feature.txt']
        # check that we have all the files we need to generate the index
        if not all(os.path.isfile(self.index + f) for f in all_files):
            for organism in self.organism:
                for file_type, name in _file_names[organism].items():
                    if isinstance(name, list):
                        for i, n in enumerate(name):
                            if not os.path.isfile(self.index + n):
                                link = _download_links[organism][file_type][i]
                                self._download_ftp_file(link, self.index)
                    else:
                        if not os.path.isfile(self.index + name):
                            link = _download_links[organism][file_type]
                            self._download_ftp_file(link, self.index)

            # generate coalignment if necessary
            coalignment_files = ['annotations.gtf', 'p_coalignment.pckl',
                                 'scid_to_feature.txt']
            if not all(os.path.isfile(self.index + f) for f in coalignment_files):
                if len(self.organism) > 1:
                    self._generate_multiorganism_coalignment()
                else:
                    self._generate_coalignment()

            # build index if necessary
            index_files = ['SA', 'SAIndex', 'Genome']
            if not all(os.path.isfile(self.index + f) for f in index_files):
                if len(self.organism) > 1:
                    self._build_multiorganism_index()
                else:
                    self._build_index()

    def align(self, fastq_records, **kwargs):

        runtime_args = self.default_alignment_args(
            fastq_records, self.n_threads, self.index, self.temp_dir)

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
        for pair in runtime_args.items():
            cmd.extend(pair)

        aln = Popen(cmd, stderr=PIPE, stdout=PIPE)
        out, err = aln.communicate()
        if err:
            raise ChildProcessError(err)

    def clean_up(self):
        rmtree(self.temp_dir)


def alignment_metadata(log_final_out, meta):
    """store a summary of the alignment run written in log_final_out"""

    with open(log_final_out, 'r') as f:
        lines = f.readlines()

        # alignment rates
        meta['mmap_rate'] = float(lines[24].strip().split('\t')[-1][:-1])
        meta['uniq_rate'] = float(lines[9].strip().split('\t')[-1][:-1])
        unmapped_rate = (float(lines[28].strip().split('\t')[-1][:-1]) +
                         float(lines[29].strip().split('\t')[-1][:-1]) +
                         float(lines[30].strip().split('\t')[-1][:-1]))
        meta['unmapped_rate'] = unmapped_rate

        # alignment numbers
        total_reads = int(lines[5].strip().split('\t')[-1])
        meta['total_reads'] = total_reads
        meta['unique_reads'] = int(lines[8].strip().split('\t')[-1])
        meta['mmapped_reads'] = int(lines[23].strip().split('\t')[-1])
        meta['unmapped_reads'] = round(unmapped_rate * total_reads / 100)

        # error rates:
        meta['mismatch_rate'] = float(lines[17].strip().split('\t')[-1][:-1])
        meta['deletion_rate'] = float(lines[18].strip().split('\t')[-1][:-1])
        meta['insertion_rate'] = float(lines[20].strip().split('\t')[-1][:-1])

        # error magnitudes:
        meta['deletion_size'] = float(lines[19].strip().split('\t')[-1])
        meta['insertion_size'] = float(lines[21].strip().split('\t')[-1])

    return meta


if __name__ == "__main__":
    import sys
    if not len(sys.argv) >= 5:
        print('usage: python3 align.py <temp_dir> <n_threads> <index> <organism1> '
              '<organism2> ... <organismN>')
        sys.exit(2)
    temp_dir = sys.argv[1]
    n_threads = int(sys.argv[2])
    index = sys.argv[3]
    organisms = sys.argv[4:]
    s = STAR(temp_dir, n_threads, index, *organisms)
    s.verify_index()