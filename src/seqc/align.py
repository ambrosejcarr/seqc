__author__ = 'ambrose'

from subprocess import Popen, PIPE, call, check_output, CalledProcessError
import os
from shutil import rmtree, copyfileobj
from collections import defaultdict
from seqc.sa_preprocess import prepare_fasta
from seqc.sa_postprocess import ordering_out_stacks, create_gtf_reduced
from seqc.log import info
import seqc
import ftplib
import gzip
import bz2


# download links for supported genomes on GEO
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
    return dict(file_names)

_file_names = _define_file_names(_download_links)


def _check_type(obj, type_):
    if not isinstance(obj, type_):
        raise TypeError('%s must be of type %s' % (obj, type_))


class STAR:

    @staticmethod
    def verify_organism(organism):
        valid_organisms = ['mm38', 'hg38', 'mmhg38', 'cs2', 'ci2']
        if organism:
            if not all(org in valid_organisms for org in organism):
                raise ValueError('Invalid organism value. Supported organisms: %r' %
                                 valid_organisms)

    @staticmethod
    def verify_temp_dir(directory):
        if not os.path.isdir(directory):
            os.makedirs(directory)
        if not directory.endswith('/'):
            directory += '/'
        return directory

    @staticmethod
    def default_alignment_args(fastq_records, n_threads, index, temp_dir):
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

    @classmethod
    def _append_phiX_to_fasta(cls, fasta, cdna=False):

        # download phiX genome
        genome_link = ('ftp://ftp.ncbi.nlm.nih.gov/genomes/Viruses/'
                       'Enterobacteria_phage_phiX174_sensu_lato_uid14015/NC_001422.fna')
        phix_genome = cls._download_ftp_file(genome_link, './')

        # add it to the fasta file
        with open(phix_genome, 'rb') as fin:
            data = fin.read()

            if cdna:  # change phix header to reflect our fake transcript id!
                data = data.decode()
                data = data.split('\n')
                data[0] = ('>ENSPXT00000000001 chromosome:NC_001422.1:1:5386:1 '
                           'gene:ENSPXG00000000001 gene_biotype:genome '
                           'transcript_biotype:genome')
                data = ('\n'.join(data)).encode()

            if fasta.endswith('.gz'):
                with gzip.open(fasta, 'ab') as fout:
                    fout.write(data)
            elif fasta.endswith('.bz2'):
                raise ValueError('appending to bz2 files is not currently supported')
            else:
                with open(fasta, 'ab') as fout:
                    fout.write(data)

        os.remove(phix_genome)

    @staticmethod
    def _append_phiX_to_gtf(gtf):
        # fake up some gtf records for a gene and transcript
        gtf_data = '\n'.join([
            '\t'.join([
                'NC_001422.1', 'RefSeq', 'gene', '1', '5386', '.', '+', '.',
                'ID "id0"; Dbxref "taxon:374840"; Is_circular "true"; gbkey "Src"; '
                'genome "genomic"; mol_type "genomic DNA"; nat-host "Escherichia coli"; '
                'gene_id "ENSPXG00000000001";'
            ]), '\t'.join([
                'NC_001422.1', 'RefSeq', 'transcript', '1', '5386', '.', '+', '.',
                'ID "id0"; Dbxref "taxon:374840"; Is_circular "true"; gbkey "Src"; '
                'genome "genomic"; mol_type "genomic DNA"; nat-host "Escherichia coli"; '
                'transcript_id "ENSPXT00000000001"; gene_id "ENSPXG00000000001";'
            ])
        ])

        if gtf.endswith('.gz'):
            with gzip.open(gtf, 'ab') as fout:
                fout.write(gtf_data.encode())
        elif gtf.endswith('.bz2'):
            raise ValueError('appending to bz2 files is not currently supported')
        else:
            with open(gtf, 'ab') as fout:
                fout.write(gtf_data.encode())

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

        return prefix + file_name

    @staticmethod
    def _gunzip_file(gzipped_file):
        # todo want to replace with gzip module, but it is too slow!
        # with gzip.open(gzipped_file, 'rb') as f:
        #     fout = gzipped_file.replace('.gz', '')
        #     data = f.read()
        #     with open(fout, 'wb') as fo:
        #         fo.write(data)
        unzipped = gzipped_file.replace('.gz', '')
        call(['gunzip -c %s > %s' % (gzipped_file, unzipped)], shell=True)

    @staticmethod
    def _gzip_files(*file_names):
        # todo want to replace with gzip module; but it is too slow!
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

    @classmethod
    def _reduce_gtf(cls, gtf):
        """eliminate records with support > 1 in ENSEMBL annotation parlance.

        this is a quick hack that invokes awk from the shell. Will not work on windows.
        """
        call("awk '{if($26~/1/){print $0}}' %s > %s" % (gtf, gtf), shell=True)

    @classmethod
    def _generate_multiorganism_coalignment(cls, index, organism, phix=True):

        # keys for coalignment relevant files
        keys = ['cdna', 'gtf']

        # set names for merged files
        names = {
            'cdna': index + '%s_cdna_all.fa' % '_'.join(organism),
            'gtf': index + '%s.gtf' % '_'.join(organism)
        }

        # merge gtf files
        files = []
        for org in organism:
            try:
                files.append(index + _file_names[org]['gtf'])
            except KeyError:
                pass  # no gtf file for this organism
        if not len(files) > 0:
            raise ValueError('user must supply at least one organism with a gtf '
                             'annotation file')
        cls._merge_files(names['gtf'], *files)

        # merge cdna files
        cdna_files = []
        for org in organism:
            for file_ in _file_names[org]['cdna']:
                cdna_files.append(index + file_)
        cls._merge_files(names['cdna'], *cdna_files)

        # add phix if requested
        if phix:
            cls._append_phiX_to_fasta(names['cdna'], cdna=True)
            cls._append_phiX_to_gtf(names['gtf'])

        # label the fasta file
        labeled_fasta = index + '%s_cdna_all_labeled.fa' % '_'.join(organism)
        prepare_fasta.standard_run(names['gtf'], names['cdna'], labeled_fasta)

        # locate and execute the julia script
        julia_script = ('/'.join(seqc.__file__.split('/')[:-2]) +
                        '/scripts/complete_sa.jl')

        # process with julia script
        info('Generating co-alignment suffix array.')
        call(['julia', julia_script, labeled_fasta, index, '50'])

        # complete processing
        scid_to_tx = index + 'scid_to_feature.txt'
        ordering_out_stacks.standard_run(
            index + 'jelly_output_stack_50.txt', scid_to_tx, labeled_fasta,
            index + 'p_coalignment_array.p')
        create_gtf_reduced.create_gtf(names['gtf'], scid_to_tx, index +
                                      'annotations.gtf')

    @classmethod
    def _generate_coalignment(cls, index, organism, phix=True):

        # organism was of length 1, convert to string
        organism = organism[0]

        # merge cDNA files
        cdna_files = [index + f for f in _file_names[organism]['cdna']]
        merged_cdna = index + organism + '_cdna_all.fa'
        cls._merge_files(merged_cdna, *cdna_files)

        # get the gtf file name
        gtf = index + _file_names[organism]['gtf']
        cls._reduce_gtf(gtf)

        # add phix if requested
        if phix:
            cls._append_phiX_to_fasta(merged_cdna, cdna=True)
            cls._append_phiX_to_gtf(gtf)

        # label the fasta with transcript ids
        labeled_fasta = index + organism + '_cdna_all_labeled.fa'
        prepare_fasta.standard_run(gtf, merged_cdna, labeled_fasta)

        # find julia script
        julia_script = ('/'.join(seqc.__file__.split('/')[:-2]) +
                        '/scripts/complete_sa.jl')

        # process with julia script
        call(['julia', julia_script, labeled_fasta, index, '50'])

        # complete processing
        scid_to_tx = index + 'scid_to_feature.txt'
        ordering_out_stacks.standard_run(
            index + 'jelly_output_stack_50.txt', scid_to_tx, labeled_fasta,
            index + 'p_coalignment_array.p')
        create_gtf_reduced.create_gtf(gtf, scid_to_tx, index + 'annotations.gtf')

    @classmethod
    def _build_multiorganism_index(cls, index, organism, n_threads, phix=True):
        """"""
        # check that genome has been unzipped
        genome_files = [index + _file_names[org]['genome'] for org in organism]
        unzipped_files = [f.replace('.gz', '') for f in genome_files]
        for i, f in enumerate(unzipped_files):
            if not os.path.isfile(f):
                cls._gunzip_file(genome_files[i])

        # merge genome files
        merged_genome = index + '%s_genome.fa' % '_'.join(organism)
        cls._merge_files(merged_genome, *unzipped_files)

        # if requested, add phiX
        if phix:
            cls._append_phiX_to_fasta(merged_genome)

        # check that coalignment files have been generated
        gtf = index + 'annotations.gtf'
        if not os.path.isfile(gtf):
            cls._generate_coalignment(index, organism)  # organism -> list of organisms

        # make index
        info('Beginning to generate STAR index.')
        star_args = [
            'STAR',
            '--genomeDir', index,
            '--runMode', 'genomeGenerate',
            '--runThreadN', str(n_threads),
            '--genomeFastaFiles', merged_genome,
            '--sjdbGTFfile', index + 'annotations.gtf',
            '--sjdbOverhang', '75']
        star = Popen(star_args, stdout=PIPE, stderr=PIPE)
        _, err = star.communicate()
        if err:
            raise ChildProcessError(err)
        info('Finished successfully. Run Complete.')

    @classmethod
    def _build_index(cls, index, organism, n_threads, phix=True):
        """"""
        # organism was singular, convert to string
        organism = organism[0]

        # check that genome has been unzipped
        genome = index + _file_names[organism]['genome']
        unzipped = genome.replace('.gz', '')
        if not os.path.isfile(unzipped):
            cls._gunzip_file(genome)

        # if requested, add phiX
        if phix:
            cls._append_phiX_to_fasta(unzipped)

        # check that coalignment files have been generated
        gtf = index + 'annotations.gtf'
        if not os.path.isfile(gtf):
            cls._generate_coalignment(index, organism, phix=phix)

        # make index
        info('Beginning to generate STAR index.')
        star_args = [
            'STAR',
            '--genomeDir', index,
            '--runMode', 'genomeGenerate',
            '--runThreadN', str(n_threads),
            '--genomeFastaFiles', unzipped,
            '--sjdbGTFfile', index + 'annotations.gtf',
            '--sjdbOverhang', '75']
        star = Popen(star_args, stdout=PIPE, stderr=PIPE)
        out, err = star.communicate()
        if err:
            raise ChildProcessError(err)
        info('Finished successfully. Run Complete.')

    @classmethod
    def build_index(cls, index, organism, n_threads, phix=True, **kwargs):

        # ensure index has a terminal '/'
        if not index.endswith('/'):
            index += '/'

        # check for expected index files. This list is not comprehensive
        all_files = ['annotations.gtf', 'p_coalignment_array.p', 'SA', 'SAIndex',
                     'Genome', 'scid_to_feature.txt']
        for_removal = []  # container for all the files we're downloading.

        info("Downloading genome files.")
        if not all(os.path.isfile(index + f) for f in all_files):
            for org in organism:
                for file_type, name in _file_names[org].items():
                    if isinstance(name, list):
                        for i, n in enumerate(name):
                            if not os.path.isfile(index + n):
                                link = _download_links[org][file_type][i]
                                dlfile = cls._download_ftp_file(link, index)
                                for_removal.extend([dlfile, dlfile.replace('.gz', '')])
                    else:
                        if not os.path.isfile(index + name):
                            link = _download_links[org][file_type]
                            dlfile = cls._download_ftp_file(link, index)
                            for_removal.extend([dlfile, dlfile.replace('.gz', '')])

            # generate coalignment if necessary
            coalignment_files = ['annotations.gtf', 'p_coalignment_array.p',
                                 'scid_to_feature.txt']
            if not all(os.path.isfile(index + f) for f in coalignment_files):
                if len(organism) > 1:
                    cls._generate_multiorganism_coalignment(index, organism, phix=phix)
                else:
                    cls._generate_coalignment(index, organism, phix=phix)

            # build index if necessary
            index_files = ['SA', 'SAIndex', 'Genome']
            if not all(os.path.isfile(index + f) for f in index_files):
                if len(organism) > 1:
                    cls._build_multiorganism_index(index, organism, n_threads, phix=phix)
                else:
                    cls._build_index(index, organism, n_threads, phix=phix)

        # these methods sometimes generate some additional files that should be removed
        for_removal.append(index + '%s_cdna_all.fa' % '_'.join(organism))
        for_removal.append(index + '%s_cdna_all_labeled.fa' % '_'.join(organism))
        for_removal.append(index + '%s_genome.fa' % '_'.join(organism))
        for_removal.append(index + '%s.gtf' % '_'.join(organism))

        for index_file in for_removal:
            if os.path.isfile(index_file):
                os.remove(index_file)

    @staticmethod
    def test_index(index, **kwargs):
        if not index.endswith('/'):
            index += '/'
        critical_files = ['annotations.gtf', 'p_coalignment_array.p', 'SA', 'SAIndex',
                          'Genome', 'scid_to_feature.txt']
        # check that we have all the files we need to generate the index
        if not all(os.path.isfile(index + f) for f in critical_files):
            print('Invalid Index')
        else:
            print('Valid Index')

    @classmethod
    def align(cls, fastq_file, index, n_threads, temp_dir, reverse_fastq_file=None,
              **kwargs):

        # check if file exists; if it does, return the filename
        if os.path.isfile(temp_dir + 'Aligned.out.sam'):
            if os.path.getsize(temp_dir + 'Aligned.out.sam') > 0:
                return temp_dir + 'Aligned.out.sam'

        runtime_args = cls.default_alignment_args(
            fastq_file, n_threads, index, temp_dir)

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

        return temp_dir + 'Aligned.out.sam'

    @classmethod
    def align_multiple_files(cls, fastq_files, index, n_threads, working_dir,
                             reverse_fastq_files=None, **kwargs):

        # use shared memory to map each of the individual cells
        kwargs['--genomeLoad'] = 'LoadAndKeep'

        # make sure forward and reverse fastq file lists match in length
        if reverse_fastq_files:
            if not len(fastq_files) == len(reverse_fastq_files):
                raise ValueError('unequal number of forward and reverse fastq files '
                                 'provided.')

        # make temporary directories for each file
        runs = []
        samfiles = []
        for fastq_file in fastq_files:
            alignment_dir = (working_dir +
                             fastq_file.split('/')[-1].replace('.fastq', '/'))

            # only add the file to runs if the samfile does not already exist
            if os.path.isfile(alignment_dir + 'Aligned.out.sam'):
                if os.path.getsize(alignment_dir + 'Aligned.out.sam') > 0:
                    samfiles.append(alignment_dir + 'Aligned.out.sam')
                    continue

            # file doesn't exist, add it to runs
            try:
                os.mkdir(alignment_dir)
            except FileExistsError:
                pass
            runs.append(cls.default_alignment_args(
                fastq_file, n_threads, index, alignment_dir)
            )
            samfiles.append(alignment_dir + 'Aligned.out.sam')

        # get output filenames

        for runtime_args in runs:
            for k, v in kwargs.items():  # overwrite or add any arguments from cmdline
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

        # construct command line arguments for STAR runs
        cmds = []
        for i, runtime_args in enumerate(runs):
            cmd = ['STAR']
            if reverse_fastq_files:
                for key, value in runtime_args.items():
                    if key == '--readFilesIn':
                        cmd.extend((key, value))
                        cmd.append(reverse_fastq_files[i])
                    else:
                        cmd.extend((key, value))
            else:
                for pair in runtime_args.items():
                    cmd.extend(pair)
            cmds.append(cmd)

        # load shared memory
        cls.load_index(index)
        try:
            for cmd in cmds:
                aln = Popen(cmd, stderr=PIPE, stdout=PIPE)
                out, err = aln.communicate()
                if err:
                    raise ChildProcessError(err.decode())
        finally:
            cls.remove_index(index)

        return samfiles

    @staticmethod
    def load_index(index):
        # set shared memory; this may bug out with larger indices; test!
        try:
            _ = check_output(['sysctl', '-w', 'kernel.shmmax=36301783210'])
            _ = check_output(['sysctl', '-w', 'kernel.shmall=36301783210'])
        except CalledProcessError:
            pass  # not available on OS X
        star_args = ['STAR', '--genomeDir', index, '--genomeLoad',
                     'LoadAndExit']
        star = Popen(star_args, stderr=PIPE, stdout=PIPE)
        return star.communicate()

    @staticmethod
    def remove_index(index):
        # remove index
        star_args = ['STAR', '--genomeDir', index, '--genomeLoad',
                     'Remove']
        star = Popen(star_args, stderr=PIPE, stdout=PIPE)
        out, err = star.communicate()
        if err:  # don't remove temp files, they have logging info
            return err

        # remove temporary files
        call(['rm', '-r', '_STARtmp/'])
        call(['rm', './Aligned.out.sam', './Log.out', './Log.progress.out'])
        return

    @staticmethod
    def clean_up(directory):
        rmtree(directory)


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
