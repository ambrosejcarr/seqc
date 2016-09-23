import re
import os
import ftplib
import gzip
import fileinput
import string
import pickle
import wikipedia
from wikipedia import PageError
from contextlib import closing
from multiprocessing import Pool
from collections import defaultdict
from copy import deepcopy
import pandas as pd
import numpy as np
from seqc import reader
from seqc.sparse_frame import SparseFrame
from collections import ChainMap
import time
import intervaltree
import seqc


#TODO: remove unused classes and files.

def first(iterable):
    return next(iter(iterable))


class Record:
    """
    Simple namespace object that makes the fields of a GTF record available. Subclassed
    to create records specific to exons, transcripts, and genes
    """

    __slots__ = ['_fields', '_attribute']

    _del_letters = string.ascii_letters.encode()
    _del_non_letters = ''.join(set(string.printable).difference(string.ascii_letters))\
        .encode()

    def __init__(self, fields: list):

        self._fields = fields
        self._attribute = {}

    def __repr__(self) -> str:
        return '<Record: %s>' % bytes(self).decode()

    def __bytes__(self) -> bytes:
        return b'\t'.join(self._fields)

    def _parse_attribute(self) -> None:
        for field in self._fields[8].rstrip(b';\n').split(b';'):
            key, *value = field.strip().split()
            self._attribute[key] = b' '.join(value).strip(b'"')

    def __hash__(self) -> int:
        """concatenate strand, start, end, and chromosome and hash the resulting bytes"""
        return hash(self._fields[6] + self._fields[3] + self._fields[4] + self._fields[0])

    @property
    def seqname(self) -> bytes:
        return self._fields[0]

    @property
    def chromosome(self) -> bytes:
        return self._fields[0]  # synonym for seqname

    @property
    def source(self) -> bytes:
        return self._fields[1]

    @property
    def feature(self) -> bytes:
        return self._fields[2]

    @property
    def start(self) -> int:
        return int(self._fields[3])

    @property
    def end(self) -> int:
        return int(self._fields[4])

    @property
    def score(self) -> bytes:
        return self._fields[5]

    @property
    def strand(self) -> bytes:
        return self._fields[6]

    @property
    def frame(self) -> bytes:
        return self._fields[7]

    @property
    def size(self) -> int:
        return self.end - self.start

    @property
    def fields(self) -> list:
        return self._fields

    def attribute(self, item):
        """
        access an item from the attribute field of a GTF file.

        :param item: item to access
        :return: value of item
        """
        try:
            return self._attribute[item]
        except KeyError:
            if not self._attribute:
                self._parse_attribute()
                return self._attribute[item]
            else:
                raise KeyError('%s is not a stored attribute of this gtf record' %
                               repr(item))

    @property
    def integer_gene_id(self) -> int:
        """ENSEMBL gene id without the organism specific prefix, encoded as an integer"""
        return int(self.attribute(b'gene_id').split(b'.')[0]
                   .translate(None, self._del_letters))

    @property
    def organism_prefix(self) -> bytes:
        """Organism prefix of ENSEMBL gene id (e.g. ENSG for human, ENSMUSG)"""
        return self.attribute(b'gene_id').translate(None, self._del_non_letters)

    @property
    def string_gene_id(self) -> bytes:
        """ENSEMBL gene id, including organism prefix."""
        return self.attribute(b'gene_id')

    @staticmethod
    def int2str_gene_id(integer_id: int, organism_prefix: bytes) -> bytes:
        """
        converts an integer gene id (suffix) to a string gene id (including organism-
        specific suffix)
        :param organism_prefix: bytes
        :param integer_id: int
        """
        bytestring = str(integer_id).encode()
        diff = 11 - len(bytestring)
        return organism_prefix + (b'0' * diff) + bytestring

    def __eq__(self, other):
        """equivalent to testing if start, end, chrom and strand are the same."""
        return hash(self) == hash(other)

    def __ne__(self, other):
        return not self.__eq__(other)


class Exon(Record):

    def __repr__(self) -> str:
        return '<Exon: %s>' % bytes(self).decode()


class Gene(Record):
    """
    Gene record object. In addition to base properties of the Record object, provides:

    :property exons: set, exon objects associated with this gene
    :method genomic_intervals: iterator, yields tuples of sequence that are covered by
      this gene.

    """

    __slots__ = ['_exons']

    def __init__(self, fields: list):
        super().__init__(fields)
        self._exons = set()

    @property
    def exons(self) -> set:
        return self._exons

    def genomic_intervals(self):
        """
        Iterator, yields tuples of genomic sequence that codes for this gene. In cases
        where multiple exons overlap, the union of the intervals is returned.
        """
        assert self.exons
        ivs = sorted(((e.start, e.end) for e in self.exons))
        ivs = sorted(ivs)
        saved = list(ivs[0])
        for st, en in ivs:
            if st <= saved[1]:
                saved[1] = max(saved[1], en)
            else:
                yield tuple(saved)
                saved[0] = st
                saved[1] = en
        yield tuple(saved)

class GeneIntervals:
    """
    Encodes genomic ranges in an Intervaltree

    :method translate: translates a genomic coordinate on a stranded chromosome into the
      gene identifier that occupies that location (if any exists)

    """

    def __init__(self, gtf: str):
        self._pos_scid_dic = self.build_pos_scid_dic(gtf)
        
    
    def build_pos_scid_dic(self, fname):
        print('using function ver 2')
        start = time.process_time()
        exons = {'-':[],'+':[]}
        res_dic = {} #The resulting dictionary holding two lists. the first is of positions sorted, the second is a list of scids relevant for that position onwards
        new_exons=False

        with open(fname) as infile:
            for line in infile:
                if line[0]=='#':
                    continue
                line_sp = line.split('\t')
                chr = self.strip_chr(line.split('\t')[0])
                if chr==-1:
                    continue
                if chr not in res_dic:
                    res_dic[chr] = {}
                    res_dic[chr]['-'] = intervaltree.IntervalTree()
                    res_dic[chr]['+'] = intervaltree.IntervalTree()
                type = line_sp[2]
                strand = line_sp[6]
                if type not in ['gene','exon','transcript']:
                    continue
                
                if type=='exon':    #Go over all exxons in a transcript and save them for later
                    new_exons = True
                    ex_start = int(line_sp[3])
                    ex_end = int(line_sp[4])                    
                    gene_id = self.strip_gene_num(line_sp[8])
                    if strand == '+':
                        exons[strand].append((ex_start, ex_end, gene_id))
                    elif strand == '-':
                        exons[strand].insert(0,(ex_start, ex_end, gene_id))
                    else:
                        raise ValueError('strand must be - or +')
                else:
                    if new_exons:  #we finished going over the tx, we need to calculate the actual 1000 bases from its end, and log the coordinates of the exons
                        chr = tx_chr  #If we're here we already read the chr and strand of the new tx or gene, these hold the values of the previous tx which are relevant to the current exons list
                        strand = tx_strand
                        tx_end = self.calc_tx_1k_end(exons[strand])
                        
                        for ex_start, ex_end, gene_id in exons[strand]:
                            if ex_end < tx_end:     #The exon ends before the relevant tx part starts - ignore it
                                continue
                            if ex_start < tx_end:   #the exon starts before the 1000 limit but ends after
                                ex_start = tx_end
                            
                            if ex_start==ex_end:    #this can happen sometimes, we don't add empty values to the tree.
                                continue
                                
                            res_dic[chr][strand].addi(ex_start, ex_end, gene_id)
                            
                        exons[strand] = []
                        new_exons = False
                    if type=='transcript':
                        #Keep track of the chr and strand
                        tx_strand = strand
                        tx_chr = chr
                        
        tot_time=time.process_time()-start
        seqc.log.info('Translator completed in {} seconds.'.format(tot_time))
        return res_dic

                    

    def calc_tx_1k_end(self, ex_list):
    # find where is the actual position that's 1k bases from the tx end ()
        remaining_len = 1000
        for (ex_start, ex_end, scid) in reversed(ex_list):
            len = ex_end - ex_start
            if len >= remaining_len:
                return ex_end - remaining_len
            remaining_len -= len
        # If we're here - tx is shorter than 1000 - return the start of the first exxon
        return ex_list[0][0]
    
    
    def translate(self, chr, strand, pos):
        if type(chr) is bytes:
            chr = chr.decode('utf-8')
        int_chr = self.strip_chr(chr)
        if int_chr == -1:   #this happens when it's a scaffold
            return -1
        if type(strand) is bytes:
            strand = strand.decode('utf-8')
            
    # return a list of all scids for the input position
        return [x.data for x in self._pos_scid_dic[int_chr][strand][pos]]
    
    @staticmethod
    def strip_chr(s):
        chr_dic = {'X':23,'Y':24,'MT':25,'M':25}
        if 'chr' in s:
            s=s[3:]
        try:
            int_chr = int(s)
        except ValueError:
            try:
                int_chr = chr_dic[s]
            except KeyError:
                return -1   #this happens when it's a scaffold
        return int_chr

    @staticmethod
    def strip_gene_num(attribute_str):
        gene_start = attribute_str.find('gene_id')
        if gene_start==-1:
            raise ValueError('Gene_id field is missing in annotations file: {}'.format(attribute_str))
        
        gene_end = attribute_str.find(';',gene_start)
        if gene_end==-1:
            raise ValueError('no ; in gene_id attribute, gtf file might be corrupted: {}'.format(attribute_str))
        
        id_start = attribute_str.find('0', gene_start)
        if id_start == -1:
            raise ValueError('Corrupt gene_id field in annotations file - {}'.format(attribute_str))
            
        # ignore the gene version if it's stored as part of the gene id after a decimal point
        id_end = attribute_str.find('.', id_start, gene_end)
        if id_end == -1:
            id_end = gene_end-1 #ignore the "
            
        return int(attribute_str[id_start:id_end])

class Reader(reader.Reader):
    """
    SubClass of reader.Reader, returns an Reader with several specialized iterator
    methods.

    :method __iter__: Iterator over all non-header records in gtf; yields Record objects.
    :method iter_genes: Iterator over all genes in gtf; yields Gene objects.
    """

    def __iter__(self):
        """return an iterator over all non-header records in gtf"""
        hook = fileinput.hook_compressed
        with fileinput.input(self._files, openhook=hook, mode='rb') as f:

            # get rid of header lines
            file_iterator = iter(f)
            first_record = next(file_iterator)
            while first_record.startswith(b'#'):
                first_record = next(file_iterator)
                continue
            yield first_record.split(b'\t')

            for record in file_iterator:  # now, run to exhaustion
                yield record.split(b'\t')

    def iter_genes(self):
        """iterate over all the records for each gene in passed gtf"""

        records = iter(self)

        # get the first gene record
        record = next(records)
        while record[2] != b'gene':
            record = next(records)

        # aggregate exons for each gene
        gene = Gene(record)
        record = next(records)
        while record:
            if record[2] == b'exon':
                gene.exons.add(Exon(record))
            elif record[2] == b'gene':
                yield gene
                gene = Gene(record)
            record = next(records)

        # yield the final gene record
        yield gene


def create_phix_annotation(phix_fasta):
    """
    Several tools in this package require ENSEMBL formatting for .gtf files. However,
    the PhiX genome provided by NCBI does not come with a .gtf file. This tool creates
    a companion gtf file for the phiX fasta file

    :param phix_fasta: str, name of the phix fasta file.
    """
    import numpy as np

    with open(phix_fasta, 'r') as f:
        header = f.readline()  # phiX has only one chromosome
        data = f.readlines()

    # concatenate data
    contig = ''
    for line in data:
        contig += line.strip()

    # get chromosome
    chromosome = header.split()[0].strip('>')
    source = 'seqc'
    score = '.'
    frame = '.'
    gene_meta = 'gene_id "PHIXG00{NUM}"; gene_name "PHIX{NAME!s}";'
    exon_meta = ('gene_id "PHIXG00{NUM}"; gene_name "PHIX{NAME!s}"; '
                 'exon_id "PHIX{NAME!s}";')

    # SEQC truncates genes at 1000b from the end of each transcript. However, phiX DNA
    # that is spiked into an experiment is not subject to library construction. Thus,
    # we will create artificial transcripts for phiX that ensure that all of the DNA is
    # correctly identified.

    length = len(contig)
    transcript_starts = np.arange(length // 1000 + 1) * 1000
    transcript_ends = np.array([min(s + 1000, length) for s in transcript_starts])

    phix_gtf = phix_fasta.replace('.fa', '.gtf')

    with open(phix_gtf, 'w') as f:
        for i, (s, e) in enumerate(zip(transcript_starts, transcript_ends)):
            # add forward strand gene
            gene = [chromosome, source, 'gene', str(s), str(e), score, '+', frame,
                    gene_meta.format(NUM=str(i + 1) * 9, NAME=i + 1)]
            f.write('\t'.join(gene) + '\n')
            exon = [chromosome, source, 'exon', str(s), str(e), score, '+', frame,
                    exon_meta.format(NUM=str(i + 1) * 9, NAME=i + 1)]
            f.write('\t'.join(exon) + '\n')
            # add reverse strand gene
            gene = [chromosome, source, 'gene', str(s), str(e), score, '-', frame,
                    gene_meta.format(NUM=str(i + 1) * 9, NAME=i + 1)]
            f.write('\t'.join(gene) + '\n')
            exon = [chromosome, source, 'exon', str(s), str(e), score, '-', frame,
                    exon_meta.format(NUM=str(i + 1) * 9, NAME=i + 1)]
            f.write('\t'.join(exon) + '\n')


def create_gene_id_to_official_gene_symbol_map(gtf: str):
    """
    create_gene_id_to_official_gene_symbol_map: map integer ENSEMBL ids to
    official gene symbols.

    :param gtf: str, filename of gtf file from which to create the map.
    """
    pattern = re.compile(
        r'(^.*?gene_id "[^0-9]*)([0-9]*)(\.?.*?gene_name ")(.*?)(".*?$)')
    gene_id_map = defaultdict(set)
    with open(gtf, 'r') as f:
        for line in f:
            match = re.match(pattern, line)
            if match:
                gene_id_map[int(match.group(2))].add(match.group(4).upper())
    return gene_id_map


def ensembl_gene_id_to_official_gene_symbol(
        data: pd.DataFrame or SparseFrame, gtf=None, gene_id_map=None):
    """convert data containing ensembl gene ids into an index of gene symbols

    :param data: cells x genes DataFrame or SparseFrame
    :param gtf: str, filename of a gtf file
    :param gene_id_map: gene_id_map constructed from
      Experiment.create_gene_id_to_official_gene_symbol_map. If converting multiple
      objects, it is much faster to only construct the map a single time.
    :return: data: pd.DataFrame or SparseFrame object
    """
    if gene_id_map is None:
        if gtf is None:
            raise ValueError('User must pass either GTF or a gene_id_map object')
        gene_id_map = create_gene_id_to_official_gene_symbol_map(gtf)
    data = deepcopy(data)
    data.columns = ['-'.join(gene_id_map[i]) for i in data.columns]

    return data


class GeneDescription:

    def __init__(self, name, short, long):
        """Simple container class for RefSeq gene descriptions

        :param str name: Gene symbol
        :param str short: Long-form gene name
        :param str long: RefSeq gene summary
        :param [str] synonyms: alternative gene names
        """
        self.name = name
        self.short = short
        self.long = long

    def __repr__(self):
        if len(self.long) > 200:
            return '{}: {}\n{}...'.format(self.name, self.short, self.long[:200])
        else:
            return '{}: {}\n{}'.format(self.name, self.short, self.long)

    @property
    def full(self):
        """return complete gene summary, regardless of length"""
        return '{}: {}\n{}'.format(self.name, self.short, self.long)


class GeneInfo:

    def __init__(self, data: dict):
        """create a GeneInfo object

        alternative constructors:
        -------------------------
        :method load: load a previously saved serialized GeneInfo object from file
        :method download: download the gene information from NCBI and create a new
          GeneInfo object

        :param data: dictionary of gene symbol: GeneDescription objects
        """
        self._data = data

    def __repr__(self):
        raise NotImplementedError  # need to implement this!

    def __getitem__(self, item):
        try:
            return self._data[item].full
        except KeyError:
            print('No description for that gene id.')
            return None

    @staticmethod
    def parse_gbff(filename):
        """fast but memory-inefficient reader for gbff files; typically ~ 1G per process
        :param filename:
        :return:
        """
        with gzip.open(filename, 'rb') as f:
            data = f.read()

        gene_to_summary = {}  # results dictionary

        # summary markers
        summary_start_pattern = b'  Summary: '  # include two spaces just in case
        summary_end_pattern = b'\nPRIMARY'
        sum_patt_len = len(summary_start_pattern)

        # description markers
        desc_start_pattern = b'\nDEFINITION'
        desc_end_pattern = b'\nACCESSION'
        desc_patt_len = len(desc_start_pattern)

        # keep track of failures for development
        no_description = 0
        no_gene_id = 0
        no_summary = 0

        # split into records; each begins with LOCUS
        records = data.strip(b'LOCUS').split(b'LOCUS')
        for record in records:
            try:  # get short gene description
                desc_start = record.index(desc_start_pattern) + desc_patt_len
                desc_end = record.index(desc_end_pattern)
                description = record[desc_start:desc_end]
                description = b' '.join(line.strip() for line in
                                        description.split(b'\n')).decode()
            except ValueError:  # description was not found
                no_description += 1
                continue

            try:  # get gene id
                gene_id_start = record.index(b'/gene="') + 7
                gene_id_end = record[gene_id_start:].index(b'"\n')
                gene_id = record[gene_id_start:][:gene_id_end].decode()
            except ValueError:
                no_gene_id += 1
                continue

            try:  # get long gene summary
                summary_start = record.index(summary_start_pattern) + sum_patt_len
                summary_end = record.index(summary_end_pattern)
                summary = record[summary_start:summary_end]
                summary = b' '.join(line.strip() for line in summary.split(b'\n')).decode()
            except ValueError:  # summary was not found
                no_summary += 1
                continue

            gene_to_summary[gene_id] = GeneDescription(gene_id, description, summary)

        return gene_to_summary

    @staticmethod
    def download_gbff(link, email):
        """map gene ids to summary information that is found in the RefSeqGene database.

        :param str link: folder on NCBI ftp containing refseq.gbff files to be downloaded
          and parsed. At the time of implementation, the human link is:
          ftp://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/RefSeqGene/
        :param email: email address; used by ncbi as a password
        """

        # parse the link for relevant ftp information
        server = link.replace('ftp://', '').split('/')[0]
        folder = '/' + '/'.join(link.replace('ftp://', '').split('/')[1:])
        username = 'anonymous'
        password = email
        download_dir = os.environ['TMPDIR']

        # start ftp; todo make contextmanaged with 'with'
        ftp = ftplib.FTP(server)
        ftp.login(username, password)
        ftp.cwd(folder)
        files = ftp.nlst()

        # parse file list, only download gbff files
        files = [f.split()[-1] for f in files if '.gbff.gz' in f]

        # download the files
        downloaded_files = [download_dir + f for f in files]
        for fin, fout in zip(files, downloaded_files):
            with open(download_dir + fin, 'wb') as open_file:
                ftp.retrbinary('RETR %s' % fin, open_file.write)
        ftp.close()

        return downloaded_files

    def save(self, filestem):
        """save a binary version of the GeneInfo object. The file will be saved with a .p
        suffix

        :param str filestem: file stem for serialized data object
        """
        with open(filestem + '.p', 'wb') as f:
            pickle.dump(self._data, f)

    @classmethod
    def load(cls, serialized_data_file):
        with open(serialized_data_file, 'rb') as f:
            data = pickle.load(f)
        return cls(data)

    @classmethod
    def download(cls, link, email):
        """download data for the GeneInfo object from NCBI

        :param str link: folder on NCBI ftp containing refseq.gbff files to be downloaded
          and parsed. At the time of implementation, the human link is:
          ftp://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/RefSeqGene/
        :param email: email address; used by ncbi as a password
        """
        files = cls.download_gbff(link, email)
        with closing(Pool()) as pool:
            dicts = pool.map(cls.parse_gbff, files)
        data = dict(ChainMap(*dicts))
        return cls(data)

    def _get_id_type(self, attr_name):
        """factory method to generate a function that gets attr_name for each id"""

        def get_id(ids):
            # user passed a list
            if len(ids) == 1 and isinstance(ids[0], (list, tuple, np.ndarray)):
                ids = ids[0]

            # user passed a single id
            if len(ids) == 1:
                try:
                    return getattr(self._data[ids[0]], attr_name)
                except KeyError:
                    return None

            # user passed a list of ids
            else:
                desc = []
                for id_ in ids:
                    try:
                        desc.append(getattr(self._data[id_], attr_name))
                    except KeyError:
                        desc.append(None)
                return pd.Series(desc, index=ids)

        return get_id

    def add_from_wikipedia(self, ids):
        added = []
        already_present = 0
        for id_ in ids:
            if id_ not in self._data.keys():
                try:
                    data = wikipedia.page(id_).summary
                    self._data[id_] = GeneDescription(id_, '', data)
                    added.append(id_)
                except PageError:
                    pass
            else:
                already_present += 1
        print('Added {!s} descriptions from wikipedia. {!s} were already present, {!s} '
              'could not be found.'.format(len(added), already_present,
                                           len(ids) - (len(added) + already_present)))
        return added

    def short(self, *ids):
        function = self._get_id_type('short')
        return function(ids)

    def long(self, *ids):
        function = self._get_id_type('long')
        return function(ids)

    def full(self, *ids):
        function = self._get_id_type('full')
        return function(ids)

    def synonynyms(self, *ids):
        function = self._get_id_type('synonyms')
        return function(ids)
