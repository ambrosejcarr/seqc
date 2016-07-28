import fileinput
import string
import bisect
import time

import seqc


def first(iterable):
    return next(iter(iterable))


class Record:

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
        return int(self.attribute(b'gene_id').split(b'.')[0]
                   .translate(None, self._del_letters))

    @property
    def organism_prefix(self) -> bytes:
        return self.attribute(b'gene_id').translate(None, self._del_non_letters)

    @property
    def string_gene_id(self) -> bytes:
        return self.attribute(b'gene_id')

    @staticmethod
    def int2str_gene_id(integer_id: int, organism_prefix: bytes) -> bytes:
        """converts an integer gene id to a string gene id
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

    __slots__ = ['_exons']

    def __init__(self, fields: list):
        super().__init__(fields)
        self._exons = set()

    @property
    def exons(self) -> set:
        return self._exons

    def genomic_intervals(self):
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

    def __init__(self, gtf: str):
        self._pos_scid_dic = self.build_pos_scid_dic(gtf)
        
    def build_pos_scid_dic(self, fname):
        start = time.process_time()
        exons = []
        res_dic = {'pos':{}, 'scid':{}} #The resulting dictionary holding two lists. the first is of positions sorted, the second is a list of scids relevant for that position onwards

        with open(fname) as infile:
            for line in infile:
                if line[0]=='#':
                    continue
                type = line.split('\t')[2]
                if type=='exon':    #Go over all exxons in a transcript and save them for later
                    ex_start = int(line.split('\t')[3])
                    ex_end = int(line.split('\t')[4])
                    scid_ind = line.split('\t')[8].find('scseq_id')
                    scid = line.split('\t')[8][scid_ind+10:-3]
                    strand = line.split('\t')[6]
                    if strand == '+':
                        exons.append((ex_start, ex_end, scid))
                    elif strand == '-':
                        exons.insert(0,(ex_start, ex_end, scid))
                    else:
                        print('strand is neither - nor + ',strand)
                else:
                    if len(exons) > 0:  #we finished going over the tx, we need to calculate thae actual 1000 bases from its end, and log the coordinates of the exons
                        tx_end = self.calc_tx_1k_end(exons)
                        for ex_start, ex_end, scid in exons:
                            if ex_end < tx_end:     #The exon ends before the relevant tx part starts
                                continue
                            if ex_start < tx_end:   #the exon starts before the 1000 limit but ends after
                                ex_start = tx_end
                        
                            try:
                                ind_start = bisect.bisect_left(res_dic['pos'][chr], ex_start)
                            except KeyError:
                                res_dic['pos'][chr] = [0]
                                res_dic['scid'][chr] = [[]]
                                ind_start = 1
                            ind_end = ind_start + 1 # This is true since exons don't overlap (or at least shouldn't)
                            if ex_start in res_dic['pos'][chr]:
                                if scid not in res_dic['scid'][chr][ind_start]:
                                    res_dic['scid'][chr][ind_start].append(scid)
                            else:
                                res_dic['pos'][chr].insert(ind_start, ex_start)
                                prev_l = list(res_dic['scid'][chr][ind_start-1])
                                if scid not in prev_l:
                                    prev_l.append(scid)
                                res_dic['scid'][chr].insert(ind_start, prev_l)
                            
                            
                            if ex_end not in res_dic['pos'][chr]:   #There is nothing to do if the position is already marked, but if it isn't we need to add it to mark end of current scid
                                res_dic['pos'][chr].insert(ind_end, ex_end)
                                prev_l = list(res_dic['scid'][chr][ind_end-1])
                                prev_l.remove(scid)
                                res_dic['scid'][chr].insert(ind_end, prev_l)
                        exons = []
                    if type=='gene':
                        try:
                            chr = int(line.split('\t')[0])
                        except ValueError:
                            #chr is 'X'
                            chr = 23
        tot_time=time.process_time()-start
        #print('time: ',tot_time)
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
    
    #TODO: Make sure what to do with the strand in translate (I take care of it when building the dic, but what is it used for in translate?)
    #TODO: when passing a position thats after the end of the chr it still returns one scid
    def translate(self, chr, strand, pos):
        chr_dic = {'X':23,'Y':24,'MT':25}
        try:
            int_chr = int(chr)
        except ValueError:
            try:
                int_chr = chr_dic[chr]
            except KeyError:
                return -1	#this happens when it's a scaffold         
            
    # return a list of all scids for the input position
        ind = bisect.bisect_right(self._pos_scid_dic['pos'][int_chr], pos) - 1
        return self._pos_scid_dic['scid'][int_chr][ind]


class Reader(seqc.reader.Reader):

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
        """iterate over all the records for each gene in passed gtf

        :yields: Gene
        """

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
