__author__ = 'ambrose'

import gzip
from collections import defaultdict, namedtuple
import pickle


class Reader:

    _Record = namedtuple('Field', ['seqname', 'source', 'feature', 'start', 'end', 'score',
                                   'strand', 'frame', 'attribute'])

    def __init__(self, gtf):
        self._gtf = gtf

    def __iter__(self):
        if self._gtf.endswith('.gz'):
            fobj = gzip.open(self._gtf, 'rt')
        else:
            fobj = open(self._gtf)
        try:
            for line in fobj:
                if line.startswith('#'):
                    continue
                line = line.strip().split('\t')
                meta = line[-1]
                dmeta = {}
                for field in meta.rstrip(';').split(';'):
                    field = field.strip().split()
                    dmeta[field[0]] = field[1].strip('"')
                line[-1] = dmeta
                yield self._Record(*line)
        finally:
            fobj.close()

    def iter_exons(self):
        for record in self:
            if record.feature == 'exon':
                yield record

    def iter_genes(self):
        for record in self:
            if record.feature == 'gene':
                yield record

    def iter_transcripts(self):
        for record in self:
            if record.feature == 'transcript':
                yield record

    def scid_to_gene(self, save=''):
        gmap = defaultdict(set)
        for record in self.iter_transcripts():
            gmap[record.attribute['scseq_id'].strip('SC')].add(
                record.attribute['gene_name'])

        # merge features that are not unique into single strings
        for k, v in gmap.items():
            gmap[k] = '-'.join(v)

        if save:
            with open(save, 'wb') as f:
                pickle.dump(gmap, f)

        return gmap
