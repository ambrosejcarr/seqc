import os
import requests
import pickle
import pandas as pd
from collections.abc import Iterable


def index_end(b: str, subsequence: str):
    """finds the position in b following the terminate of the first instance of subsequence"""
    return b.index(subsequence) + len(subsequence)


def extract_xml_field(b: str, field_name: str):
    """extract all text contained in an xml field. Calling on <Name>text</Name> will return b'text'."""
    start = index_end(b, '<%s>' % field_name)
    end = b.index('</%s>' % field_name)
    return b[start:end]


class GeneSummary:

    def __init__(self, id_, name, chromosome, synonyms, short, long):
        self._id = id_
        self._name = name
        self._chromosome = chromosome
        self._synonyms = synonyms
        self._short = short
        self._long = long

    def __repr__(self):
        if len(self.long) > 200:
            return '{} (Chr {}): {}\n{}...'.format(self.name, self._chromosome, self.short, self.long[:200])
        else:
            return '{} (Chr {}): {}\n{}'.format(self.name, self._chromosome, self.short, self.long)

    @property
    def id(self):
        return self._id

    @property
    def name(self):
        return self._name

    @property
    def chromosome(self):
        return self._chromosome

    @property
    def synonyms(self):
        return self._synonyms

    @property
    def short(self):
        return self._short

    @property
    def long(self):
        return self._long


# todo fill in way to get different description types
class NCBIGeneSummary:

    def __init__(self, id_type: str, ids: list=None, summaries: dict=None, conversion_file: str=None):
        """Container for NCBIGeneSummaries

        There are three main ways to initialize the object, depending on what kind of gene ids you are starting with:
        :method from_uids:
        :method from_symbols:
        :method from_ensembl:

        Then, there are three main ways to obtain descriptions which return Series objects when multiple ids are
        requested, and single GeneDescripion objects:
        :method describe:

        GeneDescription objects are simple containers, they can be printed in multiple ways, e.g.:
        # todo write up ways of getting different descriptions here

        To save or load a previously created NCBIGeneSummary object, use the methods:
        :method save:
        :method load:

        :param summaries: dictionary of summaries. User must either pass a dict of summaries or a list of IDs that
          should be downloaded
        :param ids: list of ids that should be downloaded. User must either pass a dict of summaries or a list of IDs
          that should be downloaded
        :param id_type: the type of ID that should be input to read summaries. default options = 'HGNC symbol',
          'Ensembl Gene ID', and 'EntrezGene ID'
        :param str conversion_file: Location of a custom conversion file containing alternative id_type(s). Default =
           None, which supports the default types only (see id_type above).
        """
        if summaries and not ids:
            self._id_type = id_type
            self._summaries = summaries
        elif ids and not summaries:
            if id_type != 'EntrezGene ID':
                converter = self._translate_to_uid(ids, id_type, file_location=conversion_file)
            else:
                converter = pd.Series(ids, index=ids)  # pass through conversion
            self._id_type = id_type
            self._summaries = self._download_summaries(ids, converter)

    @property
    def id_type(self):
        return self._id_type

    @staticmethod
    def _download_summaries(ids, converter):
        """download and parse GeneSummary objects from PubMed's EUtils framework

        :param ids: set of ids to download
        :param converter: mapping of ids to uids
        """

        def parse_summary(summary: str):
            """extract gene summary from the NCBI DocSummary for a set of (potentially multiple) genes

            :param str summary: summaries for genes
            :return dict: dictionary mapping gene UIDs to their summaries
            """
            # split by document returned, eliminating the headers (first record)
            records = summary.split('<DocumentSummary uid=')[1:]

            summaries = {}
            for r in records:

                # get uid
                uid = int(r[:r.index('>')].strip('"'))  # contained in DocumentSummary header

                # get gene name
                name = extract_xml_field(r, 'Name')
                chromosome = extract_xml_field(r, 'Chromosome')
                synonyms = extract_xml_field(r, 'OtherDesignations')
                short = extract_xml_field(r, 'Description')
                long = extract_xml_field(r, 'Summary')

                gs = GeneSummary(uid, name, chromosome, synonyms, short, long)
                summaries[inv_converter[uid]] = gs  # converter comes from outer scope

            print(summaries)
            return summaries

        # convert ids to uids
        uids = list(map(int, converter[ids].values))
        print(uids)

        # invert converter to get ids back from uids, used to create summary keys
        inv_converter = pd.Series(converter.index, converter.values)

        # download and parse over uids in sets of 500
        gene_summaries = {}
        while uids:
            current_uids = ','.join(map(str, uids[:500]))
            document_summaries = requests.get(
                'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=gene&id=%s' % current_uids).content
            gene_summaries.update(parse_summary(document_summaries.decode()))
            uids = uids[500:]  # remove processed uids from the list prior to next loop

        print(gene_summaries)
        return gene_summaries

    @staticmethod
    def _parse_gene_conversion_file(file_location=None) -> pd.DataFrame:
        """private method to parse a gene conversion table produced by BioMart using the following parameters:

        Database: Ensemble Genes
        Dataset: Homo Sapiens
        Filters:
          Gene:
            Limit to genes (external references)... with HGNC ids
        Table Attributes:
          Gene:
            Ensembl Gene ID
          External:
            EntrezGene ID
            HGNC Symbol

        This example will have columns labeled 'EntrezGene ID', 'HGNC Symbol', and 'Ensembl Gene ID'

        :return pd.DataFrame: dataframe with one column per attribute selected for the conversion file
        """
        if file_location is None:
            file_location = os.path.expanduser('~/.seqc/tools/human_gene_id_map.csv.gz')
        return pd.DataFrame.from_csv(file_location, index_col=None)

    @classmethod
    def _translate_to_uid(cls, ids, column_name, file_location=None):
        """returns a mapping to EntrezGene ids from column_name, which should describe an ID type defined in a BioMart
        conversion file. Default column_name options are 'HGNC symbol' and 'Ensembl Gene ID'. To use alternate ids,
        the user must create a custom conversion file. See cls._parse_gene_conversion_file() for details on how to
        assemble these files.

        :param list ids: a list of ids to convert from
        :param str column_name: the name of the type of id to convert from, options: ['HGNC symbol', 'Ensembl Gene ID']
        :param str file_location: (default None) an optional location of a custom conversion file
        """
        tr = cls._parse_gene_conversion_file(file_location=file_location).set_index(column_name)
        in_index = tr.index.intersection(ids)
        return tr.ix[in_index, 'EntrezGene ID']

    @classmethod
    def load(cls, filename):
        """load a stored version of this database"""
        with open(filename, 'rb') as f:
            return cls(**pickle.load(f))

    def save(self, filename):
        """save a version of this database"""
        with open(filename, 'wb') as f:
            pickle.dump({'summaries': self._summaries, 'id_type': self._id_type}, f)

    def __repr__(self):
        return 'GeneSummary object containing mappings for %d %ss' % (len(self._summaries), self._id_type)

    def __len__(self):
        return len(self._summaries)

    def __getitem__(self, identifiers):
        """return the GeneSummar(y/ies) associated with uid(s)

        :param identifiers: the identifiers to display summaries for
        :return GeneSummary or pd.Series: the summary associated with a single identifier, or a series mapping
          identifiers to summary objects
        """
        if isinstance(identifiers, Iterable) and not isinstance(identifiers, str):
            return pd.Series([self._summaries[id_] for id_ in identifiers], index=identifiers)
        else:
            return self._summaries[identifiers]

    def describe(self, identifiers):
        """return the GeneSummar(y/ies) associated with uid(s)

        :param identifiers: the identifiers to display summaries for
        :return GeneSummary or pd.Series: the summary associated with a single identifier, or a series mapping
          identifiers to summary objects
        """
        return self[identifiers]
