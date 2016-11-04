import os
import requests
import pickle
import numpy as np
import pandas as pd
from pandas.util.terminal import get_terminal_size
from contextlib import closing
from multiprocessing import Pool
from functools import partial


def index_end(b: str, subsequence: str):
    """
    finds the position in b following the terminate of the first instance of
    subsequence. If not present, raises a ValueError

    :param b:
    :param subsequence:
    :return int: position of the last character of subsequence in b
    """
    return b.index(subsequence) + len(subsequence)


def extract_xml_field(b: str, field_name: str):
    """
    extract all text contained in an xml field.

    For example, Calling on <Name>text</Name> will return b'text'.

    :param str b: string containing xml data
    :param str field_name: name of the field that is being extracted
    :return str: contents of the field
    """
    start = index_end(b, '<%s>' % field_name)
    end = b.index('</%s>' % field_name)
    return b[start:end]


class GeneSummarySeries(pd.Series):
    """simple pd.Series wrapper with a modified line width for better viewing"""

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def __unicode__(self):
        """prints a summary limited to the size of the terminal window"""
        width, _ = get_terminal_size()
        with pd.option_context('display.max_colwidth', width - 15):
            return super().__unicode__()

    def full_length(self):
        """prints a complete summary of the object"""
        with pd.option_context('display.max_colwidth', -1):
            print(super().__unicode__())

    def short(self):
        print(GeneSummarySeries(data=[v.short for v in self.values], index=self.index))

    def long(self):
        print(GeneSummarySeries(data=[v.long for v in self.values], index=self.index))


class GeneSummary:

    def __init__(self, id_, name, chromosome, synonyms, short, long):
        """simple container object for gene summaries

        :param id_:
        :param name:
        :param chromosome:
        :param synonyms:
        :param short:
        :param long:
        """
        self._id = id_
        self._name = name
        self._chromosome = chromosome
        self._synonyms = synonyms
        self._short = short
        self._long = long

    def __repr__(self):
        if len(self.long) > 200:
            return '{} (Chr {}): {}\n{}...'.format(
                self.name, self._chromosome, self.short, self.long[:200])
        else:
            return '{} (Chr {}): {}\n{}'.format(
                self.name, self._chromosome, self.short, self.long)

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

    @property
    def full(self):
        return '{} (Chr {}): {}\n{}'.format(
            self.name, self._chromosome, self.short, self.long)


def _parse_summary(uids: list, inv_converter: pd.Series):
    """obtain and parse DocSummaries for a list of uids

    :param uids: list
    :param inv_converter:
    :return dict: dictionary mapping gene UIDs to their summaries
    """
    # split by document returned, eliminating the headers (first record)
    string_uids = ','.join(map(str, uids))
    summary = requests.get(
        'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?'
        'db=gene&id=%s' % string_uids).content.decode()

    records = summary.split('<DocumentSummary uid=')[1:]

    summaries = {}
    for r in records:

        # get uid
        uid = int(r[:r.index('>')].strip('"'))  # in DocumentSummary header

        # get gene name
        name = extract_xml_field(r, 'Name')
        chromosome = extract_xml_field(r, 'Chromosome')
        synonyms = extract_xml_field(r, 'OtherDesignations')
        short = extract_xml_field(r, 'Description')
        long = extract_xml_field(r, 'Summary')
        if not long:
            long = short

        gs = GeneSummary(uid, name, chromosome, synonyms, short, long)
        try:
            summaries[inv_converter[uid]] = gs
        except TypeError:  # uid had multiple mappings
            for id_ in inv_converter[uid]:
                summaries[id_] = gs

    return summaries


class NCBIGeneSummary:

    def __init__(self, id_type: str, ids: list=None, summaries: dict=None,
                 conversion_file: str=None):
        """Container for NCBIGeneSummaries

        The NCBIGeneSummary object can be initialized from a set of ids and an id_type. If
        an unusual id_type is desired, the user must provide a custom conversion_file that
        maps those ids to EntrezGene IDs.

        Once the object has been created, gene ids can be individually described with the
        method:
        :method describe: returns descriptions of one or more ids.
        :method all: returns descriptions of all contained gene ids.

        GeneDescription objects are simple containers, they can be printed in multiple
        ways, e.g.:
        # todo write up ways of getting different descriptions here

        To save or load a previously created NCBIGeneSummary object, use the methods:
        :method save:
        :method load:

        :param summaries: dictionary of summaries. User must either pass a dict of
          summaries or a list of IDs that should be downloaded
        :param ids: list of ids that should be downloaded. User must either pass a dict
          of summaries or a list of IDs that should be downloaded
        :param id_type: the type of ID that should be input to read summaries. default
          options = 'HGNC symbol', 'Ensembl Gene ID', and 'EntrezGene ID'
        :param str conversion_file: Location of a custom conversion file containing
          alternative id_type(s). Default = None, which supports the default types only
          (see id_type above).
        """
        if summaries and not ids:
            self._id_type = id_type
            self._summaries = summaries
        elif ids and not summaries:
            if id_type != 'EntrezGene ID':
                converter = self._translate_to_uid(
                    ids, id_type, file_location=conversion_file)
            else:
                converter = pd.Series(ids, index=ids)  # pass through conversion

            filtered_ids = self._filter_ids(ids, converter)
            if len(filtered_ids) == 0:
                raise ValueError('No ids could be mapped to entrez gene ids with the '
                                 'provided conversion file.')
            self._id_type = id_type
            self._summaries = self._download_summaries(filtered_ids, converter)
            if len(filtered_ids) < len(ids):
                print('Warning: %d id(s) could not be mapped to entrez gene ids.' % (
                    len(ids) - len(filtered_ids)))
        else:
            raise RuntimeError('User passed both ids and summaries; please select only '
                               'one.')

    @property
    def id_type(self):
        return self._id_type

    @staticmethod
    def _filter_ids(ids, converter):
        """given a converter, filter any genes which lack a mapping between their existing
        ids and EntrezGene ids

        :param ids:
        :param converter:
        :return:
        """
        index_isnt_false = np.array([True if id_ else False for id_ in converter.index])
        values_arent_nan = ~np.isnan(converter)
        values_arent_false = np.array([True if id_ else False for id_ in converter])
        converter = converter[index_isnt_false & values_arent_false & values_arent_nan]
        return list(converter.index.intersection(ids))

    @staticmethod
    def _download_summaries(ids, converter):
        """download and parse GeneSummary objects from PubMed's EUtils framework

        :param ids: set of ids to download
        :param converter: mapping of ids to uids
        """

        def chunks(l, n):
            """Yield successive n-sized chunks from l."""
            for i in range(0, len(l), n):
                yield l[i:i + n]

        # convert ids to uids
        uids = list(map(int, converter[ids].values))

        # invert converter to get ids back from uids, used to create summary keys
        inv_converter = pd.Series(converter.index, converter.values)

        # parallelize the above loop
        uid_groups = list(chunks(uids, 500))

        # create parsing function
        pfunc = partial(_parse_summary, inv_converter=inv_converter)

        with closing(Pool()) as pool:
            res = pool.map(pfunc, uid_groups)
        gene_summaries = {k: v for d in res for k, v in d.items()}  # merge dictionaries

        return gene_summaries

    @staticmethod
    def _parse_gene_conversion_file(file_location=None) -> pd.DataFrame:
        """private method to parse a gene conversion table produced by BioMart using the
        following parameters:

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

        This example will have columns labeled 'EntrezGene ID', 'HGNC Symbol', and
        'Ensembl Gene ID'

        :return pd.DataFrame: dataframe with one column per attribute selected for the
          conversion file
        """
        if file_location is None:
            file_location = os.path.expanduser('~/.seqc/tools/human_gene_id_map.csv.gz')
        return pd.DataFrame.from_csv(file_location, index_col=None)

    @classmethod
    def _translate_to_uid(cls, ids, convert_from, file_location=None):
        """returns a mapping to EntrezGene ids from column_name, which should describe an
        ID type defined in a BioMart conversion file. Default column_name options are
        'HGNC symbol' and 'Ensembl Gene ID'. To use alternate ids, the user must create a
        custom conversion file. See cls._parse_gene_conversion_file() for details on how
        to assemble these files.

        :param list ids: a list of ids to convert from
        :param str convert_from: the name of the type of id to convert from,
          options: ['HGNC symbol', 'Ensembl Gene ID']
        :param str file_location: (default None) an optional location of a custom
          conversion file
        """
        tr = cls._parse_gene_conversion_file(file_location=file_location)
        tr = tr.set_index(convert_from)
        in_index = tr.index.intersection(ids)
        converted = tr.ix[in_index, 'EntrezGene ID']

        # make sure all the converted ids exist!
        converted = converted.ix[~np.isnan(converted)]
        converted = converted.ix[[True if entry else False for entry in converted]]
        return converted

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
        return 'GeneSummary object containing mappings for %d %ss' % (
            len(self._summaries), self._id_type)

    def __len__(self):
        return len(self._summaries)

    def __getitem__(self, identifier):
        """return the GeneSummar(y/ies) associated with uid(s)

        :param identifier: the identifiers to display summaries for
        :return GeneSummary for the provided identifier
        """
        return self._summaries[identifier]

    def describe(self, identifiers, attr='long'):
        """return the GeneSummar(y/ies) associated with uid(s)

        :param identifiers: the identifiers to display summaries for
        :param attr: the GeneSummary attribute to return (e.g. short, long, or chromosome)
          if None, returns the complete GeneSummary object.
        :return pd.Series: a series mapping identifiers to summary objects
        """
        if isinstance(identifiers, str):
            identifiers = [identifiers]
        gs = list(map(lambda i: self._summaries[i], identifiers))
        if attr:
            gs = list(map(lambda i: getattr(i, attr), gs))
        return GeneSummarySeries(gs, index=identifiers)

    def all(self, attr='long'):
        return self.describe(identifiers=list(self._summaries.keys()), attr=attr)
        # keys = list(self._summaries.keys())
        # return pd.Series([self[k] for k in keys], index=keys)
