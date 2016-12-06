from jinja2 import Environment, PackageLoader
from collections import OrderedDict, namedtuple
import numpy as np


ImageContent = namedtuple('ImageContent', ['image', 'caption', 'legend'])

TextContent = namedtuple('TextContent', ['text'])

DataContent = namedtuple('DataContent', ['keys', 'values'])


class Section:

    __slots__ = ['name', 'content']

    def __init__(self, name, content):
        """

        :param str name: Section name
        :param OrderedDict content: ordered dictionary containing keys corresponding to
          header information and values which are Content classes
        """
        self.name = name
        self.content = content

    def render(self, filename):
        """renders this sub-section

        :param str filename: name of rendered html file
        :return None:
        """
        env = Environment(loader=PackageLoader('seqc.summary', 'templates'))
        print(env.list_templates())
        section_template = env.get_template('section_content.html')
        rendered_section = section_template.render(
            sections=[self], section=self)  # todo this API is stupid
        with open(filename, 'w') as f:
            f.write(rendered_section)

    @classmethod
    def from_alignment_summary(cls, alignment_summary):
        """Create a summary section from an alignment summary file produced by STAR

        :param str alignment_summary: STAR alignment summary
        :return:
        """
        with open(alignment_summary, 'r') as f:
            data = f.read()
        categories = OrderedDict()

        def split_lines(block):
            keys, values = zip(*(l.split('|') for l in block.split('\n') if '|' in l))
            keys = [k.strip() for k in keys if k.strip()]
            values = [v.strip() for v in values if v.strip()]
            return DataContent(keys, values)

        # time and input read numbers
        time, _, data = data.partition('UNIQUE READS:')
        categories['Run Time'] = split_lines(time)

        # unique reads and splicing
        splicing, _, data = data.partition('MULTI-MAPPING READS:')
        categories['Unique Reads and Splicing'] = split_lines(splicing)

        # multimapping reads
        multimapping, _, unmapped = data.partition('UNMAPPED READS:')
        categories['Multimapping Reads'] = split_lines(multimapping)
        categories['Unmapped Reads'] = split_lines(unmapped)

        return cls('STAR Alignment Summary', categories)

    @classmethod
    def from_status_filters(cls, ra):
        """run after ReadArray is initialized and initial_filtering() has been run.

        :param ra: ReadArray object
        :return cls: Section containing initial filtering results
        """
        # todo replace whitespace characters with html equiv, add space b/w lines
        description = (
            'Initial filters are run over the sam file while our ReadArray database is '
            'being constructed. These filters indicate heuristic reasons why reads '
            'should be omitted from downstream operations:\n'
            'no gene: Regardless of the read\'s genomic alignment status, there was no '
            'transcriptomic alignment for this read.\n'
            'gene not unique: this indicates that more than one alignment was recovered '
            'for this read. We attempt to resolve these multi-alignments downstream. '
            'primer missing: This is an in-drop specific filter, it indices that the '
            'spacer sequence could not be identified, and thus neither a cell barcode '
            'nor an rmt were recorded for this read.\n'
            'low poly t: the primer did not display enough t-sequence in the primer '
            'tail, where these nucleotides are expected. This indicates an increased '
            'probability that this primer randomly primed, instead of hybridizing with '
            'the poly-a tail of an mRNA molecule.')
        description_section = TextContent(description)

        keys = ('no gene', 'gene not unique', 'primer missing', 'low poly t')
        values = (
            np.sum(ra['status'] & ra.filter_codes['no_gene']),
            np.sum(ra['status'] & ra.filter_codes['gene_not_unique']),
            np.sum(ra['status'] & ra.filter_codes['primer_missing']),
            np.sum(ra['status'] & ra.filter_codes['low_polyt'])
        )
        data_section = DataContent(keys, values)
        return cls(
            'Initial Filtering',
            {'Description': description_section, 'Results': data_section})

    @classmethod
    def cell_barcode_correction(cls, ra):
        """

        :param ra:
        :return:
        """
        description = ''  # todo implement
        description_section = TextContent(description)
        data_section = DataContent(
            ['cell error'],
            [np.sum(ra['status'] & ra.filter_codes['cell_error'])])
        return cls(
            'Cell Barcode Correction',
            {'Description': description_section, 'Results': data_section})

    @classmethod
    def error_correction(cls, read_array):
        raise NotImplementedError  # describe, then summarize

    @classmethod
    def resolve_multiple_alignments(cls, read_array):
        raise NotImplementedError  # describe, then summarize

    @classmethod
    def from_cell_filtering(cls, read_array):
        """

        This is the figure, but it needs a caption!

        :param seqc.read_array.ReadArray read_array: ReadArray object
        :return:
        """
        raise NotImplementedError

    @classmethod
    def from_basic_clustering_and_projection(cls):
        """What if anything do we want to do here?

        # parameterize (default = True?)
        user could specify desire for median norm, PCA + tsne.

        :return:
        """
        raise NotImplementedError

    @classmethod
    def from_overall_yield(cls):
        """

        % losses from original yield summary, can wait.

        :return:
        """
        raise NotImplementedError


class Summary:

    def __init__(self, sections):
        self.sections = sections

    def render(self):
        """This should just loop over the sections and render them all

        :return:
        """
        raise NotImplementedError


def create_seqc_summary(archive_name):
    """creates a SEQC summary from a pre-determined set of sections that are heuristically
    defined.

    # render main section
    # timing, progress, steps complete?

    # render alignment summary
    # parameters, alignment results, speed (maybe some minor figures?)
    # duration

    # render filtering (including error correction)
    # image to show relative filtering amounts
    # duration

    # render multialignment (gain, loss, etc)
    # same as above, show different categories that were recovered?

    # render cell filtering
    # this is already done, figure (add legend)

    # render some QC plots for cells
    # figures can be generate afterwards -- do some PCA, tSNE, filtering (?)

    # store the html including css and js somewhere that can be zipped and uploaded
    # for the user.


    :param str archive_name: name for the seqc archive that will be produced in .tar.gz
        format
    :return None:
    """