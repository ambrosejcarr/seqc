from jinja2 import Environment, PackageLoader
from collections import OrderedDict, namedtuple


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
    def from_filters(cls, read_array):
        """

        :param seqc.read_array.ReadArray read_array: ReadArray object
        :return:
        """

    @classmethod
    def from_read_array(cls, read_array):
        """includes error correction, filtering, multimapping

        :param seqc.read_array.ReadArray read_array:
        :return Section:
        """
        raise NotImplementedError

    def from_dense_filtration(
            self, total_mols, molecules_lost, cells_lost, cell_description):
        """ cell selection, molecule distribution, etc.

        :param total_mols:
        :param molecules_lost:
        :param cells_lost:
        :param cell_description:
        :return Section:
        """
        raise NotImplementedError

    def from_dense_matrix(self):
        raise NotImplementedError




class KVSection:

    def __init__(self, name, content):
        """

        :param name:
        :param dict content:
        """
        self._name = name
        self._content = content

    @property
    def name(self):
        return self._name

    @property
    def content(self):
        return self._content

    @property
    def keys(self):
        return list(self.content.keys())

    @property
    def values(self):
        return list(self.content.values())

    @classmethod
    def from_alignment_summary(cls, alignment_summary):
        alignment_data = OrderedDict
        with open(alignment_summary, 'r') as f:
            for line in f:
                if '|' in line:
                    k, v = [item.strip() for item in line.split('|')]
                    alignment_data[k] = v

        return cls('Alignment Results', alignment_data)


class ImageSection:

    def __init__(self, name, figure, caption, legend):
        """
        :param name:
        :param figure:
        :param caption:
        :param legend:
        """
        self._name = name
        self._figure = figure
        self._caption = caption
        self._legend = legend

    @property
    def name(self):
        return self._name

    @property
    def figure(self):
        return self._figure

    @property
    def caption(self):
        return self._caption

    @property
    def legend(self):
        return self._legend


class Summary:

    def __init__(self, sections):
        self.sections = sections

    def render(self):

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

        pass


def create_seqc_summary(archive_name):
    """creates a SEQC summary from a pre-determined set of sections that are heuristically
    defined.

    :param str archive_name: name for the seqc archive that will be produced in .tar.gz
        format
    :return None:
    """