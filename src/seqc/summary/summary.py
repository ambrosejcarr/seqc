from jinja2 import Environment, PackageLoader
import shutil
from collections import OrderedDict


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

    def from_loss_results(self):



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


def create_summary_section(sections):
    """

    :param tuple sections:
    :return: formatted html template page
    """
    env = Environment(loader=PackageLoader('seqc.summary', 'templates'))
    print(env.list_templates())
    section_template = env.get_template('kv_section.html')
    for s in sections:
        rendered_section = section_template.render(
            sections=sections,
            keys=s.keys,
            values=s.values)
        with open('./src/seqc/summary/testfile_%s.html' % s.name, 'w') as f:
            f.write(rendered_section)