import os

_path = os.path.abspath(__file__)
_dir = os.path.dirname(_path) + '/'


class TwoColumnSubSection:

    def __init__(self, name, keys, values):
        self.name = name
        self.keys = keys
        self.values = values

    def __str__(self):
        return (
            '<h2>{name}</h2>'
            '<div class="row">'
            '<div align="right" class="col-sm-4">{keys}</div>'
            '<div align="left" class="col-sm-8">{values}</div>'
            '</div>'
        ).format(name=self.name,
                 keys='<br>'.join(self.keys),
                 values='<br>'.join(self.values))

    @classmethod
    def from_alignment_summary(cls, alignment_summary):
        section_name = 'Results'
        with open(alignment_summary, 'r') as f:
            data = [l for l in f.readlines() if '|' in l]
        keys = []
        values = []
        for l in data:  # todo change to dict comprehension
            k, v = [string.strip() for string in l.split('|')]
            keys.append(k)
            values.append(v)
        return cls(section_name, keys, values)


class OneColumnSubSection:

    def __init__(self, name, contents):
        self.name = name
        self.contents = contents

    def __str__(self):
        return (
            '<h2>{name}</h2>'
            '<div align="left">{contents}</div>'
        ).format(name=self.name,
                 contents='<br>'.join(self.contents))

    @classmethod
    def from_read_array_representation(cls, ra_repr):
        return cls('Read Array Filtering Summary', ra_repr)


class Section:

    def __init__(self, name, subsections):
        """creates a Section object to be contained within a Summary Object

        :param str name: name of this section
        :param Iterable subsections: iterable of SubSection objects
        """
        self._subsections = subsections
        self._name = name

    def __str__(self):
        """return a string version of an html section for the summary"""
        raise NotImplementedError

    def sidebar_list_item(self):
        """return a hyperlinked sidebar entry for this section"""
        return '<li>\n<a href="#{section_name}">{section_name}</a>\n</li>'.format(
            section_name=self._name)

    def content(self):
        """return paragraph-wrapped content for html results"""
        return (
            '<section id="{section_name}"></section>'
            '<div class="container-fluid">'
            '<div class="row">'
            '<div class ="col-lg-12">'
            '<h1>{section_name}</h1>'
            '{subsections}'  # subsections can be any html
            '</div>'
            '</div>'
            '</div>'.format(
                section_name=self._name,
                subsections='<br>'.join(str(s) for s in self._subsections)))

    def short_summary(self):
        """small summary of this section for the global summary page"""
        raise NotImplementedError  # todo not yet used.


class Summary:

    def __init__(self, sections):
        self._sections = sections
        with open(_dir + '../summary/index.html') as f:
            self._template = f.read()

    # def __repr__(self):
    #     raise NotImplementedError

    def __str__(self):
        return self._template.format(
            sidebar_list_items=self.create_sidebar(),
            summary_sections=self.create_content())

    def add_section(self):
        raise NotImplementedError

    def create_sidebar(self):
        return '\n'.join(s.sidebar_list_item() for s in self.sections)

    def create_content(self):
        return '\n'.join(s.content() for s in self.sections)

    @property
    def sections(self):
        return self._sections


# code to create a star summary object:
# import seqc.summary
# star_summary = seqc.summary.TwoColumnSubSection.from_alignment_summary(
#     '/Users/ajc/Downloads/dj008_alignment_summary.txt')
# s1 = seqc.summary.Section('STAR Alignment', [star_summary])
# s2 = seqc.summary.OneColumnSubSection.from_read_array_representation(ra_repr)
# s2 = seqc.summary.OneColumnSubSection.from_read_array_representation(ra_repr)
# s2 = seqc.summary.Section('All ReadArray Results', [s2])
# summary = seqc.summary.Summary([s1, s2])
# with open('/Users/ajc/projects/seqc/src/summary/testfile.html', 'w') as f:
#     f.write(str(summary))
