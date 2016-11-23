from jinja2 import Environment, PackageLoader


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