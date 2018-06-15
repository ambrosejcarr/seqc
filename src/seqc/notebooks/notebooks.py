from jinja2 import Environment, FileSystemLoader
import os
import pandas as pd
import tempfile

import nbformat
from nbconvert.preprocessors import ExecutePreprocessor


class Notebook:

    def __init__(self, output_stem: str, *data):

        # strip notebook affix if user provided it; this is a common error mode
        if output_stem.endswith('.ipynb'):
            output_stem = output_stem.replace('.ipynb', '')
        self._output_stem = output_stem

        self._data = data
        self._this_dir = os.path.dirname(os.path.abspath(__file__))

    @property
    def notebook_path(self):
        return self._output_stem + '.ipynb'

    @property
    def merged_data(self):
        if isinstance(self._data, str):
            if os.path.isfile(self._data):
                return os.path.abspath(self._data)
        elif isinstance(self._data, (list, tuple)) and isinstance(self._data[0], str):
            if os.path.isfile(self._data[0]):
                return os.path.abspath(self._data[0])
        raise TypeError('Data is not a 1-length iterable or string that contains a filepath')

    def merge_data(self, merged_sample_name=None, remove_unmerged=False):
        """
        This function will merge any datasets provided as nested lists.
        Each top-level value is considered an input alias.
        Any second-level list is considered to be a group of files to be joined

        :param bool remove_unmerged: if True, this function will delete the unmerged files after
          completion
        :param str merged_sample_name: name of merged csv file
        :return None: The list of merged file names will replace the list passed to the class in
          self._datasets
        """
        dfs = [pd.read_csv(csv, index_col=0) for csv in self._data]
        df = pd.concat(
            dfs,
            keys=list(range(len(self._data))),
            names=['sample_number', 'cell_id']
        )

        if not merged_sample_name:
            merged_sample_name = self._output_stem + '_merged_data.csv'
        df.to_csv(merged_sample_name)

        # delete original files, if requested
        if remove_unmerged:
            for csv in self._data:
                os.remove(csv)

        # update file urns
        self._data = merged_sample_name

    def write_template(self):
        """write a filled ipython notebook to disk

        :return:
        """

        j2_env = Environment(loader=FileSystemLoader(self._this_dir), trim_blocks=True)
        rendered = j2_env.get_template('analysis_template.json').render(
            output_stem=self._output_stem,
            data=os.path.abspath(self.merged_data),
        )
        with open(self._output_stem + '.ipynb', 'w') as fdw:
            fdw.write(rendered)

    def run_notebook(self, notebook_filename=None):

        if not notebook_filename:
            notebook_filename = self._output_stem + '.ipynb'

        dir_ = os.getcwd()
        with open(notebook_filename) as f:
            nb = nbformat.read(f, as_version=4)

        ep = ExecutePreprocessor(timeout=600, kernel_name='python3')
        ep.preprocess(nb, {'metadata': {'path': dir_}})

        with open(notebook_filename, 'wt') as f:
            nbformat.write(nb, f)
