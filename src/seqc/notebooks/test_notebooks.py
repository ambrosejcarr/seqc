from . import notebooks
import tempfile
import pytest
import numpy as np
import pandas as pd
import uuid
import os
from seqc.core import main


@pytest.fixture()
def testing_data():
    dir_ = tempfile.mkdtemp()
    test_data = [np.random.randint(10, 110, (100, 100)) for _ in range(4)]
    test_files = []
    for f in test_data:
        filename = '{}/{}'.format(dir_, uuid.uuid4())
        pd.DataFrame(f).to_csv(filename)
        test_files.append(filename)
    return test_files


@pytest.fixture()
def merged_data(testing_data):
    output_stem = os.path.join(tempfile.mkdtemp(), 'test_notebooks')
    n = notebooks.Notebook(output_stem, *testing_data)
    n.merge_data()
    return n.merged_data


def test_template_filling(testing_data):
    output_stem = os.path.join(tempfile.mkdtemp(), 'test_notebooks')
    n = notebooks.Notebook(output_stem, *testing_data)
    n.merge_data()
    n.write_template()
    n.run_notebook()
    print(os.listdir(os.path.dirname(output_stem)))


def test_merge_api(testing_data):
    output_filename = os.path.join(tempfile.mkdtemp(), 'test_notebooks.ipynb')
    args = ['notebook', 'merge', '-o', output_filename, '-i'] + testing_data
    main.main(args)


def test_generate_api(merged_data):
    output_stem = os.path.join(tempfile.mkdtemp(), 'test_notebooks')
    args = ['notebook', 'generate', '-o', output_stem, '-i', merged_data]
    main.main(args)

