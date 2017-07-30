import subprocess
import os
import shlex
import glob
import re
import numpy as np
import pandas as pd
from scipy.special import expit


class GSEA:

    def __init__(self, correlations, output_stem=None):
        """initialize a gsea object
        :param pd.Series correlations: correlations in the range of [-1, 1] whose index
          contains gene names
        :param str output_stem: the filestem for the output data

        :method linear_scale: method to linearly scale a vector to lie on the interval
          [-1, 1]
        :method logisitc_scale: method to scale a vector by the logistic function to lie
          on the interval [-1, 1]
        :method run: run GSEA on these correlations
        """
        if not isinstance(correlations, pd.Series):
            raise TypeError('correlations must be a pandas series')
        if not ((np.min(correlations) >= -1) & (np.max(correlations) <= 1)):
            raise RuntimeError(
                'input correlations were not contained within the interval [-1, 1]. '
                'Please use JavaGSEA.linear_scale() or JavaGSEA.logistic_scale() to '
                'scale values to this interval before running.')
        self._correlations = correlations.sort_values()
        self._rnk = None
        if output_stem is None:
            self._output_stem = os.environ['TMPDIR'] + 'gsea_corr_{!s}'.format(
                np.random.randint(0, 1000000))
        elif not isinstance(output_stem, str):
            raise TypeError('output stem must be a str reference to a file prefix')
        elif output_stem.find('-') > -1:
            raise ValueError('output_stem cannot contain the dash (-) character.')
        else:
            self._output_stem = output_stem
        self._results = {}

    @property
    def correlations(self):
        return self._correlations

    @correlations.setter
    def correlations(self):
        raise RuntimeError('Please create a new object to compare different correlations')

    @property
    def results(self):
        return self._results

    @staticmethod
    def linear_scale(data: pd.Series) -> pd.Series:
        """scale input vector to interval [-1, 1] using a linear scaling
        :return correlations: pd.Series, data scaled to the interval [-1, 1]
        """
        data = data.copy()
        data -= np.min(data, axis=0)
        data /= np.max(data, axis=0) / 2
        data -= 1
        return data

    @staticmethod
    def logistic_scale(data: pd.Series) -> pd.Series:
        """scale input vector to interval [-1, 1] using a sigmoid scaling
        :return correlations: pd.Series, data scaled to the interval [-1, 1]
        """
        return pd.Series((expit(data.values) * 2) - 1, index=data.index)

    def _save_rank_file(self) -> None:
        """save the correlations to a .rnk file"""
        self._rnk = self._output_stem + '.rnk'
        df = pd.DataFrame(self._correlations).fillna(0)
        df.to_csv(self._rnk, sep='\t', header=False)

    @staticmethod
    def _gmt_options():
        """
        Private method. identifies GMT files available for mouse or human genomes
        :return: str, file options
        """

        mouse_options = os.listdir(os.path.expanduser('~/.seqc/tools/mouse'))
        human_options = os.listdir(os.path.expanduser('~/.seqc/tools/human'))
        print('Available GSEA .gmt files:\n\nmouse:\n{m}\n\nhuman:\n{h}\n'.format(
                m='\n'.join(mouse_options),
                h='\n'.join(human_options)))
        print('Please specify the gmt_file parameter as gmt_file=(organism, filename)')

    def run(self, gmt_file):
        """
        Helper function. Run GSEA on an already-ranked list of corrleations. To see
        available files, leave gmt_file parameter empty

        :param (str, str) gmt_file: organism and filename of gmt file to use
        :return (pd.DataFrame, pd.DataFrame): positive and negative GSEA enrichments
        """
        out_dir, out_prefix = os.path.split(self._output_stem)
        os.makedirs(out_dir, exist_ok=True)

        if self._rnk is None:
            self._save_rank_file()

        if not gmt_file:
            self._gmt_options()
            return
        else:
            if not len(gmt_file) == 2:
                raise ValueError('gmt_file should be a tuple of (organism, filename).')
            else:
                gmt_file = os.path.expanduser('~/.seqc/tools/{}/{}').format(*gmt_file)

        # Construct the GSEA call
        cmd = shlex.split(
            'java -cp {user}/.seqc/tools/gsea2-2.2.1.jar -Xmx1g '
            'xtools.gsea.GseaPreranked -collapse false -mode Max_probe -norm meandiv '
            '-nperm 1000 -include_only_symbols true -make_sets true -plot_top_x 0 '
            '-set_max 500 -set_min 50 -zip_report false -gui false -rnk {rnk} '
            '-rpt_label {out_prefix} -out {out_dir}/ -gmx {gmt_file}'
            ''.format(user=os.path.expanduser('~'), rnk=self._rnk, out_prefix=out_prefix,
                      out_dir=out_dir, gmt_file=gmt_file))

        # Call GSEA
        p = subprocess.Popen(cmd, stderr=subprocess.PIPE)
        _, err = p.communicate()

        # find the file that GSEA created
        if err:
            print(err.decode())
            return
        else:
            pattern = '{p}.GseaPreranked.[0-9]*'.format(p=out_prefix)
            files = os.listdir(out_dir)
            folder = None
            for f in files:
                mo = re.match(pattern, f)
                if mo:
                    folder = out_dir + '/' + mo.group(0)
        if folder is None:
            raise RuntimeError(
                'seqc.JavaGSEA was not able to recover the output of the Java '
                'executable. This likely represents a bug.')

        # recover information from run
        names = ['size', 'es', 'nes', 'p', 'fdr_q', 'fwer_p', 'rank_at_max',
                 'leading_edge']
        pos = pd.DataFrame.from_csv(glob.glob(folder + '/gsea*pos*xls')[0],
                                    sep='\t', infer_datetime_format=False, parse_dates=False).iloc[:, :-1]
        pos.drop(['GS<br> follow link to MSigDB', 'GS DETAILS'], axis=1, inplace=True)
        neg = pd.DataFrame.from_csv(glob.glob(folder + '/gsea*neg*xls')[0],
                                    sep='\t', infer_datetime_format=False, parse_dates=False).iloc[:, :-1]
        neg.drop(['GS<br> follow link to MSigDB', 'GS DETAILS'], axis=1, inplace=True)
        pos.columns, neg.columns = names, names
        self._results[gmt_file] = {'positive': pos, 'negative': neg}
        return list(self._results[gmt_file].values())
