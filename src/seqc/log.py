import json
import logging
from datetime import datetime
import pandas as pd
import seqc.stats
import seqc.exceptions
from collections import defaultdict
import os
import re


def setup_logger(filename):
    """create a simple log file in the cwd to track progress and any errors"""
    logging.basicConfig(filename=filename, level=logging.DEBUG, filemode='w')


def info(message):
    """print a timestamped update for the user.
    :param message:
    """
    logging.info(datetime.now().strftime("%Y-%m-%d %H:%M:%S") + ':' + message)


def exception():
    """log the most recent exception to an initialized logger"""
    logging.exception(datetime.now().strftime("%Y-%m-%d %H:%M:%S") + ':main:')


def notify(message):
    """print a timestamped update for the user and log it to file"""
    info(message)
    print('SEQC: ' + datetime.now().strftime("%Y-%m-%d %H:%M:%S") + ': %s' % message)


def debug(message):
    logging.debug(datetime.now().strftime("%Y-%m-%d %H:%M:%S") + ':%(module)s:%(funcName)s:' + ': %s' % message)


def print_exact_command_line(arg_line, remote=False):
    """
    log command line for easy copy pasting
    :param arg_line A string corresponding to the command that was called
    :param remote Show the command will be used remotely
    """
    remote_line = "The following command line will be run on the remote server"
    local_line = "Executing the following command line"
    line = remote_line if remote else local_line
    notify(line + ":\n\t> " + arg_line)


def args(arguments):
    """
    log namespace object from argument parser to file.

    :param arguments: namespace object, output of ArgumentParser.parse_args()
    :return: None
    """
    arguments = vars(arguments)
    info('Passed command line arguments: {}'.format(
            json.dumps(arguments, separators=(',', ': '), indent=4)))


class LogData:
    """
    Automatically parse SEQC logs

    :method parse_log: parse an individual seqc log into a pd.DataFrame with a MultiIndex
      index corresponding to the categories of the seqc log
    :method parse_multiple: parse a directory hierarchy into a pd.DataFrame with an index
      as above, but the columns recapitulate the directory structure that the logs are
      stored in.
    """

    _oldver = ('{divide}\nINPUT\n{divide}\n'
               'Total input reads:\t{n_fastq}\n'
               '{divide}\nALIGNMENT (% FROM INPUT)\n{divide}\n'
               'Total reads aligned:\t{n_sam} ({prop_al}%)\n'
               ' - Genomic alignments:\t{genomic} ({prop_gen}%)\n'
               ' - PhiX alignments:\t{phi_x} ({prop_phix}%)\n'
               ' - Transcriptome alignments:\t{trans} ({prop_trans}%)\n'
               '{divide}\nFILTERING (% FROM ALIGNMENT)\n{divide}\n'
               'Genomic alignments:\t{genomic} ({bad_gen}%)\n'
               'PhiX alignments:\t{phi_x} ({bad_phi}%)\n'
               'Incorrect barcodes:\t{wrong_cb} ({bad_cb}%)\n'
               'Missing cell barcodes:\t{no_cell} ({bad_cell}%)\n'
               'Missing RMTs (same as above):\t{no_cell} ({bad_cell}%)\n'
               'N present in RMT:\t{rmt_N} ({bad_rmtN}%)\n'
               'Insufficient poly(T):\t{poly_t} ({bad_polyt}%)\n'
               '{divide}\nCELL/MOLECULE COUNT DISTRIBUTION\n{divide}\n'
               'Total molecules:\t\t{tot_mc}\n'
               'Molecules lost:\t{mols_lost}\n'
               'Cells lost:\t{cells_lost}\n'
               'Cell description:\n{cell_desc}\n'
               '{divide}\nSUMMARY\n{divide}\n'
               'Total retained reads:\t{n_good} ({prop_good}%)\n'
               'Total reads unaligned:\t{lost_al} ({prop_un}%)\n'
               'Total reads filtered:\t{n_bad} ({prop_bad}%)\n'
               '{divide}\n'
               )

    @staticmethod
    def string_to_regex(summary: str=None) -> str:
        """
        converts the contents of seqc.stats.ExperimentalYield.output into a regex object
        that may contain duplicate definitions
        :param summary: str, optional, a summary to convert
        :return summary: str, a regex object that may contain errors
        """
        if not summary:
            summary = seqc.stats.ExperimentalYield.output
        replacements = [
            ('{divide}', '-*?'),
            ('(', '\('),
            (')', '\)'),
            ('{', '(?P<'),
            ('}', '>.*?)')
        ]
        for r in replacements:
            summary = summary.replace(r[0], r[1])
        return summary

    @staticmethod
    def identify_duplicate_patterns(regex: str) -> dict:
        """
        identifies replicated name patterns, which are not allowed in regex

        :param regex: str, pattern in which to find replicated assignments
        :return replicates: dict, contains names of replicated definition with values
          equal to the number of times each was replicated
        """

        name_pattern = '(\(\?P<)(.*?)(>\.\*\?\))'
        patterns = set()
        replicates = defaultdict(int)
        for mo in re.finditer(name_pattern, regex):
            if mo.group(2) in patterns:
                replicates[mo.group(2)] += 1
            else:
                patterns.add(mo.group(2))
        return replicates

    @staticmethod
    def replace_replicated_patterns(regex: str, duplicated_pattern: str) -> str:
        """
        replace the second definition of pattern_name with a regex-compliant reference

        :param regex: str, regex containing replicated pattern
        :param duplicated_pattern: pattern_id
        :return regex: str, pattern without duplicate group definitions
        """
        old = '(?P<{}>.*?)'.format(duplicated_pattern)
        new = '(?P={})'.format(duplicated_pattern)
        idx = regex.find(old) + len(old)
        return regex[:idx] + regex[idx:].replace(old, new)


    @classmethod
    def dictionary_to_dataframe(cls, groupdict, col_label) -> pd.DataFrame:
        """
        Warning: This function contains summary-specific information and may break or
        need to be modified when the summary is changed. This function translates
        the format parameters into interpretable columns in the dataframe

        :param groupdict: result of groupdict() call on match object generated by the
          parse_log() classmethod
        :param col_label: name of log file
        :return: pd.DataFrame containing log data
        """
        index = (
            ('total', 'input_reads'),
            ('total', 'reads_aligned'),
            ('aligned', 'genomic'),
            ('aligned', 'phi_x'),
            ('aligned', 'transcriptome'),
            ('filtered', 'genomic'),
            ('filtered', 'phi_x'),
            ('filtered', 'incorrect_barcodes'),
            ('filtered', 'no_barcodes'),
            ('filtered', 'CB_contains_N'),
            ('filtered', 'RMT_contains_N'),
            ('filtered', 'broken_capture_primer'),
            ('filtered', 'low_complexity'),
            ('summary', 'reads_retained'),
            ('summary', 'reads_not_aligned'),
            ('summary', 'reads_filtered'),
            ('summary', 'total_molecules')
        )
        data_list = ('n_fastq', 'n_sam', 'genomic', 'phi_x', 'trans', 'genomic', 'phi_x',
                     'wrong_cb', 'no_cell', 'cell_N', 'rmt_N', 'poly_t', 'dust', 'n_good',
                     'lost_al', 'n_bad', 'tot_mc')

        # account for older log version
        if groupdict['wrong_cb'] == 'NA':
            groupdict['wrong_cb'] = 0

        data = list(map(lambda x: float(groupdict[x]), data_list))

        spec_index, spec_data = cls.parse_special_fields(groupdict)

        index = pd.MultiIndex.from_tuples(index + spec_index)
        return pd.DataFrame(data + spec_data, index, columns=[col_label])

    @staticmethod
    def parse_special_fields(groupdict: dict) -> (tuple, list):
        """
        extracts information from special fields in run summary and returns
        a tuple index suitable for multiindex creation and string representations of
        data.

        :param groupdict: result of groupdict() call on match object generated by the
          parse_log() classmethod
        :returns index, data_list: (tuple, list)
        """
        lost_pattern = (
            "^\[\('low_count', (?P<low_count>[0-9]+)\), "
            "\('low_coverage', (?P<low_coverage>[0-9]+)\), "
            "\('high_mt', (?P<high_mt>[0-9]+)\), "
            "\('low_gene_detection', (?P<low_gene_detection>[0-9]+)\)\]$")

        summary_pattern = (
            "^count\s+(?P<count>[0-9]+\.[0-9]+)\s"
            "mean\s+(?P<mean>[0-9]+\.[0-9]+)\s"
            "std\s+(?P<std>[0-9]+\.[0-9]+)\s"
            "min\s+(?P<min>[0-9]+\.[0-9]+)\s"
            "25%\s+(?P<low_quartile>[0-9]+\.[0-9]+)\s"
            "50%\s+(?P<median>[0-9]+\.[0-9]+)\s"
            "75%\s+(?P<high_quartile>[0-9]+\.[0-9]+)\s"
            "max\s+(?P<max>[0-9]+\.[0-9]+)\s?")

        cell = re.match(lost_pattern, groupdict['cells_lost'], re.M).groupdict()
        mols = re.match(lost_pattern, groupdict['mols_lost'], re.M).groupdict()
        desc = re.match(summary_pattern, groupdict['cell_desc'], re.M).groupdict()

        if not all((cell, mols, desc)):
            raise ValueError('Regex failed to match log. Please check that you are using '
                             'a matched log/seqc pair.')

        index = (
            ('molecules_lost', 'low_count'),
            ('molecules_lost', 'low_coverage'),
            ('molecules_lost', 'high_mt'),
            ('molecules_lost', 'low_gene_detection'),
            ('cells_lost', 'low_count'),
            ('cells_lost', 'low_coverage'),
            ('cells_lost', 'high_mt'),
            ('cells_lost', 'low_gene_detection'),
            ('cell_summary', 'count'),
            ('cell_summary', 'mean'),
            ('cell_summary', 'std'),
            ('cell_summary', 'min'),
            ('cell_summary', '25%'),
            ('cell_summary', '50%'),
            ('cell_summary', '75%'),
            ('cell_summary', 'max')
        )
        data_list = (
            mols['low_count'], mols['low_coverage'], mols['high_mt'],
            mols['low_gene_detection'], cell['low_count'], cell['low_coverage'],
            cell['high_mt'], cell['low_gene_detection'], desc['count'], desc['mean'],
            desc['std'], desc['min'], desc['low_quartile'], desc['median'],
            desc['high_quartile'], desc['max'])
        data_list = list(map(lambda x: float(x), data_list))
        return index, data_list

    @classmethod
    def match_log(cls, log_file: str, pattern: str=None) -> dict:
        """
        create a dictionary to hold data from SEQC summary.

        :param log_file: name of the seqc log to extract information from
        :param pattern: str, optional, the value of seqc.stats.ExperimentYield.output.
          useful to parse seqc logs from older versions
        :return data: dict, argument names and values from seqc.log summary
        """
        if pattern is None:
            pattern = cls.string_to_regex()

        def get_match_object(pattern_):
            duplicates = cls.identify_duplicate_patterns(pattern_)

            for k, v in duplicates.items():
                pattern_ = cls.replace_replicated_patterns(pattern_, k)

            # add beginning and end wildcards
            pattern_ = '^.*?' + pattern_ + '.*?$'
            with open(log_file, 'r') as f:
                summary_data = f.read()
                mo = re.match(pattern_, summary_data, re.M | re.DOTALL)
                data = mo.groupdict()
            return data

        try:
            data = get_match_object(pattern)
        except AttributeError:
            data = get_match_object(cls.string_to_regex(cls._oldver))
        return data

    @classmethod
    def parse_log(cls, logfile: str) -> pd.DataFrame:
        """
        parse a SEQC log into a pd.DataFrame column with a multi-index corresponding to
        the RUN SUMMARY section of the seqc.log object.

        :param logfile: str, path to log file
        :returns df: pd.DataFrame, dataframe containing log information
        """
        mo = seqc.log.LogData.match_log(logfile)
        return seqc.log.LogData.dictionary_to_dataframe(
            mo, logfile.split('/')[-1].replace('.log', ''))

    @classmethod
    def add_log(cls, filenames, df=None):
        if not filenames:
            return df
        if df is None:
            df = cls.parse_log()


    @classmethod
    def parse_multiple(cls, directory: str, exclude: str='') -> pd.DataFrame:
        """
        parse multiple SEQC logs into a pd.DataFrame object with a multi-index
        corresponding to the RUN SUMMARY section of the seqc.log object and a column
        multi-index corresponding to the directory hierarchy containing each log.

        This function takes a root directory and parses each log within the directory and
        all sub-directories. logs matching exclude pattern are omitted.

        :param directory: str, root directory to search for logs
        :param exclude: regex pattern to exclude log names
        :returns df: pd.DataFrame, dataframe containing log information
        """
        logs = []
        for path, subdirs, files in os.walk(directory):
            for name in files:
                filepath = os.path.join(path, name)
                if filepath.endswith('.log') and re.match(exclude, filepath) is None:
                    logs.append(filepath)

        frames = [cls.parse_log(f) for f in logs]

        # create column index
        cols = pd.MultiIndex.from_tuples(list(map(lambda p: tuple(p.split('/')), logs)))
        df = pd.concat(frames, 1)
        df.columns = cols
        return df




