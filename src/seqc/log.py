import json
import logging
from datetime import datetime
import pandas as pd


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


def args(arguments):
    """
    log namespace object from argument parser to file.

    :param arguments: namespace object, output of ArgumentParser.parse_args()
    :return: None
    """
    arguments = vars(arguments)
    info('Passed command line arguments: {}'.format(
            json.dumps(arguments, separators=(',', ': '), indent=4)))


def extract_data(logfile: str) -> pd.Series:
    """
    extract data in log file into pd.Series object for parsing and comparison

    copy the output pattern from seqc.stats.ExperimentalYield and reverse-engineer a
    regex to identify the critical parameters. Put those parameters in the DataFrame

    :param logfile: str, name of a log file
    :returns data: pd.Series object containing logged information
    """

    pattern = (
        '{divide}\nINPUT\n{divide}\n'
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
        'Missing cell barcodes/RMT:\t{no_cell} ({bad_cell}%)\n'
        # 'Missing RMTs:\t\t{no_rmt} ({bad_rmt}%)\n'
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
        '{divide}\n')

    with open(logfile, 'r') as f:
        data = f.read()


