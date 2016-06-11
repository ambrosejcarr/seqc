import json
import logging
from datetime import datetime


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
