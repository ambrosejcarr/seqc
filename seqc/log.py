import json
import logging
from datetime import datetime


def setup_logger():
    """create a simple log file in the cwd to track progress and any errors"""
    logging.basicConfig(filename='seqc.log', level=logging.DEBUG, filemode='w')


def info(message):
    """print a timestamped update for the user.
    :param message:
    """
    logging.info(datetime.now().strftime("%Y-%m-%d %H:%M:%S") + ':' + message)


def exception():
    """log the most recent exception to an initialized logger"""
    logging.exception(datetime.now().strftime("%Y-%m-%d %H:%M:%S") + ':main:')


def notify(message):
    """print a timestamped update for the user and log it to file
    :param message:
    """
    info(message)
    print('SEQC: ' + datetime.now().strftime("%Y-%m-%d %H:%M:%S") + ': %s' % message)


def args(arguments):
    """
    :param arguments:
    :return:
    """
    arguments = vars(arguments)
    info('Passed command line arguments: {}'.format(
            json.dumps(arguments, separators=(',', ': '), indent=4)))
