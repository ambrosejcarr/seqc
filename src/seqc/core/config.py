import os
import configparser
from seqc.exceptions import ConfigurationError


def read_config(config_file: str=None) -> dict:
    """
    :param config_file: str, location of seqc configuration file
    :returns config: dict, dictionary of arguments and values in the configuration file
    """
    # todo add additional checks to make sure correct parameters are filled in!

    if config_file is None:
        config_file = os.path.expanduser('~/.seqc/config')
    else:
        config_file = os.path.expanduser(config_file)
    config = configparser.ConfigParser()
    if not config.read(config_file):
        raise ConfigurationError('Please run ./configure (found in the seqc '
                                 'directory) before attempting to run '
                                 'SEQC.py.')
    return config
