import seqc
import time
from contextlib import contextmanager
from functools import wraps


class EC2RuntimeError(Exception):
    pass


class VolumeCreationError(Exception):
    pass


class SpotBidError(Exception):
    pass


class BotoCallError(Exception):
    pass


class ConfigurationError(Exception):
    pass


class NoMatchError(Exception):
    pass


class ArgumentParserError(Exception):
    pass


class SparseMatrixError(Exception):
    pass

class EmptyMatrixError(Exception):
    pass



@contextmanager
def boto_errors(ident=None, time_retry=5):
    """context manager that traps and retries boto functions
    to prevent random failures -- usually during batch runs
    :param ident: name of boto call
    :param time_retry: number of second before trying again
    """

    try:
        yield
    except Exception:
        retry_str = ', retrying in {delay}s....'.format(delay=time_retry)
        if ident:
            seqc.log.notify('Error in ' + ident + retry_str)
        else:
            seqc.log.notify('Error during boto call' + retry_str)
        time.sleep(time_retry)


def retry_boto_call(func, retries=3, exception_to_raise: Exception = None, delay=5):
    """handles unexpected boto3 behavior, retries (default 3x)
    :param func: boto call to be wrapped
    :param retries: total # tries to re-call boto function
    :param exception_to_raise: raise this after repeated failure, default is BotoCallError
    :param delay: waiting time before new attempt
    """

    @wraps(func)
    def wrapper(*args, **kwargs):
        num_tries = retries
        while num_tries > 0:
            with boto_errors(func.__name__, time_retry=delay):
                return func(*args, **kwargs)
            num_tries -= 1
            if num_tries == 0:
                raise exception_to_raise if exception_to_raise else BotoCallError('Boto3 failed %i times in a row'
                                                                                  % retries)
    return wrapper
