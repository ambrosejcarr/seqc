from seqc import log
import time
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


def retry_boto_call(
        func, retries=3, exceptions_to_catch=Exception,
        exception_to_raise: Exception = None, delay_retry=5, verbose=False):
    """improved wrapper to avoid issues with Boto3 by retrying an operation until success

    :param func: boto call to be wrapped
    :param retries: total # tries to re-call boto function
    :param exceptions_to_catch: an Exception or tuple of Exception that triggers a retry
    :param exception_to_raise: raise this after repeated failures
    :param delay_retry: waiting time before new attempt
    """

    @wraps(func)
    def wrapper(*args, **kwargs):
        num_tries = retries
        while num_tries > 0:
            try:
                return func(*args, **kwargs)
            except exceptions_to_catch as e:
                if verbose:
                    log.notify('Non fatal error in function ' +
                               func.__qualname__ +
                               ' --> ' +
                               str(e) +
                               ' | retrying in {nbr}s'.format(nbr=delay_retry))
                num_tries -= 1
                if num_tries == 0:
                    real_exception = (BotoCallError('Boto3 failed repeatedly.')
                                      if not exception_to_raise else exception_to_raise)
                    log.notify('SEQC encountered repeated failures while trying to '
                               'execute {method}. Giving up after {nbr_tries} with '
                               '{sec}s spacing.'.format(
                                   method=func.__qualname__, nbr_tries=retries,
                                   sec=delay_retry))
                    raise real_exception
                time.sleep(delay_retry)
    return wrapper
