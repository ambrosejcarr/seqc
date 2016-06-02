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


@contextmanager
def boto_errors(ident=None):
    """context manager that traps and retries boto functions
    to prevent random failures -- usually during batch runs
    :param ident: name of boto call"""

    try:
        yield
    except Exception:
        if ident:
            seqc.log.notify('Error in ' + ident + ', retrying in 5s...')
        else:
            seqc.log.notify('Error during boto call, retrying in 5s...')
        time.sleep(5)


def retry_boto_call(func, retries=4):
    """handles unexpected boto3 behavior, retries (default 3x)
    :param func: boto call to be wrapped
    :param retries: total # tries to re-call boto function"""

    @wraps(func)
    def wrapper(*args, **kwargs):
        numtries = retries
        while numtries > 1:
            with boto_errors(func.__name__):
                return func(*args, **kwargs)
            numtries -= 1
            if numtries == 1:
                raise BotoCallError('Unresolvable error in boto call, exiting.')
    return wrapper
