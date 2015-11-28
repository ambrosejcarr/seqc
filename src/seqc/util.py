__author__ = 'ambrose'


def check_type(arg, type_, message):
    """utility function to raise a TypeError if an object of incorrect type is given

    args:
    -----
    arg:  argument to type-check
    type_: the type to check against
    message: the error message, if the type is found to be incorrect

    returns:
    --------
    None
    """
    if not isinstance(arg, type_):
        if message:
            raise TypeError(message)
