__author__ = 'ambrose'

import cProfile
import pstats
from memory_profiler import memory_usage
from memory_profiler import profile as mem_profile
import os


def check_type(arg, type_, arg_name):
    """utility function to raise a TypeError if an object of incorrect type is given

    args:
    -----
    arg:  argument to type-check
    type_: the type to check against
    arg_name: the name of the argument being checked

    returns:
    --------
    None
    """
    if not isinstance(arg, type_):
        raise TypeError('argument "%s" must be %s, not %s.' %
                        (arg_name, repr(type_), repr(type(arg))))



def check_dir(directory: str, arg_name: str) -> None:
    """utility function to raise a FileNotFoundError if directory is not a directory

    args:
    -----
    directory: string ostensibly pointing to a directory
    arg_name: the name of the argument that is being checked

    returns:
    --------
    None
    """
    if not os.path.isdir(directory):
        raise FileNotFoundError('argument "%s": %s, is not a directory.' %
                                arg_name, directory)


def check_file(filename: str, arg_name: str) -> None:
    """utility function to raise a FileNotFoundError if directory is not a directory

    args:
    -----
    filename: string ostensibly pointing to a file
    arg_name: the name of the argument that is being checked

    returns:
    --------
    None
    """
    if not os.path.isfile(filename):
        raise FileNotFoundError('argument "%s": %s, is not a file.' %
                                (arg_name, filename))


def time_profile(func, filename=None):
    """Decoratable time profiling function

    usage:
    ------
    >>> @time_profile(filename='time_usage.log')
    >>> def simple_sum(x, y):
    >>>    return x + y

    >>> simple_sum(1, 7)

    args:
    -----
    filename (default None): write time profiling results to this file. If None, results
      are printed to stdout

    returns:
    --------
    func, decorated by time profiling
    """

    def time_profiled_function(*args, **kwargs):

        # wrap the original function with a profiler, then run it
        profile_ = cProfile.Profile()
        profile_.enable()
        res = func(*args, **kwargs)
        profile_.disable()

        # print or dump the statistics
        stats = pstats.Stats(profile_)
        stats.strip_dirs()
        stats.sort_stats('tottime')
        if filename:
            stats.dump_stats(os.path.expanduser(filename))
        else:
            stats.print_stats(10)

        # return the normal result of the function so this can be chained
        return res

    return time_profiled_function


def mem_profile_time(func, filename=None):
    """Decoratable memory profiling function

    Note that this will slow execution because the code must be executed twice in order
    to return the normal return value due to limitations of memory_profiler

    usage:
    ------
    >>> @mem_profile(filename='memory_usage.log')
    >>> def simple_sum(x, y):
    >>>    return x + y

    >>> simple_sum(1, 7)

    args:
    -----
    filename (default None): write memory profiling results to this file. If None, results
      are printed to stdout

    returns:
    --------
    func, decorated by memory profiling
    """

    def memory_profiled_function(*args, **kwargs):

        mem_results = memory_usage((func, args, kwargs))
        if filename:
            with open(filename, 'w') as f:
                f.write(mem_results)
        else:
            print(mem_results)
        res = func(*args, **kwargs)

        return res

    return memory_profiled_function
