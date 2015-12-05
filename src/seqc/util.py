__author__ = 'ambrose'

import cProfile
import pstats
from memory_profiler import memory_usage
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
        stats.sort_stats('cumtime')
        if filename:
            stats.dump_stats(os.path.expanduser(filename))
        else:
            stats.print_stats(10)

        # return the normal result of the function so this can be chained
        return res

    return time_profiled_function


def mem_profile(func, filename=None):
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


# some old profiling code
# def memory_usage(self, sizes, directory, index):
#
#     # generate fastq files
#     for s in sizes:
#         fastq_filename = directory + str(s)
#         seqc.fastq.GenerateFastq.in_drop(s, fastq_filename)
#     forward = [directory + str(s) + '_r1.fastq' for s in sizes]
#     reverse = [directory + str(s) + '_r2.fastq' for s in sizes]
#
#     # get cell barcodes
#     with open(self.data_dir + 'in_drop/barcodes/cb_3bit.p', 'rb') as f:
#         cb = pickle.load(f)
#
#     # merge the data
#     merged_files = []
#     for f, r in zip(forward, reverse):
#         fq, _ = seqc.fastq.merge_fastq([f], [r], 'in-drop', directory, cb)
#         merged_files.append(fq)
#
#     # align the data
#     samfiles = seqc.align.STAR.align_multiple_files(
#         merged_files, index, 7, directory
#     )
#
#     ft, fp = seqc.convert_features.construct_feature_table(self.gtf, 1000)
#
#     logfile = open(directory + 'report.txt', 'w+')
#
#     h5files = [directory + '%dm.h5' % s for s in sizes]
#
#     # define and run functions
#     def test_memory_usage(idx):
#         seqc.memory_usage((seqc.arrays.ReadArray.from_samfile, (samfiles[idx], ft, fp)))
#         # arr.save_h5(h5files[index])
#
#     for i in range(len(samfiles)):
#         test_memory_usage(i)
#
#     logfile.close()
#
# @unittest.skip('')
# def test_ra_memory_usage_small(self):
#     index = self.data_dir + 'genome/mm38_chr19/'
#     working_directory = self.data_dir + 'test_ra_memory_usage/'
#     if not os.path.isdir(working_directory):
#         os.makedirs(working_directory)
#     self.memory_usage([int(1e6)], working_directory, index)
#
# @unittest.skip('')
# def test_profile_mem_usage(self):
#     samfile = ('/Users/ambrose/PycharmProjects/SEQC/src/data/test_ra_memory_usage/'
#                'merged_temp/Aligned.out.sam')
#     ft, fp = seqc.convert_features.construct_feature_table(self.gtf, 1000)
#     usage = seqc.memory_usage((seqc.arrays.ReadArray.from_samfile, (samfile, ft, fp)))
#     print(np.array(usage))
#
# @unittest.skip('')
# def test_ra_memory_usage(self):
#     data_dir = '/'.join(seqc.__file__.split('/')[:-2]) + '/data/'
#
#     # generate fastq data for testing
#     files = [data_dir + 'test_memory_usage/' + prefix for prefix in
#              ['1m', '2m', '4m', '8m', '16m']]
#     index = data_dir + 'genome/mm38_chr19/'
#     seqc.generate_in_drop_fastq_data(int(1e6), files[0])
#     seqc.generate_in_drop_fastq_data(int(2e6), files[1])
#     seqc.generate_in_drop_fastq_data(int(4e6), files[2])
#     seqc.generate_in_drop_fastq_data(int(8e6), files[3])
#     seqc.generate_in_drop_fastq_data(int(16e6), files[4])
#
#     # align the data, generating sam files
#     samfiles = seqc.align.STAR.align_multiple_files(
#         files, index, 7, data_dir + 'test_memory_usage/')
#
#     ft, fp = seqc.convert_features.construct_feature_table(self.gtf, 1000)
#
#     logfile = open(data_dir + 'test_memory_usage/report.txt', 'w+')
#
#     h5files = [f + '.h5' for f in files]
#
#     @profile(stream=logfile)
#     def test_1m():
#         arr = seqc.arrays.ReadArray.from_samfile(samfiles[0], ft, fp)
#         arr.save_h5(h5files[0])
#
#     @profile(stream=logfile)
#     def test_2m():
#         arr = seqc.arrays.ReadArray.from_samfile(samfiles[1], ft, fp)
#         arr.save_h5(h5files[1])
#
#     @profile(stream=logfile)
#     def test_4m():
#         arr = seqc.arrays.ReadArray.from_samfile(samfiles[2], ft, fp)
#         arr.save_h5(h5files[2])
#
#     @profile(stream=logfile)
#     def test_8m():
#         arr = seqc.arrays.ReadArray.from_samfile(samfiles[3], ft, fp)
#         arr.save_h5(h5files[3])
#
#     @profile(stream=logfile)
#     def test_16m():
#         arr = seqc.arrays.ReadArray.from_samfile(samfiles[4], ft, fp)
#         arr.save_h5(h5files[4])
#
#     for f in [test_1m, test_2m, test_4m, test_8m, test_16m]:
#         f()
#
#     logfile.close()
