__author__ = 'ambrose'

from multiprocessing import Process, Queue, Manager
from queue import Empty, Full
from time import sleep
import os
import seqc


def get(queue, pids=None):
    """
    iterator gets data from queue. If no data is on queue and no pid in pids is active,
    returns StopIteration
    """
    while True:
        try:
            index, data = queue.get_nowait()
            yield index, data
        except Empty:
            if pids:
                if any_alive(pids):
                    sleep(1)
                    continue
                else:
                    break
            else:
                break


def put(data, queue, pids=None):
    """
    Place data on queue. If pids is given, raises exception if all processes are dead
    """
    while True:
        try:
            queue.put_nowait(data)
            break
        except Full:
            if pids:
                if any_alive(pids):
                    sleep(1)
                    continue
                else:
                    raise RuntimeError('Attempt to place data on a queue after all '
                                       'processes to consume it have terminated')


def any_alive(pids):
    """Unix-only tool to check if any of the processes associated with pids are alive."""
    for id_ in pids:
        try:
            os.kill(id_, 0)
            return True  # at least one process alive; return True
        except ProcessLookupError:
            pass
    return False  # no process was alive; return False


def start_processes(n, target, args):
    proc = [Process(target=target, args=args) for _ in range(n)]
    for p in proc:
        p.start()
    pids = [p.pid for p in proc]
    return proc, pids


def join(processes):
    for p in processes:
        p.join()


def process_parallel(
        n_proc, h5_name, read_func, process_func, write_func, read_kwargs=None,
        process_kwargs=None, write_kwargs=None):
    """
    args:
    -----
    """

    def read(read_func_, process_queue_, kwargs, pids):
        """
        read chunks from file_iterator and places them on the processing queue.
        file iterator should take care of merging multiple files or cutting headers. It
        should also close files when iteration is complete

        Because these multithreaded processes are often i/o bound, it makes sense to
        do as little computation in this process as possible. This means that whenever
        possible, the read thread should simply take bytes chunks and dump them on the
        stack, allowing the processing threads to decompress or decode them into unicode
        strings
        """

        while not all([pids['read'], pids['process'], pids['write']]):
            print('fp view of pids: %s' % repr(pids))
            sleep(1)

        for i, chunk in enumerate(read_func_(**kwargs)):
            put((i, chunk), process_queue_, pids['process'])
            seqc.log.info('Read chunk %d.' % i)

    def process(process_func_, process_queue_, write_queue_, kwargs, pids):
        """
        Applies processing func to each chunk from the read iterator, terminating when the
        Queue is empty and the read process is dead.
        """
        while not all([pids['read'], pids['process'], pids['write']]):
            sleep(1)

        for i, data in get(process_queue_, pids['read']):
            seqc.log.info('Processing chunk %d.' % i)
            processed_data = process_func_(data, **kwargs)
            put((i, processed_data), write_queue_, pids['write'])

    def write(write_func_, write_queue_, filename, kwargs, pids):
        """
        Calls writing_func on each chunk of data retrieved from writing_queue. Please note
        that the file will likely be opened and closed, so writing func should append, not
        write.
        """
        while not all([pids['read'], pids['process'], pids['write']]):
            sleep(1)

        # open file
        h5writer = write_func_(filename)
        h5writer.create(**kwargs)
        for i, data in get(write_queue_, pids['process']):
            seqc.log.info('Writing chunk %d.' % i)
            h5writer.write(data)
        h5writer.close()

    seqc.log.setup_logger()

    # set up shared dictionary to track whether processes have started
    manager = Manager()
    pids = manager.dict()
    pids['read'] = None
    pids['write'] = None
    pids['process'] = None

    # define the number of processors
    n_read = 1
    n_write = 1
    n_process = max(n_proc - 2, 1)

    # set kwargs to empty dict if any function is not passed them.
    if not read_kwargs:
        read_kwargs = {}
    if not process_kwargs:
        process_kwargs = {}
    if not write_kwargs:
        write_kwargs = {}

    # create a read process and the process_queue it will put data onto.
    process_queue = Queue(maxsize=n_process)  # limit n chunks to # processes
    read_process, read_pids = start_processes(
        n_read, read, ([read_func, process_queue, read_kwargs, pids]))
    pids['read'] = read_pids

    # create processing processes and the write queue they will put data onto.
    write_queue = Queue(maxsize=n_process)  # limit n chunks to # processes
    processors, process_pids = start_processes(
        n_process, process,
        ([process_func, process_queue, write_queue, process_kwargs, pids]))
    pids['process'] = process_pids

    # create the write process
    write_process, write_pids = start_processes(
        n_write, write, ([write_func, write_queue, h5_name, write_kwargs, pids]))
    pids['write'] = write_pids

    # wait for each process to finish
    join(read_process)
    join(processors)
    join(write_process)
