import os
import gzip
import bz2


class Reader:
    """
    Basic reader object that seamlessly loops over multiple input files

    Can be subclassed to create readers for specific file types (fastq, gtf, etc.)
    """

    def __init__(self, files_):

        if isinstance(files_, list):
            self._files = files_
        elif isinstance(files_, str):
            self._files = [files_]

    @property
    def filenames(self):
        return self._files

    def __len__(self):
        """
        return the length of the Reader object. This depends on the implementation of
        self.__iter__(); it does not necessarily represent the length of the file in
        lines.
        """
        return sum(1 for _ in self)

    def __iter__(self):
        for f in self._files:
            if f.endswith('.gz'):
                file_input = gzip.open(f, 'rb')
            elif f.endswith('.bz2'):
                file_input = bz2.open(f, 'rb')
            else:
                file_input = open(f, 'rb')
            for record in file_input:
                yield record
            file_input.close()

    @property
    def size(self) -> int:
        """return the collective size of all files being read in bytes"""
        return sum(os.stat(f).st_size for f in self._files)
