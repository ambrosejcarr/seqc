import gzip
import bz2
from functools import lru_cache
import os


class Reader:
    """simple Reader Class, designed for inheritance across data types"""

    def __init__(self, files_):

        if isinstance(files_, list):
            self._files = files_
        elif isinstance(files_, str):
            self._files = [files_]

    @property
    def filenames(self):
        return self._files

    @lru_cache(maxsize=1)
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
        """return the colective size of all files being read in bytes"""
        return sum(os.stat(f).st_size for f in self._files)
