import os
import numpy as np
from seqc import reader


class FastqRecord:
    """Fastq record object

    Defines several properties for accessing fastq record information:
    :property name: name field
    :property sequence: sequence field
    :property name2: second name field
    :property quality: quality field

    Also defines several methods for accessing SEQC annotation fields:
    :property annotations: list of annotations
    :property metadata: dictionary of read metadata (if any present)
    :property average_quality: return the mean quality of FastqRecord
    """

    __slots__ = ['_data']

    def __init__(self, record: [bytes, bytes, bytes, bytes]):
        self._data = list(record)

    @property
    def name(self) -> bytes:
        return self._data[0]

    @name.setter
    def name(self, value: bytes):
        self._data[0] = value

    @property
    def sequence(self) -> bytes:
        return self._data[1]

    @sequence.setter
    def sequence(self, value: bytes):
        self._data[1] = value

    @property
    def name2(self) -> bytes:
        return self._data[2]

    @name2.setter
    def name2(self, value: bytes):
        self._data[2] = value

    @property
    def quality(self) -> bytes:
        return self._data[3]

    @quality.setter
    def quality(self, value: bytes):
        self._data[3] = value

    def __bytes__(self) -> bytes:
        return b''.join(self._data)

    def __str__(self) -> str:
        return bytes(self).decode()

    def __len__(self) -> int:
        return len(self.sequence)

    @property
    def annotations(self) -> list:
        """
        returns:
        --------
        list of annotations present in the fastq header
        """
        try:
            end = self.name.index(b';')
            return self.name[:end].split(b':')
        except ValueError:
            return []

    @property
    def metadata(self) -> dict:
        """
        returns:
        --------
        dictionary of annotations and fields, if any are present"""
        try:
            start = self.name.rindex(b'|')
        except ValueError:
            return {}
        fields = {}
        for field in self.name[start + 1:].split(b':'):
            k, v = field.split(b'=')
            fields[k] = v
        return fields

    def add_annotation(self, values) -> None:
        """prepends a list of annotations to the name field of self.name
        :param values:
        """
        self._data[0] = b'@' + b':'.join(values) + b';' + self.name[1:]

    def add_metadata(self, values) -> None:
        """appends a list of metadata fields to the name field of self.name
        :param values:
        """
        self.name += b'|' + b':'.join(k + '=' + v for k, v in values.items())

    def average_quality(self) -> int:
        """"""
        return np.mean(np.frombuffer(self.quality, dtype=np.int8, count=len(self)))\
            .astype(int) - 33


class Reader(reader.Reader):
    """
    Fastq Reader, defines some special methods for reading and summarizing fastq data:

    :method __iter__: Iterator over fastq Record objects
    :method __len__: return number of records in file
    :method estimate_sequence_length: estimate the length of fastq sequences in file
    """

    @staticmethod
    def record_grouper(iterable):
        args = [iter(iterable)] * 4
        return zip(*args)

    def __iter__(self):
        for record in self.record_grouper(super().__iter__()):
            yield FastqRecord(record)

    def __len__(self):
        """
        return the length of the Reader object. This depends on the implementation of
        self.__iter__(); it does not necessarily represent the length of the file in
        lines.
        """
        return sum(1 for _ in self) / 4

    def estimate_sequence_length(self):
        """
        estimate the sequence length of a fastq file from the first 10000 records of
        the file.

        :return: int mean, float standard deviation, (np.ndarray: observed lengths,
          np.ndarray: counts per length)
        """
        i = 0
        records = iter(self)
        data = np.empty(10000, dtype=int)
        while i < 10000:
            try:
                seq = next(records).sequence
            except StopIteration:  # for fastq files shorter than 10000 records
                data = data[:i]
                break
            data[i] = len(seq) - 1  # last character is a newline
            i += 1
        return np.mean(data), np.std(data), np.unique(data, return_counts=True)
    
    
def merge_paired(merge_function, fout, genomic, barcode=None) -> (str, int):
    """
    General function to annotate genomic fastq with barcode information from reverse read.
    Takes a merge_function which indicates which kind of platform was used to generate
    the data, and specifies how the merging should be done.

    :param merge_function: function from merge_functions.py
    :param fout: merged output file name
    :param genomic: fastq containing genomic data
    :param barcode: fastq containing barcode data
    :return str fout, filename of merged fastq file

    """
    directory, filename = os.path.split(fout)
    if directory and not os.path.isdir(directory):
        os.makedirs(directory, exist_ok=True)
    genomic = Reader(genomic)
    if barcode:
        barcode = Reader(barcode)
        with open(fout, 'wb') as f:
            for g, b in zip(genomic, barcode):
                r = merge_function(g, b)
                f.write(bytes(r))
    else:
        with open(fout, 'wb') as f:
            for g in genomic:
                r = merge_function(g)
                f.write(bytes(r))

    return fout


def truncate(fastq_file, lengths):
    """

    :param str fastq_file: the input fastq file
    :param [int] lengths: a list of integer lengths to truncate the input fastq file
    :return:
    """
    # get sequence length of input file
    r = Reader(fastq_file)
    length = None
    for record in r:
        length = len(record.sequence)
        break

    print('sequence length in file is %d' % length)

    # remove any lengths longer than sequence length of file
    lengths = sorted([l for l in lengths if l < length])[::-1]  # largest to smallest

    # open a bunch of files
    files = []
    for l in lengths:
        name = fastq_file.replace('.gz', '').replace('.fastq', '') + '_%d_' % l + '.fastq'
        files.append(open(name, 'wb'))

    i = 0
    indices = list(range(len(lengths)))
    for record in r:
        if i > 10e6:
            break
        for j in indices:
            record.sequence = record.sequence[:-1][:lengths[j]] + b'\n'
            record.quality = record.quality[:-1][:lengths[j]] + b'\n'
            files[j].write(bytes(record))
        i += 1

    for f in files:
        f.close()
