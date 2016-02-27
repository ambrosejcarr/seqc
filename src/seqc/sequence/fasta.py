import fileinput
from collections.abc import Iterator, Iterable
import seqc


class FastaRecord:

    __slots__ = ['_sequence', '_name', '_description']

    def __init__(self, record: list):
        self._name, *description = record[0][1:-1].split()  # remove '>' and '\n'
        self._description = b' '.join(description)
        self._sequence = b''.join(l.strip() for l in record[1:])

    @property
    def sequence(self) -> bytes:
        return self._sequence

    @property
    def name(self) -> bytes:
        return self._name

    @property
    def description(self) -> bytes:
        return self._description


class Genome:

    __slots__ = ['_chromosomes']

    def __init__(self, fasta_records: Iterable):
        chromosomes = {}
        for record in fasta_records:
            chromosomes[record.name] = record
        self._chromosomes = chromosomes

    @classmethod
    def from_file(cls, fasta_file):
        return cls(FastaReader(fasta_file))

    @property
    def chromosomes(self) -> dict:
        return self._chromosomes

    def __getitem__(self, item):
        return self.chromosomes[item]

    def __setitem__(self, key, value):
        self._chromosomes[key] = value

    def __delitem__(self, key):
        del self._chromosomes[key]

    def keys(self) -> Iterator:
        return self.chromosomes.keys()

    def values(self) -> Iterator:
        return self.chromosomes.values()

    def items(self) -> Iterator:
        return self.chromosomes.items()


class FastaReader(seqc.reader.Reader):

    def __iter__(self) -> Iterator:
        hook = fileinput.hook_compressed
        with fileinput.input(self._files, openhook=hook, mode='rb') as f:

            # get rid of file headers, get first record name
            record_iterator = iter(self)
            line = next(record_iterator)
            while not line.startswith(b'>'):
                line = next(record_iterator)

            record = [line]
            while True:
                try:
                    line = next(record_iterator)
                except StopIteration:
                    yield record
                    return

                if line.startswith(b'>'):
                    yield record
                    record = [line]
                else:
                    record.append(line)