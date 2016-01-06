import seqc
import gzip
import bz2


class Reader:

    def __init__(self, fasta):

        seqc.util.check_type(fasta, str, 'fasta')
        seqc.util.check_file(fasta, 'fasta')

        self.fasta = fasta
        try:
            fasta_iterator = iter(self)
            next(fasta_iterator)
        except:
            print('fasta is improperly formatted. Please examine file')
            raise

    def _open(self):
        """
        seamlessly open self._fasta, whether gzipped or uncompressed

        returns:
        --------
        fobj: open file object
        """
        if self.fasta.endswith('.gz'):
            fobj = gzip.open(self.fasta, 'rb')
        elif self.fasta.endswith('.bz2'):
            fobj = bz2.open(self.fasta, 'rb')
        else:
            fobj = open(self.fasta, 'rb')
        return fobj

    def __iter__(self):
        """return an iterator over all non-header records in fasta"""
        fobj = self._open()
        records = iter(fobj)

        try:
            # get rid of headers, get first header name
            line = next(records)
            while not line.startswith(b'>'):
                line = next(records)

            # set first header, get first sequence line
            name, *description = line[1:].strip().split()
            line = next(records)
            seq = line

            # yield records
            while line:
                if line.startswith(b'>'):
                    yield name, b' '.join(description), seq
                    name, *description = line[1:].strip().split()
                    seq = b''
                else:
                    seq += line[:-1]  # strip newline

                # grab the next line; catch StopIteration in case no terminal \n
                try:
                    line = next(records)
                except StopIteration:
                    break

            # yield final record
            yield name, b' '.join(description), seq

        finally:
            fobj.close()


class Fasta:

    def __init__(self, fasta_dictionary):
        self._chromosomes = fasta_dictionary

    @classmethod
    def from_file(cls, fasta_file):
        fasta_dictionary = {}
        rd = Reader(fasta_file)
        for name, _, sequence in rd:
            fasta_dictionary[name] = sequence
        return cls(fasta_dictionary)

    @property
    def chromosomes(self):
        return self._chromosomes

    def __getitem__(self, item):
        try:
            return self._chromosomes[item]
        except KeyError:
            try:  # user may have passed a string instead of a bytesobject
                return self._chromosomes[item.decode()]
            except KeyError:
                pass  # return the original error
            raise
