import pandas as pd
import os


class H5:

    def __init__(self, archive_name: str):
        """Wrapper for the pandas HDFStore class which ensures that all interactions with
        the archive result in a closed, flushed archive.

        :param archive_name: name of the h5 archive to open. If the archive does not exist
          it will be created using a blosc5 filter

        :method ls: list contents of the archive
        :method save: save an object to the h5 archive
        :method load: load an object from the archive
        """
        if os.path.isfile(archive_name):
            self._archive = pd.HDFStore(archive_name, mode='a')
            self._archive.close()
        else:
            self._archive = pd.HDFStore(
                archive_name, mode='a', complib='blosc', complevel=5)
            self._archive.close()

    def save(self, data: pd.DataFrame, location: str):
        """Save DataFrame data to the h5 archive in location.

        :param data: DataFrame object to store
        :param location: filepath to save the object in the h5 hierarchy
        """
        self._archive.open()  # must be appendable or will raise errors
        try:
            self._archive[location] = data
            self._archive.close()
        except:
            self._archive.close()
            raise

    def load(self, location: str):
        """Load and return the dataframe found at location in the archive.

        :param location: str, location of object to retrieve from h5
        :return: pd.DataFrame, object found at location
        """
        self._archive.open()
        try:
            data = self._archive[location]
            self._archive.close()
            return data
        except:
            self._archive.close()
            raise

    def ls(self):
        """list archive contents"""
        try:
            self._archive.open()
            print(self._archive)
            self._archive.close()
        except:
            self._archive.close()
            raise

    @property
    def is_open(self):
        return self._archive.is_open