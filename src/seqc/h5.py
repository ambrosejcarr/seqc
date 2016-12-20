import pandas as pd
import numpy as np
import os


class H5:

    def __init__(self, archive_name: str):
        """Wrapper for the pandas HDFStore class which ensures that all interactions with
        the archive result in a closed, flushed archive.

        In order to ensure data usability, all data must be submitted in DataFrame format.
        This decision was made to encourage users to pair metadata with sequencing data,
        and reduce the incidence of unexpected data permutation.

        :param archive_name: name of the h5 archive to open. If the archive does not exist
          it will be created using a blosc5 filter

        :method ls: list contents of the archive
        :method save: save an object to the h5 archive
        :method load: load an object from the archive
        :method remove: remove a DataFrame from the archive
        :method is_open: returns True if the h5 archive is open, else False
        """
        if os.path.isfile(archive_name):
            self._archive = pd.HDFStore(archive_name, mode='a')
            self._archive.close()
        else:
            self._archive = pd.HDFStore(
                archive_name, mode='a', complib='blosc', complevel=5)
            self._archive.close()

    def __repr__(self):
        self._archive.open()
        try:
            return repr(self._archive)
        finally:
            self._archive.close()

    def save(self, data: pd.DataFrame, location: str) -> None:
        """Save DataFrame data to the h5 archive in location.

        :param data: DataFrame object to store
        :param location: filepath to save the object in the h5 hierarchy
        """
        if not isinstance(data, pd.DataFrame):
            if isinstance(data, np.ndarray):
                res = input('np.ndarray class detected. Save as pd.DataFrame with '
                            'ascending integer indices? [y/n] ')
                if res in ['y', 'yes', 'Y', 'YES', 'True', 'true', '1']:
                    data = pd.DataFrame(data)
                else:
                    print('User elected not to save DataFrame, archive is unmodified.')
                    return
            else:
                raise TypeError('only pd.DataFrame objects can be saved using this '
                                'class. To save np.ndarray objects please see the tables '
                                'package.')
        self._archive.open()
        try:
            self._archive[location] = data
        finally:
            self._archive.close()

    def load(self, location: str) -> None:
        """Load and return the dataframe found at location in the archive.

        :param location: str, location of object to retrieve from h5
        :return: pd.DataFrame, object found at location
        """
        self._archive.open()
        try:
            return self._archive[location]
        finally:
            self._archive.close()

    def ls(self) -> None:
        """list archive contents"""
        try:
            self._archive.open()
            print(self._archive)
        finally:
            self._archive.close()

    def remove(self, location: str) -> None:
        """remove the DataFrame at location from the archive

        Note: removing a dataframe at a branch node will remove all leaves sharing this
        prefix. e.g. in an archive containing:

        /data
        /data/filtered
        /data/metadata
        /new_data/data

        removing /data would remove the first three DataFrame objects from the archive.

        :param location: location of DataFrame to remove
        :return: None
        """

        self._archive.open()
        try:
            if location not in self._archive.keys():
                raise ValueError(
                    '{} not contained in archive, nothing to remove.'.format(location))
            else:
                removed = [k for k in self._archive.keys()
                           if k.startswith(location + '/')]
                if len(removed) != 0:
                    res = input(
                        'Removing branch node {}, which is a prefix for {!a} will remove '
                        'all listed DataFrames. Continue with removal? [y/n] '.format(
                            location, removed))
                    if res not in ['y', 'yes', 'Y', 'YES', 'True', 'true', '1']:
                        print('returned without deletion.')
                        return
                self._archive.remove(location)
        finally:
            self._archive.close()

    @property
    def is_open(self) -> bool:
        return self._archive.is_open
