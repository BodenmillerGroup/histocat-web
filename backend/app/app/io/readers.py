# TmLibrary - TissueMAPS library for distibuted image analysis routines.
# Copyright (C) 2016-2019 University of Zurich.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published
# by the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
import json
import logging
import os
import re
import sys
import traceback
from abc import ABC
from abc import abstractmethod

import cv2
import h5py
import lxml.etree
import numpy as np
import pandas as pd
import yaml

from app.core.utils import same_docstring_as

logger = logging.getLogger(__name__)


class Reader(ABC):
    '''Abstract base class for reading data from files.

    Readers make use of the
    `with statement context manager <https://docs.python.org/2/reference/datamodel.html#context-managers>`_.
    and thus follow a similar syntax::

        with Reader('/path/to/file') as f:
            data = f.read()
    '''

    def __init__(self, filename: str):
        '''
        Parameters
        ----------
        filename: str
            absolute path to a file

        Raises
        ------
        OSError
            when `filename` does not exist
        '''
        if not os.path.exists(filename):
            raise OSError('File does not exist: %s' % filename)
        self.filename = filename

    def __enter__(self):
        self._stream = open(self.filename, 'r')
        return self

    def __exit__(self, except_type, except_value, except_trace):
        self._stream.close()
        if except_value:
            sys.stdout.write(
                'The following error occurred while reading from file "%s":\n%s'
                % (self.filename, str(except_value))
            )
            for tb in traceback.format_tb(except_trace):
                sys.stdout.write(tb)
            sys.exit(1)

    @abstractmethod
    def read(self):
        pass


class TextReader(Reader):
    '''Class for reading data from text files.'''

    @same_docstring_as(Reader.__init__)
    def __init__(self, filename: str):
        super(TextReader, self).__init__(filename)

    def read(self):
        '''Reads data from text file.

        Returns
        -------
        lxml.etree._Element
            xml
        '''
        logger.debug('read text from file: %s', self.filename)
        return self._stream.read()


class XmlReader(Reader):
    # TODO: implement xpath subset reading via lxml
    '''Class for reading data from files in XML format.'''

    @same_docstring_as(Reader.__init__)
    def __init__(self, filename: str):
        super(XmlReader, self).__init__(filename)

    def read(self):
        '''Reads data from XML file.

        Returns
        -------
        lxml.etree._Element
            xml
        '''
        logger.debug('read XML data from file: %s', self.filename)
        return lxml.etree.fromstring(self._stream.read())


def load_json(string: str):
    '''
    Convert JSON string to Python object.

    Parameters
    ----------
    string: str
        JSON string

    Returns
    -------
    dict or list
    '''
    return json.loads(string)


class JsonReader(Reader):
    '''Class for reading data from files in JSON format.'''

    @same_docstring_as(Reader.__init__)
    def __init__(self, filename: str):
        super(JsonReader, self).__init__(filename)

    def read(self):
        '''Reads data from JSON file.

        Returns
        -------
        dict or list
            file content
        '''
        logger.debug('read JSON data from file: %s', self.filename)
        return load_json(self._stream.read())


class YamlReader(Reader):
    '''Class for reading data from files in YAML 1.2 format.'''

    @same_docstring_as(Reader.__init__)
    def __init__(self, filename: str):
        super(YamlReader, self).__init__(filename)

    def read(self):
        '''Reads YAML file.

        Returns
        -------
        dict or list
            file content
        '''
        logger.debug('read YAML data from file: %s', self.filename)
        return yaml.load(self._stream)


class DatatableReader(Reader):
    '''Class for reading data from a HDF5 file
    using the `pytables <http://www.pytables.org/>`_ library.
    '''

    @same_docstring_as(Reader.__init__)
    def __init__(self, filename: str):
        super(DatatableReader, self).__init__(filename)

    def __enter__(self):
        logger.debug('open file: %s', self.filename)
        self._stream = pd.HDFStore(self.filename, 'r')
        return self

    def exists(self, path):
        '''Check whether a `path` exists within the file.

        Parameters
        ----------
        path: str
            absolute path to a group or dataset in the file

        Returns
        -------
        bool
            ``True`` if `path` exists and ``False`` otherwise
        '''
        if path in self._stream:
            return True
        else:
            return False

    def read(self, path: str):
        '''Reads a dataset.

        Parameters
        ----------
        path: str
            absolute path to the dataset within the file

        Returns
        -------
        pandas.DataFrame
            dataset

        Raises
        ------
        IOError
            when `path` already exists
        '''
        return self._stream.select(path)


class DatasetReader(Reader):
    '''Class for reading data from a HDF5 file
    using the `h5py <http://docs.h5py.org/en/latest/>`_ library.
    '''

    @same_docstring_as(Reader.__init__)
    def __init__(self, filename):
        super(DatasetReader, self).__init__(filename)

    def __enter__(self):
        logger.debug('open file: %s', self.filename)
        # NOTE: The file shouldn't be opened in read-only mode, because this
        # would prevent concomitant writing
        self._stream = h5py.File(self.filename, 'r+')
        return self

    def exists(self, path):
        '''Checks whether `path` exists within the file.

        Parameters
        ----------
        path: str
            absolute path to a group or dataset in the file

        Returns
        -------
        bool
            ``True`` if `path` exists and ``False`` otherwise
        '''
        if path in self._stream:
            return True
        else:
            return False

    @staticmethod
    def _is_dataset(element):
        if isinstance(element.id, h5py.h5d.DatasetID):
            return True
        else:
            return False

    @staticmethod
    def _is_group(element):
        # TODO: this doesn't work, also lists datasets
        if isinstance(element.id, h5py.h5g.GroupID):
            return True
        else:
            return True

    def list_datasets(self, path='/', pattern='.*'):
        '''Lists datasets within a given group.

        Parameters
        ----------
        path: str, optional
            absolute path to a group in the file (default: ``"/"``)
        pattern: str, optional
            regular expression pattern to filter datasets (default: ``".*"``)

        Returns
        -------
        List[str]
            names of the datasets in `path`

        Raises
        ------
        KeyError
            when `path` does not exist
        '''
        try:
            group = self._stream[path]
        except KeyError:
            raise KeyError('Group does not exist: %s' % path)
        names = list()
        r = re.compile(pattern)
        for name, value in group.iteritems():
            if self._is_dataset(value) and r.search(name):
                names.append(name)
        return names

    def list_groups(self, path, pattern='.*'):
        '''Lists groups within a given group.

        Parameters
        ----------
        path: str
            absolute path to a group in the file
        pattern: str, optional
            regular expression pattern to filter groups (default: ``".*"``)

        Returns
        -------
        List[str]
            names of the groups in `path`

        Raises
        ------
        KeyError
            when `path` does not exist
        '''
        try:
            group = self._stream[path]
        except KeyError:
            raise KeyError('Group does not exist: %s' % path)
        names = list()
        r = re.compile(pattern)
        for name, value in group.iteritems():
            if not self._is_dataset(value) and r.search(name):
                names.append(name)
        return names

    def read(self, path):
        '''Reads a dataset.

        Parameters
        ----------
        path: str
            absolute path to the dataset within the file

        Returns
        -------
        numpy.ndarray
            dataset

        Raises
        ------
        KeyError
            when `path` does not exist
        '''
        try:
            dset = self._stream[path]
        except KeyError:
            raise KeyError('Dataset does not exist: %s' % path)
        return dset[()]

    def read_subset(self, path, index=None, row_index=None, column_index=None):
        '''Reads a subset of a dataset. For *fancy-indexing* see
        `h5py docs <http://docs.h5py.org/en/latest/high/dataset.html#fancy-indexing>`_.

        Parameters
        ----------
        path: str
            absolute path to the dataset within the file
        index: int or List[int], optional
            zero-based index
        row_index: int or List[int], optional
            zero-based row index
        column_index: int or List[int], optional
            zero-based column index

        Returns
        -------
        numpy.ndarray
            dataset

        Raises
        ------
        KeyError
            when `path` does not exist
        '''
        try:
            dset = self._stream[path]
        except KeyError:
            raise KeyError('Dataset does not exist: %s' % path)
        if row_index and not column_index:
            if len(dset.shape) == 1:
                raise IndexError(
                    'Dataset dimensions do not allow row-wise indexing'
                )
            return dset[row_index, :]
        elif not row_index and column_index:
            if len(dset.shape) == 1:
                raise IndexError(
                    'Dataset dimensions do not allow column-wise indexing'
                )
            return dset[:, column_index]
        elif row_index and column_index:
            if len(dset.shape) == 1:
                raise IndexError(
                    'Dataset dimensions do not allow row/column-wise indexing'
                )
            return dset[row_index, column_index]
        elif index is not None:
            return dset[index]
        else:
            raise ValueError('No index provided.')

    def get_attribute(self, path, name):
        '''Get an attribute attached to a dataset.

        Parameters
        ----------
        path: str
            absolute path to the dataset within the file
        name: str
            name of the attribute

        Returns
        -------
        ???

        Raises
        ------
        KeyError
            when `path` does not exist
        AttributeError
            when dataset does not have an attribute called `name`
        '''
        try:
            dset = self._stream[path]
        except KeyError:
            raise KeyError('Dataset does not exist: %s' % path)
        attribute = dset.attrs.get(name)
        if not attribute:
            raise AttributeError(
                'Dataset doesn\'t have an attribute "%s": %s'
                % (name, path))
        return attribute

    def get_dims(self, path):
        '''Get the dimensions of a dataset.

        Parameters
        ----------
        path: str
            absolute path to the dataset within the file

        Returns
        -------
        Tuple[int]
            number of rows and columns of the dataset

        Raises
        ------
        KeyError
            when `path` does not exist
        '''
        try:
            dims = self._stream[path].shape
        except KeyError:
            raise KeyError('Dataset does not exist: %s' % path)
        return dims

    def get_type(self, path):
        ''''Get the data type of a dataset.

        Parameters
        ----------
        path: str
            absolute path to the dataset within the file

        Returns
        -------
        type
            data type of the dataset

        Raises
        ------
        KeyError
            when `path` does not exist
        '''
        try:
            dtype = self._stream[path].dtype
        except KeyError:
            raise KeyError('Dataset does not exist: %s' % path)
        return dtype


class ImageReader(Reader):
    '''Class for reading pixel data from standard image file formats as
    :class:`numpy.ndarray` objects.
    '''

    @same_docstring_as(Reader.__init__)
    def __init__(self, filename: str):
        super(ImageReader, self).__init__(filename)

    def __enter__(self):
        return self

    def __exit__(self, except_type, except_value, except_trace):
        if except_value:
            sys.stdout.write(
                'The following error occurred while reading from file "%s":\n%s'
                % (self.filename, str(except_value))
            )
            for tb in traceback.format_tb(except_trace):
                sys.stdout.write(tb)
            sys.exit(1)

    def read(self, dtype: np.dtype = np.uint16):
        '''Reads pixels data from image file.

        Parameters
        ----------
        dtype: type, optional
            `numpy` data type (default: ``numpy.uint16``)

        Returns
        -------
        numpy.ndarray
            pixels data
        '''
        logger.debug('read from file: %s', self.filename)
        # NOTE: The string approach fails for some PNG formats. It can be solved
        #   arr = np.fromstring(self._stream.read(), dtype)
        #   return cv2.imdecode(arr, cv2.IMREAD_UNCHANGED)
        # It can be solved by opening images with PIL and converting them to
        # numpy arrays as follows:
        #   from PIL import Image
        #   content = Image.open(self.filename)
        #   arr = np.array(content)
        #   return cv2.imdecode(arr, cv2.IMREAD_UNCHANGED)
        # However, this is way slower than reading via OpenCV directly!
        return cv2.imread(self.filename, cv2.IMREAD_UNCHANGED)
