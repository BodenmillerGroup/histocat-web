import os
from abc import abstractproperty, abstractmethod

from sqlalchemy import Column, DateTime, Integer, String
from sqlalchemy.ext.declarative import declarative_base, declared_attr
from sqlalchemy.sql.functions import now

Base = declarative_base()


class DateMixIn:
    '''
    Mixin class to automatically add columns with datetime stamps to a database table.
    '''

    # NOTE: We use the "declared_attr" property for the mixin to ensure that
    # the columns are added to the end of the columns list. This simplifies
    # table distribution.

    @declared_attr
    def created_at(cls):
        '''datetime: date and time when the row was inserted into the column'''
        return Column(DateTime, default=now())

    @declared_attr
    def updated_at(cls):
        '''datetime: date and time when the row was last updated'''
        # TODO: CREATE TRIGGER to update independent of ORM
        return Column(
            DateTime,
            default=now(),
            onupdate=now()
        )


class IdMixIn:
    '''
    Mixin class to automatically add an ID column to a database table with primary key constraint.
    '''

    id = Column(Integer, primary_key=True, autoincrement=True)


class FileSystemModel(Base, IdMixIn):
    '''
    Abstract base class for model classes, which refer to data stored on disk outside of the database.
    '''

    __abstract__ = True

    _location = Column('location', String(4096))

    @abstractproperty
    def location(self) -> str:
        '''str: location on disk

        Devired classes must implement `location` and decorate it with
        :func:`sqlalchemy.ext.hybrid.hyprid_property`.
        '''
        pass

    @location.setter
    def location(self, value: str):
        self._location = value


class DirectoryModel(FileSystemModel):
    '''
    Abstract base class for model classes, which refer to data stored in directories on disk.
    '''

    __abstract__ = True


class FileModel(FileSystemModel):
    '''
    Abstract base class for model classes, which refer to data stored in files on disk.
    '''

    __abstract__ = True

    @property
    def format(self) -> str:
        '''str: file extension, e.g. ".tif" or ".jpg"'''
        return os.path.splitext(self.location)[1]

    @abstractmethod
    def get(self):
        '''Gets the file content.'''
        pass

    @abstractmethod
    def put(self, data):
        '''Puts `data` to the file.'''
        pass
