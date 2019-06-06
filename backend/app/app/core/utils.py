import logging
import os
import time
from enum import Enum, unique
from shutil import rmtree
from typing import Tuple

import cv2
import numpy as np
import sqlalchemy

logger = logging.getLogger(__name__)


def create_directory(location):
    """Creates a directory on disk in a safe way.

    Parameters
    ----------
    location: str
        absolute path to the directory that should be created
    """
    try:
        os.makedirs(location)
    except OSError as err:
        if err.errno != 17:
            raise


class autocreate_directory_property(object):
    """Decorator class that acts like a property.
    The value represents a path to a directory on disk. The directory is
    automatically created when it doesn't exist. Once created, the value
    is cached, so that there is no reattempt to create the directory.

    Raises
    ------
    TypeError
        when the value of the property doesn't have type basestring
    ValueError
        when the value of the property is empty

    Examples
    --------
    .. code:: python

        from tmlib.utils import autocreate_directory_property

        class Foo(object):

            @autocreate_directory_property
            def my_new_directory(self):
                return '/tmp/blabla'

        foo = Foo()
        foo.my_new_directory
    """

    def __init__(self, func):
        self.__doc__ = func.__doc__
        self.func = func

    def __get__(self, obj, cls):
        if obj is None:
            return self
        value = obj.__dict__[self.func.__name__] = self.func(obj)
        if not isinstance(value, str):
            raise TypeError(
                'Value of property "%s" must have type string: %s'
                % (self.func.__name__, value)
            )
        if not value:
            raise ValueError(
                'Value of property "%s" cannot be empty.' % self.func.__name__
            )
        if not os.path.exists(value):
            logger.debug("create directory: %s", value)
            create_directory(value)
        return value


def delete_location(path: str):
    """Deletes a location on disk.

    Parameters
    ----------
    path: str
        absolute path to directory or file
    """
    if os.path.exists(path):
        logger.debug("remove location: %s", path)
        if os.path.isdir(path):
            rmtree(path)
        elif os.path.isfile(path):
            os.remove(path)


def remove_location_upon_delete(cls):
    """Decorator function for an database model class that
    automatically removes the `location` that represents an instance of the
    class on the filesystem once the corresponding row is deleted from the
    database table.

    Parameters
    ----------
    cls: tmlib.models.base.DeclarativeABCMeta
       implemenation of :class:`tmlib.models.base.FileSystemModel`

    Raises
    ------
    AttributeError
        when decorated class doesn't have a "location" attribute
    """

    def after_delete_callback(mapper, connection, target):
        delete_location(target.location)

    sqlalchemy.event.listen(cls, "after_delete", after_delete_callback)
    return cls


@unique
class FileUploadStatus(Enum):
    """Upload status of a file."""

    #: The file is registered, but upload not yet started
    WAITING = "WAITING"

    #: Upload is ongoing
    UPLOADING = "UPLOADING"

    #: Upload is complete
    COMPLETE = "COMPLETE"

    #: Upload has failed
    FAILED = "FAILED"


@unique
class Color(Enum):
    r = (0, 0, 1)  # red
    g = (0, 1, 0)  # green
    b = (1, 0, 0)  # blue
    y = (0, 1, 1)  # yellow
    c = (1, 1, 0)  # cyan
    m = (1, 0, 1)  # magenta


def colorize(image: np.ndarray, color: Color):
    image = cv2.cvtColor(image, cv2.COLOR_GRAY2RGB)
    image = image * color.value
    return image


def scale_image(image: np.ndarray, scale: float, levels: Tuple[float, float]):
    minL, maxL = levels
    if maxL <= minL:
        return image
    result = image * (scale / (maxL - minL))
    return result


def timeit(method):
    def timed(*args, **kw):
        ts = time.time()
        result = method(*args, **kw)
        te = time.time()
        if "log_time" in kw:
            name = kw.get("log_name", method.__name__.upper())
            kw["log_time"][name] = int((te - ts) * 1000)
        else:
            print("%r  %2.2f ms" % (method.__name__, (te - ts) * 1000))
        return result

    return timed
