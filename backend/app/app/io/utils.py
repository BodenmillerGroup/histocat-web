import fnmatch
import os

SCHEMA_XML_ENDING = "_schema.xml"


def locate(root: str, pattern: str) -> str:
    """
    Locate all files matching supplied filename pattern in and below supplied root directory.
    """
    for path, dirs, files in os.walk(os.path.abspath(root)):
        for filename in fnmatch.filter(files, pattern):
            yield os.path.join(path, filename)
