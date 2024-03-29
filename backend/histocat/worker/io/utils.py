import fnmatch
import os
import shutil
from pathlib import Path

IMAGE_CSV_FILENAME = "Image.csv"
CELL_CSV_FILENAME = "cell.csv"


def locate(root: str, pattern: str):
    """
    Locate all files matching supplied filename pattern in and below supplied root directory.
    """
    for path, dirs, files in os.walk(os.path.abspath(root)):
        for filename in fnmatch.filter(files, pattern):
            yield os.path.join(path, filename)


def copy_dir(src: Path, dst: str):
    return shutil.copytree(src, dst)


def copy_file(src: str, dst: str):
    return shutil.copy2(src, dst)
