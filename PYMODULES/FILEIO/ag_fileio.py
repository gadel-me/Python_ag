import os
import pdb
import fnmatch
from pathlib import Path


class FileHandler(object):
    """
    Find files, link them, etc.
    """

    def __init__(self):
        self.files = []

    def find_files(self, path, name_pattern):
        """
        Find files in folders by their ending.

        Paramters
        ---------
        path : str
            name of directory where to find the files in

        name_pattern : str
            name of file ending with '*', e.g. "*.pwscf_out"

        Sources
        -------
        https://docs.python.org/2/library/fnmatch.html#fnmatch.filter

        """
        for filename in Path(path).glob('**/*{}'.format(name_pattern)):
            self.files.append(filename)

    def files_to_strings(self):
        """
        Convert PosixPath to str.
        """
        str_filenames = []

        for filename in self.files:
            full_path = filename.resolve()
            str_filename = full_path.as_posix()
            str_filenames.append(str_filename)

        self.files = str_filenames

    def link_files(self, directory):
        """Link files to the directory.

        Parameters
        ----------
        directory : str
            name of directory to link to

        """
        for filename in self.files:
            symlink = Path(directory + filename.name)

            try:
                symlink.symlink_to(filename.resolve())
            except OSError:
                pass
