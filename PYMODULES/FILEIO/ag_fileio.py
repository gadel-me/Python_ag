import os
import fnmatch


class FileHandler(object):
    """
    Find files, link them, etc.
    """

    def __init__(self):
        self.files = []

    def find_files(self, directory, name_pattern):
        """
        Find files in folders by their ending.

        Paramters
        ---------
        directory : str
            name of directory where to find the files in

        name_pattern : str
            name of file ending with '*', e.g. "*.pwscf_out"

        Sources
        -------
        https://docs.python.org/2/library/fnmatch.html#fnmatch.filter

        """
        for filename in os.listdir(directory):
            filename = directory + filename
            if fnmatch.fnmatch(filename, name_pattern):
                self.files.append(filename)

    def link_files(self, directory):
        """Link files to the directory.

        Parameters
        ----------
        directory : str
            name of directory to link to

        """
        for filename in self.files:
            try:
                os.symlink(
                    filename, "{}/{}".format(directory, os.path.basename(filename))
                )
            except OSError:
                pass
