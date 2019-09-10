from __future__ import print_function, division
import collections

__version__ = "2017-05-05"


class LogUniverse(object):
    """
    Contains basic methods for processing log-outputs.
    """
    def __init__(self):
        """
        data = [{ikey: [value1, value2, ...], jkey: [value1, value2, ...], ...}]
        keys   str; keywords such as 'Step', 'Temp'
        value: list/tuple of floats/ints;
        """
        self.data = []

    def sort_by_key(self):
        """
        Sort data by key.
        Sources:    http://stackoverflow.com/questions/9001509/how-can-i-sort-a-dictionary-by-key
        """
        #TODO CHECK IF DICTIONARY WAS ALTERED, HAS NOT BEEN TESTED YET!
        for cdata in self.data:
            for subdata in cdata:
                cdata = collections.OrderedDict(sorted(cdata.items()))

    def cut_data(self, start, stop, keyword=None):
        """
        Shrink all data-values by start and stop values
        keyword:    str; should be ascending value, e.g. 'Step', 'Temp'
        start:      float/int; first value in list of keyword, e.g. 2500 (Step)
        stop:       float/int; last value in list of keyword, e.g. 350.6 (K)
        """
        #TODO currently not working properly
        if keyword is not None:
            for cdata in self.data:

                for ikey in cdata:
                    if ikey == keyword:
                        start_ptr = cdata[ikey].index(start)
                        stop_ptr  = cdata[ikey].index(stop)
                    break

                for jkey in cdata:
                    cdata[jkey] = cdata[jkey][start_ptr:stop_ptr]
        else:
            for cdata in self.data:
                for jkey in cdata:
                    cdata[jkey] = cdata[jkey][start:stop]

    def _int_values(self):
        """
        Convert values of keys Atoms, Bonds, Angles and Step to ints.
        """
        for cdata in self.data:
            for ckey in cdata:
                if ckey in ("Atoms", "Bonds", "Angles", "Step"):
                    cdata[ckey] = [int(i) for i in cdata[ckey]]

    def measure(self, dataframe=-1, start=0, stop=-1, keyword="PotEng", name="mean"):
        """[summary]

        [description]

        Parameters
        ----------
        dataframe : int, optional
            [description] (the default is -1, which [default_description])
        start : int, optional
            [description] (the default is 0, which [default_description])
        stop : int, optional
            [description] (the default is -1, which [default_description])
        keyword : {str}, optional
            [description] (the default is "PotEng", which [default_description])
        name : {str}, optional
            [description] (the default is "mean", which [default_description])

        Raises
        ------
        IOError
            [description]
        """
        if name == "mean":
            np.mean(self.data[dataframe][keyword][start:stop])
        elif name == "median":
            np.median(self.data[dataframe][keyword][start:stop])
        elif name == "sigma":
            np.std(self.data[dataframe][keyword][start:])
        else:
            raise IOError("Unknown name: {}".format(name))
