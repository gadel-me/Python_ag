import log_universe as logu
import collections


class Clmsv(logu.LogUniverse):
    """
    Read and write clsmv-files. Clmsv or 'column separated' files with
    '***'-filled lines as separators for different entries (e.g. Several runs
    during one simulation or input from different log-files).
    """
    def __init__(self):
        logu.LogUniverse.__init__(self)

    def read_clmsv(self, *clmsv_files):
        """
        Read a clmsv to self.data and its including dictionaries.
        """
        # /// gather data ///
        for clmsv_file in clmsv_files:

            with open(clmsv_file) as f_in:
                for line in f_in:

                    if line.startswith("#") is True:
                        continue

                    keys = line.split()
                    cdata = collections.OrderedDict(zip(keys, [[] for i in xrange(len(keys))]))

                    while 1:
                        # check if end of file, append current data
                        try:
                            line = f_in.next()
                        except StopIteration:
                                self.data.append(cdata)
                                break

                        # check if end of entry, append current data
                        if line == "\n":
                            self.data.append(cdata)
                            break
                        else:
                            values = [float(i) for i in line.split()]
                            for key, value in zip(keys, values):
                                cdata[key].append(value)

        self._int_values()

    def write_clmsv(self, clmsv_out, split=False):
        """
        Column-separate values from data-dictionary and write them into a file.
        clmsv_out       str; name of output file
        split           boolean; split several entries from self.data to different
                        files (True) or write them newline-separated to one file.
        """

        # write all data to several files
        if split is True:

            for cnum, cdata in enumerate(self.data):
                clmsv_out_file = "{}_{}.clmsv".format(clmsv_out, cnum)

                list_of_keys = cdata.keys()
                list_of_values = cdata.values()

                with open(clmsv_out_file, "w") as f_out:

                    # write header
                    for ikey in list_of_keys:
                        f_out.write("{:>20s}".format(ikey))
                    f_out.write("\n")

                    # write body
                    for row in zip(*list_of_values):
                        for number in row:
                            f_out.write("{:> 20}".format(number))
                        f_out.write("\n")

                    f_out.write("\n")

        # write all data to one file
        else:

            with open(clmsv_out, "w") as f_out:
                for cdata in self.data:
                    list_of_keys = cdata.keys()
                    list_of_values = cdata.values()

                    # write header
                    for ikey in list_of_keys:
                        f_out.write("{:>20s}".format(ikey))
                    f_out.write("\n")

                    # write body (row-wise -> zip)
                    for row in zip(*list_of_values):
                        for number in row:
                            f_out.write("{:> 20}".format(number))
                        f_out.write("\n")

                    # write tail
                    f_out.write("\n")

    def add_clmsv(self, clmsv2add):
        """
        Add another clmsv-file to already read one(s).
        """
        pass
