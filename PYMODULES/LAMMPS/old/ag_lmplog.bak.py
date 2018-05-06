#!/usr/bin/env python
from __future__ import print_function, division
import collections
import log_universe as logu

__version__ = "2017-05-05"


class LmpLog(logu.LogUniverse):
    """
    """
    def __init__(self):
        logu.LogUniverse.__init__(self)

    def read_lmplog(self, lmplog):
        keep_reading = False
        therm_keys = None
        therm_values = []

        with open(lmplog, "r") as log_in:

            # find line with thermo-keywords
            for line in log_in:
                if "Memory usage per processor" in line:
                    keep_reading = True
                    break

            # line 'Memory usage per processor' missing -> calculation did not start properly
            if keep_reading is False:
                raise UserWarning("No thermo-keywords found!" +
                                  "Are you sure your calculation started at all?")

            # /// deploy containers for the data to come
            line = log_in.next()
            therm_keys = line.split()
            for iarg in xrange(len(therm_keys)):
                therm_values.append([])
            data = collections.OrderedDict(zip(therm_keys, therm_values))
            del therm_values

            # /// fetch values for thermo-keywords
            while keep_reading:
                try:
                    line = log_in.next()
                # stop reading after end of file; that is the case when run was aborted
                except StopIteration:
                    keep_reading = False
                    break

                # Stop reading after here
                if "Loop time of" in line:
                    keep_reading = False
                    break

                # floatify values
                cdata = [float(i) for i in line.split()]

                # gather results
                for ckey, cvalue in zip(therm_keys, cdata):
                    if ckey == "Step":
                        cvalue = int(cvalue)
                    data[ckey].append(cvalue)

        self.therm_keys = therm_keys
        self.data = data

    def add_log(self, log_to_add):
        """
        Add data from another log-file to current log-file.
        """
        log = LmpLog()
        log.read_lmplog(log_to_add)

        for ikey, jkey in zip(self.data, log.data):
            self.data[jkey] += log.data[jkey]
