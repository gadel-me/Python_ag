#!/usr/bin/env python
from __future__ import print_function, division
import collections
import log_universe as logu

__version__ = "2017-05-24"


class LmpLog(logu.LogUniverse):
    """
    Class for reading log.lammps file(s).
    """
    def __init__(self):
        logu.LogUniverse.__init__(self)

    def read_lmplog(self, *lmplogs):
        """
        Read a log.lammps file with one or more thermo entries.
        """
        thermo       = None
        run          = None

        for lmplog in lmplogs:

            with open(lmplog, "r") as log_in:
                # find line with thermo-keywords
                for line in log_in:

                    if (line.startswith("thermo") and
                       not line.startswith("thermo_modify") and
                       not line.startswith("thermo_style")):

                        try:
                            thermo = int(line.split()[1])
                        except ValueError:
                            pass

                    elif line.startswith("run"):
                        # number of total timesteps
                        run = int(line.split()[1])
                    # Memory usage per processor: local machine? (or gcc?)
                    # Per MPI rank memory allocation: server? (or intel compiler?)
                    elif (line.startswith("Memory usage per processor") or
                          line.startswith("Per MPI rank memory allocation")):
                        num_frames = int((run/thermo) + 1)
                        # deploy containers for the data to come
                        line = log_in.next()
                        keys = line.split()
                        cdata = collections.OrderedDict(
                            zip(keys, [[] for i in xrange(len(keys))])
                        )

                        # cycle current entry
                        for cf in xrange(num_frames):

                            try:
                                line = log_in.next()
                            # stop reading when eof (i.e. run was aborted)
                            except StopIteration:
                                print("***Info: EOF was reached sooner than expected!")
                                self.data.append(cdata)
                                break

                            # Stop reading and append entry
                            # floatify values
                            values = [float(i) for i in line.split()]
                            for key, value in zip(keys, values):
                                cdata[key].append(value)

                        self.data.append(cdata)

        self._int_values()

    def add_log(self, log_to_add):
        """
        Add data from another log-file to current log-file.
        """
        log = LmpLog()
        log.read_lmplog(log_to_add)
        self.data.extend(log.data)
