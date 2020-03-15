import ag_clmsv as clmsv
import ag_lmplog as llog

__version__ = "2017-05-24"


class LogUnification(clmsv.Clmsv, llog.LmpLog):
    """
    Read and write several data types.
    a = LogUnification()
    a.read_lmplog(lmplog)
    a.write_clmsv("Test.clmsv", split=False)
    a.read_clmsv("Test.clmsv")
    """

    def __init__(self):
        """
        Initialize necessary classes.
        """
        clmsv.Clmsv.__init__(self)
        llog.LmpLog.__init__(self)
