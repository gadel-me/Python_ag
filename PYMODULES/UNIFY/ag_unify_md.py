
import ag_amber as agambr
import ag_gaussian as aggau
import ag_lammps as aglmps
import ag_mae as agmae
import ag_xyz as agxyz
import ag_pdb as agpdb
import ag_pw as agpw

__version__ = "2018-03-20"


class Unification(agambr.AmberStuff,
                  aggau.GauStuff,
                  aglmps.LmpStuff,
                  agmae.MaestroStuff,
                  agxyz.XYZ,
                  agpdb.PdbStuff,
                  agpw.PwStuff):
    """
    Unify all classes to read/write different file formats.
    """
    def __init__(self):
        """
        """
        agambr.AmberStuff.__init__(self)
        aggau.GauStuff.__init__(self)
        aglmps.LmpStuff.__init__(self)
        agmae.MaestroStuff.__init__(self)
        agxyz.XYZ.__init__(self)
        agpdb.PdbStuff.__init__(self)
        agpw.PwStuff.__init__(self)
