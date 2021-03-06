import os
import numpy as np
from mpi4py import MPI

#==============================================================================#
# Setup MPI
#==============================================================================#
COMM = MPI.COMM_WORLD
SIZE = COMM.Get_size()  # number of processes in communicator
RANK = COMM.Get_rank()  # process' id(s) within a communicator


class LmpSim(object):
    """
    """
    def __init__(self, tstart=None, tstop=None, pstart=None, pstop=None, logsteps=None, runsteps=None, momentum_steps=100, pc_file=None, settings_file=None, input_lmpdat=None, input_lmprst=None, inter_lmprst=None, output_lmprst=None, output_lmplog=None, output_dcd=None, output_lmpdat=None, output_name=None, gpu=False, ncores=None, dielectric=None):
        """
        ncores : list or tuple
        """
        self.thermargs = ["step", "temp", "press", "vol", "density", "cella", "cellb", "cellc", "cellalpha", "cellbeta", "cellgamma", "etotal", "pe", "evdwl", "ecoul", "ebond", "eangle", "edihed", "eimp", "enthalpy"]
        self.tstart = tstart
        self.tstop = tstop
        self.pstart = pstart
        self.pstop = pstop
        self.logsteps = logsteps
        self.runsteps = runsteps
        self.momentum_steps = momentum_steps
        self.pc_file = pc_file
        self.settings_file = settings_file
        self.input_lmpdat = input_lmpdat
        self.input_lmprst = input_lmprst
        self.inter_lmprst = inter_lmprst
        self.output_lmplog = output_lmplog
        self.output_dcd = output_dcd
        self.output_lmpdat = output_lmpdat
        self.output_lmprst = output_lmprst
        self.output_name = output_name
        self.gpu = gpu
        self.ncores = ncores
        self.dielectric = dielectric

    def load_system(self, lmp):
        """Read a lammps restart or data file."""
        if self.output_lmprst is not None and os.path.isfile(self.output_lmprst):
            lmp.command("read_restart {}".format(self.output_lmprst))
        elif self.inter_lmprst is not None and os.path.isfile(self.inter_lmprst):
            lmp.command("read_restart {}".format(self.inter_lmprst))
        elif self.input_lmprst is not None and os.path.isfile(self.input_lmprst):
            lmp.command("read_restart {}".format(self.input_lmprst))
        else:
            lmp.command("read_data {}".format(self.input_lmpdat))

    def unfix_undump(self, pylmp, lmp, used_ranks=None):
        """[summary]

        [description]

        Parameters
        ----------
        pylmp : {[type]}
            [description]
        lmp : {[type]}
            [description]
        used_ranks : {list of int}, optional
            List of ranks that were used for the lammps calculation in order to
            send the lists of fixes, dumps and groups only to those specific
            ranks via 'comm.send()'
            (the default is None, which broadcasts to all ranks)
        """
        lmp_fixes = []
        lmp_dumps = []
        lmp_groups = []

        # getting fixes, dumps and groups only works with pylammps and on rank 0
        if RANK == 0:
            for fix in pylmp.fixes:
                lmp_fixes.append(fix["name"])

            for dump in pylmp.dumps:
                lmp_dumps.append(dump["name"])

            for group in pylmp.groups:
                if group["name"] != "all":
                    lmp_groups.append(group["name"])

        # check whether to broadcast to all cores or only those involved in the
        # lammps instance
        if self.ncores is None:
            lmp_fixes = COMM.bcast(lmp_fixes, 0)
            lmp_dumps = COMM.bcast(lmp_dumps, 0)
            lmp_groups = COMM.bcast(lmp_groups, 0)
        else:
            used_ranks = list(range(self.ncores))[1:]

            if RANK == 0:
                for used_rank in used_ranks:
                    COMM.send(lmp_fixes, used_rank)
            else:
                lmp_fixes = COMM.recv(source=0)

            if RANK == 0:
                for used_rank in used_ranks:
                    COMM.send(lmp_dumps, used_rank)
            else:
                lmp_dumps = COMM.recv(source=0)

            if RANK == 0:
                for used_rank in used_ranks:
                    COMM.send(lmp_groups, used_rank)
            else:
                lmp_groups = COMM.recv(source=0)

            del used_ranks

        for fix in lmp_fixes:
            if "gpu" in fix:
                continue
            lmp.command("unfix {}".format(fix))

        for dump in lmp_dumps:
            lmp.command("undump {}".format(dump))

        for group in lmp_groups:
            try:
                lmp.command("group {} delete".format(group))
            except:
                print(("***Warning: Group {} could not be deleted in lammps!".format(group)))

    def thermo(self, lmp, hb_group="all"):
        #TODO:  hb_group is not necessary since hbonds will only be calculated
        #       for the defined atom types
        """
        Log thermodynamic data.
        """
        # compute dreiding h-bond energies, if hbond/dreiding/lj setting is set
        try:
            lmp.command("compute hb all pair hbond/dreiding/lj\n")
            lmp.command("variable n_hbond equal c_hb[1]")
            lmp.command("variable E_hbond equal c_hb[2]")

            if "v_n_hbond" not in self.thermargs:
                self.thermargs.append("v_n_hbond")

            if "v_E_hbond" not in self.thermargs:
                self.thermargs.append("v_E_hbond")

        except:
            pass

        lmp.command("thermo_style custom " + " ".join(self.thermargs))
        lmp.command("thermo_modify lost warn flush yes")
        #lmp.command("thermo_modify line multi format float %g")
        lmp.command("thermo {}".format(self.logsteps))

    def dump(self, lmp, ndcd=None, nrst=None, unwrap=False):
        """
        Dump dcd and lammps restart files.
        """
        # default values
        if ndcd is None:
            ndcd = self.logsteps

        if nrst is None:
            nrst = self.logsteps * 10

        # trajectory
        lmp.command("dump trajectory all dcd {} {}".format(ndcd, self.output_dcd))
        lmp.command("restart {} {} {}".format(nrst, self.inter_lmprst, self.inter_lmprst))

        if unwrap is True:
            lmp.command("dump_modify trajectory unwrap yes")

    def read_dump(self, lmp, filename, nframes):
        """The read dump command from lammps.

        Shortcut for the read dump command - xyz coordinates always included

        Arguments:
            lmp {[type]} -- [description]
            file {[type]} -- [description]
        """
        lmp.command(f"read_dump {filename} {nframes} x y z box no format xyz")

    def fix_berendsen(self, lmp, group, ensemble, keyword, integrator="nve"):
        """
        """
        # nve
        lmp.command("fix integrator {} {}".format(group, integrator))

        # nvt
        if (self.tstart is not None and self.tstop is not None) and (ensemble == "npt" or ensemble == "nvt"):
            lmp.command("fix thermostat {} temp/berendsen {} {} 0.5".format(group, self.tstart, self.tstop))

        # npt
        if (self.pstart is not None and self.pstop is not None) and ensemble == "npt":
            lmp.command("fix barostat {} press/berendsen {} {} {} 50".format(group, keyword, self.pstart, self.pstop))

    def fix_langevin(self, lmp, group, ensemble, keyword, integrator="nve"):
        """
        """
        # nve
        lmp.command("fix integrator {} {}".format(group, integrator))

        # nvt
        if (self.tstart is not None and self.tstop is not None) and (ensemble == "npt" or ensemble == "nvt"):
            lmp.command("fix thermostat {} langevin {} {} 100 {}".format(group, self.tstart, self.tstop, np.random.randint(0, 10000000)))

        # npt
        if (self.pstart is not None and self.pstop is not None) and ensemble == "npt":
            lmp.command("fix barostat {} press/berendsen {} {} {} 50".format(group, keyword, self.pstart, self.pstop))

    def fix_hoover(self, lmp, group, ensemble, keyword=None):
        nvt_section = "fix integrator {} {} temp {} {} 0.5".format(group, ensemble, self.tstart, self.tstop)

        if (self.pstart is not None and self.pstop is not None) and ensemble == "npt" and keyword is not None:
            npt_section = " {} {} {} 50".format(keyword, self.pstart, self.pstop)
            lmp.command(nvt_section + npt_section)
        elif ensemble == "nph":
            raise Warning("Not implemented yet!")
        else:
            lmp.command(nvt_section)

    def minimize(self, lmp, style="cg", keyword=None):
        """

        Sources
        ---
        https://lammps.sandia.gov/doc/min_modify.html

        Parameters
        ----------
        lmp : Lammps-object
        style : str
            minimization style to use: sd, cg, fire
        keyword : str
            relaxation keyword if box should be scaled as well: iso, aniso, tri

        """
        if keyword is not None:
            lmp.command("fix box_relax all box/relax {} {}".format(keyword, self.pstart))

        lmp.command("min_style {}".format(style))
        lmp.command("min_modify dmax 2.0")
        lmp.command("minimize 1.0e-9 1.0e-12 100000 1000000")

        # tidy up
        if keyword is not None:
            lmp.command("unfix box_relax")

    def use_gpu(self, lmp, neigh=True):
        """
        """
        if neigh is True:
            lmp.command("package gpu 1 neigh yes")
        else:
            lmp.command("package gpu 1 neigh no")

        lmp.command("suffix gpu")
