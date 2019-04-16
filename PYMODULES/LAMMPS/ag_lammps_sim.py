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
    def __init__(self, tstart=None, tstop=None, pstart=None, pstop=None, logsteps=None, runsteps=None, pc_file=None, settings_file=None, input_lmpdat=None, input_lmprst=None, inter_lmprst=None, output_lmprst=None, output_lmplog=None, output_dcd=None, output_lmpdat=None, output_name=None, gpu=False):
        """
        """
        self.thermargs = ["step", "temp", "press", "vol", "density", "cella", "cellb", "cellc", "cellalpha", "cellbeta", "cellgamma", "etotal", "pe", "evdwl", "ecoul", "ebond", "eangle", "edihed", "eimp", "enthalpy"]
        self.tstart = tstart
        self.tstop = tstop
        self.pstart = pstart
        self.pstop = pstop
        self.logsteps = logsteps
        self.runsteps = runsteps
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

    def unfix_undump(self, pylmp, lmp):
        """
        Remove all fixes and dumps.
        """
        lmp_fixes = []
        lmp_dumps = []

        if RANK == 0:
            for fix in pylmp.fixes:
                lmp_fixes.append(fix["name"])
            for dump in pylmp.dumps:
                lmp_dumps.append(dump["name"])

        lmp_fixes = COMM.bcast(lmp_fixes, 0)
        lmp_dumps = COMM.bcast(lmp_dumps, 0)

        for fix in lmp_fixes:
            if "gpu" in fix:
                continue
            lmp.command("unfix {}".format(fix))

        for dump in lmp_dumps:
            lmp.command("undump {}".format(dump))

    def thermo(self, lmp, hb_group="all"):
        """
        Log thermodynamic data.
        """
        # compute dreiding h-bond energies, if hbond/dreiding/lj setting is set
        try:
            lmp.command("compute hb all pair hbond/dreiding/lj\n")
            lmp.command("variable n_hbond equal c_hb[1]")
            lmp.command("variable E_hbond equal c_hb[2]")
            self.thermargs.append("v_n_hbond")
            self.thermargs.append("v_E_hbond")
        except:
            pass

        lmp.command("thermo_style custom " + " ".join(self.thermargs))
        lmp.command("thermo_modify lost warn flush yes")
        #lmp.command("thermo_modify line multi format float %g")
        lmp.command("thermo {}".format(self.logsteps))

    def dump(self, lmp, ndcd=None, nrst=None, unwrap=True):
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

    def fix_berendsen(self, lmp, group, ensemble, keyword):
        """
        """
        # nve
        lmp.command("fix integrator {} nve".format(group))

        # nvt
        if (self.tstart is not None and self.tstop is not None) and (ensemble == "npt" or ensemble == "nvt"):
            lmp.command("fix thermostat {} temp/berendsen {} {} 0.5".format(group, self.tstart, self.tstop))

        # npt
        if (self.pstart is not None and self.pstop is not None) and ensemble == "npt":
            lmp.command("fix barostat {} press/berendsen {} {} {} 50".format(group, keyword, self.pstart, self.pstop))

    def fix_hoover(self, lmp, group, ensemble, keyword):
        nvt_section = "fix integrator {} {} temp {} {} 0.5".format(group, ensemble, self.tstart, self.tstop)
        npt_section = " {} {} {} 50".format(keyword, self.pstart, self.pstop)

        if (self.pstart is not None and self.pstop is not None) and ensemble == "npt":
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
        if keyword:
            lmp.command("fix box_relax all box/relax {} {}".format(keyword, self.pstart))

        lmp.command("min_style {}".format(style))
        lmp.command("min_modify dmax 2.0")
        lmp.command("minimize 1.0e-9 1.0e-12 100000 1000000")

        # tidy up
        if keyword:
            lmp.command("unfix box_relax")

    def use_gpu(self, lmp, neigh=True):
        """
        """
        if neigh is True:
            lmp.command("package gpu 1 neigh yes")
        else:
            lmp.command("package gpu 1 neigh no")

        lmp.command("suffix gpu")
