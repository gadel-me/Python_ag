from __future__ import print_function, division
import os
import copy
import shutil as sl
import re
import argparse
import math
import pdb
import time
import numpy as np
from natsort import natsorted
#import itertools as it
import scipy.stats
from mpi4py import MPI
from lammps import lammps, PyLammps
import Transformations as cgt
import md_elements as mde
import md_box as mdb
#import ag_unify_md as agum
import ag_geometry as agm
import ag_lammps as aglmp
import ag_unify_log as agul
import ag_vectalg as agv
import ag_statistics as ags

#==============================================================================#
# Setup MPI
#==============================================================================#

comm = MPI.COMM_WORLD
size = comm.Get_size()  # number of processes in communicator
rank = comm.Get_rank()  # process' id(s) within a communicator

#==============================================================================#
# Helper functions and variables
#==============================================================================#
_sqrt_3 = math.sqrt(3)


class LmpShortcuts(object):
    """
    """
    def __init__(self, tstart=None, tstop=None, pstart=None, pstop=None, logsteps=None, runsteps=None, pc_file=None, settings_file=None, input_lmpdat=None, input_lmprst=None, inter_lmprst=None, output_lmprst=None, output_lmplog=None, output_dcd=None, output_lmpdat=None, output_name=None, gpu=False):
        """
        """
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

    def read_system(self, lmp):
        """
        Read a restart file from another stage, previous run or a data file instead.
        """
        # reading order: output restart -> input restart -> data
        if self.output_lmprst is not None and os.path.isfile(self.output_lmprst):
            lmp.command("read_restart {}".format(self.output_lmprst))
        elif self.input_lmprst is not None and os.path.isfile(self.input_lmprst):
            lmp.command("read_restart {}".format(self.input_lmprst))
        else:
            lmp.command("read_data {}".format(self.input_lmpdat))
            lmp.command("velocity all create {} 4928459 mom yes rot yes dist gaussian".format(self.tstart))

    def unfix_undump(self, pylmp, lmp):
        """
        Remove all fixes and dumps.
        """
        lmp_fixes = []
        lmp_dumps = []

        if rank == 0:
            for fix in pylmp.fixes:
                lmp_fixes.append(fix["name"])
            for dump in pylmp.dumps:
                lmp_dumps.append(dump["name"])

        lmp_fixes = comm.bcast(lmp_fixes, 0)
        lmp_dumps = comm.bcast(lmp_dumps, 0)

        for fix in lmp_fixes:
            if "gpu" in fix:
                continue
            lmp.command("unfix {}".format(fix))

        for dump in lmp_dumps:
            lmp.command("undump {}".format(dump))

    def thermo(self, lmp):
        """
        Log thermodynamic data.
        """
        thermargs = ["step", "temp", "press", "vol", "density", "cella", "cellb", "cellc", "cellalpha", "cellbeta", "cellgamma", "etotal", "pe", "evdwl", "ecoul", "ebond", "eangle", "edihed", "eimp", "enthalpy"]
        lmp.command("thermo_style custom " + " ".join(thermargs))
        lmp.command("thermo_modify lost warn flush yes")
        #lmp.command("thermo_modify line multi format float %g")
        lmp.command("thermo {}".format(self.logsteps))

    def dump(self, lmp, unwrap=True):
        """
        """
        # trajectory
        lmp.command("dump trajectory all dcd {} {}".format(self.logsteps, self.output_dcd))
        lmp.command("restart {} {} {}".format(self.logsteps*50, self.inter_lmprst, self.inter_lmprst))

        if unwrap is True:
            lmp.command("dump_modify trajectory unwrap yes")

    def berendsen(self, lmp, group="all"):
        """
        """
        lmp.command("fix integrator {} nve".format(group))
        lmp.command("fix thermostat {} temp/berendsen {} {} 0.5".format(group, self.tstart, self.tstop))
        lmp.command("fix barostat {} press/berendsen iso {} {} 50".format(group, self.pstart, self.pstop))

    def use_gpu(self, lmp, neigh=True):
        """
        """
        if neigh:
            lmp.command("package gpu 1 neigh yes")
        else:
            lmp.command("package gpu 1 neigh no")

        lmp.command("suffix gpu")

    def minimize(self, lmp, style="cg", box_relax=False):
        """
        """
        if box_relax:
            lmp.command("fix box_relax all box/relax aniso {}".format(self.pstart))

        lmp.command("min_style {}".format(style))
        lmp.command("min_modify dmax 2.0")
        lmp.command("minimize 1.0e-9 1.0e-12 100000 1000000")


#def check_energy_convergence(logfiles,
#                             keyword="PotEng",
#                             percentage=100,
#                             debug=False):
#    """
#    TODO No matter the outcome, save an image of the qq plot with the
#    TODO respective data after several attempts to achieve normality.
#    TODO The bigger the sample size the more probable it is that the normal
#    TODO tests state to wrongly reject H0.
#
#    Calculate the skewness, p-value and z-score, to check if the values
#    of 'keyword' are normally distributed (i.e. the MD-Simulation run ended
#    successfully). "Generally speaking, the p-value is the probability of an
#    outcome different that what was expected given the null hypothesis - in this
#    case the probability of getting a skewness different from that of a normal
#    distribution (which is 0 because of symmetry).
#
#    Normality tests:
#
#    >   The Shapiro-Wilk test evaluates a data sample and quantifies how likely it is that the data
#        was drawn from a Gaussian distribution, named for Samuel Shapiro and Martin Wilk.
#
#    >   The D'Agostino's K^2 test calculates summary statistics from the data, namely kurtosis and
#        skewness, to determine if the data distribution departs from the normal distribution,
#        named for Ralph D'Agostino.
#        - Skew is a quantification of how much a distribution is pushed left or right,
#          a measure of asymmetry in the distribution.
#        - Kurtosis quantifies how much of the distribution is in the tail. It is a simple and commonly
#          used statistical test for normality.
#
#    >   Anderson-Darling Test is a statistical test that can be used to evaluate whether a data
#        sample comes from one of among many known data samples, named for Theodore Anderson and
#        Donald Darling.
#
#    A quick word on normality tests in general:
#        The theory for this test is based on the probability of getting a rational number from a
#        truly continuous distribution defined on the reals. The main goal of this test is to quickly
#        give a p-value for those that feel it necessary to test the uninteresting and uninformative
#        null hypothesis that the data represents an exact normal, and allows the user to then move
#        on to much more important questions, like "is the data close enough to the normal to use
#        normal theory inference?". After running this test (or better instead of running this and
#        any other test of normality) you should ask yourself what it means to test for normality
#        and why you would want to do so. Then plot the data and explore the interesting/useful
#        questions.
#
#    Source: https://machinelearningmastery.com/a-gentle-introduction-to-normality-tests-in-python/
#            https://www.rdocumentation.org/packages/TeachingDemos/versions/2.10/topics/SnowsPenultimateNormalityTest
#            https://stackoverflow.com/questions/7781798/seeing-if-data-is-normally-distributed-in-r/7788452#7788452
#
#    Input:
#        > logfiles      list; all written lammps-logfiles
#        > keyword       str; keywords from thermo-output
#        > percentage   int; number of last values in %,
#                        e.g. 80 means last 80 % of all values
#        > debug         boolean; enable debug messaging
#        > stage         str; which stage of the kwz-approach (relevant for
#                        debugging)
#    Output:
#        > skewness, pvalue, zscore
#        > data          list; all data for keyword from all files given
#    """
#    data = []
#
#    # gather all values from all logfiles given
#    log_data = agul.LogUnification()
#    log_data.read_lmplog(*logfiles)
#    data = [value for data_index in xrange(len(log_data.data)) for value in
#            log_data.data[data_index][keyword]]
#    num_values = len(data)
#    percentage /= 100  # percentage to per cent
#    testdata = data[-int(percentage * num_values):]
#    len_testdata = len(testdata)
#    min_value_testdata = min(testdata)
#
#    if debug is True:
#        print(("***Info: Smallest value of "
#               "{} % of all measured data is: {}").format(percentage * 100,
#                                                          min_value_testdata))
#
#    # normality tests
#    p_alpha = 0.05  # alpha helps interpreting the p-value from the normality test at hand
#
#    # Shapiro-Wilk Test (only for ~ 2000 samples)
#    normal_shapiro = False
#
#    if len_testdata <= 2000:
#        write_to_log("Shapiro-Wilk Test:\n")
#        stat_shapiro, p_shapiro = scipy.stats.shapiro(testdata)
#        normal_shapiro = p_shapiro > p_alpha
#
#        if normal_shapiro:
#            write_to_log("{} > {}: Sample looks Gaussian (fail to reject H0)\n".format(p_shapiro, p_alpha))
#        else:
#            write_to_log("{} < {}: Sample does not look Gaussian (reject H0)\n".format(p_shapiro, p_alpha))
#
#        write_to_log("\n")
#
#    # D'Agostino's K^2 Test
#    write_to_log("D'Agostino's K^2 Test:\n")
#    stat_agostino, p_agostino = scipy.stats.normaltest(testdata)
#    normal_agostino = p_agostino > p_alpha
#
#    if normal_agostino:
#        write_to_log("Sample looks Gaussian (fail to reject H0)\n")
#    else:
#        write_to_log("Sample does not look Gaussian (reject H0)\n")
#
#    write_to_log("\n")
#
#    # Anderson-Darling Test
#    write_to_log("Anderson-Darling Test:\n")
#    result_anderson = scipy.stats.anderson(testdata, dist="norm")
#    write_to_log("\t> Statistic: {}\n".format(result_anderson.statistic))
#
#    for idx in range(len(result_anderson.critical_values)):
#        sl_anderson = result_anderson.significance_level[idx]
#        cv_anderson = result_anderson.critical_values[idx]
#
#        # check if the null hypothesis can be rejected (H0: normal distributed)
#        normal_anderson = result_anderson.statistic < cv_anderson
#
#        if normal_anderson:
#            write_to_log("\t> {:> .3f}: {:> .3f}, data looks normal (fail to reject H0)\n".format(sl_anderson, cv_anderson))
#        else:
#            write_to_log("\t> {:> .3f}: {:> .3f}, data does not look normal (reject H0)\n".format(sl_anderson, cv_anderson))
#
#        # if H0 is ok at any confidence level, stop further testing
#        if normal_anderson is True:
#            break
#
#    write_to_log("\n")
#
#    # Compute the skewness of a data set, For normally distributed data,
#    # the skewness should be about 0
#    skewness = scipy.stats.skew(testdata)
#    write_to_log("Skewness (should be -0.3 < skewness < 0.3): {:> .12f}\n".format(skewness))
#
#    # determine if the skewness is close enough to 0 (statistically speaking)
#    stat_skewness, p_skewness = scipy.stats.skewtest(testdata)
#    write_to_log("P-Value (skewness, should be > 0.05): {:> .12f}\n".format(p_skewness))
#    write_to_log("\n")
#    normal_skewness = (-0.3 < skewness < 0.3) and (p_skewness > 0.05)
#
#    # determine if the kurtosis is close enough to 0 (statistically speaking)
#    normal_kurtosis = False
#
#    if len_testdata > 20:
#        kurtosis = scipy.stats.kurtosis(testdata)
#        write_to_log("Kurtosis (should be -0.3 < kurtosis < 0.3): {:> .12f}\n".format(kurtosis))
#        stat_kurtosis, p_kurtosis = scipy.stats.kurtosistest(testdata)
#        normal_kurtosis = (-0.3 < kurtosis < 0.3) and (p_kurtosis > 0.05)
#        write_to_log("Statistic (Kurtosis): {:> 4.12f}\n".format(stat_kurtosis))
#        write_to_log("P-Value (Kurtosis, should be > 0.05): {:> 4.12f}\n".format(p_kurtosis))
#        write_to_log("\n")
#    else:
#        write_to_log("Need more than 20 values for Kurtosis-Test!\n")
#
#    write_to_log("\n")
#
#    # chi squared test
#    num_bins = int(np.sqrt(len(data)))
#    histo, bin_edges = np.histogram(data, bins=num_bins, normed=False)
#    a1, b1 = scipy.stats.norm.fit(data)
#    cdf = scipy.stats.norm.cdf(bin_edges, a1, b1)
#    #scaling_factor = len(data) * (x_max - x_min) / num_bins
#    scaling_factor = len(data)
#    # expected frequencies (haufigkeiten)
#    expected_values = scaling_factor * np.diff(cdf)
#    chisquare_results = scipy.stats.chisquare(histo, f_exp=expected_values, ddof=2)
#    #write_to_log(chisquare_results[1], "\n")
#
#    # only use tests that are allowed for the amount of given values for the shape and normality
#    if len_testdata > 20:
#        normal_shape = normal_skewness or normal_kurtosis
#    else:
#        normal_shape = normal_skewness
#
#    # accept normal distribution if any test states is is normally distributed
#    if len_testdata <= 2000:
#        normal_distributed = normal_shape and (normal_shapiro or normal_agostino or normal_anderson)
#    else:
#        normal_distributed = normal_shape and (normal_agostino or normal_anderson)
#
#    return (normal_distributed, min_value_testdata)


def _molecules_radii(lmp_sys):
    """
    """
    # get center of geometry and sphere of each solvate molecule which is placed
    # inside the solvent box
    radii_sphere = []
    cogs_solvate = []

    for molecule in lmp_sys.molecules:
        sphere, cog = lmp_sys.get_mol_radius(-1, *molecule)
        radii_sphere.append(sphere)
        cogs_solvate.append(cog)

    return (radii_sphere, cogs_solvate)


def read_system(lmpdat, dcd=None, frame_idx_start=-1, frame_idx_stop=-1):
    """
    Read a lammps data file and a dcd file on top (optional).

    This function simplifies the process of reading a lammps data file and loading
    a dcd file on top of it for further coordinates. It returns a Universe object
    which has many methods to manipulate it with.

    Parameters
    ----------
    lmpdat : str
        Name of the lammps data file

    dcd : str, optional
        Name of the dcd file with the same amounts of atoms as the lammps data
        file

    frame_idx : int, optional
        Index of the frame to use from the dcd file.

    Returns
    -------
    md_sys : LmpStuff object
        An object of LmpStuff which can be further processed.

    """
    # read output from quenching
    md_sys = aglmp.LmpStuff()
    md_sys.read_lmpdat(lmpdat)

    if dcd is not None:
        md_sys.import_dcd(dcd)
        # since we are only interested in one frame, delete all others
        md_sys.ts_coords = []
        md_sys.ts_boxes = []

        # enable reading the last frame with negative indexing
        if frame_idx_stop == -1:
            md_sys.read_frames(frame=frame_idx_start - 1, to_frame=frame_idx_stop)
        else:
            md_sys.read_frames(frame=frame_idx_start, to_frame=frame_idx_stop + 1)

        md_sys.close_dcd()

    return md_sys


def merge_sys(sys_a, sys_b, frame_idx_a=-1, frame_idx_b=-1, pair_coeffs=False):
    """
    """
    sys_both = copy.deepcopy(sys_a)
    sys_both.extend_universe(sys_b, u1_frame_id=frame_idx_a, u2_frame_id=frame_idx_b, mode="merge")

    if pair_coeffs is True:
        sys_both.mix_pair_types(mode="ii", mix_style="arithmetic")

    sys_both.fetch_molecules_by_bonds()
    sys_both.mols_to_grps()
    return sys_both


def _write_data(lmpdat_out, lmpdat_a, lmpdat_b=None, dcd_a=None, dcd_b=None, frame_idx_a=-1, frame_idx_b=-1, pair_coeffs=False):
    """
    """
    sys_a = read_system(lmpdat_a, dcd_a, frame_idx_a)

    if lmpdat_b is not None:
        sys_b = read_system(lmpdat_b, dcd_b, frame_idx_b)
        # since we only read one frame, only this frame combined will be written
        sys_ab = merge_sys(sys_a, sys_b, pair_coeffs=pair_coeffs)
    else:
        sys_ab = sys_a

    sys_ab.change_indices(incr=1, mode="increase")
    sys_ab.write_lmpdat(lmpdat_out, cgcmm=True)


def check_aggregate(mdsys, frame_id=-1, atm_atm_dist=4, unwrap=False, debug=False):
    """
    Check if several molecules form an aggregate.

    #TODO Buggy, does not recognize aggregates properly when they form a chain?

    At the moment this works only for orthogonal cells (reason: sub-cell length)
    Check if the aggregate did not get dissolved in the process. This is achieved
    by calculating the center of geometry for each molecule and subsequent com-
    parison of the distance between each center of geometry. The smallest distance
    must not be larger than four times the radius of the biggest molecule in the
    system.

    Input:
        > mdsys             ag_unify_md.Unification; system output from 'SYS-PREPARATION'-step
        > frame_id          int; frame to process
        > atm_atm_dist      float; radius around an atom in which another atom
                            should be positioned
        > debug             boolean; True if further output should be given,
                            default=False

    Output:
        > aggregate_ok      boolean; True, if aggregate is still intact, False if
                            a part of it drifted away
    """
    aggregate_ok = False

    # lammps may internally have the coordinates wrapped -> unwrap cell when
    # necessary (e.g. coordinates were taken from restart)
    if unwrap is True:
        print("***Check Aggregate Info: Unwrapping cell")
        mdsys.unwrap_cell(frame_id)

    # maximal possible distance between two atoms
    cube_side = atm_atm_dist / _sqrt_3  # -> atm_atm_dist/(2*math.sqrt(3)) * 2

    if debug is True:
        print("***Check Aggregate Info: Cube side {}".format(cube_side))

    mdsys.create_linked_cells(frame_id, rcut_a=cube_side, rcut_b=cube_side,
                              rcut_c=cube_side)

    # include same molecule or it gets missing in our aggregates set
    close_atoms, aggregates = mdsys.chk_atm_dist(frame_id,
                                                 min_dist=atm_atm_dist,
                                                 exclude_same_molecule=False,
                                                 get_aggregates=True,
                                                 debug=False)

    if debug is True:
        print("***Check Aggregate Info: Number of cubes in direction a {}".format(mdsys.ts_lnk_cls[frame_id].ra))
        print("***Check Aggregate Info: Number of cubes in direction b {}".format(mdsys.ts_lnk_cls[frame_id].rb))
        print("***Check Aggregate Info: Number of cubes in direction c {}".format(mdsys.ts_lnk_cls[frame_id].rc))
        print("***Group IDs of all aggregates: {}".format(aggregates))
        print("***IDs of close atoms: {}".format(" ".join([str(i) for i in close_atoms])))

    # aggregate is only o.k. if all molecules are part of it
    if len(aggregates) == 1 and len(aggregates[0]) == len(mdsys.molecules):
        aggregate_ok = True

        if debug is True:
            print("***Check Aggregate Info: Aggregate looks fine!")

    return aggregate_ok


def cut_box(lmpdat_out, lmpdat, box, dcd=None, frame_idx=-1):
    """
    Cut a smaller solvent box from a bigger one that will be used during the simulation.

    Given a (in the best case) larger solvent box, a smaller one will be cut.
    Having only as many solvent molecules as absolutely necessary reduces the
    calculation time. This is possible since the potential energy for a group
    of atoms may be calculated with lammps. The center of the box will be at
    the origin.
    Only reads one frame.

    > lmpdat            str; lammps data file of the system to cut from
    > box               Box; size of the box to cut out
    > lmpdat_out       str; name of new lammps data file with cut coordinates
    """
    md_sys = read_system(lmpdat, dcd, frame_idx_start=frame_idx - 1, frame_idx_stop=frame_idx)

    # generate planes that describe the box
    cart_box = copy.deepcopy(box)
    cart_box.box_lmp2cart()

    # define planes which form the box to cut (plane vectors always need the same
    # origin or they will be shifted)
    plane_ab = agv.get_plane(cart_box.crt_a, cart_box.crt_b)
    plane_ca = agv.get_plane(cart_box.crt_c, cart_box.crt_a)
    plane_bc = agv.get_plane(cart_box.crt_b, cart_box.crt_c)

    # opposite planes to ab, ca and bc
    plane_ab_c = [-1 * i for i in agv.get_plane(cart_box.crt_a, cart_box.crt_b, cart_box.crt_c)]
    plane_ca_b = [-1 * i for i in agv.get_plane(cart_box.crt_c, cart_box.crt_a, cart_box.crt_b)]
    plane_bc_a = [-1 * i for i in agv.get_plane(cart_box.crt_b, cart_box.crt_c, cart_box.crt_a)]

    # cut plane
    md_sys.cut_shape(-1, True, plane_ab, plane_ab_c, plane_ca, plane_ca_b, plane_bc, plane_bc_a)

    # enlarge box a little so atoms at the edge will not clash due to pbc
    cart_box.box_cart2lat()
    cart_box.ltc_a += 4
    cart_box.ltc_b += 4
    cart_box.ltc_c += 4
    cart_box.box_lat2lmp()

    # def new box vectors by given box angles
    md_sys.ts_boxes[0] = cart_box

    md_sys.ts_boxes[0].lmp_xy = None
    md_sys.ts_boxes[0].lmp_xz = None
    md_sys.ts_boxes[0].lmp_yz = None

    md_sys.mols_to_grps()
    md_sys.change_indices(incr=1, mode="increase")
    md_sys.write_lmpdat(lmpdat_out, cgcmm=True)


# System Preparation ==========================================================#

def _rot_matrix():
    """
    Create an arbitrary rotational matrix.

    Returns
    -------
    m_arb_xyz : array of arrays
        The arbitrary rotational matrix

    """
    # rotation matrix, consisting of 3 arbitrary rotations A, B and C
    m_arb_rot_x = agm.arb_rot_matrix([1, 0, 0])  # A
    m_arb_rot_y = agm.arb_rot_matrix([0, 1, 0])  # B
    m_arb_rot_z = agm.arb_rot_matrix([0, 0, 1])  # C
    # multiply all three matrices; ABC=A(BC)=(AB)C
    m_arb_xy = np.matmul(m_arb_rot_x, m_arb_rot_y)  # AB
    m_arb_xyz = np.matmul(m_arb_xy, m_arb_rot_z)  # ABC
    return m_arb_xyz


def _rotate_sys(md_sys):
    """
    Rotate system arbitrarily.
    """
    m_arb_xyz = _rot_matrix()
    sys_add_all_atms = range(len(md_sys.atoms))
    md_sys.mm_atm_coords(-1, m_arb_xyz, False, *sys_add_all_atms)


def _shift_sys(md_sys, radius, radius_buffer=1):
    """
    """
    # enlarge radius so atoms do not clash
    radius += radius_buffer
    rn_pos = agm.points_on_sphere(npoints=1, ndim=3, radius=radius)[0]
    mx_trans = cgt.translation_matrix(rn_pos)
    sys_add_all_atms = range(len(md_sys.atoms))
    md_sys.mm_atm_coords(-1, mx_trans, False, *sys_add_all_atms)


def _create_new_box(md_sys):
    """
    """
    # add new orthogonal box
    md_sys.ts_boxes = []  # delete all previous unneeded boxes
    box_diameter = md_sys.get_system_radius(-1) + 30
    pi_2 = math.pi / 2
    new_box = mdb.Box(boxtype="lattice", ltc_a=box_diameter, ltc_b=box_diameter,
                      ltc_c=box_diameter, ltc_alpha=pi_2, ltc_beta=pi_2,
                      ltc_gamma=pi_2)
    md_sys.ts_boxes.append(new_box)
    md_sys.ts_boxes[0].box_lat2lmp(triclinic=False)


def _check_sys(md_sys, atm_idx_max):
    """
    """
    # check interatomic distances
    md_sys.create_linked_cells(-1)
    close_atms = md_sys.chk_atm_dist(frame_id=-1, min_dist=1.2)
    close_atms = sorted(close_atms)

    # check if new atoms were placed too near to the main system by comparing
    # the atom indices (atom indices lower than those from the original main
    # system are not counted)
    try:
        return close_atms[-1] > atm_idx_max
    except IndexError:
        return True


def sysprep(lmpdat_out, lmpdat_main, lmpdat_add, dcd_add=None, frame_idx=-1):
    """
    Prepare the system for the next docking step.

    Since the kawska-zahn approach adds an agglomerate of molecules each stage,
    the system must be prepared for lammps beforehand. Therefor a sphere with
    radius r1 is created around the main system as well as a sphere with radius
    r2 around the agglomerate to add. Both systems are combined, checked for clashes
    and eventually a new lammps data file is written which can further be
    utilized. If

    Parameters
    ----------
    lmpdat_out : str
        Name of merged lammps data file to write

    lmpdat_main : str
        Name of the lammps data file to add another system to

    lmpdat_add : str
        Name of the lammps data file which is added to the main molecular
        system

    Returns
    -------
    success : bool
        True if successful, False otherwise.

    Writes a new lammps data file called sysprep_out_index

    """
    # read and transpose the main sys to the origin
    main_sys = read_system(lmpdat_main)
    main_sys.transpose_by_cog(-1, [0, 0, 0], copy=False)
    _natoms = len(main_sys.atoms)

    # read and transpose the add sys to the origin
    add_sys = read_system(lmpdat_add, dcd_add, frame_idx)
    add_sys.transpose_by_cog(-1, [0, 0, 0], copy=False)

    # rotate add sys
    _rotate_sys(add_sys)

    # shift add sys to sphere around main sys
    main_sys_radius = main_sys.get_system_radius(-1)
    add_sys_radius = add_sys.get_system_radius(-1)
    kwz_radius = main_sys_radius + add_sys_radius

    _shift_sys(add_sys, kwz_radius)

    # merge both systems
    main_sys.extend_universe(add_sys, mode="merge")
    _create_new_box(main_sys)

    # group atoms by bonds
    main_sys.fetch_molecules_by_bonds()
    main_sys.mols_to_grps()

    success = _check_sys(main_sys, _natoms)

    # write an output lammps data only if everything worked out
    if success is True:
        # write new data file
        main_sys.change_indices(incr=1, mode="increase")
        main_sys.write_lmpdat(lmpdat_out, frame_id=0, title="System ready" +
                              "for docking", cgcmm=True)
    return success


# Quenching ===================================================================#
def quench(lmpcuts, lmpdat_main):
    """
    """
    main_sys = read_system(lmpdat_main)
    natoms_main_sys = len(main_sys.atoms)
    del main_sys

    lmp = lammps()
    pylmp = PyLammps(ptr=lmp)
    lmp.command("log {} append".format(lmpcuts.output_lmplog))

    if lmpcuts.gpu is True:
        lmpcuts.use_gpu(lmp, neigh=True)

    lmp.file(lmpcuts.settings_file)
    # change box type to not be periodic
    #lmp.command("boundary f f f")
    lmpcuts.read_system(lmp)
    lmp.command("fix ic_prevention all momentum 100 linear 1 1 1 angular rescale")
    lmpcuts.dump(lmp)

    if lmpcuts.pc_file is not None:
        lmp.file(lmpcuts.pc_file)

    # intermediate restart files
    lmp.command("restart {0} {1} {1}".format(lmpcuts.logsteps * 10, lmpcuts.output_lmprst))

    # define the atoms that may move during the simulation
    lmp.command("group grp_add_sys id > {}".format(natoms_main_sys))
    lmp.command("group grp_main_sys id <= {}".format(natoms_main_sys))
    #lmp.command("fix freeze grp_main_sys setforce {0} {0} {0}".format(0.0))

    # pre-optimization
    lmp.command("min_style cg")
    lmp.command("min_modify dmax 0.5")
    lmp.command("minimize 1.0e-5 1.0e-8 10000 100000")

    # set an additional push to the added atoms that should be docked
    # towards the origin at (0/0/0)
    # (prevents losing atoms due to being localized outside the cutoff)
    if rank == 0:
        prep_sys = read_system(lmpcuts.input_lmpdat)
        cog = agm.get_cog(prep_sys.ts_coords[-1][natoms_main_sys + 1:])
        cog /= np.linalg.norm(cog, axis=0)  # unit vector
        cog *= -1
    else:
        cog = None

    cog = comm.bcast(cog, 0)
    lmp.command(("fix force grp_add_sys addforce {0} {1} {2} every 10000").format(*cog))

    # barostatting, thermostatting
    lmp.command("fix integrator grp_add_sys nvt temp {0} {1} 0.1".format(lmpcuts.tstart, lmpcuts.tstop))
    lmp.command("run {}".format(lmpcuts.runsteps))
    # remove pushing force
    lmp.command("unfix force")

    # post-optimization
    lmp.command("min_style quickmin")
    lmp.command("minimize 1.0e-5 1.0e-8 10000 100000")

    # write restart file
    lmpcuts.unfix_undump(pylmp, lmp)
    lmp.command("reset_timestep 0")
    lmp.command("write_restart {}".format(lmpcuts.output_lmprst))
    lmp.command("clear")
    lmp.close()

    # check aggregate
    if rank == 0:
        quench_sys = read_system(lmpcuts.input_lmpdat, dcd=lmpcuts.output_dcd)
        quench_success = check_aggregate(quench_sys)
    else:
        quench_success = False

    quench_success = comm.bcast(quench_success, 0)
    return quench_success


# Annealing ===================================================================#
def relax_box(lmpcuts):
    """
    """
    lmp = lammps()
    pylmp = PyLammps(ptr=lmp)
    lmp.command("log {} append".format(lmpcuts.output_lmplog))

    if lmpcuts.gpu is True:
        lmpcuts.use_gpu(lmp, neigh=True)

    lmp.file(lmpcuts.settings_file)
    lmpcuts.read_system(lmp)
    lmpcuts.thermo(lmp)
    lmp.command("fix ic_prevention all momentum 100 linear 1 1 1 angular rescale")
    lmpcuts.dump(lmp)

    if lmpcuts.pc_file is not None:
        lmp.file(lmpcuts.pc_file)

    # minimize cut box if not done already
    if lmpcuts.input_lmprst is None or os.path.isfile(lmpcuts.input_lmprst):
        lmpcuts.minimize(lmp, style="cg")

    # barostatting, thermostatting
    lmpcuts.berendsen(lmp)

    #lmp.command("fix nose_hoover grp_add_sys nvt temp {0} {1} 0.1".format(lmpcuts.tstart, lmpcuts.tstop))

    lmp.command("run {}".format(lmpcuts.runsteps))
    lmpcuts.unfix_undump(pylmp, lmp)
    lmp.command("reset_timestep 0")
    lmp.command("write_restart {}".format(lmpcuts.output_lmprst))
    lmp.command("clear")
    lmp.close()


def _fix_indent_ids(radii, cogs, group_name, scale_start=1, scale_stop=12):
    """
    Create the strings that are used as IDs for the lammps fix indent command.

    Parameters
    ----------
    group_name : str
        Name for the group, e.g. molecule or atom
    radii : list of floats
        Radii of all molecules to create a sphere around
    cogs : list of floats
        Centers of geometry of all molecules, i.e. center of each sphere

    Returns
    -------
    indent_fixes : list of lists
        Contains all ready up fixes for the lammps indent fix

    """
    if rank == 0:
        print(radii)

    fix_str = "fix {0} all indent 1 sphere {c[0]} {c[1]} {c[2]} {1} side out"
    all_indent_fixes = []

    for scaling_factor in xrange(scale_start, scale_stop):
        indent_fixes = []
        # set fixes for void_indent
        for idx, (radius, cog) in enumerate(zip(radii, cogs)):
            fix_name = "{}_{}_fix_indent".format(group_name, idx)
            grown_radius = radius * 0.1 * scaling_factor
            cur_fix_str = fix_str.format(fix_name, grown_radius, c=cog)
            indent_fixes.append(cur_fix_str)
        all_indent_fixes.append(indent_fixes)

    return all_indent_fixes


def _lmp_indent(lmp, indents, runsteps, keep_last_fixes=False):
    """
    Perform indentation runs to grow spheres around certain molecules.

    Parameters
    ----------
    lmp : lammps-object
    indents : list of lists of str
    runsteps : int
    keep_last_fixes : bool

    """
    # make sphere grow (i.e. create and delete fixes for fix indent)
    #pdb.set_trace()

    for indent_fixes in indents[:-1]:
        for indent_fix in indent_fixes:
            lmp.command(indent_fix)
        lmp.command("run {}".format(runsteps))

        for indent_fix in indent_fixes:
            lmp.command("unfix {}".format(indent_fix.split()[1]))

    # last run for sphere growth
    for indent_fix in indents[-1]:
        lmp.command(indent_fix)
    lmp.command("run {}".format(runsteps))

    # unfix last fixes
    if keep_last_fixes is False:
        for indent_fix in indents[-1]:
            lmp.command("unfix {}".format(indent_fix.split()[1]))


def _check_imaginary_clashes(sys_both, sys_a, sys_b, dcd_b):
    """
    Check for imaginary clashes between atoms of two systems before combining them.

    Parameters
    ----------
    sys_both : Universe
    sys_a : Universe
    sys_b : Universe

    Returns
    -------
    close_atms_a : list of int
        Atom ids of atoms that need additional spheres because they clash with
        solvent atoms otherwise

    """
    #sys_both.reset_cells()
    # read latest solvent coordinates and boxes
    sys_b.import_dcd(dcd_b)
    sys_b.read_frames(frame=-2, to_frame=-1)
    sys_b.close_dcd()

    # concatenate solute and latest solvent coordinates
    sys_both.ts_coords.append(np.concatenate((sys_a.ts_coords[-1], sys_b.ts_coords[-1])))
    sys_both.ts_boxes = sys_b.ts_boxes
    sys_both.create_linked_cells(-1, rcut_a=2, rcut_b=2, rcut_c=2)
    close_atms_ab = sys_both.chk_atm_dist(-1, min_dist=1.0, exclude_same_molecule=True)
    # only get close atoms from system a
    close_atms_a = [i for i in close_atms_ab if i <= len(sys_a.atoms)]

    sys_b.reset_cells()
    sys_both.reset_cells()
    return close_atms_a


def create_voids(lmpcuts, lmpdat_solvate, dcd_solvate):
    """
    """
    # solvate with molecule radii and cogs
    solvate_sys = read_system(lmpdat_solvate, dcd_solvate)
    radii_mol, cogs_mol = _molecules_radii(solvate_sys)
    indent_strs = _fix_indent_ids(radii_mol, cogs_mol, "molecule", scale_start=10, scale_stop=15)

    # load lammps and lammps settings
    lmp = lammps()
    pylmp = PyLammps(ptr=lmp)
    lmp.command("log {} append".format(lmpcuts.output_lmplog))

    if lmpcuts.gpu is True:
        lmpcuts.use_gpu(lmp)

    lmp.file(lmpcuts.settings_file)
    lmpcuts.read_system(lmp)
    lmpcuts.thermo(lmp)
    lmp.command("fix ic_prevention all momentum 100 linear 1 1 1 angular rescale")
    lmpcuts.dump(lmp)

    if lmpcuts.pc_file is not None:
        lmp.file(lmpcuts.pc_file)

    lmpcuts.berendsen(lmp)
    #lmp.command("unfix integrator")
    #lmp.command("fix integrator all nve/limit 1.0")
    _lmp_indent(lmp, indent_strs, lmpcuts.runsteps, keep_last_fixes=True)

    # solution system
    solvent_sys = read_system(lmpcuts.input_lmpdat)
    solution_sys = merge_sys(solvate_sys, solvent_sys)
    solution_sys.reset_cells()

    # check solvate atoms with too close contacts to solvent atoms
    close_atoms = _check_imaginary_clashes(solution_sys, solvate_sys, solvent_sys, lmpcuts.output_dcd)
    cogs_atoms = [solvate_sys.ts_coords[-1][i] for i in close_atoms]
    radii_atoms = [mde.elements_mass_radii[round(solvate_sys.atm_types[solvate_sys.atoms[i].atm_key].weigh, 1)] for i in close_atoms]

    # gather close atoms and remember which were close
    all_close_atoms = []
    all_radii_atoms = []
    all_cogs_atoms = []

    all_close_atoms.extend(close_atoms)
    all_radii_atoms.extend(radii_atoms)
    all_cogs_atoms.extend(cogs_atoms)

    factor_start = 10
    factor_stop = factor_start + 1

    if close_atoms != []:
        # move solvent molecules away from close solvate atoms
        for _ in xrange(5):
            indent_strs = _fix_indent_ids(all_radii_atoms, all_cogs_atoms, "atom", scale_start=factor_start, scale_stop=factor_stop)
            _lmp_indent(lmp, indent_strs, lmpcuts.runsteps, keep_last_fixes=False)
            close_atoms = _check_imaginary_clashes(solution_sys, solvate_sys, solvent_sys, lmpcuts.output_dcd)

            # add new close atoms to present ones or stop indenting
            if close_atoms == []:
                break

            # add further close atoms to present ones
            for atm_idx in close_atoms:
                print(close_atoms)
                if atm_idx not in all_close_atoms:
                    all_close_atoms.append(atm_idx)
                    all_radii_atoms.append(mde.elements_mass_radii[round(solvate_sys.atm_types[solvate_sys.atoms[i].atm_key].weigh, 1)])
                    all_cogs_atoms.append(solvate_sys.ts_coords[-1][atm_idx])

            # dynamically grow sphere around atoms
            factor_start += 1
            factor_stop = factor_start + 1

    lmpcuts.unfix_undump(pylmp, lmp)
    lmp.command("write_restart {}".format(lmpcuts.output_lmprst))
    lmp.command("clear")
    lmp.close()

    del (solvate_sys, solvent_sys, solution_sys)

    # combine solvent and solute and write a lammps data file
    if rank == 0:
        _write_data(lmpcuts.output_lmpdat, lmpdat_solvate, lmpdat_b=lmpcuts.input_lmpdat, dcd_a=dcd_solvate, dcd_b=lmpcuts.output_dcd)

    return close_atoms == []


def anneal_productive(lmpcuts):
    """
    Use nose-hoover baro- and thermostat for productive annealing run.
    """
    lmp = lammps()
    pylmp = PyLammps(ptr=lmp)
    lmp.command("log {} append".format(lmpcuts.output_lmplog))

    if lmpcuts.gpu is True:
        lmpcuts.use_gpu(lmp, neigh=True)

    lmp.file(lmpcuts.settings_file)
    lmpcuts.read_system(lmp)
    lmpcuts.thermo(lmp)
    lmp.command("fix ic_prevention all momentum 100 linear 1 1 1 angular rescale")
    lmpcuts.dump(lmp)

    if lmpcuts.pc_file is not None:
        lmp.file(lmpcuts.pc_file)

    lmp.command("fix barostat all npt temp {} {} 0.5 aniso {} {} 50".format(lmpcuts.tstart, lmpcuts.tstop, lmpcuts.pstart, lmpcuts.pstop))
    lmp.command("run {}".format(lmpcuts.runsteps))
    lmpcuts.unfix_undump(pylmp, lmp)
    lmp.command("reset_timestep 0")
    lmp.command("write_restart {}".format(lmpcuts.output_lmprst))
    lmp.command("clear")
    lmp.close()


def requench(lmpcuts):
    """
    """
    lmp = lammps()
    pylmp = PyLammps(ptr=lammps)
    lmp.file(lmpcuts.settings_file)
    lmpcuts.read_system(lmp)
    lmpcuts.thermo(lmp)
    lmp.command("fix ic_prevention all momentum 100 linear 1 1 1 angular rescale")
    lmpcuts.dump(lmp)

    if lmpcuts.pc_file is not None:
        lmp.file(lmpcuts.pc_file)

    lmpcuts.minimize(lmp, style="cg", box_relax=False)
    lmpcuts.unfix_undump(pylmp, lmp)
    lmp.command("reset_timestep 0")
    lmp.command("write_restart {}".format(lmpcuts.output_lmprst))
    lmp.command("clear")
    lmp.close()
