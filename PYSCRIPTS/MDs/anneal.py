from lammps import lammps
from mpi4py import MPI
import argparse

#==============================================================================#
# Setup MPI
#==============================================================================#

comm = MPI.COMM_WORLD
size = comm.Get_size()  # number of processes in communicator
rank = comm.Get_rank()  # process' id(s) within a communicator


#==============================================================================#
# Helping functions
#==============================================================================#
class bcolors:
    header = '\033[95m'
    blue = '\033[94m'
    green = '\033[92m'
    yellow = '\033[93m'
    red = '\033[91m'
    endc = '\033[0m'
    bold = '\033[1m'
    underline = '\033[4m'


def anneal(lmpdat, tmax, steps, save_step, cycles, output, settings_file=None):
    """
    """

    thermargs = ["step", "temp", "pe", "eangle", "edihed", "eimp", "evdwl",
                 "ecoul", "ebond", "enthalpy"]

    lmp = lammps(cmdargs=["-screen", "lmp_out.txt"])
    lmp.command("log {}.lmplog".format(output))

    # load default settings if no file with settings was supplied
    if settings_file is None:
        lmp.command("units metal")
        lmp.command("boundary p p p")
        lmp.command("dimension 3")
        lmp.command("atom_style full")
        lmp.command("pair_style lj/cut/coul/cut 30.0")
        lmp.command("bond_style harmonic")
        lmp.command("angle_style harmonic")
        lmp.command("dihedral_style charmm")
        lmp.command("improper_style cvff")
        lmp.command("special_bonds amber")
        lmp.command("timestep 0.001")

    lmp.command("read_data {}".format(lmpdat))
    # make lammps calculate the value of the entity (bond, angle, dihedral)
    lmp.command("thermo_style custom " + " ".join(thermargs))
    lmp.command("thermo_modify lost warn flush yes")
    lmp.command("thermo {}".format(save_step))
    lmp.command("fix MOMENT all momentum 100 linear 1 1 1 angular")
    lmp.command("dump DUMP all dcd {} {}.dcd".format(save_step, output))
    lmp.command("fix NVE all nve")

    # heat up
    #lmp.command("velocity all create {} 8675309 mom yes rot yes dist gaussian".format(300))
    for _ in xrange(cycles):
        print(bcolors.red + "Annealing" + bcolors.endc)
        lmp.command("fix TFIX all langevin {} {} 100 24601".format(0.0, tmax))

        try:
            lmp.command("run {}".format(steps))
        except:
            print("***Error: Simulation crashed (annealing)! Force constants too high?")
            MPI.COMM_WORLD.Abort()

        print(bcolors.yellow + "Quenching" + bcolors.endc)
        lmp.command("fix TFIX all langevin {} {} 100 24601".format(tmax, 0.0))

        try:
            lmp.command("run {}".format(steps*2))
        except:
            print("***Error: Simulation crashed (annealing)! Force constants too high?")
            MPI.COMM_WORLD.Abort()

        print(bcolors.green + "Minimization" + bcolors.endc)
        try:
            lmp.command("minimize 1e-6 1e-9 2000000 100000")
        except:
            print("***Error:  Simulation crashed (minimization)!")
            MPI.COMM_WORLD.Abort()


if __name__ == "__main__":
    PARSER = argparse.ArgumentParser()
    PARSER.add_argument("lmpdat", help="Lammps' data-file")
    PARSER.add_argument("-tmax", type=float, default=900)
    PARSER.add_argument("-steps", type=int, default=1000000)
    PARSER.add_argument("-save_step", type=int, default=1000)
    PARSER.add_argument("-cycles", type=int, default=5)
    PARSER.add_argument("-set", default=None)
    PARSER.add_argument("-out", default="DEFAULTNAME", help="Name of energy files.")
    ARGS = PARSER.parse_args()
    anneal(ARGS.lmpdat, ARGS.tmax, ARGS.steps, ARGS.save_step, ARGS.cycles, ARGS.out)
