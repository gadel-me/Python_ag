#!/usr/bin/env python

# from mpi4py import MPI
from lammps import lammps, PyLammps

# FOR FURTHER READING:
#       http://lammps.sandia.gov/doc/tutorial_pylammps.html#pylammps-tutorial

lmp = lammps()
L = PyLammps(ptr=lmp)

# process one line at a time
with open("OUT.in", "r") as f_in:
    for line in f_in:
        lmp.command(line)

# lmp.file("OUT.in")
# me = MPI.COMM_WORLD.Get_rank()
# nprocs = MPI.COMM_WORLD.Get_size()
# print "Proc %d out of %d procs has" % (me,nprocs),lmp
# MPI.Finalize()  # shutdown mpi properly
