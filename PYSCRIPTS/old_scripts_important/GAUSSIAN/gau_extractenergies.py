from __future__ import print_function
import sys
import os
#import re
#import subprocess32 as sp32
#import shutil as sl
#import psutil
#import math
from natsort import natsorted

script, dir_logfiles = sys.argv

all_files = os.listdir(dir_logfiles)
all_files = natsorted(all_files)

# remove files ending with gau
for index, cur_file in enumerate(all_files):
    if ".gau" in cur_file:
        del all_files[index]

# file with energies
energy_file = "{}_energies.dat".format(dir_logfiles)  # "quantum.dat"

with open(energy_file, "w") as en_out:
    for cur_file in all_files:
        # make sure we reade the right files
        if cur_file.endswith(".log") or cur_file.endswith(".out"):
            with open(dir_logfiles+"/"+cur_file) as en_in:
                for line in en_in:

                    # get angle from simulation title
                    if "Rotated urea-group by angle:" in line:
                        line = line.split()

                        if "----------" in en_in.next():
                            angle = line[-1]
                            en_out.write("{:>5s}".format(angle))
                        else:
                            pass

                    # get line with energy
                    elif "SCF Done" in line:
                        line = line.split()
                        energy = line[4]
                        en_out.write("{:>20}\n".format(energy))
                    else:
                        pass
        else:
            pass
