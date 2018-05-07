#!/usr/bin/env python
from __future__ import print_function
import sys
import os
import subprocess32 as sp32
import re
import shutil as sl

script, g09_start, gaussian_files_dir = sys.argv

start_g09_temp = g09_start + ".temp"
# run script "s_g09.py" for every gaussian-input-file
for cur_file in os.listdir(gaussian_files_dir):

    # only process gaussian-input files
    if cur_file.endswith(".gau"):
        gaussian_in = "{}/{}".format(gaussian_files_dir, cur_file)

        # change stuff in the s_g09-script
        with open(g09_start, "r") as s_g09, open(start_g09_temp, "w") as s_g09_temp:
            for line in s_g09:

                # change job-name
                if "#SBATCH --job-name" in line:
                    job_name = re.sub('[!.gau]', '', gaussian_in)
                    line = "#SBATCH --job-name={}\n".format(job_name)

                s_g09_temp.write(line)

        # change script name to previous
        sl.move(start_g09_temp, g09_start)
        # run s_g09-script
        sp32.call(["sbatch", g09_start, gaussian_in], bufsize=-1)
        #sp32.call(["python", g09_start, gaussian_in], bufsize=-1)
