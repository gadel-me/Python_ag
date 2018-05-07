#!/usr/bin/env python
from __future__ import print_function
import sys
import os

script, gaussian_files_dir = sys.argv

for cur_file in os.listdir(gaussian_files_dir):
    success = False
    # check only output-files
    if cur_file.endswith(".log"):
        # full/relative path of file needed
        cur_file = gaussian_files_dir + cur_file
        with open(cur_file, "r") as s_g09_out:
            for line in s_g09_out:
                # check if phrase in file
                if "Normal termination of Gaussian 09" in line:
                    success = True

        if not success:
            print("Something went wrong with {}".format(cur_file))
