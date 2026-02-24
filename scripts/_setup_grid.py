#!/usr/local/bin/python
# Authors:
#          Mathieu Renzo <mrenzo@arizona.edu>
#
# Keywords: files

# Copyright (C) 2026 Mathieu Renzo

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or (at
# your option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see http://www.gnu.org/licenses/.

__author__ = ["Mathieu Renzo <mrenzo@arizona.edu>"]

import os
from MESAreader import get_src_col
from short_and_iso import prepare_SNEC_input
# from utilsLib import checkFolder, gitPush


def get_num_cells(INPUT_PROFILE):
    """ get number of mesh points for grid """
    src, col = get_src_col(INPUT_PROFILE)
    IMAX = len(src[:, col.index("zones")])
    return int(IMAX)


def setup_one_model(INPUT_SHORT, INPUT_COMP_FILE,
                    parameters_template,
                    OUTDIR="/tmp/Data",
                    FINAL_ENERGY="0.0d0",
                    NI_MASS="0.0d0",
                    ):
    """ Assumes OUTDIR exists"""
    # IMAX = get_num_cells(INPUT_PROFILE)

    with open(parameters_template, "r") as p:
        parameters_file = p.read()

    parameters_file = parameters_file.replace("INPUT_PROFILE", INPUT_SHORT)
    parameters_file = parameters_file.replace("INPUT_COMP_FILE", INPUT_COMP_FILE)
    parameters_file = parameters_file.replace("OUTDIR", OUTDIR)
    parameters_file = parameters_file.replace("FINAL_ENERGY", FINAL_ENERGY)
    parameters_file = parameters_file.replace("IMAX", f"{1000:d}")
    parameters_file = parameters_file.replace("NI_MASS", NI_MASS)

    with open(OUTDIR+'/parameters', 'w') as p:
        p.write(parameters_file)

    print("")
    print("Done, setup in "+OUTDIR)
    print("")


if __name__ == "__main__":
    print("")
    print("Setting up SNEC grid!")
    print("")

    # SNEC_ROOT="/home/u20/mrenzo/codes/SNEC/" # cluster_ua
    SNEC_ROOT = "/home/mrenzo/Documents/Research/codes/SNEC-1.01/" # ua_w
    OUTDIR_ROOT='/home/mrenzo/Runs/SNEC_grid/long_runs/15Msun/' # s30VdJNL_0.33/'
    # if below is None use SNEC provided model
    INPUT_MESA_FILE = None # "/home/mrenzo/Runs/LMXRB/CCSN_progenitors/s30VdJNL_0.33_onset_cc.data" # None #
    # final_energies need to be strings including a decimal point and d for exponential notation, or SNEC will complain
    final_energies = ["1.0d51", "0.0d0"]


    # check if folder exists and user wants to erase it
    if os.path.isdir(OUTDIR_ROOT):
        print("Folder "+OUTDIR_ROOT+" exists!")
        go_on = input("Erase it and re-create it? [Y,y]")
        if go_on.lower() != "y":
            sys.exit()
        else:
            print("Ok, cleaning up...")
            os.system("rm -rf "+OUTDIR_ROOT)
            print("...removed pre-existing "+OUTDIR_ROOT)

    for FINAL_ENERGY in final_energies:
        print("working on "+f"{FINAL_ENERGY}")

        OUTDIR = OUTDIR_ROOT+"/bomb_final_e_"+str(FINAL_ENERGY)+"/"

        # make directory and SNEC output directory
        os.system("mkdir -p "+OUTDIR+"/Data/")


        if INPUT_MESA_FILE:
            INPUT_SHORT, INPUT_COMP_FILE = prepare_SNEC_input(INPUT_MESA_FILE, OUTDIR,
                                                              SNEC_ROOT+'/scripts/')
        else:
            # Use SNEC provided progenitor
            INPUT_SHORT = SNEC_ROOT+"/profiles/15Msol_RSG.short"
            INPUT_COMP_FILE = SNEC_ROOT+"/profiles/15Msol_RSG.iso.dat"

        setup_one_model(INPUT_SHORT, INPUT_COMP_FILE,
                        parameters_template=SNEC_ROOT+"/parameters_template",
                        OUTDIR=OUTDIR,
                        FINAL_ENERGY=FINAL_ENERGY,
                        )

        # copy SNEC executable -- assumes code is already compiled
        os.system("cp "+SNEC_ROOT+"snec "+OUTDIR)

        # make submission script
        with open(OUTDIR_ROOT+"run_all.sh", "w") as f:
            run_all_bash = """#!/bin/bash

for dir in */; do
    cd \"$dir\"
    ./snec >> ./Data/output.txt 2>&1 &
    cd ..
done
pgrep -a snec
"""
            # oneliner -- no terminal ouput
            # "#!/bin/bash\n\nfor dir in */; do\n\tcd \"$dir\"\n\t./snec &\n\tcd ..\ndone"
            f.write(run_all_bash)
        # give -x permissions
        os.chmod(OUTDIR_ROOT+"run_all.sh", 0o755)

        print ("Done grid setup in")
        print("")
        print(OUTDIR_ROOT)
        print("")
