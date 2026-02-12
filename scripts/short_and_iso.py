#!/usr/bin/env python

import os
import sys


def prepare_SNEC_input(MESA_core_collapse_pfile,
                       where_to_save="/tmp/",
                       where_the_scripts_are="/home/mrenzo/Documents/Research/codes/SNEC-1.01/scripts/"):
    os.system(
        "python "
        + where_the_scripts_are
        + "mesa_to_GR1D.py "
        + MESA_core_collapse_pfile
        + " "
        + where_to_save+'/'
        + MESA_core_collapse_pfile
        + ".short"
    )
    os.system(
        "python "
        + where_the_scripts_are
        + "MESA_isotopes.py "
        + MESA_core_collapse_pfile
        + " "
        + where_to_save+'/'
        + MESA_core_collapse_pfile
        + ".iso.dat"
    )
    # os.system("ls "+where_to_save)
    return (where_to_save+'/' + MESA_core_collapse_pfile + ".short",
            where_to_save+'/' + MESA_core_collapse_pfile + ".iso.dat")


if __name__ == "__main__":
    pfile, compfile = prepare_SNEC_input(sys.argv[1], sys.argv[2])
