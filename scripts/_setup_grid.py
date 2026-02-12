# script to setup grid of SNEC runs
import os
from MESAreader import get_src_col
from short_and_iso import prepare_SNEC_input


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
    OUTDIR_ROOT='/tmp/SNEC_test/'
    INPUT_MESA_FILE = "/home/mrenzo/Runs/LMXRB/CCSN_progenitors/s30VdJNL_0.33_onset_cc.data"
    final_energies = ["1.0d51", # 1.0d50, 1.0d49, 1.0d48, 1d47, 1d46, 1d45, 0.0d0, -1d45,
                      "-1d46"]

    # check if folder exists and user wants to erase it
    if os.path.isdir(OUTDIR_ROOT):
        print("Folder "+OUTDIR_ROOT+" exists!")
        go_on = input("Erase it and re-create it? [Y,y]")
        if go_on.lower() != "y":
            exit()
        else:
            print("Ok, cleaning up...")
            os.system("rm -rf "+OUTDIR_ROOT)
            print("...removed pre-existing "+OUTDIR_ROOT)

    for FINAL_ENERGY in final_energies:
        print("working on "+f"{FINAL_ENERGY}")

        OUTDIR = OUTDIR_ROOT+"/bomb_final_e_"+str(FINAL_ENERGY)+"/"

        # make directory and SNEC output directory
        os.system("mkdir -p "+OUTDIR+"/Data/")

        INPUT_SHORT, INPUT_COMP_FILE = prepare_SNEC_input(INPUT_MESA_FILE, OUTDIR,
                                                          SNEC_ROOT+'/scripts/')

        setup_one_model(INPUT_SHORT, INPUT_COMP_FILE,
                        parameters_template=SNEC_ROOT+"/parameters",
                        OUTDIR=OUTDIR,
                        FINAL_ENERGY=FINAL_ENERGY,
                        )
        print ("Done grid setup in")
        print("")
        print(OUTDIR_ROOT)
        print("")
