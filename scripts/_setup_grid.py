# script to setup grid of SNEC runs
import os
from MESAreader import get_src_col
from short_and_iso import prepare_SNEC_input

global SNEC_ROOT
SNEC_ROOT="/home/u20/mrenzo/codes/SNEC/" # cluster_ua


def get_num_cells(INPUT_PROFILE):
    """ get number of mesh points for grid """
    src, col = get_src_col(INPUT_PROFILE)
    IMAX = len(src[:, col.index("zones")])
    return int(IMAX)


def setup_one_model(INPUT_PROFILE, INPUT_COMP_FILE,
                    OUTDIR=SNEC_ROOT+"Data",
                    FINAL_ENERGY=0,
                    parameters_template=SNEC_ROOT+"parameters"
                    ):
    os.system("mkdir -p "+OUTDIR)

    IMAX = get_num_cells(INPUT_PROFILE)

    with open(parameters_template, "r") as p:
        parameters_file = p.read()

    parameters_file = parameters_file.replace("OUTDIR", OUTDIR)
    parameters_file = parameters_file.replace("INPUT_PROFILE", INPUT_PROFILE)
    parameters_file = parameters_file.replace("INPUT_COMP_FILE", INPUT_COMP_FILE)
    parameters_file = parameters_file.replace("FINAL_ENERGY", FINAL_ENERGY)
    parameters_file = parameters_file.replace("IMAX", IMAX)
    parameters_file = parameters_file.replace("NI_MASS", NI_MASS)

    with open(OUTDIR+'/parameters', 'w') as p:
        p.write(parameters_file)

    print("")
    print("Done, setup in "+OUTDIR)


if __name__ == "__main":
    print("")
    print("Setting up SNEC grid!")
    print("")

    SNEC_ROOT = "/home/mrenzo/Documents/Research/codes/SNEC-1.01/"
    OUTDIR_ROOT='/tmp/'
    INPUT_MESA_FILE = "/home/mrenzo/Runs/LMXRB/CCSN_progenitors/s30VdJNL_0.33_onset_cc.data"


    final_energies = [1.0d51, # 1.0d50, 1.0d49, 1.0d48, 1d47, 1d46, 1d45, 0.0d0, -1d45,
                      -1d46]

    for FINAL_ENERGY in final_energies:
        OUTDIR = OUTDIR_ROOT+"/bomb_final_e_"+str(FINAL_ENERGY)+"/"
        INPUT_PROFILE, INPUT_COMP_FILE = prepare_SNEC_input(INPUT_MESA_FILE, OUTDIR,
                                                            SNEC_ROOT+'/scripts/')
        setup_one_model(INPUT_PROFILE, INPUT_COMP_FILE,
                        OUTDIR=OUTDIR,
                        FINAL_ENERGY=FINAL_ENERGY,
                        parameters_template=SNEC_ROOT+"/parameters"
                        )
