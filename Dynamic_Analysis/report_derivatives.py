from shutil import ReadError
import numpy as np
import json
import matplotlib.pyplot as plt
from scipy.linalg import eig, block_diag
from control import ctrb

from os import mkdir, rmdir, walk, remove, listdir
from os.path import exists as path_exists

from eigensolver import Solver

if __name__ == "__main__":
    
    # change plot text parameters
    plt.rcParams["font.family"] = "Serif"
    plt.rcParams["font.size"] = 8.0
    plt.rcParams["axes.labelsize"] = 8.0
    # plt.rcParams['lines.linewidth'] = 1.0
    # plt.rcParams["xtick.minor.visible"] = True
    # plt.rcParams["ytick.minor.visible"] = True
    # plt.rcParams["xtick.direction"] = plt.rcParams["ytick.direction"] = "in"
    # plt.rcParams["xtick.bottom"] = plt.rcParams["xtick.top"] = True
    # plt.rcParams["ytick.left"] = plt.rcParams["ytick.right"] = True
    # plt.rcParams["xtick.major.width"] = plt.rcParams["ytick.major.width"] = 1.0
    # plt.rcParams["xtick.minor.width"] = plt.rcParams["ytick.minor.width"] = 1.0
    # plt.rcParams["xtick.major.size"] = plt.rcParams["ytick.major.size"] = 5.0
    # plt.rcParams["xtick.minor.size"] = plt.rcParams["ytick.minor.size"] = 2.5
    plt.rcParams["mathtext.fontset"] = "dejavuserif"
    plt.rcParams['figure.dpi'] = 300.0
    plt.rcParams['figure.constrained_layout.use'] = True

    # Convair 990 https://www.nasa.gov/centers/dryden/pdf/87810main_H-693.pdf
    # some aircraft (plots, not values) in "" by Edward Seckel

    # Order:
    # Cessna book
    # Our F16
    # Teper
    # Heffley
    # McRuer
    # Blakelock
    # Phillips
    run_files = [
    ### Class I -- small light
    "Cessna_172.json",
    "navion.json",
    "F_2B.json", #"VZ_4.json",
    ### Class II -- medium-weight, low-to-medium maneuverability
    "X_15.json", "HL_10.json", "Lockheed_Jetstar.json", "Convair_880M.json", 
    "F_105B.json",
    "C_47.json", #"XC_142.json",
    ### Class III -- large, heavy, low-to-medium maneuverability
    "boeing_747.json", "C_5A.json", "XB_70A.json",
    "DC_8.json",
    # "9_8_2.json",
    ### Class IV -- high-maneuverability
    "F_16.json",
    "NT_33A.json", "F_104A.json", "F_4C.json",
    "A_7A.json", "A_4D.json",
    "F_94A.json", "F_15.json",
    ]
    run_files = ["aircraft_database/" + i for i in run_files]
    num_craft = len(run_files)
    eigensolved = [0.0] * num_craft
    eigensolvedV = [0.0] * num_craft
    # check directory for plots, remove if there
    delta_folder = "phasor_plots/phasor_delta"
    if path_exists(delta_folder):
        for filename in listdir(delta_folder):
            remove(delta_folder + "/" + filename)
    else:
        mkdir(delta_folder)
    for i in range(num_craft):
        eigensolved[i] = Solver(run_files[i],report=False)

    print("reporting CW vs CYrbar")
    for i in range(num_craft):
        ei = eigensolved[i]
        CY_rbreve = ei.CY_rbar*ei.g*ei.bw/2.0/ei.vo**2.0
        print(ei.aircraft_name, ei.CW, CY_rbreve)

    print("\n\n")
    print("reporting CW vs CLqbar")
    for i in range(num_craft):
        ei = eigensolved[i]
        CL_qbreve = ei.CL_qbar*ei.g*ei.cwbar/2.0/ei.vo**2.0
        print(ei.aircraft_name, ei.CW, CL_qbreve)