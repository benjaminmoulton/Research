import eigensolver as es
import numpy as np
from matplotlib import pyplot as plt

def print_hq(hq):
    if type(hq) != str:
        print(" \& {:> 8.4f}".format(hq),end="")
    else:
        print(" \&  {:^7}".format(hq),end="")


if __name__ == "__main__":
    # plt.rcParams['text.usetex'] = True
    # plt.rcParams["font.family"] = "Times New Roman"
    # plt.rcParams["font.size"] = 16
    # import matplotlib.colors.ColorConverter().to_rgba as rgba

    # es.Solver("aircraft_database/F16_bolander.json",report=True)
    # es.Solver("aircraft_database/NT_33A.json",report=True)


    # Convair 990 https://www.nasa.gov/centers/dryden/pdf/87810main_H-693.pdf

    run_files = [
    # Class I
    "Cessna_172.json",
    "navion.json",
    # Class II
    "X_15.json", "HL_10.json", "Lockheed_Jetstar.json", "Convair_880M.json", 
    "F_105B.json", "DC_8.json",
    # Class III
    "boeing_747.json", "C_5A.json", "XB_70A.json",
    # Class IV
    "F16_bolander.json",
    "NT_33A.json", "F_104A.json", "F_4C.json",
    "A_7A.json", "A_4D.json"
    ]
    run_files = ["aircraft_database/" + i for i in run_files]
    num_craft = len(run_files)
    eigensolved = [0.0] * num_craft
    for i in range(num_craft):
        eigensolved[i] = es.Solver(run_files[i],report=False)

    # create dictionaries for pretty titling of print out
    modenames = {
        "sp" : "Short Period",
        "ph" : "Phugoid",
        "sl" : "Spiral",
        "ro" : "Roll",
        "dr" : "Dutch Roll"
    }
    propertynames = {
        "zt" : "Damping Ratio",
        "wn" : "Natural Frequency, rad/s",
        "Sg" : "Damping Rate, 1/s"
    }

    # report handling qualities
    max_name = max([len(plane.aircraft_name) for plane in eigensolved]) + 2
    side = ["lon"] * 4 + ["lat"] * 4
    mode = ["sp"] * 2 + ["ph"] * 2 + ["ro","sl"] + ["dr"] * 2
    haqu = ["wn","zt"] * 2 + ["Sg"] * 2 + ["wn","zt"]
    column_header = ["wt/ft","linear","approx","buckham"]
    for i in range(len(side)):
        mo = mode[i]
        si = side[i]
        hq = haqu[i]
        print(modenames[mo],propertynames[hq])

        print("{:^{}}".format("aircraft",max_name),end="")
        for i in range(4):
            print_hq(column_header[i])
        print(" \\\\")


        for j in range(num_craft):
            print("{:^{}}".format(eigensolved[j].aircraft_name,max_name),\
                end="")

            # print actual value
            if mo=="sp" and eigensolved[j].split_short_period: print_hq("---")
            elif mo=="ph" and eigensolved[j].split_phugoid: print_hq("---")
            else: print_hq(eigensolved[j].HQ[hq+"_"+mo])
            if mo=="ro" or mo=="sl":
                print_hq(eigensolved[j].b[si][hq][eigensolved[j].b[si][mo]])
            else:
                print_hq(eigensolved[j].b[si][hq][eigensolved[j].b[si][mo]][0])
            print_hq(eigensolved[j].b[mo+hq])
            print_hq(eigensolved[j].b[mo+hq+"hq"])

            print(" \\\\")
        print()