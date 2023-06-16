import eigensolver as es
import numpy as np
from matplotlib import pyplot as plt
import json

if __name__ == "__main__":
    plt.rcParams['text.usetex'] = True
    plt.rcParams["font.family"] = "Times New Roman"
    plt.rcParams["font.size"] = 14
    plt.rcParams['lines.linewidth'] = 0.85
    plt.rcParams["xtick.minor.visible"] = True
    plt.rcParams["ytick.minor.visible"] = True
    plt.rcParams["xtick.direction"] = plt.rcParams["ytick.direction"] = "in"
    plt.rcParams["xtick.bottom"] = plt.rcParams["xtick.top"] = True
    plt.rcParams["ytick.left"] = plt.rcParams["ytick.right"] = True
    plt.rcParams["xtick.major.width"] = plt.rcParams["ytick.major.width"] = 1.0
    plt.rcParams["xtick.minor.width"] = plt.rcParams["ytick.minor.width"] = 1.0
    plt.rcParams["xtick.major.size"] = plt.rcParams["ytick.major.size"] = 10.0
    plt.rcParams["xtick.minor.size"] = plt.rcParams["ytick.minor.size"] = 5.0
    plt.rcParams["axes.labelpad"] = "10"


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
    "F_16.json",
    "NT_33A.json", "F_104A.json", "F_4C.json",
    "A_7A.json", "A_4D.json"
    ]
    craft_class = [
        "I","I",
        "II","II","II","II","II","II",
        "III","III","III",
        "IV","IV","IV","IV","IV","IV",
    ]
    run_files = ["aircraft_database/" + i for i in run_files]
    save_directory = "eigfigs/"
    num_craft = len(run_files)
    eigensolved = [0.0] * num_craft
    for i in range(num_craft):
        # get json file
        input_dict = json.loads(open(run_files[i]).read())

        # change CL, rho
        input_dict["aerodynamics"]["CL"]["0"] = 0.25
        input_dict["analysis"]["density[slugs/ft^3]"] = 2.37689261479e-03
        eigensolved[i] = es.Solver(input_dict,report=False)

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



    # plot eigenvalues
    klik = ["trad","buck","diml"]
    craft_marker = {"I":".","II":"v","III":"2","IV":"x"}
    # short period eigenvalues
    for k in range(3):
        fig, ax = plt.subplots(tight_layout=True)
        for i in range(num_craft):
            if k == 0:
                dicti = eigensolved[i].p
                evals = "evals"
            if k == 1:
                dicti = eigensolved[i].b
                evals = "evals"
            if k == 2:
                dicti = eigensolved[i].b
                evals = "dimevals"
            i_sp = dicti["lon"]["sp"]
            for j in i_sp:
                eig = dicti["lon"][evals][j]
                ax.plot(np.real(eig),np.imag(eig),\
                    marker=craft_marker[craft_class[i]],c="k")
                # ax.text(np.real(eig),np.imag(eig) + 0.0025,str(i+1),\
                #     ha="center",size=10.0)
        # get limits
        xlim = ax.get_xlim(); ax.set_xlim(xlim)
        ylim = ax.get_ylim(); ax.set_ylim(ylim)
        # plot axes
        ax.plot([xlim[0],xlim[1]],[0.0,0.0],"k")
        ax.plot([0.0,0.0],[ylim[0],ylim[1]],"k")

        plt.savefig(save_directory + "eigs_sameCL_" + klik[k] + ".png")
        if False:
            plt.show()
        else:
            plt.close()
    
    # # run boeing 747
    # file_boeing = run_files[8]

    # # initialize weights, altitudes, velocities
    # n_W = 5; n_rho = 5; n_CL = 5
    # # W = np.linspace(400000.0,636600.0,n_W)
    # W = np.linspace(400000.0,636600.0,n_W)
    # rho = np.linspace(2.37689261479e-03,1.26725882157e-03,n_rho) # SL, 20Kft
    # # V = np.linspace(221.0,850.0,n_V)
    # CL = np.linspace(0.1,0.8,n_CL)

    # for i in range(n_W):
    #     for j in range(n_rho):
    #         for k in range(n_CL):
    #             # get json file
    #             input_dict = json.loads(open(run_files[i]).read())

    #             # change W, CL, rho
    #             Sw = input_dict["aircraft"]["wing_area[ft^2]"]
    #             input_dict["aircraft"]["weight[lbf]"] = W[i]
    #             input_dict["analysis"]["density[slugs/ft^3]"] = rho[j]
    #             V = (W[i]/0.5/rho[j]/CL[k]/Sw)**0.5
    #             print(i,j,k,V)
    #             input_dict["aerodynamics"]["CL"]["0"] = CL[k]
    #             eigensolved[i] = es.Solver(input_dict,report=False)

