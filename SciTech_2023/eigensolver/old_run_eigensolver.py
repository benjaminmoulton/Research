import eigensolver as es
import numpy as np
from matplotlib import pyplot as plt
plt.rcParams['text.usetex'] = True
plt.rcParams["font.family"] = "Times New Roman"
plt.rcParams["font.size"] = 16
# import matplotlib.colors.ColorConverter().to_rgba as rgba

# es.Solver("aircraft_database/F16_bolander.json",report=True)
es.Solver("aircraft_database/NT_33A.json",report=True)
quit()

print("Initializing...")
info = es.Solver("F16.json",report=False).input_dict

# initialize arrays
W_empty = 18900.0
W_max = 42300.0
n_W = 5
W = np.linspace(W_empty,W_max,n_W)
CL_lo = 0.1
CL_max = 0.8
n_CL = 5
CL = np.linspace(CL_lo,CL_max,n_CL)
rho_lo  = 0.00028652 # 55,000 ft altitude
rho_max = 0.0023769 # sea level altitude
n_rho = 5
rho = np.geomspace(rho_lo,rho_max,n_rho)

# initialize eigenvalue 3D arrays
eg = {}
modes = ["sp","ph","ro","sl","dr"]
side = ["lon","lon","lat","lat","lat"]
for mo in modes:
    eg[mo] = np.zeros((3,n_W,n_CL,n_rho),dtype=complex)
spCAP = np.zeros((n_W,n_CL,n_rho))
spzeta = np.zeros((n_W,n_CL,n_rho))
maxv = 0.0

# # initialize koefficients dictionary
# knam = ["Kz,a","Kx,a","Km,a","Kz,ahat","Kx,ahat","Km,ahat","Kz,muhat",
# "Kx,muhat","Km,muhat","Ky,b","Kl,b","Kn,b","Ky,pbreve","Kl,pbreve","Kn,pbreve",
# "Kz,qbreve","Kx,qbreve","Km,qbreve","Ky,rbreve","Kl,rbreve","Kn,rbreve"]
# koef = {}
# for name in knam:
#     koef[name] = np.zeros((n_W,n_CL,n_rho))

print("Running solver...")
for i in range(n_W):
    # set weight in dictionary
    info["aircraft"]["weight[lbf]"] = W[i]

    for j in range(n_CL):
        # set CL in dictionary
        info["aerodynamics"]["CL"]["0"] = CL[j]

        for k in range(n_rho):
            # set rho in dictionary
            info["analysis"]["density[slugs/ft^3]"] = rho[k]

            # run analysis
            run = es.Solver(info,report=False)
            bd = run.g / run.vo
            if maxv < run.vo:
                maxv = run.vo
            
            # pull out short period handling qualities info
            spCAP[i,j,k] = run.b["lon"]["wn"][run.b["lon"]["sp"][0]]**2. / \
                run.CL_a * run.CW
            spzeta[i,j,k] = run.b["lon"]["zeta"][run.b["lon"]["sp"][0]]
            
            # if si == "lon":
            #     td = 2. * run.vo / run.cwbar
            # else:
            #     td = 2. * run.vo / run.bw

            # # pull out koefs
            # for name in knam:
            #     koef[name][i,j,k] = run.b["K"][name]

            for mo,si in zip(modes,side):
                # store eigenvalues
                if type(run.h[si][mo]) == np.ndarray:
                    eg[mo][0,i,j,k] = run.h[si]["evals"][run.h[si][mo][0]]
                    eg[mo][1,i,j,k] = run.b[si]["evals"][run.b[si][mo][0]]
                    eg[mo][2,i,j,k] = run.b[si]["evals"][run.b[si][mo][0]]*bd
                    # eg[mo][3,i,j,k] = run.h[si]["evals"][run.h[si][mo][0]]*td
                else:
                    eg[mo][0,i,j,k] = run.h[si]["evals"][run.h[si][mo]]
                    eg[mo][1,i,j,k] = run.b[si]["evals"][run.b[si][mo]]
                    eg[mo][2,i,j,k] = run.b[si]["evals"][run.b[si][mo]]*bd
                    # eg[mo][3,i,j,k] = run.h[si]["evals"][run.h[si][mo]]*td
# print("max v= ",maxv,"ft/s")

# plot eigenvalues

# report
print("Plotting...")


mkstl = dict(marker="o",fillstyle="top",ms=6,mec=(0,0,0,1),mfc=(0,0,0,1),
mfcalt=(0,0,0,1),label="",ls="none")
ms_small = 4
ms_big = 8
i_step = (ms_big-ms_small)/(n_W-1)
a_small = 0.0
a_big = 1.0
j_step = (a_big-a_small)/(n_CL-1)
k_step = (a_big-a_small)/(n_rho-1)


# #########################################
# for name in knam:
#     for i in range(n_W):
#             for j in range(n_CL):
#                 for k in range(n_rho):
#                     mkstl["ms"] = i * i_step + ms_small
#                     mkstl["mfc"] = (1,0,0, j * j_step + a_small)
#                     mkstl["mfcalt"] = (0,0,1, k * k_step + a_small)
#                     plt.plot(koef[name][i,j,k],koef[name][i,j,k],**mkstl)
#     plt.title(name)
#     plt.show()
# quit()
# #########################################

# create dictionaries for legend
# min weight
loW = dict(mkstl)
loW["label"] = "W = {:<,.0f} lbf".format(W_empty)
loW["mfc"] = (0,0,0,0)
loW["ms"] = ms_small
loW["fillstyle"] = "none"
# min weight
mdW = dict(loW)
mdW["label"] = "W = {:<,.0f} lbf".format((W_empty+W_max)/2)
mdW["ms"] = (ms_small+ms_big)/2
# max weight
hiW = dict(loW)
hiW["label"] = "W = {:<,.0f} lbf".format(W_max)
hiW["ms"] = ms_big
# min CL
loCL = dict(mkstl)
loCL["label"] = "$C_L$ = {:<.1f}".format(CL_lo)
loCL["mfc"] = (1,0,0,0)
loCL["mfcalt"] = (0,0,0,0)
loCL["mec"] = (0,0,0,0)
loCL["ms"] = ms_big
loCL["fillstyle"] = "top"
# mid CL
mdCL = dict(loCL)
mdCL["label"] = "$C_L$ = {:<.1f}".format((CL_lo+CL_max)/2)
mdCL["mfc"] = (1,0,0,0.5)
# max CL
hiCL = dict(loCL)
hiCL["label"] = "$C_L$ = {:<.1f}".format(CL_max)
hiCL["mfc"] = (1,0,0,1)
# min rho
lorho = dict(loCL)
lorho["label"] = r"$\rho$"+" = {:<.3e} slugs/ft$^3$".format(rho_lo)
lorho["mfc"] = (0,0,0,0)
lorho["mfcalt"] = (0,0,1,0)
# mid rho
mdrho = dict(lorho)
mdrho["label"] = r"$\rho$"+" = {:<.3e} slugs/ft$^3$".format((rho_lo+rho_max)/2)
mdrho["mfcalt"] = (0,0,1,0.5)
# max rho
hirho = dict(lorho)
hirho["label"] = r"$\rho$"+" = {:<.3e} slugs/ft$^3$".format(rho_max)
hirho["mfcalt"] = (0,0,1,1)

# increasing size = increasing weight
# increasing red opacity = increasing CL
# increasing blue opacity = increasing density

show_plots = False
show_dim = False
lname = ["T","A","D","t"]
print("\tEigenvalues...")
# plot eigenvalues
for mo,si in zip(modes,side):
    print("\t\t" + mo + "...")
    for l in range(3):
        for i in range(n_W):
            for j in range(n_CL):
                for k in range(n_rho):
                    mkstl["ms"] = i * i_step + ms_small
                    mkstl["mfc"] = (1,0,0, j * j_step + a_small)
                    mkstl["mfcalt"] = (0,0,1, k * k_step + a_small)
                    plt.plot(eg[mo][l,i,j,k].real,eg[mo][l,i,j,k].imag,**mkstl)

        # titles and such
        # plt.legend()
        units = ""
        if lname[l] in ["D","t"]:
            units = ", 1/s"
        plt.xlabel("real" + units)
        plt.ylabel("imaginary" + units)
        # if mo not in ["ro","sl"]:
        #     plt.xscale("symlog")
        #     plt.yscale("symlog")
        # plt.tight_layout()
        fig = plt.gcf()
        fig.savefig("eigfigs/" + mo + lname[l] + ".png",bbox_inches="tight")
        if show_plots:
            plt.show()
        elif show_dim and lname[l] == "D":
            plt.show()
        else:
            plt.close()

        if lname[l] == "D":
            # legend info
            # min weight
            plt.plot(-1,-1,**loW)
            # mid weight
            plt.plot(-1,-1,**mdW)
            # max weight
            plt.plot(-1,-1,**hiW)
            # min CL
            plt.plot(-1,-1,**loCL)
            # mid CL
            plt.plot(-1,-1,**mdCL)
            # max CL
            plt.plot(-1,-1,**hiCL)
            # min rho
            plt.plot(-1,-1,**lorho)
            # mid rho
            plt.plot(-1,-1,**mdrho)
            # max rho
            plt.plot(-1,-1,**hirho)

            # axis limits
            plt.xlim(0,1)
            plt.ylim(0,1)
            plt.legend(loc="center")
            plt.axis("off")
            plt.tight_layout()
            plt.savefig("eigfigs/" + mo + "L.png",bbox_inches="tight")
            plt.close()


# handling qualities
print("\tHandling qualities...")

print("\t\tShort Period...")

# plot level lines
def lgmid(a,b):
    return np.exp((np.log(a)+np.log(b))/2.)
txtd = dict(ha="center",va="center")
ld = dict(
A1 = np.array([
    [1.30,0.28],
    [0.35,0.28],
    [0.35,3.60],
    [1.30,3.60],
    [1.30,0.28]
]),
A2 = np.array([
    [2.00,0.15],
    [0.25,0.15],
    [0.25,10.0],
    [2.00,10.0],
    [2.00,0.15]
]),
A3 = np.array([
    [5.00,0.15],
    [0.13,0.15],
    [0.13,30.0]
]),
B1 = np.array([
    [2.00,0.085],
    [0.30,0.085],
    [0.30,3.600],
    [2.00,3.600],
    [2.00,0.085]
]),
B2 = np.array([
    [2.00,0.038],
    [0.20,0.038],
    [0.20,10.00],
    [2.00,10.00],
    [2.00,0.038]
]),
B3 = np.array([
    [5.00,0.038],
    [0.13,0.038],
    [0.13,30.0]
]),
C1 = np.array([
    [1.30,0.15],
    [0.35,0.15],
    [0.35,3.60],
    [1.30,3.60],
    [1.30,0.15]
]),
C2 = np.array([
    [2.00,0.096],
    [0.25,0.096],
    [0.25,10.00],
    [2.00,10.00],
    [2.00,0.096]
]),
C3 = np.array([
    [5.00,0.096],
    [0.13,0.096],
    [0.13,30.00]
])
)

phase = ["A","B","C"]

for p in phase:
    print("\t\t\t" + p + " flight phases...")
    plt.plot(ld[p+"1"][:,0],ld[p+"1"][:,1],"k",lw=2)
    plt.plot(ld[p+"2"][:,0],ld[p+"2"][:,1],"k",lw=2)
    plt.plot(ld[p+"3"][:,0],ld[p+"3"][:,1],"k",lw=2)

    for i in range(n_W):
        for j in range(n_CL):
            for k in range(n_rho):
                mkstl["ms"] = i * i_step + ms_small
                mkstl["mfc"] = (1,0,0, j * j_step + a_small)
                mkstl["mfcalt"] = (0,0,1, k * k_step + a_small)
                plt.plot(spzeta[i,j,k],spCAP[i,j,k],**mkstl)

    plt.xscale("log")
    plt.yscale("log")
    plt.xlabel("Short-Period Damping Ratio," + r"$\zeta_{sp}$")
    plt.ylabel("CAP, 1/s$^2$")
    plt.xlim(0.05,ld[p+"3"][0,0])
    plt.ylim(0.01,ld[p+"3"][2,1])
    plt.text(lgmid(ld[p+"1"][0,0],ld[p+"1"][1,0]),\
        lgmid(ld[p+"1"][0,1],ld[p+"1"][3,1]),"Level 1",**txtd)
    plt.text(lgmid(ld[p+"2"][0,0],ld[p+"2"][1,0]),\
        lgmid(ld[p+"2"][3,1],ld[p+"1"][3,1]),"Level 2",**txtd)
    plt.text(lgmid(ld[p+"2"][0,0],ld[p+"2"][1,0]),\
        lgmid(ld[p+"2"][2,1],ld[p+"3"][2,1]),"Level 3",**txtd)
    plt.text(lgmid(ld[p+"2"][0,0],ld[p+"2"][1,0]),\
        lgmid(0.01,ld[p+"3"][0,1]),"Level 4",**txtd)
    fig = plt.gcf()
    fig.savefig("hqfigs/spHQ" + p + ".png",bbox_inches="tight")
    if show_plots:
        plt.show()
    else:
        plt.close()