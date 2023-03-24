from shutil import ReadError
import numpy as np
import json
import matplotlib.pyplot as plt
from scipy.linalg import eig
import numpy as np
from matplotlib import pyplot as plt

class Solver:
    """A class which solves an eigenproblem based on an aircraft 
    derivatives input file.

    Parameters
    ----------
    input_vars : string or dict , optional
        Must be a .json file or python dictionary.

    Raises
    ------
    TypeError
        If the input_vars type is not a dictionary or the file path to 
        a .json file
    """
    def __init__(self,input_vars={},report=True,show_evector=False):

        # report
        self.report = report
        if self.report:
            print("Running Aircraft Dynamics EigenSolver written by Ben Moulton\n")

        # get info or raise error
        self.input = input_vars
        self._get_input_vars(input_vars)

        # retrieve info
        self._retrieve_info()

        # create A and B matrices, 
        self._phillips_matrices_and_approximations()
        self._hunsaker_matrices_and_approximations()
        self._buckingham_matrices_and_approximations()

        # determine longitudinal dimensional properties
        self._longitudinal_dimensional_properties(self.p)
        self._longitudinal_dimensional_properties(self.h)
        self._longitudinal_dimensional_properties(self.b,is_alternate=True)

        # determine lateral dimensional properties
        self._lateral_dimensional_properties(self.p)
        self._lateral_dimensional_properties(self.h)
        self._lateral_dimensional_properties(self.b,is_alternate=True)

        # report properties
        if report:
            self._full_report(show_evector)


    def _get_input_vars(self,input_vars):
        # get info or raise error

        # determine if the input_vars is a file or a dictionary
        input_vars_type = type(input_vars)

        # dictionary
        if input_vars_type == dict:
            self.input_dict = input_vars
        
        # json file
        elif input_vars_type == str and input_vars.split(".")[-1] == "json":
            # import json file from file path
            json_string = open(input_vars).read()

            # save to vals dictionary
            self.input_dict = json.loads(json_string)

        # raise error
        else:
            raise TypeError("input_vars must be json file path, or " + \
                "dictionary, not {0}".format(input_vars_type))


    def _retrieve_info(self):
        """Method which retrieves the information and stores it globally.
        """
        
        # store variables from file input dictionary
        
        # store aircraft input values
        aircraft = self.input_dict.get("aircraft",{})
        self.aircraft_name = aircraft.get("name","Jimmy")
        self.EL = aircraft.get("launch_kinetic_energy[ft-lbf]",0.0)
        self.Sw = aircraft.get("wing_area[ft^2]")
        self.bw = aircraft.get("wing_span[ft]")
        self.cwbar = self.Sw / self.bw
        self.W = aircraft.get("weight[lbf]")
        self.Ixx = aircraft.get("Ixx[slug-ft^2]")
        self.Iyy = aircraft.get("Iyy[slug-ft^2]")
        self.Izz = aircraft.get("Izz[slug-ft^2]")
        self.Ixy = aircraft.get("Ixy[slug-ft^2]")
        self.Ixz = aircraft.get("Ixz[slug-ft^2]")
        self.Iyz = aircraft.get("Iyz[slug-ft^2]")
        
        # store analysis input values
        analysis = self.input_dict.get("analysis",{})
        self.rho = analysis.get("density[slugs/ft^3]")
        
        # store aerodynamics input values
        aero = self.input_dict.get("aerodynamics",{})
        CL = aero.get("CL")
        self.CL0        = CL.get("0")
        self.CL_a       = CL.get("alpha")
        self.CL_qbar    = CL.get("qbar")
        self.CL_ahat    = CL.get("alpha_hat")
        self.CL_uhat    = CL.get("mu_hat",0.0)
        CY = aero.get("CS")
        self.CY_b       = CY.get("beta")
        self.CY_pbar    = CY.get("Lpbar") * self.CL0 + CY.get("pbar")
        self.CY_rbar    = CY.get("rbar")
        self.CY_da      = CY.get("da",0.0)
        self.CY_dr      = CY.get("dr",0.0)
        CD = aero.get("CD")
        self.CD0        = CD.get("L0")
        self.CD1        = CD.get("L")
        self.CD2        = CD.get("L2")
        self.CD_qbar    = CD.get("L2qbar") * self.CL0**2. + \
            CD.get("Lqbar") * self.CL0 + CD.get("qbar")
        self.CD_ahat    = CD.get("alpha_hat",0.0)
        self.CD_uhat    = CD.get("mu_hat",0.0)
        Cl = aero.get("Cl")
        self.Cl_b       = Cl.get("beta")
        self.Cl_pbar    = Cl.get("pbar")
        self.Cl_rbar    = Cl.get("Lrbar") * self.CL0 + Cl.get("rbar")
        self.Cl_da      = Cl.get("da",0.0)
        self.Cl_dr      = Cl.get("dr",0.0)
        Cm = aero.get("Cm")
        self.Cmo        = Cm.get("0")
        self.Cm_a       = Cm.get("alpha")
        self.Cm_qbar    = Cm.get("qbar")
        self.Cm_ahat    = Cm.get("alpha_hat")
        self.Cm_uhat    = Cm.get("mu_hat",0.0)
        Cn = aero.get("Cn")
        self.Cn_b       = Cn.get("beta")
        self.Cn_pbar    = Cn.get("Lpbar") * self.CL0 + Cn.get("pbar")
        self.Cn_rbar    = Cn.get("rbar")
        self.Cn_da      = Cn.get("da",0.0)
        self.Cn_dr      = Cn.get("dr",0.0)

        # store handling quality input values
        self.HQ = self.input_dict.get("handling_qualities",{})
        if "zt_sp" in self.HQ and "wn_sp" in self.HQ:
            self.split_short_period = False
        elif "S1_sp" in self.HQ and "S2_sp" in self.HQ:
            self.split_short_period = True
        else:
            raise ReadError("missing short-period handling-qualities.")
        if "zt_ph" in self.HQ and "wn_ph" in self.HQ:
            self.split_phugoid = False
        elif "S1_ph" in self.HQ and "S2_ph" in self.HQ:
            self.split_phugoid = True
        else:
            raise ReadError("missing phugoid handling-qualities.")

        # set gravity
        self.g = 32.17#4

        # set angles
        self.thetao = 0.0
        self.phio = 0.0
        self.ct = np.cos(self.thetao)
        self.st = np.sin(self.thetao)
        self.cp = np.cos(self.phio)

        # calculate velocity
        self.vo = (self.W / 0.5 / self.rho / self.CL0 / self.Sw)**0.5

        # calculate lift and drag
        self.CW = self.CL0
        self.CLo = self.CL0 * self.ct / self.cp
        self.CDo = self.CD0 + self.CD1 * self.CLo + self.CD2 * self.CLo**2.
        self.CD_a = self.CD1 * self.CL_a + 2. * self.CD2 * self.CLo * self.CL_a
        if self.report:
            print("Running solver for {}".format(self.aircraft_name))
            print("Vo   = ", self.vo)
            print("CLo  = ", self.CLo)
            print("CDo  = ", self.CDo)
            print("CD,a = ", self.CD_a)
            print()

        # initialize information dictionaries
        self.h = {}
        self.h["name"] = "Hunsaker"
        self.p = {}
        self.p["name"] = "Phillips"
        self.b = {}
        self.b["name"] = "Buckingham"


    def _phillips_matrices_and_approximations(self):
        """Method which creates the matrices to be solved.
        """

        # initialized constants
        Rpx = self.rho * self.Sw * self.cwbar / 4.0 / (self.W / self.g)
        Rgx = self.g * self.cwbar / 2. / self.vo**2.
        Ryy = self.rho * self.Sw * self.cwbar**3. / 8. / self.Iyy
        Rpy = self.rho * self.Sw * self.bw / 4.0 / (self.W / self.g)
        Rgy = self.g * self.bw / 2. / self.vo**2.
        Rxx = self.rho * self.Sw * self.bw**3. / 8. / self.Ixx
        # Rxz = 8. * self.Ixz / self.rho / self.Sw / self.bw**3.
        Rzz = self.rho * self.Sw * self.bw**3. / 8. / self.Izz
        # LONGITUDINAL
        # values
        Rz_ahat = Rpx * - self.CL_ahat
        Rm_ahat = Ryy * self.Cm_ahat
        Rx_mu = Rpx * -2. * self.CDo
        Rz_mu = Rpx * -2. * self.CLo
        Rm_mu = Ryy * 2. * self.Cmo
        Rx_a = Rpx * (self.CLo - self.CD_a)
        Rz_a = Rpx * (-self.CL_a - self.CDo)
        Rm_a = Ryy * self.Cm_a
        Rx_qbar = Rpx * - self.CD_qbar
        Rz_qbar = Rpx * - self.CL_qbar
        Rm_qbar = Ryy * self.Cm_qbar
        # initialize A and B matrices
        A = np.zeros((6,6))
        A[0,0] = Rx_mu
        A[0,1] = Rx_a
        A[0,2] = Rx_qbar
        A[0,5] = - Rgx * self.ct
        A[1,0] = Rz_mu
        A[1,1] = Rz_a
        A[1,2] = 1. + Rz_qbar
        A[1,5] = - Rgx * self.st
        A[2,0] = Rm_mu
        A[2,1] = Rm_a
        A[2,2] = Rm_qbar
        A[3,0] = A[4,1] = self.ct
        A[3,1] = self.st
        A[4,5] = - self.ct
        A[3,5] = A[4,0] = - self.st
        A[5,2] = 1.

        B = np.identity(6)
        B[1,1] = 1. - Rz_ahat
        B[2,1] = - Rm_ahat

        # calculate C matrix
        # self.p["A"] = A * 1.
        # self.p["B"] = B * 1.
        C = np.matmul(np.linalg.inv(B),A)

        # print("A matrix")
        # for i in range(6):
        #     for j in range(6):
        #         print("\t{:>17.12f}".format(A[i,j]),end="")
        #     print()

        # print("B matrix")
        # for i in range(6):
        #     for j in range(6):
        #         print("\t{:>17.12f}".format(B[i,j]),end="")
        #     print()

        # calculate eigenvalues and vectors
        self.p["lon"] = {}
        self.p["lon"]["A"] = A
        self.p["lon"]["B"] = B
        self.p["lon"]["evals"],self.p["lon"]["evecs"] = eig(C)

        # calculate amplitude and phase angles
        self.p["lon"]["amp"] = np.sqrt(np.real(self.p["lon"]["evecs"])**2. + \
            np.imag(self.p["lon"]["evecs"])**2.)
        self.p["lon"]["phase"] = np.rad2deg(\
            np.arctan2(np.imag(self.p["lon"]["evecs"]),\
            np.real(self.p["lon"]["evecs"])))

        # LATERAL
        # values
        ixz = self.Ixz / self.Ixx
        izx = self.Ixz / self.Izz
        Ry_b = Rpy * self.CY_b
        Rl_b = Rxx * self.Cl_b
        Rn_b = Rzz * self.Cn_b
        Ry_pbar = Rpy * self.CY_pbar
        Rl_pbar = Rxx * self.Cl_pbar
        Rn_pbar = Rzz * self.Cn_pbar
        Ry_rbar = Rpy * self.CY_rbar
        Rl_rbar = Rxx * self.Cl_rbar
        Rn_rbar = Rzz * self.Cn_rbar
        # initialize A and B matrices
        A *= 0.
        A[0,0] = Ry_b
        A[0,1] = Ry_pbar
        A[0,2] = Ry_rbar - 1.
        A[0,4] = Rgy * self.ct
        A[1,0] = Rl_b
        A[1,1] = Rl_pbar
        A[1,2] = Rl_rbar
        A[2,0] = Rn_b
        A[2,1] = Rn_pbar
        A[2,2] = Rn_rbar
        A[3,0] = A[4,1] = 1.
        A[3,5] = self.ct
        A[4,2] = self.st / self.ct
        A[5,2] = 1. / self.ct

        B = np.identity(6)
        B[1,2] = - ixz
        B[2,1] = - izx

        # calculate C matrix
        C = np.matmul(np.linalg.inv(B),A)

        # calculate eigenvalues and vectors
        self.p["lat"] = {}
        self.p["lat"]["A"] = A
        self.p["lat"]["B"] = B
        self.p["lat"]["evals"],self.p["lat"]["evecs"] = eig(C)

        # calculate amplitude and phase angles
        self.p["lat"]["amp"] = np.sqrt(np.real(self.p["lat"]["evecs"])**2. + \
            np.imag(self.p["lat"]["evecs"])**2.)
        self.p["lat"]["phase"] = np.rad2deg(\
            np.arctan2(np.imag(self.p["lat"]["evecs"]),np.real(self.p["lat"]["evecs"])))

        
        # calculate Short Period approximation
        self.p["spSg"] = - self.vo / self.cwbar * (Rz_a + Rm_qbar + Rm_ahat)
        self.p["spt99"] = np.log(0.01) / - self.p["spSg"]
        self.p["spt2x"] = np.log(2.0) / - self.p["spSg"]
        self.p["spWD"] = 2. *self.vo / self.cwbar * np.sqrt(np.abs(\
            (Rz_a*Rm_qbar + Rm_a) - 0.25*(Rz_a + Rm_qbar + Rm_ahat)**2.))
        self.p["spT"] = 2. * np.pi / self.p["spWD"]

        # calculate Phugoid approximation
        Rs = Rm_a / (Rm_a - Rz_a * Rm_qbar)
        Rd = Rx_a * Rm_qbar / (Rm_a - Rz_a * Rm_qbar)
        Rp = Rgx * Rs * (Rz_a + Rm_qbar) / (Rm_a - Rz_a * Rm_qbar)
        self.p["phSg"] = -2. * self.vo / self.cwbar * (Rz_mu/2. * \
            (Rx_mu/Rz_mu + Rd - Rp))
        self.p["pht99"] = np.log(0.01) / - self.p["phSg"]
        self.p["pht2x"] = np.log(2.0) / - self.p["phSg"]
        self.p["phWD"] = 2. * self.vo / self.cwbar * Rz_mu/2. * \
            np.sqrt( -4.*Rgx/Rz_mu*Rs - (Rx_mu/Rz_mu+Rd)**2.)
        self.p["phT"] = 2. * np.pi / self.p["phWD"]

        # calculate Roll approximation
        self.p["roSg"] = - self.rho * self.Sw * self.bw**2. * self.vo / 4. / \
            self.Ixx * self.Cl_pbar
        self.p["rot99"] = np.log(0.01) / - self.p["roSg"]

        # calculate Spiral approximation
        self.p["slSg"] = self.g / self.vo * \
            (self.Cl_b*self.Cn_rbar - self.Cl_rbar*self.Cn_b) / \
            (self.Cl_b*self.Cn_pbar - self.Cl_pbar*self.Cn_b)
        self.p["slt99"] = np.log(0.01) / - self.p["slSg"]
        self.p["slt2x"] = np.log(2.0) / - self.p["slSg"]

        # calculate Dutch Roll approximation
        RDRs = (Rl_b*(Rgy - (1.-Ry_rbar)*Rn_pbar) - Ry_b*Rl_rbar*Rn_pbar) / \
            Rl_pbar
        RDRc = Rl_rbar * Rn_pbar / Rl_pbar
        RDRp = Rgy*(Rl_rbar*Rn_b-Rl_b*Rn_rbar)/Rl_pbar/(Rn_b+Ry_b*Rn_rbar) - \
            RDRs/Rl_pbar
        self.p["drSg"] = - self.vo / self.bw * (Ry_b + Rn_rbar - RDRc + RDRp)
        self.p["drt99"] = np.log(0.01) / - self.p["drSg"]
        self.p["drt2x"] = np.log(2.0) / - self.p["drSg"]

        a = (1.-Ry_rbar)*Rn_b 
        b = Ry_b * Rn_rbar 
        c = RDRs 
        d = - 0.25*( Ry_b + Rn_rbar)**2.
        e = complex(a + b + c + d)
        self.p["drWD"] = np.abs(2. * self.vo / self.bw * np.sqrt( e ))
        self.p["drT"] = 2. * np.pi / self.p["drWD"]

        # print("A matrix")
        # for i in range(6):
        #     for j in range(6):
        #         print("\t{:>17.12f}".format(A[i,j]),end="")
        #     print()

        # print("B matrix")
        # for i in range(6):
        #     for j in range(6):
        #         print("\t{:>17.12f}".format(B[i,j]),end="")
        #     print()

        # print("C matrix")
        # for i in range(6):
        #     for j in range(6):
        #         print("\t{:>17.12f}".format(C[i,j]),end="")
        #     print()

        # print("\nLONGITUDINAL")
        # print("Eigenvalues")
        # for i in range(6):
        #     print("\t{:>17.12f}".format(self.h["lon"]["evals"][i]))

        # print("\nEigenvectors")
        # for i in range(6):
        #     for j in range(6):
        #         print("\t{:>17.12f}".format(self.h["lon"]["evecs"][i,j]),end="")
        #     print()

        # print("\nAmplitudes")
        # for i in range(6):
        #     for j in range(6):
        #         print("\t{:>17.12f}".format(self.h["lon"]["amp"][i,j]),end="")
        #     print()

        # print("\nPhase Angle (deg)")
        # for i in range(6):
        #     for j in range(6):
        #         print("\t{:>17.12f}".format(self.h["lon"]["phase"][i,j]),end="")
        #     print()

        # print("\nLATERAL")
        # print("Eigenvalues")
        # for i in range(6):
        #     print("\t{:>17.12f}".format(self.h["lat"]["evals"][i]))

        # print("Eigenvectors")
        # for i in range(6):
        #     for j in range(6):
        #         print("\t{:>17.12f}".format(self.h["lat"]["evecs"][i,j]),end="")
        #     print()

        # print("\nAmplitudes")
        # for i in range(6):
        #     for j in range(6):
        #         print("\t{:>17.12f}".format(self.h["lat"]["amp"][i,j]),end="")
        #     print()

        # print("\nPhase Angle (deg)")
        # for i in range(6):
        #     for j in range(6):
        #         print("\t{:>17.12f}".format(self.h["lat"]["phase"][i,j]),end="")
        #     print()
        # print()

        a = 10


    def _hunsaker_matrices_and_approximations(self):
        """Method which creates the matrices to be solved.
        """

        # initialized constants
        Rpx = 4.0 * self.W / self.g / self.rho / self.Sw / self.cwbar
        Rgx = self.g * self.cwbar / 2. / self.vo**2.
        Ryy = 8. * self.Iyy / self.rho / self.Sw / self.cwbar**3.
        Rpy = 4.0 * self.W / self.g / self.rho / self.Sw / self.bw
        Rgy = self.g * self.bw / 2. / self.vo**2.
        Rxx = 8. * self.Ixx / self.rho / self.Sw / self.bw**3.
        Rxz = 8. * self.Ixz / self.rho / self.Sw / self.bw**3.
        Rzz = 8. * self.Izz / self.rho / self.Sw / self.bw**3.

        # LONGITUDINAL
        # initialize A and B matrices
        A = np.zeros((6,6))
        A[0,0] = - 2. * self.CDo
        A[0,1] = self.CLo - self.CD_a
        A[0,2] = - self.CD_qbar
        A[0,5] = - Rpx * Rgx * self.ct
        A[1,0] = - 2. * self.CLo
        A[1,1] = - self.CL_a - self.CDo
        A[1,2] = - self.CL_qbar + Rpx
        A[1,5] = - Rpx * Rgx * self.st
        A[2,0] = 2. * self.Cmo
        A[2,1] = self.Cm_a
        A[2,2] = self.Cm_qbar
        A[3,0] = A[4,1] = self.ct
        A[3,1] = self.st
        A[4,5] = - self.ct
        A[3,5] = A[4,0] = - self.st
        A[5,2] = 1.

        B = np.zeros((6,6))
        B[0,0] = Rpx + self.CD_uhat
        B[0,1] = self.CD_ahat
        B[1,0] = self.CL_uhat
        B[1,1] = Rpx + self.CL_ahat
        B[2,0] = - self.Cm_uhat
        B[2,1] = - self.Cm_ahat
        B[2,2] = Ryy
        B[3,3] = B[4,4] = B[5,5] = 1.

        # calculate C matrix
        C = np.matmul(np.linalg.inv(B),A)

        # print("A matrix")
        # for i in range(6):
        #     for j in range(6):
        #         print("\t{:>17.12f}".format(A[i,j]),end="")
        #     print()

        # print("B matrix")
        # for i in range(6):
        #     for j in range(6):
        #         print("\t{:>17.12f}".format(B[i,j]),end="")
        #     print()

        # calculate eigenvalues and vectors
        self.h["lon"] = {}
        self.h["lon"]["A"] = A
        self.h["lon"]["B"] = B
        self.h["lon"]["evals"],self.h["lon"]["evecs"] = eig(C)

        # calculate amplitude and phase angles
        self.h["lon"]["amp"] = np.sqrt(np.real(self.h["lon"]["evecs"])**2. + \
            np.imag(self.h["lon"]["evecs"])**2.)
        self.h["lon"]["phase"] = np.rad2deg(\
            np.arctan2(np.imag(self.h["lon"]["evecs"]),\
            np.real(self.h["lon"]["evecs"])))

        # LATERAL
        # initialize A and B matrices
        A *= 0.
        A[0,0] = self.CY_b
        A[0,1] = self.CY_pbar
        A[0,2] = self.CY_rbar - Rpy
        A[0,4] = Rpy * Rgy * self.ct
        A[1,0] = self.Cl_b
        A[1,1] = self.Cl_pbar
        A[1,2] = self.Cl_rbar
        A[2,0] = self.Cn_b
        A[2,1] = self.Cn_pbar
        A[2,2] = self.Cn_rbar
        A[3,0] = A[4,1] = 1.
        A[3,5] = self.ct
        A[4,2] = self.st / self.ct
        A[5,2] = 1. / self.ct

        B *= 0.
        B[0,0] = Rpy
        B[1,1] = Rxx
        B[1,2] = B[2,1] = - Rxz
        B[2,2] = Rzz
        B[3,3] = B[4,4] = B[5,5] = 1.

        # calculate C matrix
        C = np.matmul(np.linalg.inv(B),A)

        # calculate eigenvalues and vectors
        self.h["lat"] = {}
        self.h["lat"]["A"] = A
        self.h["lat"]["B"] = B
        self.h["lat"]["evals"],self.h["lat"]["evecs"] = eig(C)

        # calculate amplitude and phase angles
        self.h["lat"]["amp"] = np.sqrt(np.real(self.h["lat"]["evecs"])**2. + \
            np.imag(self.h["lat"]["evecs"])**2.)
        self.h["lat"]["phase"] = np.rad2deg(\
            np.arctan2(np.imag(self.h["lat"]["evecs"]),np.real(self.h["lat"]["evecs"])))

        
        # calculate Short Period approximation
        Asp = Ryy * (Rpx + self.CL_ahat)
        Bsp = Ryy * (self.CL_a + self.CDo) - self.Cm_qbar * (Rpx + self.CL_ahat) - \
            self.Cm_ahat * (Rpx - self.CL_qbar)
        Csp = - self.Cm_qbar * (self.CL_a + self.CDo) - \
            self.Cm_a * (Rpx - self.CL_qbar)
        self.h["spSg"] = self.vo / self.cwbar * Bsp / Asp
        self.h["spt99"] = np.log(0.01) / - self.h["spSg"]
        self.h["spt2x"] = np.log(2.0) / - self.h["spSg"]
        self.h["spWD"] = self.vo / self.cwbar * np.sqrt(np.abs(Bsp**2. - \
            4.*Asp*Csp)/Asp)
        self.h["spT"] = 2. * np.pi / self.h["spWD"]

        # calculate Phugoid approximation
        SD = self.g / self.vo * self.CDo / self.CLo
        Sq = self.g / self.vo * (self.CLo-self.CD_a)*self.Cm_qbar/(Rpx*self.Cm_a + \
            (self.CDo+self.CL_a)*self.Cm_qbar)
        Rps = Rpx*self.Cm_a / (Rpx*self.Cm_a + (self.CDo+self.CL_a)*self.Cm_qbar)
        Sp = - self.g / self.vo * Rgx * Rps * \
            (Rpx*self.Cm_qbar - Ryy*(self.CDo+self.CL_a)) / \
            (Rpx*self.Cm_a + (self.CDo+self.CL_a)*self.Cm_qbar)
        self.h["phSg"] = SD + Sq + Sp
        self.h["pht99"] = np.log(0.01) / - self.h["phSg"]
        self.h["pht2x"] = np.log(2.0) / - self.h["phSg"]
        self.h["phWD"] = np.sqrt(2.*(self.g/self.vo)**2. * Rps - (SD+Sq)**2.)
        self.h["phT"] = 2. * np.pi / self.h["phWD"]

        # calculate Roll approximation
        self.h["roSg"] = - self.rho * self.Sw * self.bw**2. * self.vo / 4. / \
            self.Ixx * self.Cl_pbar
        self.h["rot99"] = np.log(0.01) / - self.h["roSg"]

        # calculate Spiral approximation
        self.h["slSg"] = self.g / self.vo * \
            (self.Cl_b*self.Cn_rbar - self.Cl_rbar*self.Cn_b) / \
            (self.Cl_b*self.Cn_pbar - self.Cl_pbar*self.Cn_b)
        self.h["slt99"] = np.log(0.01) / - self.h["slSg"]
        self.h["slt2x"] = np.log(2.0) / - self.h["slSg"]

        # calculate Dutch Roll approximation
        RDRs = (self.Cl_b*(Rgy*Rpy*Rzz - (Rpy-self.CY_rbar)*self.Cn_pbar) - \
            self.CY_b*self.Cl_rbar*self.Cn_pbar) / Rpy / Rzz / self.Cl_pbar
        self.h["drSg"] = - self.vo / self.bw * (self.CY_b/Rpy + self.Cn_rbar/Rzz - \
            self.Cl_rbar*self.Cn_pbar/self.Cl_pbar/Rzz + Rgy*\
            (self.Cl_rbar*self.Cn_b-self.Cl_b*self.Cn_rbar) / self.Cl_pbar / \
            (self.Cn_b + self.CY_b*self.Cn_rbar/Rpy) - Rxx*RDRs/self.Cl_pbar)
        self.h["drt99"] = np.log(0.01) / - self.h["drSg"]
        self.h["drt2x"] = np.log(2.0) / - self.h["drSg"]
        self.h["drWD"] = 2. * self.vo / self.bw * np.sqrt( (1. - self.CY_rbar/Rpy)*\
            self.Cn_b/Rzz + self.CY_b*self.Cn_rbar/Rpy/Rzz + RDRs - \
            0.25 * ( self.CY_b/Rpy + self.Cn_rbar/Rzz )**2. )
        self.h["drT"] = 2. * np.pi / self.h["drWD"]



        # print("A matrix")
        # for i in range(6):
        #     for j in range(6):
        #         print("\t{:>17.12f}".format(A[i,j]),end="")
        #     print()

        # print("B matrix")
        # for i in range(6):
        #     for j in range(6):
        #         print("\t{:>17.12f}".format(B[i,j]),end="")
        #     print()

        # print("C matrix")
        # for i in range(6):
        #     for j in range(6):
        #         print("\t{:>17.12f}".format(C[i,j]),end="")
        #     print()

        # print("\nLONGITUDINAL")
        # print("Eigenvalues")
        # for i in range(6):
        #     print("\t{:>17.12f}".format(self.h["lon"]["evals"][i]))

        # print("\nEigenvectors")
        # for i in range(6):
        #     for j in range(6):
        #         print("\t{:>17.12f}".format(self.h["lon"]["evecs"][i,j]),end="")
        #     print()

        # print("\nAmplitudes")
        # for i in range(6):
        #     for j in range(6):
        #         print("\t{:>17.12f}".format(self.h["lon"]["amp"][i,j]),end="")
        #     print()

        # print("\nPhase Angle (deg)")
        # for i in range(6):
        #     for j in range(6):
        #         print("\t{:>17.12f}".format(self.h["lon"]["phase"][i,j]),end="")
        #     print()

        # print("\nLATERAL")
        # print("Eigenvalues")
        # for i in range(6):
        #     print("\t{:>17.12f}".format(self.h["lat"]["evals"][i]))

        # print("Eigenvectors")
        # for i in range(6):
        #     for j in range(6):
        #         print("\t{:>17.12f}".format(self.h["lat"]["evecs"][i,j]),end="")
        #     print()

        # print("\nAmplitudes")
        # for i in range(6):
        #     for j in range(6):
        #         print("\t{:>17.12f}".format(self.h["lat"]["amp"][i,j]),end="")
        #     print()

        # print("\nPhase Angle (deg)")
        # for i in range(6):
        #     for j in range(6):
        #         print("\t{:>17.12f}".format(self.h["lat"]["phase"][i,j]),end="")
        #     print()
        # print()

        a = 10


    def _buckingham_matrices_and_approximations(self):
        """Method which creates the matrices to be solved.
        """

        # initialized constants        
        Rpx = self.rho * self.Sw * self.cwbar * self.g / 4. / self.W
        Rpy = self.rho * self.Sw * self.bw    * self.g / 4. / self.W
        CWinv = 1. / self.CW
        Ryy2 = self.rho * self.vo**4.*self.Sw*self.cwbar/2./self.g**2./self.Iyy
        Rxx2 = self.rho * self.vo**4.*self.Sw*self.bw   /2./self.g**2./self.Ixx
        Rzz2 = self.rho * self.vo**4.*self.Sw*self.bw   /2./self.g**2./self.Izz
        Ryy = self.rho*self.vo**2.*self.Sw*self.cwbar**2./4./self.g/self.Iyy
        Rxx = self.rho*self.vo**2.*self.Sw*self.bw**2./4./self.g/self.Ixx
        Rzz = self.rho*self.vo**2.*self.Sw*self.bw**2./4./self.g/self.Izz

        # LONGITUDINAL
        # values
        ixz = self.Ixz / self.Ixx
        izx = self.Ixz / self.Izz
        Kx_muhat = Rpx * self.CD_uhat
        Kz_muhat = Rpx * self.CL_uhat
        Km_muhat = Ryy * self.Cm_uhat
        Kx_ahat = Rpx * self.CD_ahat
        Kz_ahat = Rpx * self.CL_ahat
        Km_ahat = Ryy * self.Cm_ahat
        Kx_mu = CWinv * -2. * self.CDo
        Kz_mu = CWinv *  2. * self.CLo
        Km_mu = Ryy2 * 2. * self.Cmo
        Kx_a = CWinv * ( self.CLo -      self.CD_a)
        Kz_a = CWinv * (-self.CDo - self.CL_a)
        Km_a = Ryy2 * self.Cm_a
        Ky_b = CWinv * self.CY_b
        Kl_b = Rxx2 * self.Cl_b
        Kn_b = Rzz2 * self.Cn_b
        Ky_pbreve = Rpy * self.CY_pbar
        Kl_pbreve = Rxx * self.Cl_pbar
        Kn_pbreve = Rzz * self.Cn_pbar
        Kx_qbreve = Rpx * self.CD_qbar
        Kz_qbreve = Rpx * self.CL_qbar
        Km_qbreve = Ryy * self.Cm_qbar
        Ky_rbreve = Rpy * self.CY_rbar
        Kl_rbreve = Rxx * self.Cl_rbar
        Kn_rbreve = Rzz * self.Cn_rbar
        # initialize A and B matrices
        A = np.zeros((6,6))
        A[0,0] = Kx_mu
        A[0,1] = Kx_a
        A[0,2] = - Kx_qbreve
        A[0,5] = - self.ct
        A[1,0] = - Kz_mu
        A[1,1] = Kz_a
        A[1,2] = 1. - Kz_qbreve
        A[1,5] = - self.st
        A[2,0] = Km_mu
        A[2,1] = Km_a
        A[2,2] = Km_qbreve
        A[3,0] = A[4,1] = self.ct
        A[3,1] = self.st
        A[4,5] = - self.ct
        A[3,5] = A[4,0] = - self.st
        A[5,2] = 1.

        B = np.identity(6)
        B[0,0] = 1. + Kx_muhat
        B[0,1] = Kx_ahat
        B[1,0] = Kz_muhat
        B[1,1] = 1. + Kz_ahat
        B[2,0] = - Km_muhat
        B[2,1] = - Km_ahat

        # calculate C matrix
        C = np.matmul(np.linalg.inv(B),A)

        # print("A matrix")
        # for i in range(6):
        #     for j in range(6):
        #         print("\t{:>17.12f}".format(A[i,j]),end="")
        #     print()

        # print("B matrix")
        # for i in range(6):
        #     for j in range(6):
        #         print("\t{:>17.12f}".format(B[i,j]),end="")
        #     print()

        # calculate eigenvalues and vectors
        self.b["lon"] = {}
        self.b["lon"]["A"] = A
        self.b["lon"]["B"] = B
        self.b["lon"]["evals"],self.b["lon"]["evecs"] = eig(C)

        # calculate amplitude and phase angles
        self.b["lon"]["amp"] = np.sqrt(np.real(self.b["lon"]["evecs"])**2. + \
            np.imag(self.b["lon"]["evecs"])**2.)
        self.b["lon"]["phase"] = np.rad2deg(\
            np.arctan2(np.imag(self.b["lon"]["evecs"]),\
            np.real(self.b["lon"]["evecs"])))

        # LATERAL
        # initialize A and B matrices
        A *= 0.
        A[0,0] = Ky_b
        A[0,1] = Ky_pbreve
        A[0,2] = Ky_rbreve - 1.
        A[0,4] = self.ct
        A[1,0] = Kl_b
        A[1,1] = Kl_pbreve
        A[1,2] = Kl_rbreve
        A[2,0] = Kn_b
        A[2,1] = Kn_pbreve
        A[2,2] = Kn_rbreve
        A[3,0] = A[4,1] = 1.
        A[3,5] = self.ct
        A[4,2] = self.st / self.ct
        A[5,2] = 1. / self.ct

        B = np.identity(6)
        B[1,2] = - ixz
        B[2,1] = - izx

        # calculate C matrix
        C = np.matmul(np.linalg.inv(B),A)

        # calculate eigenvalues and vectors
        self.b["lat"] = {}
        self.b["lat"]["A"] = A
        self.b["lat"]["B"] = B
        self.b["lat"]["evals"],self.b["lat"]["evecs"] = eig(C)

        # calculate amplitude and phase angles
        self.b["lat"]["amp"] = np.sqrt(np.real(self.b["lat"]["evecs"])**2. + \
            np.imag(self.b["lat"]["evecs"])**2.)
        self.b["lat"]["phase"] = np.rad2deg(\
            np.arctan2(np.imag(self.b["lat"]["evecs"]),np.real(self.b["lat"]["evecs"])))

        
        # report for scitech paper table

        # print("$C_{{L,\\alpha}}$ & {:<.5f} & $K_{{z,\\alpha}}$ & {:<.5f} \\\\".format(self.CL_a,Kz_a))
        # print("$C_{{D,\\alpha}}$ & {:<.5f} & $K_{{x,\\alpha}}$ & {:<.5f}  \\\\".format(self.CD_a,Kx_a))
        # print("$C_{{m,\\alpha}}$ & {:<.5f} & $K_{{m,\\alpha}}$ & {:<.5f} \\\\".format(self.Cm_a,Km_a))
        # print("$C_{{L,\\hat \\alpha}}$ & {:<.5f} & $K_{{z,\\hat \\alpha}}$ & {:<.5f}  \\\\".format(self.CL_ahat,Kz_ahat))
        # print("$C_{{D,\\hat \\alpha}}$ & {:<.5f} & $K_{{x,\\hat \\alpha}}$ & {:<.5f} \\\\".format(self.CD_ahat,Kx_ahat))
        # print("$C_{{m,\\hat \\alpha}}$ & {:<.5f} & $K_{{m,\\hat \\alpha}}$ & {:<.5f} \\\\".format(self.Cm_ahat,Km_ahat))
        # print("$C_{{L,\\hat \\mu}}$ & {:<.5f} & $K_{{z,\\hat \\mu}}$ & {:<.5f} \\\\".format(self.CL_uhat,Kz_muhat))
        # print("$C_{{D,\\hat \\mu}}$ & {:<.5f} & $K_{{x,\\hat \\mu}}$ & {:<.5f} \\\\".format(self.CD_uhat,Kx_muhat))
        # print("$C_{{m,\\hat \\mu}}$ & {:<.5f} & $K_{{m,\\hat \\mu}}$ & {:<.5f} \\\\".format(self.Cm_uhat,Km_muhat))
        # print("$C_{{Y,\\beta}}$ & {:<.5f} & $K_{{y,\\beta}}$ & {:<.5f} \\\\".format(self.CY_b,Ky_b))
        # print("$C_{{\\ell,\\beta}}$ & {:<.5f} & $K_{{\\ell,\\beta}}$ & {:<.5f} \\\\".format(self.Cl_b,Kl_b))
        # print("$C_{{n,\\beta}}$ & {:<.5f} & $K_{{n,\\beta}}$ & {:<.5f} \\\\".format(self.Cn_b,Kn_b))
        # print("$C_{{Y,\\bar{{p}}}}$ & {:<.5f} & $K_{{y,\\breve p}}$ & {:<.5f}  \\\\".format(self.CY_pbar,Ky_pbreve))
        # print("$C_{{\\ell,\\bar{{p}}}}$ & {:<.5f} & $K_{{\\ell,\\breve p}}$ & {:<.5f} \\\\".format(self.Cl_pbar,Kl_pbreve))
        # print("$C_{{n,\\bar{{p}}}}$ & {:<.5f} & $K_{{n,\\breve p}}$ & {:<.5f} \\\\".format(self.Cn_pbar,Kn_pbreve))
        # print("$C_{{L,\\bar{{q}}}}$ & {:<.5f} & $K_{{z,\\breve q}}$ & {:<.5f}  \\\\".format(self.CL_qbar,Kz_qbreve))
        # print("$C_{{D,\\bar{{q}}}}$ & {:<.5f} & $K_{{x,\\breve q}}$ & {:<.5f} \\\\".format(self.CD_qbar,Kx_qbreve))
        # print("$C_{{m,\\bar{{q}}}}$ & {:<.5f} & $K_{{m,\\breve q}}$ & {:<.5f} \\\\".format(self.Cm_qbar,Km_qbreve))
        # print("$C_{{Y,\\bar{{r}}}}$ & {:<.5f} & $K_{{y,\\breve r}}$ & {:<.5f} \\\\".format(self.CY_rbar,Ky_rbreve))
        # print("$C_{{\\ell,\\bar{{r}}}}$ & {:<.5f} & $K_{{\\ell,\\breve r}}$ & {:<.5f} \\\\".format(self.Cl_rbar,Kl_rbreve))
        # print("$C_{{n,\\bar{{r}}}}$ & {:<.5f} & $K_{{n,\\breve r}}$ & {:<.5f} \\\\".format(self.Cn_rbar,Kn_rbreve))

        # calculate Short Period approximation
        self.b["spSg"] = - 0.5 * self.g / self.vo * (Kz_a + Km_qbreve + Km_ahat)
        self.b["spt99"] = np.log(0.01) / - self.b["spSg"]
        self.b["spt2x"] = np.log(2.0) / - self.b["spSg"]
        self.b["spWD"] = self.g / self.vo * np.sqrt(np.abs(\
            (Kz_a*Km_qbreve + Km_a) - 0.25*(Kz_a + Km_qbreve + Km_ahat)**2.))
        self.b["spT"] = 2. * np.pi / self.b["spWD"]
        self.b["spwn"] = ( self.b["spSg"]**2. + self.b["spWD"]**2. )**0.5
        self.b["spzt"] = self.b["spSg"] / self.b["spwn"]

        # calculate Phugoid approximation
        RPHs = Km_a / (Km_a - Kz_a * Km_qbreve)
        RPHd = Kx_a * Km_qbreve / (Km_a - Kz_a * Km_qbreve)
        RPHp = RPHs * (Kz_a + Km_qbreve) / (Km_a - Kz_a * Km_qbreve)
        self.b["phSg"] = -self.g / self.vo * (-Kz_mu/2. *(Kx_mu/-Kz_mu+RPHd-RPHp))
        self.b["pht99"] = np.log(0.01) / - self.b["phSg"]
        self.b["phWD"] = self.g / self.vo * -Kz_mu/2. * \
            np.sqrt( -4./-Kz_mu*RPHs - (Kx_mu/-Kz_mu+RPHd)**2.)
        self.b["phT"] = 2. * np.pi / self.b["phWD"]
        self.b["phwn"] = ( self.b["phSg"]**2. + self.b["phWD"]**2. )**0.5
        self.b["phzt"] = self.b["phSg"] / self.b["phwn"]

        # calculate Roll approximation
        self.b["roSg"] = - self.g / self.vo * Kl_pbreve
        self.b["rot99"] = np.log(0.01) / - self.b["roSg"]

        # calculate Spiral approximation
        self.b["slSg"] = self.g / self.vo * \
            (Kl_b*Kn_rbreve - Kl_rbreve*Kn_b) / \
            (Kl_b*Kn_pbreve - Kl_pbreve*Kn_b)
        self.b["slt99"] = np.log(0.01) / - self.b["slSg"]
        self.b["slt2x"] = np.log(2.0) / - self.b["slSg"]

        # calculate Dutch Roll approximation
        RDRs = (Kl_b*(1. - (1.-Ky_rbreve)*Kn_pbreve) - \
            Ky_b*Kl_rbreve*Kn_pbreve) / Kl_pbreve
        # RDRs = (Kl_b*1.0 - \
        #     Ky_b*Kl_rbreve*Kn_pbreve) / Kl_pbreve # # test error buildup
        RDRc = Kl_rbreve * Kn_pbreve / Kl_pbreve
        RDRp = (Kl_rbreve * Kn_b - Kl_b * Kn_rbreve) / Kl_pbreve / \
            (Kn_b + Ky_b * Kn_rbreve) - RDRs / Kl_pbreve
        # RDRp = 0.0*(Kl_rbreve * Kn_b - Kl_b * Kn_rbreve) / Kl_pbreve / \
        #     (Kn_b + Ky_b * Kn_rbreve) - RDRs / Kl_pbreve # # test error buildup
        self.b["drSg"] = - 0.5*self.g/self.vo * (Ky_b + Kn_rbreve - RDRc + RDRp)
        self.b["drt99"] = np.log(0.01) / - self.b["drSg"]
        self.b["drt2x"] = np.log(2.0) / - self.b["drSg"]
        self.b["drWD"] = self.g/self.vo * np.sqrt( (1.-Ky_rbreve)*Kn_b +\
            Ky_b * Kn_rbreve + RDRs - 0.25*( Ky_b + Kn_rbreve)**2. )
        self.b["drT"] = 2. * np.pi / self.b["drWD"]
        self.b["drwn"] = ( self.b["drSg"]**2. + self.b["drWD"]**2. )**0.5
        self.b["drzt"] = self.b["drSg"] / self.b["drwn"]

        # # print((1. - (1.-Ky_rbreve)*Kn_pbreve))
        # # print(self.Cn_rbar*self.CY_b - self.Cn_b*self.CY_rbar)
        # # print(self.Cl_rbar*self.Cn_b - self.Cl_b*self.Cn_rbar)
        # # print((self.p["drSg"]-self.b["drSg"])/self.p["drSg"])
        # # print((self.p["drWD"]-self.b["drWD"])/self.p["drWD"])
        # pWn = ( self.p["drSg"]**2. + self.p["drWD"]**2. ) ** 0.5
        # bWn = ( self.b["drSg"]**2. + self.b["drWD"]**2. ) ** 0.5
        # print((pWn-bWn)/pWn)
        # pz = self.p["drSg"] / pWn
        # bz = self.b["drSg"] / bWn
        # print((pz-bz)/pz)
        # print()


        # ratios of import
        pSwcw_4 = self.rho * self.Sw * self.cwbar / 4.
        pSwbw_4 = self.rho * self.Sw * self.bw    / 4.
        dynF = 0.5 * self.rho * self.vo**2. * self.Sw

        ### Longitudinal derivatives
        self.D_udot = pSwcw_4 * self.CD_uhat
        self.L_udot = pSwcw_4 * self.CL_uhat
        self.m_udot = pSwcw_4 * self.cwbar * self.Cm_uhat
        self.D_adot = pSwcw_4 * self.vo * self.CD_ahat
        self.L_adot = pSwcw_4 * self.vo * self.CL_ahat
        self.m_adot = pSwcw_4 * self.vo * self.cwbar * self.Cm_ahat
        self.Do = dynF * self.CDo
        self.Lo = dynF * self.CLo
        self.mo = dynF * self.cwbar * self.Cmo
        self.D_a = dynF * self.CD_a
        self.L_a = dynF * self.CL_a
        self.m_a = dynF * self.cwbar * self.Cm_a
        self.D_q = pSwcw_4 * self.vo * self.CD_qbar
        self.L_q = pSwcw_4 * self.vo * self.CL_qbar
        self.m_q = pSwcw_4 * self.vo * self.cwbar * self.Cm_qbar
        ### Lateral Derivatives
        self.Y_b = dynF * self.CY_b
        self.l_b = dynF * self.bw * self.Cl_b
        self.n_b = dynF * self.bw * self.Cn_b
        self.Y_p = pSwbw_4 * self.vo * self.CY_pbar
        self.l_p = pSwbw_4 * self.vo * self.bw * self.Cl_pbar
        self.n_p = pSwbw_4 * self.vo * self.bw * self.Cn_pbar
        self.Y_r = pSwbw_4 * self.vo * self.CY_rbar
        self.l_r = pSwbw_4 * self.vo * self.bw * self.Cl_rbar
        self.n_r = pSwbw_4 * self.vo * self.bw * self.Cn_rbar

        # other important terms
        self.l_mp = -self.m_a/self.L_a - self.m_q/self.W * self.g/self.vo
        self.l_np = -self.m_a/self.L_a
        self.ryyb = ( self.g * self.Iyy / self.W )**0.5
        self.rxxb = ( self.g * self.Ixx / self.W )**0.5
        self.RG0 = self.Lo / self.Do
        
        # simplified handling qualities
        self.b["spwnhq"] = np.abs(( self.g*self.L_a/self.W*self.l_mp )**0.5 / \
            self.ryyb)
        num = - self.ryyb/self.Iyy * (-self.ryyb**2./self.vo*self.L_a + \
            self.m_q + self.m_adot)
        den = 2. * ( self.g*self.L_a/self.W*self.l_mp )**0.5
        self.b["spzthq"] = np.abs(num / den)
        self.b["phwnhq"] = self.g/self.vo * (2.*self.l_np/self.l_mp)**0.5
        a = 1./self.RG0 * (self.l_mp/2./self.l_np)**0.5
        c = -self.ryyb**2./self.vo + self.m_q / self.L_a
        b = self.g/self.vo * (2.*self.l_np/self.l_mp**3.)**0.5 * c
        self.b["phzthq"] = a + b
        self.b["roSghq"] = -self.l_p / self.Ixx
        num = self.l_b * self.n_r - self.l_r * self.n_b
        den = self.l_b * self.n_p - self.l_p * self.n_b
        self.b["slSghq"] = self.g/self.vo * num / den
        e = self.Ixx * self.l_b * self.n_p / self.l_p**2.
        f = - self.l_r * self.n_p / self.l_p
        a = 1./self.Izz*(self.n_r + e + f)
        g = - self.Ixx * self.l_b / self.l_p**2.
        b = self.g/self.vo * ( 1./self.W*self.Y_b + g )
        Sg = -0.5*( a + b )
        i = - self.rxxb**2. * self.Y_b / 2. / self.l_p
        j = - self.Ixx / self.Izz * self.n_r / 2. / self.l_p
        h = 1. + i + j
        c = 1./self.Izz * ( self.n_b - self.l_b * self.n_p / self.l_p * h )
        k = 1. / self.W / self.Izz * self.Y_b * self.n_r
        d = self.g/self.vo * ( k + self.l_b / self.l_p * h )
        Wn = ( c + d )**0.5
        self.b["drwnhq"] = np.abs(Wn)
        self.b["drzthq"] = Sg/np.abs(Wn)

        a = (Ky_b + Kn_rbreve) / Kl_pbreve
        b = Kl_b * (1.-Kn_pbreve)/Kl_pbreve
        c = Kn_b + Ky_b * Kn_rbreve
        self.b["drwnhq"] = self.g/self.vo*( -0.5*a*b+c+b )**0.5
        a = Ky_b + Kn_rbreve
        b = Kl_rbreve * Kn_pbreve / Kl_pbreve
        c = Kl_b * (1.-Kn_pbreve)/Kl_pbreve**2.
        Sg = -0.5*self.g/self.vo*(a-b-c)
        self.b["drzthq"] = Sg/self.b["drwnhq"]

        # a = (self.Y_b/self.W + self.vo*self.n_r/self.g/self.Izz) / (self.vo*self.l_p/self.g/self.Ixx)
        # b = self.vo**2.*self.l_b/self.g**2./self.Ixx * (1.-self.vo*self.n_p/self.g/self.Izz)/(self.vo*self.l_p/self.g/self.Ixx)
        # c = self.vo**2.*self.n_b/self.g**2./self.Izz + self.Y_b/self.W * self.vo*self.n_r/self.g/self.Izz
        # self.b["drwn"] = self.g/self.vo*( -0.5*a*b+c+b )**0.5
        
        sm_l = self.Cl_b/self.CY_b
        sm_m = -self.Cm_a/self.CL_a
        sm_n = -self.Cn_b/self.CY_b
        # print(self.aircraft_name,sm_n>sm_m)
        # print(self.aircraft_name, sm_m*self.cwbar,sm_n*self.bw)
        # cy = "{:> 6.2f}{:> 6.2f}{:> 6.2f}{:> 6.2f}{:> 6.2f}".format(self.CY_b,self.CY_rbar,self.CY_pbar,self.CY_da,self.CY_dr)
        # cl = "{:> 6.2f}{:> 6.2f}{:> 6.2f}{:> 6.2f}{:> 6.2f}".format(self.Cl_b,self.Cl_rbar,self.Cl_pbar,self.Cl_da,self.Cl_dr)
        # cn = "{:> 6.2f}{:> 6.2f}{:> 6.2f}{:> 6.2f}{:> 6.2f}".format(self.Cn_b,self.Cn_rbar,self.Cn_pbar,self.Cn_da,self.Cn_dr)
        # print(cy,"  ",cl,"  ",cn)

        a = - Kl_rbreve * Kn_pbreve / Kl_pbreve
        b_num = Kl_rbreve * Kn_b - Kl_b * Kn_rbreve
        b_den = Kn_b + Ky_b * Kn_rbreve
        b = b_num / b_den / Kl_pbreve
        c = -(Kl_b*(1.-(1.-Ky_rbreve)*Kn_pbreve) - Ky_b*Kl_rbreve*Kn_pbreve)/Kl_pbreve**2.
        Sdr = -0.5 * self.g / self.vo * ( Ky_b + Kn_rbreve + a + b + c )
        d = (Kl_rbreve * Kn_b - Kl_b * Kn_rbreve)/(Kn_b + Ky_b * Kn_rbreve)
        e = -(Kl_b*(1.-(1.-Ky_rbreve)*Kn_pbreve) - Ky_b*Kl_rbreve*Kn_pbreve)/Kl_pbreve
        f = - Kl_rbreve * Kn_pbreve
        deff = d + e + f
        g = 0.5 * (Ky_b + Kn_rbreve) / Kl_pbreve * deff
        h = 0.25 / Kl_pbreve**2. * deff**2.
        i = (1.-Ky_rbreve)*Kn_b + Ky_b*Kn_rbreve - e
        wndr = self.g / self.vo * (g + h + i)**0.5
        self.b["drwnhq"] = wndr
        self.b["drzthq"] = Sdr/wndr # these are correct



        # print("A matrix")
        # for i in range(6):
        #     for j in range(6):
        #         print("\t{:>17.12f}".format(A[i,j]),end="")
        #     print()

        # print("B matrix")
        # for i in range(6):
        #     for j in range(6):
        #         print("\t{:>17.12f}".format(B[i,j]),end="")
        #     print()

        # print("C matrix")
        # for i in range(6):
        #     for j in range(6):
        #         print("\t{:>17.12f}".format(C[i,j]),end="")
        #     print()

        # print("\nLONGITUDINAL")
        # print("Eigenvalues")
        # for i in range(6):
        #     print("\t{:>17.12f}".format(self.h["lon"]["evals"][i]))

        # print("\nEigenvectors")
        # for i in range(6):
        #     for j in range(6):
        #         print("\t{:>17.12f}".format(self.h["lon"]["evecs"][i,j]),end="")
        #     print()

        # print("\nAmplitudes")
        # for i in range(6):
        #     for j in range(6):
        #         print("\t{:>17.12f}".format(self.h["lon"]["amp"][i,j]),end="")
        #     print()

        # print("\nPhase Angle (deg)")
        # for i in range(6):
        #     for j in range(6):
        #         print("\t{:>17.12f}".format(self.h["lon"]["phase"][i,j]),end="")
        #     print()

        # print("\nLATERAL")
        # print("Eigenvalues")
        # for i in range(6):
        #     print("\t{:>17.12f}".format(self.h["lat"]["evals"][i]))

        # print("Eigenvectors")
        # for i in range(6):
        #     for j in range(6):
        #         print("\t{:>17.12f}".format(self.h["lat"]["evecs"][i,j]),end="")
        #     print()

        # print("\nAmplitudes")
        # for i in range(6):
        #     for j in range(6):
        #         print("\t{:>17.12f}".format(self.h["lat"]["amp"][i,j]),end="")
        #     print()

        # print("\nPhase Angle (deg)")
        # for i in range(6):
        #     for j in range(6):
        #         print("\t{:>17.12f}".format(self.h["lat"]["phase"][i,j]),end="")
        #     print()
        # print()

        self.b["K"] = {}
        self.b["K"]["Kz,a"] = Kz_a
        self.b["K"]["Kx,a"] = Kx_a
        self.b["K"]["Km,a"] = Km_a
        self.b["K"]["Kz,ahat"] = Kz_ahat
        self.b["K"]["Kx,ahat"] = Kx_ahat
        self.b["K"]["Km,ahat"] = Km_ahat
        self.b["K"]["Kz,muhat"] = Kz_muhat
        self.b["K"]["Kx,muhat"] = Kx_muhat
        self.b["K"]["Km,muhat"] = Km_muhat
        self.b["K"]["Ky,b"] = Ky_b
        self.b["K"]["Kl,b"] = Kl_b
        self.b["K"]["Kn,b"] = Kn_b
        self.b["K"]["Ky,pbreve"] = Ky_pbreve
        self.b["K"]["Kl,pbreve"] = Kl_pbreve
        self.b["K"]["Kn,pbreve"] = Kn_pbreve
        self.b["K"]["Kz,qbreve"] = Kz_qbreve
        self.b["K"]["Kx,qbreve"] = Kx_qbreve
        self.b["K"]["Km,qbreve"] = Km_qbreve
        self.b["K"]["Ky,rbreve"] = Ky_rbreve
        self.b["K"]["Kl,rbreve"] = Kl_rbreve
        self.b["K"]["Kn,rbreve"] = Kn_rbreve


    def _eigenvalue_properties(self,method,is_longitudinal,is_alternate=False):

        # reference length
        if is_longitudinal:
            ref_len = self.cwbar
            side = "lon"
        else:
            ref_len = self.bw
            side = "lat"
        
        if is_alternate:
            redim = self.g / self.vo
        else:
            redim = 2. * self.vo / ref_len

        
        # initialize properties arrays
        method[side]["Sg"] = np.zeros(6,)
        method[side]["99damp"] = np.zeros(6,)
        method[side]["double"] = np.zeros(6,)
        method[side]["wd"] = np.zeros(6,)
        method[side]["T"] = np.zeros(6,)
        method[side]["wn"] = np.zeros(6,)
        method[side]["zt"] = np.zeros(6,)
        method[side]["dimevals"] = method[side]["evals"] * 0.0
        
        # cycle through eigenvalues
        for i in range(method[side]["evals"].shape[0]):

            # dimensionalize eigenvalue
            method[side]["dimevals"][i] = method[side]["evals"][i] * redim

            # calculate damping rate
            method[side]["Sg"][i] = - np.real(method[side]["evals"][i]) * \
                redim
            
            if np.real(method[side]["evals"][i]) != 0.0:
                # calculate 99% damping time / doubling time
                method[side]["99damp"][i] = np.log(0.01) / - \
                    method[side]["Sg"][i]
                method[side]["double"][i] = np.log(2.0) / - \
                    method[side]["Sg"][i]
            
            if np.imag(method[side]["evals"][i]) != 0.0:
                # calculate damped natural frequency, period
                method[side]["wd"][i] = \
                    np.abs(np.imag(method[side]["evals"][i])) * redim 
                method[side]["T"][i] = 2.*np.pi / method[side]["wd"][i]
                
                # calculate natural frequency and damping ratio
                e1 = method[side]["evals"][i]
                e2 = np.conjugate(e1)
                method[side]["wn"][i] = np.real(np.sqrt(e1*e2) * redim)
                method[side]["zt"][i] = np.real(-(e1+e2) /2./np.sqrt(e1*e2))


    def _longitudinal_dimensional_properties(self,method,is_alternate=False):
        """Method which determines the longitudinal dimensional properties of 
        the eigenvalues"""
        
        # determine indexes for rigid body modes
        i_rb = np.argwhere(np.abs(method["lon"]["evals"]) == 0.0).T[0]

        # determine indexes for short period mode
        i_sp = i_rb * 0
        i_sp[0] = np.abs(method["lon"]["evals"]).argmax()
        # determine "second" index
        sp_not_complex = np.imag(method["lon"]["evals"][i_sp[0]]) == 0.0
        i_sp[1] = np.abs(np.delete(method["lon"]["evals"],i_sp[0])).argmax()
        if i_sp[1] >= i_sp[0]:
            i_sp[1] += 1
        
        # determine indexes for phugoid mode
        i_ph = np.delete(np.array(range(6)),np.concatenate((i_rb,i_sp)))

        # calculate eigenvalue properties
        self._eigenvalue_properties(method,True,is_alternate)

        # save indices
        method["lon"]["rb"] = i_rb
        method["lon"]["sp"] = i_sp
        method["lon"]["ph"] = i_ph

        # # create list of modes
        # method["lon"]["modes"] = [""] * 6
        # for i in range(6):
        #     if i in i_rb:
        #         method["lon"]["modes"][i] = "rb"
        #     if i in i_sp:
        #         method["lon"]["modes"][i] = "sp"
        #     if i in i_ph:
        #         method["lon"]["modes"][i] = "ph"

        a = 10


    def _lateral_dimensional_properties(self,method,is_alternate=False):
        """Method which determines the lateral dimensional properties of 
        the eigenvalues"""

        # determine the rigid body modes
        i_rb = np.argwhere(np.abs(method["lat"]["evals"]) == 0.0).T[0]
        
        # determine roll index
        i_ro = np.abs(method["lat"]["evals"].real).argmax()

        # determine the other two
        if method["lat"]["evals"][i_ro].imag == 0.0:
            # dutch roll must be the only complex evals
            i_dr = np.argwhere(np.abs(method["lat"]["evals"].imag) != 0.0).T[0]

            # spiral remains
            i_sl = np.delete(np.array(range(6)),\
                np.concatenate((i_rb,i_dr,np.array([i_ro]))))[0]
        else:
            # spiral is complex conjugate
            i_sl = np.argwhere(method["lat"]["evals"] == \
                np.conj(method["lat"]["evals"][i_ro])).T[0,0]

            # dutch roll remains
            i_dr = np.delete(np.array(range(6)),\
                np.concatenate((i_rb,np.array([i_ro,i_sl]))))
        
        # raise error if modes not properly identified
        if len(i_dr) > 2:
            raise IndexError("Dutch Roll mode corresponds to two indices, " + \
                "not ",i_dr.tolist())

        # calculate eigenvalue properties
        self._eigenvalue_properties(method,False,is_alternate)

        # save indices
        method["lat"]["rb"] = i_rb
        method["lat"]["ro"] = i_ro
        method["lat"]["sl"] = i_sl
        method["lat"]["dr"] = i_dr
        

    def _report(self,method1,method2,method3,is_longitudinal,show_evector=False):

        if is_longitudinal:
            typename = "LONGITUDINAL"
            si = "lon"
            evec_sym = ["\u0394\u03BC","\u0394\u03B1","\u0394q",
            "\u0394\u03BEx","\u0394\u03BEz","\u0394\u03B8"]
            approx_modes = ["sp","ph"]
        else:
            typename = "LATERAL"
            si = "lat"
            evec_sym = ["\u0394\u03B2","\u0394p","\u0394r",
            "\u0394\u03BEy","\u0394\u03C6","\u0394\u03C8"]
            approx_modes = ["sl","ro","dr"]

        # print type
        num = 133#58
        print("{:*^{}}".format(" " + typename + " ",num+2))
        print("{:*^{}}".format("",num+2))

        dimeval = "Dimensionless Eigenvalue = "
        dimevec = "Dimensionless Eigenvector:"
        colwidth = 35
        decwidth = 12
        rform = " {:> " + "{}".format(colwidth) +"."+"{}".format(decwidth)+"f}"

        # string titles
        ti_sigma = "Damping Rate [1/sec] ="
        ti_99dam = "99% Damping Time [sec] ="
        ti_doubl = "Doubling Time [sec] ="
        ti_wd    = "Damp Nat Freq [rad/sec] ="
        ti_T     = "Period [sec] ="
        ti_wn    = "Natural Freq [rad/sec] ="
        ti_zeta  = "Damping Ratio ="
        modenames = {
            "sp" : "Short Period",
            "ph" : "Phugoid",
            "sl" : "Spiral",
            "ro" : "Roll",
            "dr" : "Dutch Roll"
        }

        # cycle through each eigenvalue, report
        for i in range(6):
            print(("{:=^{}}"+" {:-^{}}"*3).format("",len(dimeval),\
                " " + method1["name"] + " ",colwidth,\
                " " + method2["name"] + " ",colwidth,\
                " " + method3["name"] + " ",colwidth))

            # report eigenvalue
            print((dimeval + rform*3 + "\n").format(method1[si]["evals"][i],\
                method2[si]["evals"][i],method3[si]["evals"][i]))
            
            # bool statements
            m1real = method1[si]["evals"][i].real != 0.0
            m1imag = method1[si]["evals"][i].imag != 0.0
            m1conv = method1[si]["evals"][i].real <  0.0
            m2real = method2[si]["evals"][i].real != 0.0
            m2imag = method2[si]["evals"][i].imag != 0.0
            m2conv = method2[si]["evals"][i].real <  0.0
            m3real = method3[si]["evals"][i].real != 0.0
            m3imag = method3[si]["evals"][i].imag != 0.0
            m3conv = method3[si]["evals"][i].real <  0.0
            anyreal = m1real or m2real or m3real
            anyimag = m1imag or m2imag or m3imag
            anyconv = m1conv or m2conv or m3conv
            anydive = not m1conv or not m2conv or not m3conv
            
            # report properties if necessary
            if anyreal:
                # damping rate
                print("{:<{}s}".format(ti_sigma,len(dimeval)),end="")
                print((rform*3).format(method1[si]["Sg"][i],\
                    method2[si]["Sg"][i],method3[si]["Sg"][i]))

                # 99% damping time
                if anyconv:
                    print("{:<{}s}".format(ti_99dam,len(dimeval)),end="")
                    print((rform*3).format(method1[si]["99damp"][i],\
                        method2[si]["99damp"][i],method3[si]["99damp"][i]))
                # doubling time
                if anydive:
                    print("{:<{}s}".format(ti_doubl,len(dimeval)),end="")
                    print((rform*3).format(method1[si]["double"][i],\
                        method2[si]["double"][i],method3[si]["double"][i]))

                # complex
                if anyimag:
                    # wd
                    print("{:<{}s}".format(ti_wd,len(dimeval)),end="")
                    print((rform*3).format(method1[si]["wd"][i],\
                        method2[si]["wd"][i],method3[si]["wd"][i]))

                    # T
                    print("{:<{}s}".format(ti_T,len(dimeval)),end="")
                    print((rform*3).format(method1[si]["T"][i],\
                        method2[si]["T"][i],method3[si]["T"][i]))

                    # wn
                    print("{:<{}s}".format(ti_wn,len(dimeval)),end="")
                    print((rform*3).format(method1[si]["wn"][i],\
                        method2[si]["wn"][i],method3[si]["wn"][i]))
                    
                    # zeta
                    print("{:<{}s}".format(ti_zeta,len(dimeval)),end="")
                    print((rform*3).format(method1[si]["zt"][i],\
                        method2[si]["zt"][i],method3[si]["zt"][i]))
                
                print()
                
            # report eigenvector
            if show_evector:
                print(dimevec,end="")
                print(3*" {:^{}s} {:^{}s}".format("Real",17,"Imaginary",17))
                
                for j in range(6):
                    print(("{:^{}s}" + " {:> 17.12f}"*6).format(evec_sym[j],\
                        len(dimevec),method1[si]["evecs"][j,i].real,\
                        method1[si]["evecs"][j,i].imag,\
                        method2[si]["evecs"][j,i].real,\
                        method2[si]["evecs"][j,i].imag,\
                        method3[si]["evecs"][j,i].real,\
                        method3[si]["evecs"][j,i].imag))
                
                # report eigenvector
                print(len(dimevec)*" " + 3*" {:^{}s} {:^{}s}".format(\
                    "Amplitude",17,"Phase[deg]",17))
                
                for j in range(6):
                    print(("{:^{}s}" + " {:> 17.12f}"*6).format(evec_sym[j],\
                        len(dimevec),method1[si]["amp"][j,i],\
                        method1[si]["phase"][j,i],method2[si]["amp"][j,i],\
                        method2[si]["phase"][j,i],method3[si]["amp"][j,i],\
                        method3[si]["phase"][j,i]))

                # new line at end
                print()
        
        # report approximations
        title = typename + " APPROXIMATIONS"
        print("{:=^{}}".format(" " + title + " ",num+2))

        for mo in approx_modes:
            # header
            print(("{:=^{}}"+" {:-^{}}"*3).format(" " + modenames[mo] + " ",\
                len(dimeval)," " + method1["name"] + " ",colwidth,\
                " " + method2["name"] + " ",colwidth,\
                " " + method3["name"] + " ",colwidth))
            # mode name
            # print("{:-^{}}".format(,num+2))

            # damping rate
            print("{:<{}s}".format(ti_sigma,len(dimeval)),end="")
            app = "Sg"
            print((rform*3).format(method1[mo + app],\
                method2[mo + app],method3[mo + app]))

            # 99% damping time
            app = "t99"
            if method1[mo + app] > 0.:
                print("{:<{}s}".format(ti_99dam,len(dimeval)),end="")
                print((rform*3).format(method1[mo + app],\
                    method2[mo + app],method3[mo + app]))

            # doubling time
            app = "t2x"
            if mo != "ro" and method1[mo + app] > 0.:
                print("{:<{}s}".format(ti_doubl,len(dimeval)),end="")
                print((rform*3).format(method1[mo + app],\
                    method2[mo + app],method3[mo + app]))

            # damped natural frequency
            app = "WD"
            if mo + app in method1 and mo + app in method2 and mo + app in method3:
                print("{:<{}s}".format(ti_wd,len(dimeval)),end="")
                print((rform*3).format(method1[mo + app],\
                    method2[mo + app],method3[mo + app]))

            # period
            app = "T"
            if mo + app in method1 and mo + app in method2 and mo + app in method3:
                print("{:<{}s}".format(ti_T,len(dimeval)),end="")
                print((rform*3).format(method1[mo + app],\
                    method2[mo + app],method3[mo + app]))
        
        print()

        # report eigenvalues both dimensional and nondimensional for easy use in latex
        print("{:+^{}}".format(" " + typename + " EIGENVALUES" + " ",num+2))
        dec = " {:> 16.12f} "
        comp = dec+"$\\pm$"+dec+"$i$"
        real = dec+"  +  "+dec+"$i$"

        print("{:=^{}}".format(" NonDimensional ",93))
        print((" {:-^{}}"*2).format(" " + method2["name"] + " ",45,\
            " " + method3["name"] + " ",46))

        for i in range(6):
            if method2[si]["evals"][i].imag == 0.0:
                string = real
            else:
                string = comp
            
            print((string + " &").format(method2[si]["evals"][i].real,\
                method2[si]["evals"][i].imag),end="")
            
            if method3[si]["evals"][i].imag == 0.0:
                string = real
            else:
                string = comp
            
            print((string + " \\\\").format(method3[si]["evals"][i].real,\
                method3[si]["evals"][i].imag))

        print("{:=^{}}".format(" Dimensional ",93))
        print((" {:-^{}}"*2).format(" " + method2["name"] + " ",45,\
            " " + method3["name"] + " ",46))

        for i in range(6):
            if method2[si]["dimevals"][i].imag == 0.0:
                string = real
            else:
                string = comp
            
            print((string + " &").format(method2[si]["dimevals"][i].real,\
                method2[si]["dimevals"][i].imag),end="")
            
            if method3[si]["dimevals"][i].imag == 0.0:
                string = real
            else:
                string = comp
            
            print((string + " \\\\").format(method3[si]["dimevals"][i].real,\
                method3[si]["dimevals"][i].imag))
        
        print("\n")


    def _full_report(self,show_evector=False):
        """Method used to report eigenpair properties to the screen.
        """

        # report longitudinal, lateral
        is_longitudinal = [True,False]
        for typepos in is_longitudinal:
            self._report(self.p,self.h,self.b,typepos,show_evector)
            
        # # C matrix
        # C = self.h["lon"]["A"] / self.b["lon"]["A"]
        # D = self.h["lon"]["B"] / self.b["lon"]["B"]
        # E = self.h["lat"]["A"] / self.b["lat"]["A"]
        # F = self.h["lat"]["B"] / self.b["lat"]["B"]

        # print("C matrix")
        # for i in range(6):
        #     for j in range(6):
        #         print("\t{:>17.12f}".format(C[i,j]),end="")
        #     print()

        # print("D matrix")
        # for i in range(6):
        #     for j in range(6):
        #         print("\t{:>17.12f}".format(D[i,j]),end="")
        #     print()

        # print("E matrix")
        # for i in range(6):
        #     for j in range(6):
        #         print("\t{:>17.12f}".format(E[i,j]),end="")
        #     print()

        # print("F matrix")
        # for i in range(6):
        #     for j in range(6):
        #         print("\t{:>17.12f}".format(F[i,j]),end="")
        #     print()

        a = 10



def print_hq(hq,print_bold=False):
    if type(hq) != str:
        string = " \& {:> 8.4f}".format(hq)
    else:
        string = " \&  {:^7s}".format(hq)
    if print_bold:
        print("\033[1m" + string + "\033[0m", end="")
    else:
        print(string, end="")
    


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
    # plt.rcParams["mathtext.fontset"] = "dejavuserif"

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
    # "Cessna_172.json",
    # "navion.json",
    "VZ_4.json", "F_2B.json",
    ### Class II -- medium-weight, low-to-medium maneuverability
    # "X_15.json", "HL_10.json", "Lockheed_Jetstar.json", "Convair_880M.json", 
    # "F_105B.json",
    "C_47.json", "XC_142.json",
    ### Class III -- large, heavy, low-to-medium maneuverability
    # "boeing_747.json", "C_5A.json", "XB_70A.json",
    # "DC_8.json",
    # "9_8_2.json"
    ### Class IV -- high-maneuverability
    # "F16_bolander.json",
    # "NT_33A.json", "F_104A.json", "F_4C.json",
    # "A_7A.json", "A_4D.json",
    "F_94A.json", "F_15.json",
    ]
    run_files = ["aircraft_database/" + i for i in run_files]
    num_craft = len(run_files)
    eigensolved = [0.0] * num_craft
    for i in range(num_craft):
        eigensolved[i] = Solver(run_files[i],report=False)

    # create dictionaries for pretty title-ing of print out
    modes = ["sp","ph","sl","ro","dr"]
    sides = ["lon","lon","lat","lat","lat"]
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
    evecnames = {
        "lon" : [r"$\Delta \mu$", r"$\Delta \alpha$", r"$\Delta \breve{q}$", 
            r"$\Delta \breve{x}$", r"$\Delta \breve{z}$", r"$\Delta \theta$"],
        "lat" : [r"$\Delta \beta$", r"$\Delta \breve{p}$", r"$\Delta \breve{r}$", 
            r"$\Delta \breve{y}$", r"$\Delta \phi$", r"$\Delta \psi$"]
    }
    cs = ["#F5793A","#A95AA1","#85C0F9","#0F2080"]
    clrs = {
        "lon" : [cs[0],cs[0],cs[1],cs[2],cs[2],cs[3]],
        "lat" : [cs[0],cs[1],cs[1],cs[2],cs[3],cs[3]]
    }
    ms = ["o","^","s"]
    mrks = {
        "lon" : [ms[0],ms[2],ms[1],ms[0],ms[2],ms[1]],
        "lat" : [ms[1],ms[0],ms[2],ms[1],ms[0],ms[2]]
    }
    # base directory
    bd = "phasor_plots/"
    # show plots
    show = False
    print("plotting...")

    # for plotting    
    num = 100
    t = np.linspace(0.0,4.*np.pi,num=num)
    r = np.linspace(0.0,2.0,num=num)

    for i in range(num_craft):
        dyn = eigensolved[i]
        print(" "*11 + dyn.aircraft_name + "...")
        for j in range(len(modes)):
            method = dyn.b
            si = sides[j]
            mo = modes[j]
            direc = bd + "phasor_" + modenames[mo].lower().replace(" ","_")+"/"
            file = direc + mo + "_"+ dyn.aircraft_name.lower().replace(" ","_")
            file = file.replace("-","_") + ".png"
            # print(file)
            mbind = dyn.b[si][mo]
            if mo in ["sp","ph","dr"]:
                mbind = mbind[0]
            mpind = dyn.p[si][mo]
            if mo in ["sp","ph","dr"]:
                mpind = mpind[0]

            # print(dyn.b["lon"]["evals"])
            A_b = dyn.b[si][ "amp" ][:, mbind]
            P_b = np.deg2rad(dyn.b[si]["phase"][:, mbind])
            A_p = dyn.p[si][ "amp" ][:, mpind]
            P_p = np.deg2rad(dyn.p[si]["phase"][:, mpind])
            # print(A_b)
            # print(np.rad2deg(P_b))

            # determine phase difference between methods
            b_p_diff = P_b[0] - P_p[0]
            P_b = P_b - b_p_diff

            # plot
            fig, ax = plt.subplots(1,2,figsize=(10,5),
                subplot_kw={"projection" : "polar"})
            for k in range(len(A_b)):
                ax[1].plot([0,P_b[k]],[0,A_b[k]],c=clrs[si][k],
                    marker=mrks[si][k],label=evecnames[si][j])
                ax[0].plot([0,P_p[k]],[0,A_p[k]],c=clrs[si][k],
                    marker=mrks[si][k],mfc="none")
            # ax.set_rmax(2)
            # ax.set_rticks([0.5, 1, 1.5, 2])  # Less radial ticks
            ax[0].set_rlabel_position(60.0)
            ax[0].grid(True)
            ax[0].set_xticks(np.pi/180. * np.linspace(0.0, 360.0, 16, endpoint=False))
            ax[0].set_xlabel("Phase Angle [deg]")
            ax[0].text(np.deg2rad(80.0),0.4*max(A_p),"Amplitude",rotation=60.0)
            ax[0].set_title("Traditional Dimensionless Form")
            ax[1].set_rlabel_position(60.0)
            ax[1].grid(True)
            ax[1].set_xticks(np.pi/180. * np.linspace(0.0, 360.0, 16, endpoint=False))
            ax[1].set_xlabel("Phase Angle [deg]")
            ax[1].text(np.deg2rad(80.0),0.4*max(A_b),"Amplitude",rotation=60.0)
            ax[1].set_title("Alternate Dimensionless Form")
            ax[1].legend(bbox_to_anchor=(-0.05, 1.0), loc='upper right')
            fig.suptitle(dyn.aircraft_name + " " + modenames[mo])
            plt.tight_layout()
            plt.savefig(file)#,dpi=300.0)
            if show:
                plt.show()
            else:
                plt.close()














    quit()

    # report handling qualities
    max_name = max([len(plane.aircraft_name) for plane in eigensolved]) + 2
    side = ["lon"] * 4 + ["lat"] * 4
    mode = ["sp"] * 2 + ["ph"] * 2 + ["ro","sl"] + ["dr"] * 2
    haqu = ["wn","zt"] * 2 + ["Sg"] * 2 + ["wn","zt"]
    column_header = ["wt/ft","linear","approx","redimmd","pd-a&b"]
    for i in range(len(side)):
        mo = mode[i]
        si = side[i]
        hq = haqu[i]
        print(modenames[mo],propertynames[hq])

        print("{:^{}}".format("aircraft",max_name),end="")
        for i in range(len(column_header)):
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
            
            app = eigensolved[j].b[mo+hq]
            buck = eigensolved[j].b[mo+hq+"hq"]
            if app != 0.0:
                pd = (app-buck)/app
                if abs(pd) >= 0.2:
                    print_bold = True
                else:
                    print_bold = False
            else:
                pd = "---"
                print_bold=False
            print_hq(pd,print_bold)

            print(" \\\\")
        print()




