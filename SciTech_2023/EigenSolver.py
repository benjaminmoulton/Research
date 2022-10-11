import numpy as np
import json
import matplotlib.pyplot as plt
from scipy.linalg import eig

class Solver:
    """A class which uses solves an eigenproblem based on an aircraft 
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
    def __init__(self,input_vars={}):

        # report
        print("Running EigenSolver written by Ben Moulton\n")

        # get info or raise error
        self.input = input_vars
        self._get_input_vars(input_vars)

        # retrieve info
        self._retrieve_info()

        # calculate eigenvalues traditional way
        self._eigs_traditional()

        # calculate eigenvalues alternate way
        self._eigs_alternate()

        # report to user
        self._report()


    def _get_input_vars(self,input_vars):
        # get info or raise error

        # determine if the input_vars is a file or a dictionary
        input_vars_type = type(input_vars)

        # dictionary
        if input_vars_type == dict:
            self.input_dict = input_vars
        
        # json file
        elif input_vars_type == str and input_vars.split(".")[-1] == "json":
            self.input_dict = self._get_json(input_vars)

        # raise error
        else:
            raise TypeError("input_vars must be json file path, or " + \
                "dictionary, not {0}".format(input_vars_type))
    

    def _get_json(self,file_path):
        # import json file from file path
        json_string = open(file_path).read()

        # save to vals dictionary
        input_dict = json.loads(json_string)
        
        return input_dict


    def _retrieve_info(self):
        """Method which retrieves the information and stores it globally.
        """
        
        # store variables from file input dictionary
        
        # store aircraft input values
        aircraft = self.input_dict.get("aircraft",{})
        self.Sw = aircraft.get("wing_area")
        self.b = aircraft.get("wing_span")
        self.cbar = self.Sw / self.b
        
        # store flight input values
        flight = self.input_dict.get("flight",{})
        self.fl_V     = flight.get("air_speed")
        self.fl_rho   = flight.get("density")
        self.fl_W     = flight.get("weight")
        self.fl_theta = flight.get("climb")
        self.fl_phi   = flight.get("bank")
        
        # store reference input values
        ref = self.input_dict.get("reference",{})
        self.l       = ref.get("length")
        self.V       = ref.get("airspeed")
        self.rho     = ref.get("density")
        self.W       = ref.get("weight")
        self.Ixx     = ref.get("Ixx")
        self.Iyy     = ref.get("Iyy")
        self.Izz     = ref.get("Izz")
        self.Ixy     = ref.get("Ixy")
        self.Ixz     = ref.get("Ixz")
        self.Iyz     = ref.get("Iyz")
        self.hx      = ref.get("hx")
        self.hy      = ref.get("hy")
        self.hz      = ref.get("hz")
        self.CD0     = ref.get("CD0")
        self.CD1     = ref.get("CD1")
        self.CD2     = ref.get("CD2")
        self.CL_M    = ref.get("CL_M")
        self.CD_M    = ref.get("CD_M")
        self.Cm_M    = ref.get("Cm_M")
        self.CL_a    = ref.get("CL_a")
        self.Cl_a    = ref.get("Cl_a")
        self.Cm_a    = ref.get("Cm_a")
        self.Cn_a    = ref.get("Cn_a")
        self.CL_ahat = ref.get("CL_ahat")
        self.CD_ahat = ref.get("CD_ahat")
        self.Cm_ahat = ref.get("Cm_ahat")
        self.CL_uhat = ref.get("CL_uhat")
        self.CD_uhat = ref.get("CD_uhat")
        self.Cm_uhat = ref.get("Cm_uhat")
        self.CY_b    = ref.get("CY_b")
        self.Cl_b    = ref.get("Cl_b")
        self.Cm_b    = ref.get("Cm_b")
        self.Cn_b    = ref.get("Cn_b")
        self.CY_p    = ref.get("CY_p")
        self.Cl_p    = ref.get("Cl_p")
        self.Cn_p    = ref.get("Cn_p")
        self.CL_q    = ref.get("CL_q")
        self.CD_q    = ref.get("CD_q")
        self.Cm_q    = ref.get("Cm_q")
        self.CY_r    = ref.get("CY_r")
        self.Cl_r    = ref.get("Cl_r")
        self.Cn_r    = ref.get("Cn_r")

        # calculate new values from reference values (for alt nondim)
        F = 0.5 * self.rho * self.V**2. * self.Sw
        Mc = F * self.cbar
        Mb = F * self.b
        self.D0 = self.CD0 * F
        self.D1 = self.CD1 * F
        self.D2 = self.CD2 * F
        self.L_M = self.CL_M * F
        self.D_M = self.CD_M * F
        self.m_M = self.Cm_M * Mc
        self.L_a = self.CL_a * F
        self.l_a = self.Cl_a * Mb
        self.m_a = self.Cm_a * Mc
        self.n_a = self.Cn_a * Mb
        self.L_ahat = self.CL_ahat * F
        self.D_ahat = self.CD_ahat * F
        self.m_ahat = self.Cm_ahat * Mc
        self.L_uhat = self.CL_uhat * F
        self.D_uhat = self.CD_uhat * F
        self.m_uhat = self.Cm_uhat * Mc
        self.Y_b = self.CY_b * F
        self.l_b = self.Cl_b * Mb
        self.m_b = self.Cm_b * Mc
        self.n_b = self.Cn_b * Mb
        self.Y_p = self.CY_p * F
        self.l_p = self.Cl_p * Mb
        self.n_p = self.Cn_p * Mb
        self.L_q = self.CL_q * F
        self.D_q = self.CD_q * F
        self.m_q = self.Cm_q * Mc
        self.Y_r = self.CY_r * F
        self.l_r = self.Cl_r * Mb
        self.n_r = self.Cn_r * Mb

        # set gravity
        self.g = 32.17#4


    def _eigs_traditional(self):
        """Method which calculates the eigenvectors and values for the
        traditional setup of equations, both longitudinal and lateral.
        """

        # initial calcs to be used in components
        wg = self.W / self.g
        A = self.rho * self.Sw * self.cbar / 4. / self.W * self.g
        B = self.rho * self.Sw * (self.cbar**3.) / 8. / self.Iyy
        CL0 = self.W * 2. / self.rho / (self.V**2.) / self.Sw
        CD0 = self.CD0 + self.CD1 * CL0 + self.CD2 * CL0**2.
        CD_a = self.CD1 * self.CL_a  +  2. * self.CD2 * CL0 * self.CL_a
        # print(CD_a)
        # # print(CD0,self.CD0)
        CD0 = 0.05; CD_a = 0.35
        th = 0.0
        st = np.sin(th)
        ct = np.cos(th)

        # calculate components
        Rgx = self.g * self.cbar / 2. / (self.V**2.)
        Rz_ah = A * -self.CL_ahat
        Rm_ah = B * self.Cm_ahat
        Rx_mu = A * -2. * CD0
        Rz_mu = A * -2. * CL0
        Rm_mu = 0.0
        Rx_a = A * (CL0 - CD_a) 
        Rz_a = A * (-self.CL_a - CD0)
        Rm_a = B * self.Cm_a
        Rx_qbar = A * (-self.CD_q)
        Rz_qbar = A * (-self.CL_q)
        Rm_qbar = B * self.Cm_q
        ixz = self.Ixz / self.Ixx
        izx = self.Ixz / self.Izz
        Av = self.rho * self.Sw * self.b / 4. / self.W * self.g
        Bv = self.rho * self.Sw * (self.b**3.) / 8. / self.Ixx
        Cv = self.rho * self.Sw * (self.b**3.) / 8. / self.Izz
        Ry_B = Av * self.CY_b
        Rl_B = Bv * self.Cl_b
        Rn_B = Cv * self.Cn_b
        Ry_pbar = Av * self.CY_p
        Rl_pbar = Bv * self.Cl_p
        Rn_pbar = Cv * self.Cn_p
        Ry_rbar = Av * self.CY_r
        Rl_rbar = Bv * self.Cl_r
        Rn_rbar = Cv * self.Cn_r
        Rgy = self.g * self.b / 2. / (self.V**2.)

        # report short period approximate eigenvalue
        rp = (Rz_a + Rm_qbar + Rm_ah) / 2.
        ip = ( (Rz_a * Rm_qbar - Rm_a) - rp**2. )**0.5
        self.l_sp = np.complex(rp,ip)
        # print(ip)
        # print(l_sp)
        # print("sp {:>14.12f}\u00B1{:>14.12f}".format(np.real(l_sp),\
        # np.abs(np.imag(l_sp))))

        # report phugoid approximate eigenvalue
        Rs = Rm_a / (Rm_a - Rz_a * Rm_qbar)
        Rd = Rx_a * Rm_qbar / (Rm_a - Rz_a * Rm_qbar)
        Rp = Rgx * Rs * (Rz_a + Rm_qbar) / (Rm_a - Rz_a * Rm_qbar)
        rp = Rz_mu / 2. * (Rx_mu / Rz_mu + Rd - Rp)
        ip = Rz_mu / 2. * ( -4. *Rgx / Rz_mu *Rs - ( Rx_mu / Rz_mu + Rd )**2. )**0.5
        self.l_ph = np.complex(rp,ip)
        # print(l_ph)
        # print("ph {:>14.12f}\u00B1{:>14.12f}".format(np.real(l_ph),\
        # np.abs(np.imag(l_ph))))

        # calculate roll approx eigenvalue
        self.l_r = np.complex(Rl_pbar,0.0)

        # calculate spiral approx eigenvalue
        rp = - Rgy * ( Rl_B * Rn_rbar - Rl_rbar * Rn_B ) / \
            ( Rl_B * Rn_pbar - Rl_pbar * Rn_B )
        self.l_s = np.complex(rp,0.0)

        # calculate dutch roll approx eigenvalue
        RDs = (Rl_B * (Rgy - (1. - Ry_rbar) * Rn_pbar) - Ry_B * Rl_rbar * \
            Rn_pbar) / Rl_pbar
        RDc = Rl_rbar * Rn_pbar / Rl_pbar
        RDp = Rgy * (Rl_rbar * Rn_B - Rl_B * Rn_rbar) / Rl_pbar / \
            (Rn_B + Ry_B * Rn_rbar) - RDs / Rl_pbar
        rp = (Ry_B + Rn_rbar - RDc + RDp) / 2.
        ip = ( (1.-Ry_rbar)*Rn_B + Ry_B*Rn_rbar + RDs - \
            ( (Ry_B+Rn_rbar)/2. )**2. )**0.5
        self.l_dr = np.complex(rp,ip)

        # print(wg)
        # print(self.cbar)
        # print(A)
        # print(B)
        # print()
        # print(Rgx)
        # print(CL0)
        # print(Rz_ah)
        # print(Rm_ah)
        # print()
        # print(Rx_mu)
        # print(Rz_mu)
        # print(Rm_mu)
        # print(Rx_a)
        # print()
        # print(Rz_a)
        # print(Rm_a)
        # print(Rx_qbar)
        # print(Rz_qbar)
        # print(Rm_qbar)
        # print()
        # print()

        # calculate A,B,C and D matrices
        A = [[Rx_mu, Rx_a, Rx_qbar,  0, 0, -Rgx*ct],
            [Rz_mu,   Rz_a, 1+Rz_qbar,  0, 0, -Rgx*st],
            [Rm_mu, Rm_a, Rm_qbar, 0, 0, 0],
            [ct, st, 0, 0, 0, -st],
            [-st, ct, 0, 0, 0, -ct],
            [0, 0, 1, 0, 0, 0]]
        A = np.array(A)

        B = [[1, 0, 0,  0, 0, 0],
            [0, 1-Rz_ah, 0,  0, 0, 0],
            [0, -Rm_ah, 1, 0, 0, 0],
            [0, 0, 0, 1, 0, 0],
            [0, 0, 0, 0, 1, 0],
            [0, 0, 0, 0, 0, 1]]
        B = np.array(B)
        
        # # calculate A,B,C and D matrices
        # # Example 8.2.1
        # A = [[-0.00071, 0.0003, 0,  0, 0, -0.00278],
        #     [-0.00556,   -0.03151, .97309,  0, 0, -Rgx*st],
        #     [0, -0.0022, -0.03212, 0, 0, 0],
        #     [1, 0, 0, 0, 0,  0],
        #     [0, 1, 0, 0, 0, -1],
        #     [0, 0, 1, 0, 0,  0]]
        # A = np.array(A)

        # B = [[1, 0, 0,  0, 0, 0],
        #     [0, 1.01133, 0,  0, 0, 0],
        #     [0, 0.01404, 1, 0, 0, 0],
        #     [0, 0, 0, 1, 0, 0],
        #     [0, 0, 0, 0, 1, 0],
        #     [0, 0, 0, 0, 0, 1]]
        # B = np.array(B)

        # print(A)
        # print(B)
        
        C = [[Ry_B, Ry_pbar, Ry_rbar-1,  0, Rgy*ct, 0],
            [Rl_B,   Rl_pbar, Rl_rbar,  0, 0, 0],
            [Rn_B, Rn_pbar, Rn_rbar, 0, 0, 0],
            [1, 0, 0, 0, 0, ct],
            [0, 1, np.tan(th), 0, 0, 0],
            [0, 0, 1./ct, 0, 0, 0]]
        C = np.array(C)

        D = [[1, 0, 0,  0, 0, 0],
            [0, 1, -ixz,  0, 0, 0],
            [0, -izx, 1, 0, 0, 0],
            [0, 0, 0, 1, 0, 0],
            [0, 0, 0, 0, 1, 0],
            [0, 0, 0, 0, 0, 1]]
        D = np.array(D)

        # make eigenvalue calculations
        lonval, lonvec = eig(A,B)
        latval, latvec = eig(C,D)

        # save globally
        self.trad_lon_val = lonval
        self.trad_lon_vec = lonvec
        self.trad_lat_val = latval
        self.trad_lat_vec = latvec


    def _eigs_alternate(self):
        """Method which calculates the eigenvectors and values for the
        alternate setup of equations, both longitudinal and lateral.
        """

        # initial calcs to be used in components
        wg = self.W / self.g
        A = self.g / self.V / self.W
        B = self.V / self.g / self.Iyy
        CL0 = self.W * 2. / self.rho / (self.V**2.) / self.Sw
        L0 = self.W
        D0 = (self.CD0 + self.CD1 * CL0 + self.CD2 * CL0**2.) * \
            0.5 * self.rho * self.V**2. * self.Sw
        D_a = (self.CD1 * self.CL_a  +  2. * self.CD2 * CL0 * self.CL_a) * \
            0.5 * self.rho * self.V**2. * self.Sw
        L_a = self.CL_a * 0.5 * self.rho * self.V**2. * self.Sw
        # # print(CD0,self.CD0)
        # CD0 = 0.05; CD_a = 0.35
        th = 0.0
        st = np.sin(th)
        ct = np.cos(th)

        # calculate components
        Kz_ah = A * self.L_ahat
        Km_ah = B * self.m_ahat
        Kx_mu = -2. * D0 / self.W
        Kz_mu =  2. * L0 / self.W 
        Km_mu = 0.0
        Kx_a = (L0 - D_a) / self.W
        Kz_a = -(D0 + L_a) / self.W
        Km_a = B * self.V / self.g * self.m_a
        Kx_qbar = A * self.D_q
        Kz_qbar = A * self.L_q
        Km_qbar = B * self.m_q
        txz = self.Ixz / self.Ixx
        tzx = self.Ixz / self.Izz
        
        Av = A
        Bv = self.V**2. / self.g**2. / self.Ixx
        Cv = self.V**2. / self.g**2. / self.Izz

        Ky_B = self.Y_b / self.W
        Kl_B = Bv * self.l_b
        Kn_B = Cv * self.n_b
        Ky_pbar = Av * self.Y_p
        Kl_pbar = self.V / self.g / self.Ixx * self.l_p
        Kn_pbar = self.V / self.g / self.Izz * self.n_p
        Ky_rbar = Av * self.Y_r
        Kl_rbar = self.V / self.g / self.Ixx * self.l_r
        Kn_rbar = self.V / self.g / self.Izz * self.n_r

        # # report short period approximate eigenvalue
        # ap = (Kz_a + Km_qbar + Km_ah) / 2.
        # bp = ( ap**2. - (Kz_a * Km_qbar - Km_a) )**0.5
        # l_sp = self.g * self.cbar / 2. / self.V**2. * np.complex(ap,bp)
        # print(bp)
        # print(l_sp)
        a = 1. + Kz_ah
        b = -Kz_a - Km_qbar * (1.+Kz_ah) + Km_ah*Kz_qbar
        c = Kz_a*Km_qbar - Kz_qbar*Km_a
        l_sp1 = (-b + (b**2. - 4.*a*c)**0.5) / (2.*a)
        l_sp2 = (-b - (b**2. - 4.*a*c)**0.5) / (2.*a)
        print(l_sp1)
        print(l_sp2)
        

        # # report for SciTech 2022 paper
        # print("& $K_{x,\\alpha}$ & "+"{:<8.5f} \\\\".format(Kx_a))
        # print("& $\imath_{xz}$ & "+"{:<8.5f} \\\\".format(txz))
        # print("& $\imath_{zx}$ & "+"{:<8.5f} \\\\".format(tzx))
        # print("& $K_{z,\\alpha}$ & "+"{:<8.5f} \\\\".format(Kz_a))
        # print("& $K_{m,\\alpha}$ & "+"{:<8.5f} \\\\".format(Km_a))
        # print()
        # print("& $K_{z,\\hat \\alpha}$ & "+"{:<8.5f} \\\\".format(Kz_ah))
        # print("& $K_{x,\\hat \\alpha}$ & "+"0 \\\\")
        # print("& $K_{m,\\hat \\alpha}$ & "+"{:<8.5f} \\\\".format(Km_ah))
        # print()
        # print("& $K_{z,\\hat \\mu}$ & "+"0 \\\\")
        # print("& $K_{x,\\hat \\mu}$ & "+"0 \\\\")
        # print("& $K_{m,\\hat \\mu}$ & "+"0 \\\\")
        # print()
        # print("& $K_{y,\\beta}$ & "+"{:<8.5f} \\\\".format(Ky_B))
        # print("& $K_{\\ell,\\beta}$ & "+"{:<8.5f} \\\\".format(Kl_B))
        # print("& $K_{n,\\beta}$ & "+"{:<8.5f} \\\\".format(Kn_B))
        # print()
        # print("& $K_{y,\\breve p}$ & "+"{:<8.5f} \\\\".format(Ky_pbar))
        # print("& $K_{\\ell,\\breve p}$ & "+"{:<8.5f} \\\\".format(Kl_pbar))
        # print("& $K_{n,\\breve p}$ & "+"{:<8.5f} \\\\".format(Kn_pbar))
        # print()
        # print("& $K_{z,\\breve q}$ & "+"{:<8.5f} \\\\".format(Kz_qbar))
        # print("& $K_{x,\\breve q}$ & "+"{:<8.5f} \\\\".format(Kx_qbar))
        # print("& $K_{m,\\breve q}$ & "+"{:<8.5f} \\\\".format(Km_qbar))
        # print()
        # print("& $K_{y,\\breve r}$ & "+"{:<8.5f} \\\\".format(Ky_rbar))
        # print("& $K_{\\ell,\\breve r}$ & "+"{:<8.5f} \\\\".format(Kl_rbar))
        # print("& $K_{n,\\breve r}$ & "+"{:<8.5f} \\\\".format(Kn_rbar))
        # print()

        # calculate A,B,C and D matrices
        A = [[Kx_mu,  Kx_a, -Kx_qbar,  0, 0, -ct],
            [-Kz_mu, Kz_a, Kz_qbar,  0, 0, -st],
            [Km_mu, Km_a, Km_qbar, 0, 0, 0],
            [ct, st, 0, 0, 0, -st],
            [-st, ct, 0, 0, 0, -ct],
            [0, 0, 1, 0, 0, 0]]
        A = np.array(A)

        B = [[1, 0, 0,  0, 0, 0],
            [0, 1+Kz_ah, 0,  0, 0, 0],
            [0, -Km_ah, 1, 0, 0, 0],
            [0, 0, 0, 1, 0, 0],
            [0, 0, 0, 0, 1, 0],
            [0, 0, 0, 0, 0, 1]]
        B = np.array(B)

        C = [[Ky_B, Ky_pbar, Ky_rbar-1,  0, ct, 0],
            [Kl_B,   Kl_pbar, Kl_rbar,  0, 0, 0],
            [Kn_B, Kn_pbar, Kn_rbar, 0, 0, 0],
            [1, 0, 0, 0, 0, ct],
            [0, 1, np.tan(th), 0, 0, 0],
            [0, 0, 1./ct, 0, 0, 0]]
        C = np.array(C)

        D = [[1, 0, 0,  0, 0, 0],
            [0, 1, -txz,  0, 0, 0],
            [0, -tzx, 1, 0, 0, 0],
            [0, 0, 0, 1, 0, 0],
            [0, 0, 0, 0, 1, 0],
            [0, 0, 0, 0, 0, 1]]
        D = np.array(D)

        # make eigenvalue calculations
        lonval, lonvec = eig(A,B)
        latval, latvec = eig(C,D)

        # save globally
        self.altr_lon_val = lonval
        self.altr_lon_vec = lonvec
        self.altr_lat_val = latval
        self.altr_lat_vec = latvec


    def _report(self):
        """Method which reports the calculated values to the user.
        """

        # whether to report dimensional or not
        dimensional = False

        # dimensional ratio
        if dimensional:
            drg = 2. * self.V / self.cbar
            drt = 2. * self.V / self.b
            dagt = self.g / self.V # * self.cbar / 2. / self.V**2.
            dsr = "Dimensional "
        else:
            drg = 1.
            drt = 1.
            dagt = 1.
            dsr = ""

        # report longitudinal eigenvalues
        t = 22
        d = 12
        n = t + d

        # report approx values
        print("Approximate Values")
        print("sp {:>+14.12f}\u00B1{:>14.12f}".format(np.real(self.l_sp),\
        np.abs(np.imag(self.l_sp))))
        print("ph {:>+14.12f}\u00B1{:>14.12f}".format(np.real(self.l_ph),\
        np.abs(np.imag(self.l_ph))))
        print("r  {:>+14.12f}+{:>14.12f}".format(np.real(self.l_r),\
        np.abs(np.imag(self.l_r))))
        print("s  {:>+14.12f}+{:>14.12f}".format(np.real(self.l_s),\
        np.abs(np.imag(self.l_s))))
        print("dr {:>+14.12f}\u00B1{:>14.12f}".format(np.real(self.l_dr),\
        np.abs(np.imag(self.l_dr))))

        print("{:^{}s}".format(self.input.split(".")[0],2*n))
        print("{:^{}s}".format("Longitudinal",2*n))
        print("-"*2*n)
        print("{:^{}s} {:^{}s}".format(dsr + "Traditional",n,\
            dsr + "Alternate",n))
        print("-"*2*n)
        for i in range(self.trad_lon_val.shape[0]):
            print("{:>{}.{}f} {:>{}.{}f}".format(self.trad_lon_val[i]*drg,n,d,\
                self.altr_lon_val[i]*dagt,n,d))
        print()

        print("{:^{}s}".format("Lateral",2*n))
        print("-"*2*n)
        print("{:^{}s} {:^{}s}".format(dsr + "Traditional",n,\
            dsr + "Alternate",n))
        print("-"*2*n)
        for i in range(self.trad_lon_val.shape[0]):
            print("{:>{}.{}f} {:>{}.{}f}".format(self.trad_lat_val[i]*drt,n,d,\
                self.altr_lat_val[i]*dagt,n,d))



