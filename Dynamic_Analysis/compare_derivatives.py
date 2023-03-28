import numpy as np
import json
import matplotlib.pyplot as plt
from scipy.linalg import eig

class Aircraft:
    def __init__(self):
        return


class Comparison:
    """A class which compares several aircraft 
    derivatives from multiple input files.

    Parameters
    ----------
    input_vars : list of strings or dicts, optional
        Each item must be a .json file or python dictionary.

    Raises
    ------
    TypeError
        If the input_vars type is not a dictionary or the file path to 
        a .json file
    """
    def __init__(self,input_vars=[{}],report=True):

        # report
        self.report = report
        if self.report:
            print("Running Comparison written by Ben Moulton\n")

        # get info or raise error
        self.num_craft = len(input_vars)
        self.input = input_vars
        self._get_input_vars(input_vars)

        # retrieve info
        self._retrieve_info()

        # calculate dimensional derivatives
        self._get_derivatives()

        # calculate ratios
        self._get_ratios()

        # report ratios
        if self.report:
            self._report_ratios()


    def _get_input_vars(self,input_vars):
        # get info or raise error
        # run through items in list
        self.Craft = [Aircraft() for item in input_vars]
        for index, item in enumerate(input_vars):

            # determine if the input_vars is a file or a dictionary
            input_vars_type = type(item)

            # dictionary
            if input_vars_type == dict:
                self.Craft[index].input_dict = item
            
            # json file
            elif input_vars_type == str and item.split(".")[-1] =="json":
                # import json file from file path
                json_string = open(item).read()

                # save to vals dictionary
                self.Craft[index].input_dict = json.loads(json_string)

            # raise error
            else:
                raise TypeError("input_vars item {0} ".format(index) + \
                    "must be json file path, or " + \
                    "dictionary, not {0}".format(input_vars_type))


    def _retrieve_info(self):
        """Method which retrieves the information and stores it globally.
        """
        
        # store variables from file input dictionary for each aircraft
        for i in range(self.num_craft):
            zelf = self.Craft[i]
        
            # store aircraft input values
            aircraft = zelf.input_dict.get("aircraft",{})
            zelf.aircraft_name = aircraft.get("name","Jimmy")
            zelf.name = aircraft.get("short_name","Jimmy")
            zelf.EL = aircraft.get("launch_kinetic_energy[ft-lbf]",0.0)
            zelf.Sw = aircraft.get("wing_area[ft^2]")
            zelf.bw = aircraft.get("wing_span[ft]")
            zelf.cwbar = zelf.Sw / zelf.bw
            zelf.W = aircraft.get("weight[lbf]")
            zelf.Ixx = aircraft.get("Ixx[slug-ft^2]")
            zelf.Iyy = aircraft.get("Iyy[slug-ft^2]")
            zelf.Izz = aircraft.get("Izz[slug-ft^2]")
            zelf.Ixy = aircraft.get("Ixy[slug-ft^2]")
            zelf.Ixz = aircraft.get("Ixz[slug-ft^2]")
            zelf.Iyz = aircraft.get("Iyz[slug-ft^2]")
            
            # store analysis input values
            analysis = zelf.input_dict.get("analysis",{})
            zelf.rho = analysis.get("density[slugs/ft^3]")
            
            # store aerodynamics input values
            aero = zelf.input_dict.get("aerodynamics",{})
            CL = aero.get("CL")
            zelf.CL0        = CL.get("0")
            zelf.CL_a       = CL.get("alpha")
            zelf.CL_qbar    = CL.get("qbar")
            zelf.CL_ahat    = CL.get("alpha_hat")
            zelf.CL_uhat    = CL.get("mu_hat",0.0)
            CY = aero.get("CS")
            zelf.CY_b       = CY.get("beta")
            zelf.CY_pbar    = CY.get("Lpbar") * zelf.CL0 + \
                CY.get("pbar")
            zelf.CY_rbar    = CY.get("rbar")
            zelf.CY_da      = CY.get("da",0.0)
            CD = aero.get("CD")
            zelf.CD0        = CD.get("L0")
            zelf.CD1        = CD.get("L")
            zelf.CD2        = CD.get("L2")
            zelf.CD_qbar    = CD.get("L2qbar")*zelf.CL0**2.+\
                CD.get("Lqbar") * zelf.CL0 + CD.get("qbar")
            zelf.CD_ahat    = CD.get("alpha_hat",0.0)
            zelf.CD_uhat    = CD.get("mu_hat",0.0)
            Cl = aero.get("Cl")
            zelf.Cl_b       = Cl.get("beta")
            zelf.Cl_pbar    = Cl.get("pbar")
            zelf.Cl_rbar    = Cl.get("Lrbar") * zelf.CL0 + \
                Cl.get("rbar")
            Cm = aero.get("Cm")
            zelf.Cm0        = Cm.get("0")
            zelf.Cm_a       = Cm.get("alpha")
            zelf.Cm_qbar    = Cm.get("qbar")
            zelf.Cm_ahat    = Cm.get("alpha_hat")
            zelf.Cm_uhat    = Cm.get("mu_hat",0.0)
            Cn = aero.get("Cn")
            zelf.Cn_b       = Cn.get("beta")
            zelf.Cn_pbar    = Cn.get("Lpbar") * zelf.CL0 + \
                Cn.get("pbar")
            zelf.Cn_rbar    = Cn.get("rbar")

            # set gravity
            zelf.g = 32.17#4

            # set angles
            zelf.thetao = 0.0
            zelf.phio = 0.0
            zelf.ct = np.cos(zelf.thetao)
            zelf.st = np.sin(zelf.thetao)
            zelf.cp = np.cos(zelf.phio)

            # calculate velocity
            zelf.vo = (zelf.W / 0.5 / zelf.rho / \
                zelf.CL0 / zelf.Sw)**0.5

            # calculate lift and drag
            zelf.CW = zelf.CL0
            zelf.CLo = zelf.CL0 * zelf.ct / \
                zelf.cp
            zelf.CDo = zelf.CD0 + zelf.CD1 * \
                zelf.CLo + zelf.CD2 * zelf.CLo**2.
            zelf.CD_a = zelf.CD1 * zelf.CL_a + \
                2. * zelf.CD2 * zelf.CLo * zelf.CL_a
            zelf.Cmo = zelf.Cm0
            # if self.report:
            #     print("Running solver for {}".format(\
            #         zelf.aircraft_name))
            #     print("Vo  = ", zelf.vo)
            #     print("CDo  = ", zelf.CDo)
            #     print("CD,a = ", zelf.CD_a)
            #     print()

            "this is here"


    def _get_derivatives(self):
        # redimensionalize derivatives
        for i in range(self.num_craft):
            zelf = self.Craft[i]

            # ratios of import
            pSwcw_4 = zelf.rho * zelf.Sw * zelf.cwbar / 4.
            pSwbw_4 = zelf.rho * zelf.Sw * zelf.bw    / 4.
            dynF = 0.5 * zelf.rho * zelf.vo**2. * zelf.Sw

            ### Longitudinal derivatives
            # D,udot
            zelf.D_udot = pSwcw_4 * zelf.CD_uhat
            # L,udot
            zelf.L_udot = pSwcw_4 * zelf.CL_uhat
            # m,udot
            zelf.m_udot = pSwcw_4 * zelf.cwbar * zelf.Cm_uhat
            # D,adot
            zelf.D_adot = pSwcw_4 * zelf.vo * zelf.CD_ahat
            # L,adot
            zelf.L_adot = pSwcw_4 * zelf.vo * zelf.CL_ahat
            # m,adot
            zelf.m_adot = pSwcw_4 * zelf.vo * zelf.cwbar * zelf.Cm_ahat
            # Do
            zelf.Do = dynF * zelf.CDo
            # Lo
            zelf.Lo = dynF * zelf.CLo
            # m
            zelf.mo = dynF * zelf.cwbar * zelf.Cmo
            # D,a
            zelf.D_a = dynF * zelf.CD_a
            # L,a
            zelf.L_a = dynF * zelf.CL_a
            # m,a
            zelf.m_a = dynF * zelf.cwbar * zelf.Cm_a
            # D,q
            zelf.D_q = pSwcw_4 * zelf.vo * zelf.CD_qbar
            # L,q
            zelf.L_q = pSwcw_4 * zelf.vo * zelf.CL_qbar
            # m,q
            zelf.m_q = pSwcw_4 * zelf.vo * zelf.cwbar * zelf.Cm_qbar

            ### Lateral Derivatives
            # Y,b
            zelf.Y_b = dynF * zelf.CY_b
            # l,b
            zelf.l_b = dynF * zelf.bw * zelf.Cl_b
            # n,b
            zelf.n_b = dynF * zelf.bw * zelf.Cn_b
            # Y,p
            zelf.Y_p = pSwbw_4 * zelf.vo * zelf.CY_pbar
            # l,p
            zelf.l_p = pSwbw_4 * zelf.vo * zelf.bw * zelf.Cl_pbar
            # n,p
            zelf.n_p = pSwbw_4 * zelf.vo * zelf.bw * zelf.Cn_pbar
            # Y,r
            zelf.Y_r = pSwbw_4 * zelf.vo * zelf.CY_rbar
            # l,r
            zelf.l_r = pSwbw_4 * zelf.vo * zelf.bw * zelf.Cl_rbar
            # n,r
            zelf.n_r = pSwbw_4 * zelf.vo * zelf.bw * zelf.Cn_rbar

            # also save new nondimensionalizations
            Rpx = zelf.rho * zelf.Sw * zelf.cwbar * zelf.g / 4. / zelf.W
            Rpy = zelf.rho * zelf.Sw * zelf.bw    * zelf.g / 4. / zelf.W
            CWinv = 1. / zelf.CW
            Ryy2 = zelf.rho * zelf.vo**4.*zelf.Sw*zelf.cwbar/2./zelf.g**2./zelf.Iyy
            Rxx2 = zelf.rho * zelf.vo**4.*zelf.Sw*zelf.bw   /2./zelf.g**2./zelf.Ixx
            Rzz2 = zelf.rho * zelf.vo**4.*zelf.Sw*zelf.bw   /2./zelf.g**2./zelf.Izz
            Ryy = zelf.rho*zelf.vo**2.*zelf.Sw*zelf.cwbar**2./4./zelf.g/zelf.Iyy
            Rxx = zelf.rho*zelf.vo**2.*zelf.Sw*zelf.bw**2./4./zelf.g/zelf.Ixx
            Rzz = zelf.rho*zelf.vo**2.*zelf.Sw*zelf.bw**2./4./zelf.g/zelf.Izz

            # LONGITUDINAL
            # values
            zelf.ixz = zelf.Ixz / zelf.Ixx
            zelf.izx = zelf.Ixz / zelf.Izz
            zelf.Kx_muhat = Rpx * zelf.CD_uhat
            zelf.Kz_muhat = Rpx * zelf.CL_uhat
            zelf.Km_muhat = Ryy * zelf.Cm_uhat
            zelf.Kx_ahat = Rpx * zelf.CD_ahat
            zelf.Kz_ahat = Rpx * zelf.CL_ahat
            zelf.Km_ahat = Ryy * zelf.Cm_ahat
            zelf.Kx_mu = CWinv * -2. * zelf.CDo
            zelf.Kz_mu = CWinv *  2. * zelf.CLo
            zelf.Km_mu = Ryy2 * 2. * zelf.Cm0
            zelf.Kx_a = CWinv * ( zelf.CLo -      zelf.CD_a)
            zelf.Kz_a = CWinv * (-zelf.CDo - zelf.CL_a)
            zelf.Km_a = Ryy2 * zelf.Cm_a
            zelf.Ky_b = CWinv * zelf.CY_b
            zelf.Kl_b = Rxx2 * zelf.Cl_b
            zelf.Kn_b = Rzz2 * zelf.Cn_b
            zelf.Ky_pbreve = Rpy * zelf.CY_pbar
            zelf.Kl_pbreve = Rxx * zelf.Cl_pbar
            zelf.Kn_pbreve = Rzz * zelf.Cn_pbar
            zelf.Kx_qbreve = Rpx * zelf.CD_qbar
            zelf.Kz_qbreve = Rpx * zelf.CL_qbar
            zelf.Km_qbreve = Ryy * zelf.Cm_qbar
            zelf.Ky_rbreve = Rpy * zelf.CY_rbar
            zelf.Kl_rbreve = Rxx * zelf.Cl_rbar
            zelf.Kn_rbreve = Rzz * zelf.Cn_rbar

            # other important terms
            zelf.l_mp = -zelf.m_a/zelf.L_a - zelf.m_q/zelf.W * zelf.g/zelf.vo
            zelf.l_np = -zelf.m_a/zelf.L_a


    def _get_ratios(self):
        # create ratios names list
        self.ratios_names = [
            ["Kl,b","-Kl,b*Kn,p","Kl,b*Ky,r*Kn,p","-KY,b*Kl,r*Kn,p"],
            ["Kn,b","Ky,r*Kn,b","Ky,b*Kn,r"],
            ["Kl,r*Kn,b","-Kl,b*Kn,r"],
            ["Kn,b","Ky,b*Kn,r"],
            ["-Kl,b*(1-Kn,p)/Kl,p","(Kl,r*Kn,b-Kl,b*Kn,r)/(Kn,b+Ky,b*Kn,r)","-Kl,r*Kn,p"],
            ["Ky,b","Kn,r"],
            ["a*a","2*a*b","2*a*c","b*b","2*b*c","c*c",],
            ["A","B","C","D","E","F",],
            ["n,b","l,b*n,p/l,p"],
            ["1/W/Izz*Y,b*n,r","l,b/l,p"],
            ["1","2","3","4","5"],
            ["Km,a","Kz,a*Km,q"],
            ["Km,a*Kz,a","Km,a*Km,q"],
            ["A","B","C"],
            ["-Kz,mu","Frac"],
            ["a","b","c"],
            ["Do","L,a"],
            ["W*m,a","g/vo*L,a*m,q"],
            ["Kz,a","Km,q","Km,adot"],
            ["l,b*n,p","-l,p*n,b","l,b*n,r","-l,r*n,b"],
            ["CW","CL,qbreve"],
            ["m,q","L,a"],
            ["-ryyb^2/Vo","m,q/L,a"],
            ["-ryyb^2/Vo*L,a","m,q","m,adot"],
            ["Aa","Ab","Ac","B","C","D","E","F"],
            ["1","-rxx^2*Y,b/2/l,p","-Ixx/Izz*n,r/2/l,p"],
            ["Ixx","Izz"],
            ["n,p","n,r"],
            ["-l,p","-l,r"],
            ["l,b*n,p","l,b*n,r","-l,p*n,b","-l,r*n,b"],
            ["-Kz,mu","Frac"],
            ["2*Lo/W*lnp/lmp","-Do**2/W**2"],
            ["1","Ky,r"],
            ["CS,pbar","CS,rbar","CS,da"],
            ["CW","CS,rbreve"],
            ["Kn,pbreve","Kl,pbreve"],
            ["hnp_l/bw","lnp_m/cwbar","lnp_n/bw"],
            ["hnp_l","lnp_m","lnp_n"],
        ]
        self.num_ratios = len(self.ratios_names)
        self.ratlen = [len(self.ratios_names[i]) for i in range(self.num_ratios)]
        self.maxratlen = max(self.ratlen)
        self.ratios_terms = ["Dutch Roll wn"] * self.num_ratios
        for i in range(10,11): self.ratios_terms[i] = "Dutch Roll Sn"
        for i in range(11,14): self.ratios_terms[i] = "Phugoid Sn"
        for i in range(14,18): self.ratios_terms[i] = "Phugoid wn"
        for i in range(18,19): self.ratios_terms[i] = "Short Period Sn"
        for i in range(19,20): self.ratios_terms[i] = "Spiral Sn"
        for i in range(20,21): self.ratios_terms[i] = "Phillips Term Drop"
        for i in range(21,23): self.ratios_terms[i] = "Phugoid Sn"
        for i in range(23,24): self.ratios_terms[i] = "Short Period Sn"
        for i in range(30,32): self.ratios_terms[i] = "Phugoid wn"

        ratios = np.zeros((self.num_ratios,self.num_craft,self.maxratlen))
        sortcomps = np.zeros((self.num_ratios,self.maxratlen),dtype=int)

        for i in range(self.num_ratios):
            for j in range(self.num_craft):
                zelf = self.Craft[j]
                g = -1

                # new ratio
                if i == 0:
                    g+=1;ratios[i,j,g] = zelf.Kl_b
                    g+=1;ratios[i,j,g] = -zelf.Kl_b * zelf.Kn_pbreve
                    g+=1;ratios[i,j,g] = zelf.Kl_b*zelf.Ky_rbreve*zelf.Kn_pbreve
                    g+=1;ratios[i,j,g] = -zelf.Ky_b*zelf.Kl_rbreve*zelf.Kn_pbreve
                elif i == 1:
                    g+=1;ratios[i,j,g] = zelf.Kn_b
                    g+=1;ratios[i,j,g] = zelf.Ky_rbreve * zelf.Kn_b
                    g+=1;ratios[i,j,g] = zelf.Ky_b * zelf.Kn_rbreve
                elif i == 2:
                    g+=1;ratios[i,j,g] = zelf.Kl_rbreve * zelf.Kn_b
                    g+=1;ratios[i,j,g] = -zelf.Kl_b * zelf.Kn_rbreve
                elif i == 3:
                    g+=1;ratios[i,j,g] = zelf.Kn_b
                    g+=1;ratios[i,j,g] = zelf.Ky_b * zelf.Kn_rbreve
                elif i == 4:
                    a_num = zelf.Kl_rbreve*zelf.Kn_b - zelf.Kl_b*zelf.Kn_rbreve
                    a_den = zelf.Kn_b + zelf.Ky_b * zelf.Kn_rbreve
                    a = a_num / a_den
                    b = -zelf.Kl_b * (1. - zelf.Kn_pbreve) / zelf.Kl_pbreve
                    c = -zelf.Kl_rbreve * zelf.Kn_pbreve
                    g+=1;ratios[i,j,g] = a
                    g+=1;ratios[i,j,g] = b
                    g+=1;ratios[i,j,g] = c
                elif i == 5:
                    g+=1;ratios[i,j,g] = zelf.Ky_b
                    g+=1;ratios[i,j,g] = zelf.Kn_rbreve
                elif i == 6:
                    a_num = zelf.Kl_rbreve*zelf.Kn_b - zelf.Kl_b*zelf.Kn_rbreve
                    a_den = zelf.Kn_b + zelf.Ky_b * zelf.Kn_rbreve
                    a = a_num / a_den
                    b = -zelf.Kl_b * (1. - zelf.Kn_pbreve) / zelf.Kl_pbreve
                    c = -zelf.Kl_rbreve * zelf.Kn_pbreve
                    g+=1;ratios[i,j,g] = a*a
                    g+=1;ratios[i,j,g] = 2*a*b
                    g+=1;ratios[i,j,g] = 2*a*c
                    # g+=1;ratios[i,j,g] = b*a
                    g+=1;ratios[i,j,g] = b*b
                    g+=1;ratios[i,j,g] = 2*b*c
                    # g+=1;ratios[i,j,g] = c*a
                    # g+=1;ratios[i,j,g] = c*b
                    g+=1;ratios[i,j,g] = c*c
                elif i == 7:
                    a_num = zelf.Kl_rbreve*zelf.Kn_b - zelf.Kl_b*zelf.Kn_rbreve
                    a_den = zelf.Kn_b + zelf.Ky_b * zelf.Kn_rbreve
                    a = a_num / a_den
                    b = -zelf.Kl_b * (1. - zelf.Kn_pbreve) / zelf.Kl_pbreve
                    c = -zelf.Kl_rbreve * zelf.Kn_pbreve
                    A = 0.5 * (zelf.Ky_b+zelf.Kn_rbreve)/zelf.Kl_pbreve*(a+b+c)
                    B = 0.25/zelf.Kl_pbreve**2. * (a+b+c)**2.
                    C = zelf.Kn_b
                    D = zelf.Ky_b * zelf.Kn_rbreve
                    E =  zelf.Kl_b / zelf.Kl_pbreve
                    F = -zelf.Kl_b * zelf.Kn_pbreve / zelf.Kl_pbreve
                    g+=1;ratios[i,j,g] = A
                    g+=1;ratios[i,j,g] = B
                    g+=1;ratios[i,j,g] = C
                    g+=1;ratios[i,j,g] = D
                    g+=1;ratios[i,j,g] = E
                    g+=1;ratios[i,j,g] = F
                elif i == 8:
                    g+=1;ratios[i,j,g] =  zelf.n_b 
                    g+=1;ratios[i,j,g] = -zelf.l_b*zelf.n_p/zelf.l_p 
                elif i == 9:
                    g+=1;ratios[i,j,g] =  1./zelf.W/zelf.Izz*zelf.Y_b*zelf.n_r
                    g+=1;ratios[i,j,g] =  zelf.l_b/zelf.l_p
                elif i == 10:
                    a_num = zelf.Kl_rbreve*zelf.Kn_b - zelf.Kl_b*zelf.Kn_rbreve
                    a_den = zelf.Kn_b + zelf.Ky_b * zelf.Kn_rbreve
                    b = -zelf.Kl_b * (1. - zelf.Kn_pbreve) / zelf.Kl_pbreve
                    c = -zelf.Kl_rbreve * zelf.Kn_pbreve
                    V1 = zelf.Ky_b
                    V2 = zelf.Kn_rbreve
                    V3 = - zelf.Kl_rbreve * zelf.Kn_pbreve / zelf.Kl_pbreve
                    V4 = a_num / a_den / zelf.Kl_pbreve
                    V5 = b / zelf.Kl_pbreve                
                    g+=1;ratios[i,j,g] = V1
                    g+=1;ratios[i,j,g] = V2
                    g+=1;ratios[i,j,g] = V3
                    g+=1;ratios[i,j,g] = V4
                    g+=1;ratios[i,j,g] = V5
                elif i == 11:
                    g+=1;ratios[i,j,g] = zelf.Km_a
                    g+=1;ratios[i,j,g] = -zelf.Kz_a * zelf.Km_qbreve
                elif i == 12:
                    g+=1;ratios[i,j,g] = zelf.Km_a * zelf.Kz_a
                    g+=1;ratios[i,j,g] = zelf.Km_a * zelf.Km_qbreve
                elif i == 13:
                    A = - zelf.Kx_mu / zelf.Kz_mu
                    den = zelf.Km_a - zelf.Kz_a * zelf.Km_qbreve
                    B = zelf.Kx_a * zelf.Km_qbreve / den
                    C = - zelf.Km_a * (zelf.Kz_a + zelf.Km_qbreve)/den**2.
                    g+=1;ratios[i,j,g] = A
                    g+=1;ratios[i,j,g] = B
                    g+=1;ratios[i,j,g] = C
                elif i == 14:
                    A = - zelf.Kz_mu
                    den = zelf.Km_a - zelf.Kz_a * zelf.Km_qbreve
                    B = zelf.Kx_a * zelf.Km_qbreve / den
                    g+=1;ratios[i,j,g] = A
                    g+=1;ratios[i,j,g] = B
                elif i == 15:
                    A = - zelf.Kz_mu
                    den = zelf.Km_a - zelf.Kz_a * zelf.Km_qbreve
                    B = zelf.Kx_a * zelf.Km_qbreve / den
                    C = - zelf.Km_a * (zelf.Kz_a + zelf.Km_qbreve)/den**2.
                    a = ( -zelf.Kx_mu/zelf.Kz_mu + B + C )**2.
                    b = 4./zelf.Kz_mu * zelf.Km_a / den
                    c = -( zelf.Kx_mu / zelf.Kz_mu )**2.
                    g+=1;ratios[i,j,g] = a
                    g+=1;ratios[i,j,g] = b
                    g+=1;ratios[i,j,g] = c
                elif i == 16:
                    g+=1;ratios[i,j,g] = zelf.Do
                    g+=1;ratios[i,j,g] = zelf.L_a
                elif i == 17:
                    g+=1;ratios[i,j,g] = zelf.W*zelf.m_a
                    g+=1;ratios[i,j,g] = zelf.g/zelf.vo*zelf.L_a*zelf.m_q
                elif i == 18:
                    g+=1;ratios[i,j,g] = zelf.Kz_a
                    g+=1;ratios[i,j,g] = zelf.Km_qbreve
                    g+=1;ratios[i,j,g] = zelf.Km_ahat
                elif i == 19:
                    g+=1;ratios[i,j,g] =  zelf.l_b * zelf.n_p
                    g+=1;ratios[i,j,g] = -zelf.l_p * zelf.n_b
                    g+=1;ratios[i,j,g] =  zelf.l_b * zelf.n_r
                    g+=1;ratios[i,j,g] = -zelf.l_r * zelf.n_b
                elif i == 20:
                    g+=1;ratios[i,j,g] =  zelf.CW
                    g+=1;ratios[i,j,g] =  zelf.CL_qbar * zelf.g*zelf.cwbar/2./zelf.vo**2.
                elif i == 21:
                    g+=1;ratios[i,j,g] = zelf.m_q
                    g+=1;ratios[i,j,g] = zelf.L_a
                elif i == 22:
                    g+=1;ratios[i,j,g] = -(zelf.g*zelf.Iyy/zelf.W)/zelf.vo
                    g+=1;ratios[i,j,g] = zelf.m_q/zelf.L_a
                elif i == 23:
                    g+=1;ratios[i,j,g] = -(zelf.g*zelf.Iyy/zelf.W)/zelf.vo*zelf.L_a
                    g+=1;ratios[i,j,g] = zelf.m_q
                    g+=1;ratios[i,j,g] = zelf.m_adot
                elif i == 24:
                    a_num = zelf.Kl_rbreve*zelf.Kn_b - zelf.Kl_b*zelf.Kn_rbreve
                    a_den = zelf.Kn_b + zelf.Ky_b * zelf.Kn_rbreve
                    a = a_num / a_den
                    b = -zelf.Kl_b * (1. - zelf.Kn_pbreve) / zelf.Kl_pbreve
                    c = -zelf.Kl_rbreve * zelf.Kn_pbreve
                    Aa = 0.5 * (zelf.Ky_b+zelf.Kn_rbreve)/zelf.Kl_pbreve*(a)
                    Ab = 0.5 * (zelf.Ky_b+zelf.Kn_rbreve)/zelf.Kl_pbreve*(b)
                    Ac = 0.5 * (zelf.Ky_b+zelf.Kn_rbreve)/zelf.Kl_pbreve*(c)
                    B = 0.25/zelf.Kl_pbreve**2. * (a+b+c)**2.
                    C = zelf.Kn_b
                    D = zelf.Ky_b * zelf.Kn_rbreve
                    E =  zelf.Kl_b / zelf.Kl_pbreve
                    F = -zelf.Kl_b * zelf.Kn_pbreve / zelf.Kl_pbreve
                    g+=1;ratios[i,j,g] = Aa
                    g+=1;ratios[i,j,g] = Ab
                    g+=1;ratios[i,j,g] = Ac
                    g+=1;ratios[i,j,g] = B
                    g+=1;ratios[i,j,g] = C
                    g+=1;ratios[i,j,g] = D
                    g+=1;ratios[i,j,g] = E
                    g+=1;ratios[i,j,g] = F
                elif i == 25:
                    rxx2 = zelf.g*zelf.Ixx/zelf.W
                    g+=1;ratios[i,j,g] = 1.0
                    g+=1;ratios[i,j,g] = -rxx2*zelf.Y_b/2./zelf.l_p
                    g+=1;ratios[i,j,g] = -zelf.Ixx/zelf.Izz*zelf.n_r/2./zelf.l_p
                elif i == 26:
                    g+=1;ratios[i,j,g] = zelf.Ixx
                    g+=1;ratios[i,j,g] = zelf.Izz
                elif i == 27:
                    g+=1;ratios[i,j,g] = zelf.n_p
                    g+=1;ratios[i,j,g] = zelf.n_r
                elif i == 28:
                    g+=1;ratios[i,j,g] = -zelf.l_p
                    g+=1;ratios[i,j,g] = -zelf.l_r
                elif i == 29:
                    g+=1;ratios[i,j,g] = zelf.n_p  * zelf.l_b
                    g+=1;ratios[i,j,g] = zelf.n_r  * zelf.l_b
                    g+=1;ratios[i,j,g] = -zelf.l_p * zelf.n_b
                    g+=1;ratios[i,j,g] = -zelf.l_r * zelf.n_b
                elif i == 30:
                    A = - zelf.Kz_mu
                    den = zelf.Km_a - zelf.Kz_a * zelf.Km_qbreve
                    B = zelf.Kx_a * zelf.Km_qbreve / den
                    g+=1;ratios[i,j,g] = A
                    g+=1;ratios[i,j,g] = B
                elif i == 31:
                    g+=1;ratios[i,j,g] = 2.*zelf.Lo/zelf.W*zelf.l_np/zelf.l_mp
                    g+=1;ratios[i,j,g] = -zelf.Do**2./zelf.W**2.
                elif i == 32:
                    g+=1;ratios[i,j,g] = 1.0
                    g+=1;ratios[i,j,g] = zelf.Ky_rbreve
                elif i == 33:
                    g+=1;ratios[i,j,g] = zelf.CY_pbar
                    g+=1;ratios[i,j,g] = zelf.CY_rbar
                    g+=1;ratios[i,j,g] = zelf.CY_da
                elif i == 34:
                    g+=1;ratios[i,j,g] = zelf.CW
                    g+=1;ratios[i,j,g] = zelf.CY_rbar * zelf.g*zelf.bw/2./zelf.vo**2.
                elif i == 35:
                    g+=1;ratios[i,j,g] = zelf.Kn_pbreve
                    g+=1;ratios[i,j,g] = zelf.Kl_pbreve
                elif i == 36:
                    g+=1;ratios[i,j,g] =   zelf.l_b / zelf.bw    / zelf.Y_b
                    g+=1;ratios[i,j,g] = - zelf.m_a / zelf.cwbar / zelf.L_a
                    g+=1;ratios[i,j,g] = - zelf.n_b / zelf.bw    / zelf.Y_b
                elif i == 37:
                    g+=1;ratios[i,j,g] =   zelf.l_b / zelf.Y_b
                    g+=1;ratios[i,j,g] = - zelf.m_a / zelf.L_a
                    g+=1;ratios[i,j,g] = - zelf.n_b / zelf.Y_b
                
        
        # take absolute value
        # ratios = np.abs(ratios)
        powers = ratios * 0.0
        avgpwr = np.zeros((self.num_ratios,self.maxratlen))
        # get power
        for i in range(self.num_ratios):
            for j in range(self.ratlen[i]):
                for k in range(self.num_craft):
                    if ratios[i,k,j] != 0.0:
                        powers[i,k,j] = np.floor(np.log10(np.abs(ratios[i,k,j])))
                
                # determine average power
                avgpwr[i,j] = np.average(powers[i,:,j])
                
            # reorder by largest magnitude
            avgpwr[i,:self.ratlen[i]] -= np.max(avgpwr[i,:self.ratlen[i]])

            sortcomps[i,:self.ratlen[i]] = \
                np.flip(np.argsort(avgpwr[i,:self.ratlen[i]]))



        self.avgpwr = avgpwr
        self.ratios = ratios
        self.skip_up_to = 0 # 24 # 
        self.sortcomps = sortcomps


    def _report_ratios(self):
        # report
        print("Reporting...")

        # run through ratios
        for i in range(self.skip_up_to,self.num_ratios):
            print(self.ratios_terms[i])
            for h in range(len(self.ratios_names[i])):
                hsort = self.sortcomps[i,h]
                print(str(h+1)+".",self.ratios_names[i][hsort],end="")
                if h != len(self.ratios_names[i]) - 1:
                    print(end=",  ")
            print()
            for j in range(-1,self.ratlen[i]):
                    for k in range(self.num_craft):
                        if j==-1:
                            if k == 0:
                                print("{:^3s}\t".format("num"),end="")
                            print("{:^6s}".format(self.Craft[k].name),end="")
                            if k != self.num_craft - 1:
                                print(" \t",end="")
                            else:
                                print(" \t{:^6s}".format("avgpwr"),end="")
                        else:
                            jsort = self.sortcomps[i,j]
                            if k == 0:
                                print("{:^3s}\t".format(str(j+1)),end="")
                            print("{:<+6.0e}\t".format(self.ratios[i][k][jsort]),end="")
                            if k == self.num_craft -1:
                                print("{:<+5.3f}".format(self.avgpwr[i,jsort]),end="")
                    print()

            # organize results
            def_cuts = np.argwhere(self.avgpwr[i]<-2.0)[:,0].tolist()
            all_cuts = np.argwhere(self.avgpwr[i]<-1.0)[:,0].tolist()
            may_cuts = [i for i in all_cuts if i not in def_cuts]
            print("Results:")
            if def_cuts:
                print("\tCut:   ",end="")
                for d in range(len(def_cuts)):
                    print(self.ratios_names[i][def_cuts[d]],end="")
                    if d != len(def_cuts)-1: print(", ",end="")
                    else: print()
            if may_cuts:
                print("\tMaybe: ",end="")
                for d in range(len(may_cuts)):
                    print(self.ratios_names[i][may_cuts[d]],end="")
                    if d != len(may_cuts)-1: print(", ",end="")
                    else: print()
            if not def_cuts and not may_cuts:
                print("\tno change")
            print()

if __name__ == "__main__":

    # compare general aviation, RC glider, and fighter aircraft
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
    "F16_bolander.json",
    "NT_33A.json", "F_104A.json", "F_4C.json",
    "A_7A.json", "A_4D.json",
    "F_94A.json", "F_15.json",
    ]
    run_files = ["aircraft_database/" + i for i in run_files]
    Comparison(run_files)

