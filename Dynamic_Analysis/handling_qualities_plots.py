import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

from eigensolver import Solver

def analyze_aircraft(file_name):
    # run eigensolver
    sys = Solver(file_name,report=False)

    ## Phugoid
    lnpm_0 = - sys.m_a/sys.L_a
    lmpm_0 = lnpm_0 - sys.m_q*sys.g/sys.vo/sys.W
    RG = sys.Lo/sys.Do
    ryyb = ( sys.g*sys.Iyy/sys.W )**0.5
    ayyb = sys.m_q/sys.L_a - ryyb**2./sys.vo

    # shift by 0.005 lnpm/cwbar
    sm_final = 0.30
    step = 0.0025
    num = int(sm_final/step) + 1
    sm = np.linspace(0.,sm_final,num=num)
    sm[0] = 0.0001
    sm2 = np.linspace(0.,sm_final,num=int(num/2))
    sm2[0] = 0.0001
    lnpm = sm*sys.cwbar
    lmpm = lnpm - sys.m_q*sys.g/sys.vo/sys.W
    
    # calculate phugoid damping ratio
    z_ph = 1./RG*(lmpm/2./lnpm)**0.5 \
        + sys.g/sys.vo*(lnpm/2./lmpm**3.)**0.5*ayyb

    # for i in range(5):
    #     lmpm_n = lmpm + (i+1)*0.1*ryyb
    #     ayyb_n = (lnpm - lmpm_n)/(sys.g/sys.vo/sys.W)/sys.L_a - ryyb**2./sys.vo

    #     z_ph_n = 1./RG*(lmpm_n/2./lnpm)**0.5 \
    #         + sys.g/sys.vo*(lnpm/2./lmpm_n**3.)**0.5*ayyb_n
        
    #     plt.plot(sm,z_ph_n,c=str((i+1)*0.1))
    
    # calculate exact
    Cm_a_0 = sys.Cm_a*1.
    z_ph_exact = sm2*0.
    for i in range(sm2.shape[0]):
        sys.Cm_a = -sm2[i]*sys.CL_a
        sys._buckingham_matrices_and_approximations()
        sys._longitudinal_dimensional_properties(sys.b,is_alternate=True)
        z_ph_exact[i] = sys.b["lon"]["zt"][sys.b["lon"]["ph"][0]]*1.
    # return to normal
    sys.Cm_a = Cm_a_0*1.
    sys._buckingham_matrices_and_approximations()
    sys._longitudinal_dimensional_properties(sys.b,is_alternate=True)

    # calculate bad approx
    z_ph_other = z_ph*0. + 1./2.**0.5/RG
    
    subdict = {
        "figsize" : (3.25,2.4375),
        "constrained_layout" : True,
        "sharex" : True
    }

    # plot grid
    plt.grid(which="major",lw=0.6,ls="-",c="0.5")
    plt.grid(which="minor",lw=0.5,ls="dotted",c="0.5")
    
    # plot
    plt.plot(sm2[1:],z_ph_exact[1:],c="k",marker="o",mfc="none",ms=4.,\
        ls="none",label="exact")
    plt.plot(sm,z_ph,c="k",label="Eq. (111)")
    plt.plot(sm,z_ph_other,c="k",ls="--",label="Eq. (113)")
    # plt.plot([lnpm_0 - lmpm_0,lnpm_0 - lmpm_0],[0,1],c="k",ls="-.")
    # plt.plot(lnpm_0/sys.cwbar,sys.b["lon"]["zt"][sys.b["lon"]["ph"][0]]*1.,\
    #     marker="o",c="b",mfc="none",ms=4.,ls="none")
    plt.xlabel(r"Pitching static margin, $l_{np_m}/\bar{c}_w$")
    plt.ylabel(r"Phugoid damping ratio, $\zeta_{ph}$")
    plt.legend()#bbox_to_anchor=(1.05, 1.0), loc='upper right')
    folder = "hand_qual_plots/"
    savedict = dict(transparent=True,format="png",dpi=300.0)
    name = "z_ph_"
    name += sys.aircraft_name.lower().replace(" ","_").replace("-","_")
    plt.savefig(folder + "phugoid/" + name + ".png",**savedict)
    if False:
        plt.show()
    else:
        plt.close()
    

    ## Spiral & Dutch-roll
    hnpl_0 =   sys.l_b/sys.Y_b
    lnpn_0 = - sys.n_b/sys.Y_b
    # print(sys.aircraft_name,hnpl_0/sys.bw,lnpn_0/sys.bw)
    # shift by 0.0025 lnpn/bw
    nsm_final = 0.30
    step = 0.0025
    num = int(nsm_final/step) + 1
    nsm = np.linspace(0.,nsm_final,num=num)
    nsm[0] = 0.001
    nsm2 = np.linspace(0.,nsm_final,num=int(num/3))
    nsm2[0] = 0.001
    lnpn = nsm*sys.bw

    # run for various rolling static margins, shift by 0.05 lnpn/bw
    lsm_final = 0.20
    step = 0.05
    num = int(lsm_final/step) + 1
    lsm = np.linspace(0.,lsm_final,num=num)
    lsm[0] = 0.0001
    hnpl_array = lsm*sys.bw

    # initialize figures
    fsl, asl = plt.subplots()
    fdr, adr = plt.subplots()

    for i,hnpl in enumerate(hnpl_array):
        # calculate spiral time to double
        Ep = hnpl*sys.n_p + lnpn*sys.l_p
        Er = hnpl*sys.n_r + lnpn*sys.l_r
        t_sl = - np.log(2.)*sys.vo/sys.g*Ep/Er
        t_sl_other = - np.log(2.)*sys.Izz*hnpl/Er

        # calculate dutch roll damping ratio
        rzzb = ( sys.g*sys.Izz/sys.W )**0.5
        azzb = sys.n_r/sys.Y_b + rzzb**2./sys.vo
        hmpl = hnpl + sys.l_r*sys.g/sys.vo/sys.W
        lmpn = lnpn - sys.n_r*sys.g/sys.vo/sys.W

        dr_a = azzb + sys.Ixx*(sys.n_p/sys.l_p**2.*hmpl)
        dr_b = - sys.l_r*sys.n_p/sys.Y_b/sys.l_p
        dr_c = - sys.g/sys.vo/sys.l_p*Er/lmpn
        s_dr = -0.5*( sys.Y_b/sys.Izz*(dr_a + dr_b) + dr_c )
        dr_d = hmpl*sys.n_p/sys.l_p - sys.W*rzzb**2.*hnpl/sys.vo/sys.l_p
        dr_e = - sys.Y_b/sys.Izz*(lmpn + dr_d + 0.25*sys.Y_b/sys.Izz*azzb**2.)
        wn_dr = ( s_dr**2. + dr_e )**0.5
        z_dr = s_dr / wn_dr
        z_dr_other = 0.5*(-sys.Y_b/sys.Izz/lmpn)**0.5*azzb

        # phillips
        Cn_b = -nsm*sys.CY_b
        Cl_b = lsm[i]*sys.CY_b
        Rpy = sys.rho * sys.Sw * sys.bw / 4.0 / (sys.W / sys.g)
        Rgy = sys.g * sys.bw / 2. / sys.vo**2.
        Rxx = sys.rho * sys.Sw * sys.bw**3. / 8. / sys.Ixx
        Rzz = sys.rho * sys.Sw * sys.bw**3. / 8. / sys.Izz
        Ry_b = Rpy * sys.CY_b
        Rl_b = Rxx * Cl_b
        Rn_b = Rzz * Cn_b
        Rl_pbar = Rxx * sys.Cl_pbar
        Rn_pbar = Rzz * sys.Cn_pbar
        Ry_rbar = Rpy * sys.CY_rbar
        Rl_rbar = Rxx * sys.Cl_rbar
        Rn_rbar = Rzz * sys.Cn_rbar
        RDRs = (Rl_b*(Rgy - (1.-Ry_rbar)*Rn_pbar) - Ry_b*Rl_rbar*Rn_pbar) / \
            Rl_pbar
        RDRc = Rl_rbar * Rn_pbar / Rl_pbar
        RDRp = Rgy*(Rl_rbar*Rn_b-Rl_b*Rn_rbar)/Rl_pbar/(Rn_b+Ry_b*Rn_rbar) - \
            RDRs/Rl_pbar
        drSg = - sys.vo / sys.bw * (Ry_b + Rn_rbar - RDRc + RDRp)
        a = (1.-Ry_rbar)*Rn_b 
        b = Ry_b * Rn_rbar 
        c = RDRs 
        d = - 0.25*( Ry_b + Rn_rbar)**2.
        e = abs(a + b + c + d)
        drWD = np.abs(2. * sys.vo / sys.bw * np.sqrt( e ))
        drZ = drSg/drWD
        z_dr = drZ


        # get positive indices for spiral
        try:
            ips = np.argwhere(t_sl >= 0.0)[0,0]
        except:
            ips = None
        try:
            ipo = np.argwhere(t_sl_other >= 0.0)[-1,0]
        except:
            ipo = None
        
        # calculate exact
        Cn_b_0 = sys.Cn_b*1.
        Cl_b_0 = sys.Cl_b*1.
        t_sl_exact = nsm2*0.
        z_dr_exact = nsm2*0.
        for j in range(2,nsm2.shape[0]):
            sys.Cn_b = -nsm2[j]*sys.CY_b
            sys.Cl_b =  lsm[i]*sys.CY_b
            sys._buckingham_matrices_and_approximations()
            sys._lateral_dimensional_properties(sys.b,is_alternate=True)
            t_sl_exact[j] = sys.b["lat"]["double"][sys.b["lat"]["sl"]]*1.
            try:
                z_dr_exact[j] = sys.b["lat"]["zt"][sys.b["lat"]["dr"][0]]*1.
            except:
                z_dr_exact[j] = sys.b["lat"]["zt"][sys.b["lat"]["sl"][0]]*1.
        # return to normal
        sys.Cl_b = Cl_b_0*1.
        sys.Cn_b = Cn_b_0*1.
        sys._buckingham_matrices_and_approximations()
        sys._lateral_dimensional_properties(sys.b,is_alternate=True)

        # plot spiral
        asl.plot(nsm2[2:],t_sl_exact[2:],c=str(i/num*0.95),marker="o",\
            mfc="none",ms=4.,ls="none")
        asl.plot(nsm[ips:],t_sl[ips:],c=str(i/num*0.95))
        asl.plot(nsm[:ipo],t_sl_other[:ipo],ls="--",c=str(i/num*0.95))

        # plot spiral
        adr.plot(nsm2[2:],z_dr_exact[2:],c=str(i/num*0.95),marker="o",\
            mfc="none",ms=4.,ls="none")
        adr.plot(nsm,z_dr,c=str(i/num*0.95))
        adr.plot(nsm,z_dr_other,ls="--",c=str(i/num*0.95))
    
    # legend elements
    lbl_name = r"$h_{np_\ell}/b_w =$ "
    lgnd_elms_sp = [
        Line2D([0], [0], c='k', ls='none',lw=1,marker="o",\
            mfc="none",ms=4., label='exact'),
        Line2D([0], [0], c='k', ls='-',lw=1, label='Eq. (115)'),
        Line2D([0], [0], c='k', ls='--',lw=1, label='Eq. (118)'),
        Line2D([0], [0], c=str(0/num*0.95), ls='-',lw=1, label=lbl_name + "0.0"),
        Line2D([0], [0], c=str(2/num*0.95), ls='-',lw=1, label=lbl_name + "0.1"),
        Line2D([0], [0], c=str(4/num*0.95), ls='-',lw=1, label=lbl_name + "0.2")
    ]
    lgnd_elms_dr = [
        Line2D([0], [0], c='k', ls='none',lw=1,marker="o",\
            mfc="none",ms=4., label='exact'),
        Line2D([0], [0], c='k', ls='-',lw=1, label='Eq. (122)'),
        Line2D([0], [0], c='k', ls='--',lw=1, label='Eq. (125)'),
        Line2D([0], [0], c=str(0/num*0.95), ls='-',lw=1, label=lbl_name + "0.0"),
        Line2D([0], [0], c=str(2/num*0.95), ls='-',lw=1, label=lbl_name + "0.1"),
        Line2D([0], [0], c=str(4/num*0.95), ls='-',lw=1, label=lbl_name + "0.2")
    ]

    t_2m = 240.
    asl.set_ylim((0.,t_2m))
    asl.set_xlabel(r"Yawing static margin, $l_{np_n}/b_w$")
    asl.set_ylabel(r"Spiral time to double, $\tau_{sl}$")
    asl.legend(handles=lgnd_elms_sp)
    name = "t_sl_"
    name += sys.aircraft_name.lower().replace(" ","_").replace("-","_")
    fsl.savefig(folder + "spiral/" + name + ".png",**savedict)
    if False:
        plt.show()
    else:
        plt.close()

    limit = 1.
    adr.set_ylim((-limit,limit))
    adr.set_xlabel(r"Yawing static margin, $l_{np_n}/b_w$")
    adr.set_ylabel(r"Dutch-roll damping ratio, $\zeta_{dr}$")
    adr.legend(handles=lgnd_elms_dr)
    name = "z_dr_"
    name += sys.aircraft_name.lower().replace(" ","_").replace("-","_")
    fdr.savefig(folder + "dutch_roll/" + name + ".png",**savedict)
    if False:
        plt.show()
    else:
        plt.close()





if __name__ == "__main__":


    # change plot text parameters
    plt.rcParams["font.family"] = "Serif"
    plt.rcParams["font.size"] = 12.0
    plt.rcParams["axes.labelsize"] = 12.0
    plt.rcParams['lines.linewidth'] = 1.0
    plt.rcParams["xtick.minor.visible"] = True
    plt.rcParams["ytick.minor.visible"] = True
    plt.rcParams["xtick.direction"] = plt.rcParams["ytick.direction"] = "in"
    plt.rcParams["xtick.bottom"] = plt.rcParams["xtick.top"] = True
    plt.rcParams["ytick.left"] = plt.rcParams["ytick.right"] = True
    plt.rcParams["xtick.major.width"] = plt.rcParams["ytick.major.width"] = 0.75
    plt.rcParams["xtick.minor.width"] = plt.rcParams["ytick.minor.width"] = 0.75
    plt.rcParams["xtick.major.size"] = plt.rcParams["ytick.major.size"] = 5.0
    plt.rcParams["xtick.minor.size"] = plt.rcParams["ytick.minor.size"] = 2.5
    # plt.rcParams["axes.labelpad"] = 2.0
    # plt.rcParams["font.weight"] = "bold"
    # plt.rcParams["figure.constrained_layout.hspace"] = 0.0
    # plt.rcParams["figure.subplot.hspace"] = 0.0
    plt.rcParams["mathtext.fontset"] = "dejavuserif"
    plt.rcParams['figure.dpi'] = 200.0
    plt.rcParams['figure.constrained_layout.use'] = True

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
    num_craft = len(run_files)
    for i in range(num_craft):
        analyze_aircraft(run_files[i])

    # folder = "aircraft_database/"
    # filename = "Cessna_172.json"
    # analyze_aircraft(file_name=folder+filename)