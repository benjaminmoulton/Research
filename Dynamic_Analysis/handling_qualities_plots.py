import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import ScalarFormatter

from eigensolver import Solver

def analyze_aircraft(file_name):
    # run eigensolver
    sys = Solver(file_name,report=False)
    print("running {}".format(file_name))

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

    # calculate short-period damping ratio
    z_sp = -0.5*(sys.L_a/sys.Iyy/lmpm)**0.5*(sys.m_adot/sys.L_a + ayyb)
    CAP_sp = sys.g*lmpm/ryyb**2.
    
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
    z_sp_exact = sm2*0.
    CAP_sp_exact = sm2*0.
    craft_presented = ["Navion","F-94A","Lockheed Jetstar","Boeing 747"]
    for i in range(sm2.shape[0]):
        sys.Cm_a = -sm2[i]*sys.CL_a
        sys._buckingham_matrices_and_approximations()
        sys._longitudinal_dimensional_properties(sys.b,is_alternate=True)
        z_sp_exact[i] = sys.b["lon"]["zt"][sys.b["lon"]["sp"][0]]*1.
        CAP_sp_exact[i] = sys.b["lon"]["wn"][sys.b["lon"]["sp"][0]]**2.\
            /sys.CL_a*sys.CW
        z_ph_exact[i] = sys.b["lon"]["zt"][sys.b["lon"]["ph"][0]]*1.
        if z_sp_exact[i] == 0. and sys.aircraft_name in craft_presented:
            e1 = sys.b["lon"]["evals"][sys.b["lon"]["sp"][0]]
            e2 = sys.b["lon"]["evals"][sys.b["lon"]["sp"][1]]
            z_sp_exact[i] = np.real(-(e1+e2) /2./np.sqrt(e1*e2))
            e1 = sys.b["lon"]["dimevals"][sys.b["lon"]["sp"][0]]
            e2 = sys.b["lon"]["dimevals"][sys.b["lon"]["sp"][1]]
            wn = np.real(np.sqrt(e1*e2))
            CAP_sp_exact[i] = wn**2./sys.CL_a*sys.CW
            # print(sm2[i])
            # print(sys.b["lon"]["sp"])
            # print(sys.b["lon"]["ph"])
            # print(sys.b["lon"]["zt"][sys.b["lon"]["ph"][0]])
            # print(sys.b["lon"]["zt"][sys.b["lon"]["sp"][0]])
            # print(sys.b["lon"]["evals"])
            # print()
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
    aircraft_name = sys.aircraft_name.lower().replace(" ","_").replace("-","_")

    # initialize figures
    fsp, asp = plt.subplots()
    fsC, asC = plt.subplots()
    fph, aph = plt.subplots()

    # plot grid
    asp.grid(which="major",lw=0.6,ls="-",c="0.5")
    asp.grid(which="minor",lw=0.5,ls="dotted",c="0.5")
    aph.grid(which="major",lw=0.6,ls="-",c="0.5")
    aph.grid(which="minor",lw=0.5,ls="dotted",c="0.5")
    asC.grid(which="major",lw=0.6,ls="-",c="0.5")
    asC.grid(which="minor",lw=0.5,ls="dotted",c="0.5")
    
    # plot
    asp_y2 = asp.twinx()
    asp_y2.plot(sm2[1:],CAP_sp_exact[1:],c="0.25",marker="o",mfc="none",ms=4.,\
        ls="none")
    asp_y2.plot(sm,CAP_sp,c="0.25",label="Eq. (105)")
    asp.plot(sm2[1:],z_sp_exact[1:],c="k",marker="o",mfc="none",ms=4.,\
        ls="none",label="exact")
    asp.plot(sm,z_sp,c="k",label="Eq. (106)")
    # plot asC limits based on flight phase
    if sys.aircraft_name in ["Lockheed Jetstar","Boeing 747"]:
        # landing lines
        asC.loglog([5.,.15,.15],[.038,.038,20.],c="k")
        asC.loglog([.2,.2,2.,2.],[.038,10.,10.,.038],c="k")
        asC.loglog([2.,.3,.3,2.],[.085,.085,3.6,3.6],c="k")
        l1pos = (0.7,0.4)
        l4y = 0.03
    elif sys.aircraft_name in ["Navion","F-94A"]:
        # cruise lines
        asC.loglog([5.,.15,.15],[.096,.096,20.],c="k")
        asC.loglog([.25,.25,2.,2.],[.096,10.,10.,.096],c="k")
        asC.loglog([1.3,.35,.35,1.3,1.3],[.15,.15,3.6,3.6,.15],c="k")
        l1pos = (0.7,1.3)
        l4y = 0.04
    else:
        # plot cruise lines
        asC.loglog([5.,.15,.15],[.096,.096,20.],c="k")
        asC.loglog([.25,.25,2.,2.],[.096,10.,10.,.096],c="k")
        asC.loglog([1.3,.35,.35,1.3,1.3],[.15,.15,3.6,3.6,.15],c="k")
        l1pos = (0.7,1.3)
        l4y = 0.04
    asC.loglog(z_sp,CAP_sp,c="0.25")
    CAPlist = np.logspace(np.log10(CAP_sp[int(num/20)]),np.log10(CAP_sp[-1]),
        10,base=10.,endpoint=False)
    ind_list = [np.argwhere(CAP_sp <= CAPlist_it)[-1,0] for CAPlist_it in CAPlist]
    for i in ind_list:
        asC.annotate('',
        xytext=(z_sp[i],CAP_sp[i]),
        xy=(z_sp[i+1],CAP_sp[i+1]),
        arrowprops=dict(arrowstyle="->", color="0.25"),
        size=10.
        )
    # plot level labels
    alfa = 0.8
    text = asC.text(l1pos[0],l1pos[1],"Level 1",va="center",ha="center",
        bbox=dict(facecolor="w",linewidth=0,alpha=alfa,
        boxstyle="Square, pad=0.0"))
    text = asC.text(0.7,6.,"Level 2",va="center",ha="center",
        bbox=dict(facecolor="w",linewidth=0,alpha=alfa,
        boxstyle="Square, pad=0.0"))
    text = asC.text(0.7,13.,"Level 3",va="center",ha="center",
        bbox=dict(facecolor="w",linewidth=0,alpha=alfa,
        boxstyle="Square, pad=0.0"))
    text = asC.text(0.7,l4y,"Level 4",va="center",ha="center",
        bbox=dict(facecolor="w",linewidth=0,alpha=alfa,
        boxstyle="Square, pad=0.0"))
    # asp.plot(sm,z_sp_other,c="k",ls="--",label="Eq. (113)")
    #
    aph.plot(sm2[1:],z_ph_exact[1:],c="k",marker="o",mfc="none",ms=4.,\
        ls="none",label="exact")
    aph.plot(sm,z_ph,c="k",label="Eq. (111)")
    aph.plot(sm,z_ph_other,c="k",ls="--",label="Eq. (113)")
    # aph.plot([lnpm_0 - lmpm_0,lnpm_0 - lmpm_0],[0,1],c="k",ls="-.")
    # aph.plot(lnpm_0/sys.cwbar,sys.b["lon"]["zt"][sys.b["lon"]["ph"][0]]*1.,\
    #     marker="o",c="b",mfc="none",ms=4.,ls="none")
    if sys.aircraft_name == "Navion":
        ph_limit = 0.5
        sp_limit = None
        y2_limit = None
    elif sys.aircraft_name == "F-94A":
        ph_limit = 0.1
        sp_limit = 0.2
        y2_limit = 0.05
    elif sys.aircraft_name == "Lockheed Jetstar":
        ph_limit = 0.25
        sp_limit = 0.35
        y2_limit = 0.09
    elif sys.aircraft_name == "Boeing 747":
        ph_limit = 0.45
        sp_limit = None
        y2_limit = None
    else:
        ph_limit = 1.
        sp_limit = None
        y2_limit = None
    asp.set_xlim((0.,sm_final))
    asp.set_ylim(bottom=sp_limit)
    asp_y2.set_ylim(bottom=y2_limit)
    aph.set_xlim((0.,sm_final))
    aph.set_ylim((0.,ph_limit))
    asC.set_xlim((0.1,5.))
    asC.set_ylim((0.02,20.))
    asC.xaxis.set_major_formatter(ScalarFormatter())
    asC.yaxis.set_major_formatter(ScalarFormatter())
    # asC.semilogy()
    # asC.semilogx()
    asp.set_xlabel(r"Pitch static margin, $l_{np_m}/\bar{c}_w$")
    asp.set_ylabel(r"Short-period damping ratio, $\zeta_{sp}$")
    asp_y2.set_ylabel(r"Short-period CAP")
    aph.set_xlabel(r"Pitch static margin, $l_{np_m}/\bar{c}_w$")
    aph.set_ylabel(r"Phugoid damping ratio, $\zeta_{ph}$")
    asC.set_xlabel(r"Short-period damping ratio, $\zeta_{sp}$")
    asC.set_ylabel(r"Short-period CAP")
    aph.legend()
    lgnd_elms_sp = [
        Line2D([0], [0], c="k",marker="o",mfc="none",ms=4.,\
        ls="none",label="exact"),
        Line2D([0], [0], c='k', label='Eq. (106)'),
        Line2D([0], [0], c='0.25', label='Eq. (105)')
    ]
    asp.legend(handles=lgnd_elms_sp,loc="upper center")
    folder = "hand_qual_plots/"
    savedict = dict(transparent=True,format="png",dpi=300.0)
    name = "z_sp_" + aircraft_name
    fsp.savefig(folder + "short_period/" + name + ".png",**savedict)
    name = "z_ph_" + aircraft_name
    fph.savefig(folder + "phugoid/" + name + ".png",**savedict)
    name = "CAP_v_z_sp_" + aircraft_name
    fsC.savefig(folder + "short_period_CAP_v_z/" + name + ".png",**savedict)
    if False:
        plt.show()
    else:
        plt.close('all')
    

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

    # run for various roll static margins, shift by 0.05 lnpn/bw
    lsm_final = 0.20
    step = 0.05
    num = int(lsm_final/step) + 1
    lsm = np.linspace(0.,lsm_final,num=num)
    lsm[0] = 0.0001
    hnpl_array = lsm*sys.bw

    # initialize figures
    fsl, asl = plt.subplots()
    fdr, adr = plt.subplots()
    fdC, adC = plt.subplots()

    # plot grid
    asl.grid(which="major",lw=0.6,ls="-",c="0.5")
    asl.grid(which="minor",lw=0.5,ls="dotted",c="0.5")
    adr.grid(which="major",lw=0.6,ls="-",c="0.5")
    adr.grid(which="minor",lw=0.5,ls="dotted",c="0.5")
    adC.grid(which="major",lw=0.6,ls="-",c="0.5")
    adC.grid(which="minor",lw=0.5,ls="dotted",c="0.5")

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

        # calculate dutch roll CAP
        CAP_dr = sys.g*lmpn/rzzb**2.


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
        CAP_dr_exact = nsm2*0.
        for j in range(2,nsm2.shape[0]):
            sys.Cn_b = -nsm2[j]*sys.CY_b
            sys.Cl_b =  lsm[i]*sys.CY_b
            sys._buckingham_matrices_and_approximations()
            sys._lateral_dimensional_properties(sys.b,is_alternate=True)
            t_sl_exact[j] = sys.b["lat"]["double"][sys.b["lat"]["sl"]]*1.
            z_dr_exact[j] = sys.b["lat"]["zt"][sys.b["lat"]["dr"][0]]*1.
            CAP_dr_exact[j] = -sys.b["lat"]["wn"][sys.b["lat"]["dr"][0]]**2./\
                sys.CY_b*sys.CW
        # return to normal
        sys.Cl_b = Cl_b_0*1.
        sys.Cn_b = Cn_b_0*1.
        sys._buckingham_matrices_and_approximations()
        sys._lateral_dimensional_properties(sys.b,is_alternate=True)

        # plot spiral
        asl.plot(nsm[:ipo],t_sl_other[:ipo],ls="--",c=str(i/num*0.95))
        asl.plot(nsm[ips:],t_sl[ips:],c=str(i/num*0.95))
        asl.plot(nsm2[2:],t_sl_exact[2:],c=str(i/num*0.95),marker="o",\
            mfc="none",ms=4.,ls="none")

        # plot Dutch-roll damping ratio
        adr.plot(nsm,z_dr_other,ls="--",c=str(i/num*0.95))
        adr.plot(nsm,z_dr,c=str(i/num*0.95))
        adr.plot(nsm2[2:],z_dr_exact[2:],c=str(i/num*0.95),marker="o",\
            mfc="none",ms=4.,ls="none")
        
        # plot Dutch-roll CAP
        if i == 0:
            adC.plot(nsm,CAP_dr,c=str(i/num*0.95))
            adC.plot(nsm2[2:],CAP_dr_exact[2:],c=str(i/num*0.95),marker="o",\
                mfc="none",ms=4.,ls="none")
    
    # legend elements
    lbl_name = r"$h_{np_\ell}/b_w =$ "
    lgnd_elms_sl = [
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
    lgnd_elms_CP = [
        Line2D([0], [0], c='k', ls='none',lw=1,marker="o",\
            mfc="none",ms=4., label='exact'),
        Line2D([0], [0], c='k', ls='-',lw=1, label='Eq. (126)')
    ]

    t_2m = 240.
    asl.set_xlim((0.,sm_final))
    asl.set_ylim((0.,t_2m))
    asl.set_xlabel(r"Yaw static margin, $l_{np_n}/b_w$")
    asl.set_ylabel(r"Spiral time to double, $\tau_{sl}$")
    asl.legend(handles=lgnd_elms_sl)
    name = "t_sl_" + aircraft_name
    fsl.savefig(folder + "spiral/" + name + ".png",**savedict)

    if sys.aircraft_name == "Navion":
        limit = 0.75
    elif sys.aircraft_name == "F-94A":
        limit = 0.4
    elif sys.aircraft_name == "Lockheed Jetstar":
        limit = 0.6
    elif sys.aircraft_name == "Boeing 747":
        limit = 0.6
    else:
        limit = 1.
    adr.set_xlim((0.,sm_final))
    adr.set_ylim((0.,limit))
    adr.set_xlabel(r"Yaw static margin, $l_{np_n}/b_w$")
    adr.set_ylabel(r"Dutch-roll damping ratio, $\zeta_{dr}$")
    adr.legend(handles=lgnd_elms_dr)
    name = "z_dr_" + aircraft_name
    fdr.savefig(folder + "dutch_roll/" + name + ".png",**savedict)

    adC.set_xlim((0.,sm_final))
    # adC.set_ylim((0.,limit))
    adC.set_xlabel(r"Yaw static margin, $l_{np_n}/b_w$")
    adC.set_ylabel(r"Dutch-roll CAP")
    adC.legend(handles=lgnd_elms_CP)
    name = "CAP_dr_" + aircraft_name
    fdC.savefig(folder + "dutch_roll_CAP/" + name + ".png",**savedict)
    if False:
        plt.show()
    else:
        plt.close('all')





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
    "F_16.json",
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