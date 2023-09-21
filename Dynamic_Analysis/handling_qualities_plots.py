import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import ScalarFormatter
from matplotlib.legend_handler import HandlerTuple

from eigensolver import Solver

def analyze_aircraft(file_name):
    # run eigensolver
    print("running {}".format(file_name))
    sys = Solver(file_name,report=False)

    ## Phugoid
    lnpm_0 = - sys.m_a/sys.L_a
    lmpm_0 = lnpm_0 - sys.m_q*sys.g/sys.vo/sys.W
    RG = sys.Lo/sys.Do
    ryyb = ( sys.g*sys.Iyy/sys.W )**0.5
    ayyb = sys.m_q/sys.L_a - ryyb**2./sys.vo

    # shift by 0.005 lnpm/cwbar
    dm_final = 0.30
    step = 0.0025
    num = int(dm_final/step) + 1
    dm = np.linspace(0.,dm_final,num=num)
    dm[0] = 0.0001
    dm2 = np.linspace(0.,dm_final,num=int((num-1)/2)+1)
    dm2[0] = 0.0001
    lmpm = dm*ryyb
    lnpm = lmpm + sys.m_q*sys.g/sys.vo/sys.W

    # calculate short-period damping ratio
    z_sp = -0.5*np.sqrt(sys.L_a/sys.Iyy/lmpm)*(sys.m_adot/sys.L_a + ayyb)
    CAP_sp = sys.g*lmpm/ryyb**2.
    
    # calculate phugoid damping ratio
    z_ph = 1./RG*np.sqrt(np.abs(lmpm/2./lnpm)) \
        + sys.g/sys.vo*np.sqrt(np.abs(lnpm/2./lmpm**3.))*ayyb
    np.set_printoptions(linewidth=np.inf)
    
    # calculate exact
    Cm_a_0 = sys.Cm_a*1.
    z_ph_exact = dm2*0.
    z_sp_exact = dm2*0.
    CAP_sp_exact = dm2*0.
    craft_presented = ["Navion","F-94A","Lockheed Jetstar","Boeing 747"]
    # evals = np.zeros((6,dm2.shape[0]),dtype=complex)
    plot_sp = True
    plot_ph = True
    for i in reversed(range(dm2.shape[0])):
        dx = -(lmpm_0 - dm2[i]*ryyb)
        # print(lmpm_0/ryyb,dm2[i],dx)
        sys._buckingham_matrices_and_approximations(cg_shift=[dx,0.,0.])
        sys._longitudinal_dimensional_properties(sys.b,is_alternate=True)
        if plot_sp and np.imag(sys.b["lon"]["evals"][sys.b["lon"]["sp"][0]]):
            z_sp_exact[i] = sys.b["lon"]["zt"][sys.b["lon"]["sp"][0]]*1.
            CAP_sp_exact[i] = sys.b["lon"]["wn"][sys.b["lon"]["sp"][0]]**2.\
                /sys.CL_a*sys.CW
        else:
            z_sp_exact[i] = None
            CAP_sp_exact[i] = None
            plot_sp = False
        if plot_ph and np.imag(sys.b["lon"]["evals"][sys.b["lon"]["ph"][0]]):
            z_ph_exact[i] = sys.b["lon"]["zt"][sys.b["lon"]["ph"][0]]*1.
        else:
            z_ph_exact[i] = None
            plot_ph = False
        # ###   ###   ###   ###   ###
        # evals[0:2,i] = sys.b["lon"]["evals"][sys.b["lon"]["rb"]]
        # evals[2:4,i] = sys.b["lon"]["evals"][sys.b["lon"]["sp"]]
        # evals[4:6,i] = sys.b["lon"]["evals"][sys.b["lon"]["ph"]]
        # ###   ###   ###   ###   ###
    # cs = ["k","k","b","r","g","m"]
    # flg, alg = plt.subplots()
    # for i in range(6):
    #     if cs[i] not in []:#"m","k","g"]:
    #         for k in range(len(evals[i])):
    #             alg.plot(np.real(evals[i,k]),np.imag(evals[i,k]),".",c=cs[i],\
    #                 ms=float(k/dm2.shape[0])*5. + 1.0)#,c="k")
    #         alg.plot(np.real(evals[i,:]),np.imag(evals[i,:]),linewidth=0.2,c=cs[i])
    #         alg.set_title(sys.aircraft_name)
    # folder = "hand_qual_plots/"
    # aircraft_name = sys.aircraft_name.lower().replace(" ","_").replace("-","_")
    # name = "eigs_" + aircraft_name
    # savedict = dict(transparent=True,format="png",dpi=300.0)
    # flg.savefig(folder + "lon_eigvals/" + name + ".png",**savedict)
            
    # return to normal
    sys._buckingham_matrices_and_approximations()
    sys._longitudinal_dimensional_properties(sys.b,is_alternate=True)

    # calculate bad approx
    z_ph_other = dm*0. + 1./2.**0.5/RG
    
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
    asp_y2.plot(dm2[1:],CAP_sp_exact[1:],c="0.25",marker="x",mfc="none",ms=4.,\
        ls="none")
    asp_y2.plot(dm,CAP_sp,c="0.25",label="Eq. (78)")
    # asp_y2.plot(lmpm_0/ryyb,sys.b["lon"]["wn"][sys.b["lon"]["sp"][0]]**2.\
    #     /sys.CL_a*sys.CW,\
    #     c="b",marker="x",mfc="none",ms=4.)
    asp.plot(dm,z_sp,c="k",label="Eq. (79)")
    asp.plot(dm2[1:],z_sp_exact[1:],c="k",marker="o",mfc="none",ms=4.,\
        ls="none",label="exact")
    # asp.plot(lmpm_0/ryyb,sys.b["lon"]["zt"][sys.b["lon"]["sp"][0]]*1.,\
    #     c="b",marker="o",mfc="none",ms=4.)
    # plot asC limits based on flight phase
    if sys.aircraft_name in ["F-94A"]:
        # maneuvering lines -- A flight phases
        asC.loglog([5.,.15,.15],[.15,.15,20.],c="k")
        asC.loglog([.25,.25,2.,2.],[.15,10.,10.,.15],c="k")
        asC.loglog([1.3,.35,.35,1.3,1.3],[.28,.28,3.6,3.6,.28],c="k")
        l1pos = (0.7,1.3)
        l4y = 0.04
    elif sys.aircraft_name in ["Navion"]:
        # cruise lines -- B flight phases
        asC.loglog([5.,.15,.15],[.038,.038,20.],c="k")
        asC.loglog([.2,.2,2.,2.],[.038,10.,10.,.038],c="k")
        asC.loglog([2.,.3,.3,2.,2.],[.085,.085,3.6,3.6,0.085],c="k")
        l1pos = (0.7,0.4)
        l4y = 0.03
    elif sys.aircraft_name in ["Lockheed Jetstar","Boeing 747"]:
        # landing lines -- C flight phases
        asC.loglog([5.,.15,.15],[.096,.096,20.],c="k")
        asC.loglog([.25,.25,2.,2.],[.096,10.,10.,.096],c="k")
        asC.loglog([1.3,.35,.35,1.3,1.3],[.15,.15,3.6,3.6,.15],c="k")
        l1pos = (0.7,1.3)
        l4y = 0.04
    else:
        # plot cruise lines
        asC.loglog([5.,.15,.15],[.038,.038,20.],c="k")
        asC.loglog([.2,.2,2.,2.],[.038,10.,10.,.038],c="k")
        asC.loglog([2.,.3,.3,2.,2.],[.085,.085,3.6,3.6,0.085],c="k")
        l1pos = (0.7,0.4)
        l4y = 0.03
    asC.loglog(z_sp,CAP_sp,c="0.25")
    # CAPlist = np.logspace(np.log10(CAP_sp[int(num/20)]),np.log10(CAP_sp[-1]),
    #     10,base=10.,endpoint=False)
    # ind_list = [np.argwhere(CAP_sp <= CAPlist_it)[-1,0] for CAPlist_it in CAPlist]
    for i in range(0,len(CAP_sp)-1,int((num-1)/6)): # ind_list: #
        asC.annotate('',
        xytext=(z_sp[i],CAP_sp[i]),
        xy=(z_sp[i+1],CAP_sp[i+1]),
        arrowprops=dict(arrowstyle="->", color="0.25"),
        size=10.
        )
    # asC.plot(sys.b["lon"]["zt"][sys.b["lon"]["sp"][0]]*1.,\
    #     sys.b["lon"]["wn"][sys.b["lon"]["sp"][0]]**2.\
    #     /sys.CL_a*sys.CW,\
    #     c="b",marker="x",mfc="none",ms=4.)
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
    # asp.plot(dm,z_sp_other,c="k",ls="--",label="Eq. (86)")
    #
    aph.plot(dm[1:],z_ph[1:],c="k",label="Eq. (84)")
    aph.plot(dm,z_ph_other,c="k",ls="--",label="Eq. (86)")
    aph.plot(dm2[1:],z_ph_exact[1:],c="k",marker="o",mfc="none",ms=4.,\
        ls="none",label="exact")
    # aph.plot(lmpm_0/ryyb,sys.b["lon"]["zt"][sys.b["lon"]["ph"][0]]*1.,\
    #     c="b",marker="o",mfc="none",ms=4.)
    # aph.plot([lnpm_0 - lmpm_0,lnpm_0 - lmpm_0],[0,1],c="k",ls="-.")
    # aph.plot(lnpm_0/sys.cwbar,sys.b["lon"]["zt"][sys.b["lon"]["ph"][0]]*1.,\
    #     marker="o",c="b",mfc="none",ms=4.,ls="none")
    # set plot limits
    if sys.aircraft_name == "Navion":
        ph_limit = 0.5
        sp_limit = 3.0
        y2_limit = 1.75
    elif sys.aircraft_name == "F-94A":
        ph_limit = 0.1
        sp_limit = 2.0
        y2_limit = 1.2
    elif sys.aircraft_name == "Lockheed Jetstar":
        ph_limit = 0.25
        sp_limit = 2.0
        y2_limit = 1.0
    elif sys.aircraft_name == "Boeing 747":
        ph_limit = 0.45
        sp_limit = 2.5
        y2_limit = 0.25
    else:
        ph_limit = 1.
        sp_limit = None
        y2_limit = None
    # plot HQ limits
    #                for phugoid
    # aph.plot([0.0,dm_final],[0.00,0.00],c="0.25")
    # aph.plot([0.0,dm_final],[0.04,0.04],c="0.25")
    # aph.plot(2*[(lmpm_0-lnpm_0)/ryyb],[0.0,ph_limit],c="0.2",ls="-.")
    aph.fill_between([0.0,dm_final],2*[0.04],2*[0.00],color="0.5",alpha=0.4)
    aph.text(dm_final/2.,0.0,"Level 2",va="bottom",ha="center",
        bbox=dict(facecolor="w",linewidth=0,alpha=alfa,
        boxstyle="Square, pad=0.0"))
    aph.text(dm_final/2.,(0.04+ph_limit)/2.,"Level 1",va="bottom",ha="center",
        bbox=dict(facecolor="w",linewidth=0,alpha=0.6,
        boxstyle="Square, pad=0.0"))
    # implement plot limits
    asp.set_xlim((0.,dm_final))
    asp.set_ylim(top=sp_limit)
    asp_y2.set_ylim(top=y2_limit)
    aph.set_xlim((0.,dm_final))
    aph.set_ylim((0.,ph_limit))
    asC.set_xlim((0.1,5.))
    asC.set_ylim((0.02,20.))
    asC.xaxis.set_major_formatter(ScalarFormatter())
    asC.yaxis.set_major_formatter(ScalarFormatter())
    # asC.semilogy()
    # asC.semilogx()
    asp.set_xlabel(r"Pitch dynamic margin, $l_{mp_m}/r_{yy_b}$")
    asp.set_ylabel(r"Short-period damping ratio, $\zeta_{sp}$")
    asp_y2.set_ylabel(r"Short-period CAP [s$^{-2}$]")
    aph.set_xlabel(r"Pitch dynamic margin, $l_{mp_m}/r_{yy_b}$")
    aph.set_ylabel(r"Phugoid damping ratio, $\zeta_{ph}$")
    asC.set_xlabel(r"Short-period damping ratio, $\zeta_{sp}$")
    asC.set_ylabel(r"Short-period CAP [s$^{-2}$]")
    aph.legend()
    han0,lab0 = asp.get_legend_handles_labels()
    han1,lab1 = asp_y2.get_legend_handles_labels()
    han0.insert(1,han1[0])
    lab0.insert(1,lab1[0])
    han0[-1] = (
        Line2D([0], [0], c="k",marker="o",mfc="none",ms=4.,\
            ls="none"),
        Line2D([0], [0], c="0.25",marker="x",mfc="none",ms=4.,\
            ls="none")
    )
    # print(han0)
    # quit()
    # lgnd_elms_sp = [
    #     Line2D([0], [0], c='k', label='Eq. (79)'),
    #     Line2D([0], [0], c='0.25', label='Eq. (78)'),
    #     Line2D([0], [0], c="k",marker="o",mfc="none",ms=4.,\
    #     ls="none",label="exact")
    # ]
    asp.legend(handles=han0,labels=lab0,loc="upper center",\
        handler_map={tuple: HandlerTuple(ndivide=None)})
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
    # quit()
    

    ## Spiral & Dutch-roll
    hnpl_0 =   sys.l_b/sys.Y_b
    hmpl_0 = hnpl_0 + sys.l_r*sys.g/sys.vo/sys.W
    rxxb = ( sys.g*sys.Ixx/sys.W )**0.5
    lnpn_0 = - sys.n_b/sys.Y_b
    lmpn_0 = lnpn_0 - sys.n_r*sys.g/sys.vo/sys.W
    lmpn_sm_0 = - sys.n_r*sys.g/sys.vo/sys.W
    rzzb = ( sys.g*sys.Izz/sys.W )**0.5
    azzb = sys.n_r/sys.Y_b + rzzb**2./sys.vo
    # print("axis  {:^12s} {:^12s} {:^12s} {:^12s} {:^12s} {:^12s}".format(\
    #     "np","np/cref","mp","mp/rref","ratecdot","Inertia"))
    # print("roll  {:> 12.8f} {:> 12.8f} {:> 12.8f} {:> 12.8f} {:> 12.8f} {:> 12.3f}".format(\
    #     hnpl_0,hnpl_0/sys.bw,hmpl_0,hmpl_0/rxxb,hmpl_0*sys.g/rxxb**2.,sys.Ixx))
    # print("pitch {:> 12.8f} {:> 12.8f} {:> 12.8f} {:> 12.8f} {:> 12.8f} {:> 12.3f}".format(\
    #     lnpm_0,lnpm_0/sys.cwbar,lmpm_0,lmpm_0/ryyb,lmpm_0*sys.g/ryyb**2.,sys.Iyy))
    # print("yaw   {:> 12.8f} {:> 12.8f} {:> 12.8f} {:> 12.8f} {:> 12.8f} {:> 12.3f}".format(\
    #     lnpn_0,lnpn_0/sys.bw,lmpn_0,lmpn_0/rzzb,lmpn_0*sys.g/rzzb**2.,sys.Izz))
    # print(sys.aircraft_name,hnpl_0/sys.bw,lnpn_0/sys.bw)
    ndm_final = 1.5
    step = 0.0025 # 0.00625 # 0.0125 # 
    num = int(ndm_final/step) + 1
    num_n = num + 0
    ndm = np.linspace(0.,ndm_final,num=num)
    ndm[0] = 0.001
    ndm2 = np.linspace(0.,ndm_final,num=int(num/15))
    ndm2[0] = 0.001
    lmpn = ndm*rzzb
    lnpn = lmpn + sys.n_r*sys.g/sys.vo/sys.W

    # run for various roll static margins, shift by 0.05 lnpn/bw
    ldm_final = 1.50
    step = 0.375
    num = int(ldm_final/step) + 1
    ldm = np.linspace(0.,ldm_final,num=num)
    ldm[0] = 0.0001
    hmpl = ldm*rxxb
    hnpl = hmpl - sys.l_r*sys.g/sys.vo/sys.W

    # initialize figures
    fsl, asl = plt.subplots()
    fdz, adz = plt.subplots()
    fds, ads = plt.subplots()
    fdw, adw = plt.subplots()
    fdC, adC = plt.subplots()

    # plot grid
    asl.grid(which="major",lw=0.6,ls="-",c="0.5")
    asl.grid(which="minor",lw=0.5,ls="dotted",c="0.5")
    adz.grid(which="major",lw=0.6,ls="-",c="0.5")
    adz.grid(which="minor",lw=0.5,ls="dotted",c="0.5")
    ads.grid(which="major",lw=0.6,ls="-",c="0.5")
    ads.grid(which="minor",lw=0.5,ls="dotted",c="0.5")
    adw.grid(which="major",lw=0.6,ls="-",c="0.5")
    adw.grid(which="minor",lw=0.5,ls="dotted",c="0.5")
    adC.grid(which="major",lw=0.6,ls="-",c="0.5")
    adC.grid(which="minor",lw=0.5,ls="dotted",c="0.5")
    adC_z_lim = (0.01,5.)
    adC_C_lim = (0.02,20.)

    # evals = np.zeros((6,ndm2.shape[0]),dtype=complex)
    # flg, alg = plt.subplots()

    for i,hmpli in enumerate(hmpl):
        hnpli = hmpli - sys.l_r*sys.g/sys.vo/sys.W
        # calculate spiral time to double
        Ep = hnpli*sys.n_p + lnpn*sys.l_p
        Er = hnpli*sys.n_r + lnpn*sys.l_r
        t_sl = - np.log(2.)/sys.W*(sys.l_p/hnpli + sys.n_p/lnpn)/\
            (hmpli/hnpli - lmpn/lnpn)
        t_sl_other = - np.log(2.)/sys.W*sys.g/sys.vo*(sys.Izz/lnpn)/\
            (hmpli/hnpli - lmpn/lnpn)

        dr_a = azzb + sys.Ixx*(sys.n_p/sys.l_p**2.*hmpli)
        dr_b = - sys.l_r*sys.n_p/sys.Y_b/sys.l_p
        dr_c = sys.g/sys.vo/sys.l_p*Er/lmpn
        s_dr = -0.5*( sys.Y_b/sys.Izz*(dr_a + dr_b) + dr_c )
        dr_d = hmpli*sys.n_p/sys.l_p - sys.W*rzzb**2.*hnpli/sys.vo/sys.l_p
        dr_e = - sys.Y_b/sys.Izz*(lmpn + dr_d + 0.25*sys.Y_b/sys.Izz*azzb**2.)
        wn_dr = ( s_dr**2. + dr_e )**0.5
        z_dr = s_dr / wn_dr
        s_dr_other = -0.5*sys.Y_b/sys.Izz*azzb + 0.*lmpn
        wn_dr_other = (-sys.Y_b/sys.Izz*lmpn)**0.5
        z_dr_other = s_dr_other / wn_dr_other

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
        t_sl_exact = ndm2*0.
        z_dr_exact = ndm2*0.
        wn_dr_exact = ndm2*0.
        s_dr_exact = ndm2*0.
        CAP_dr_exact = ndm2*0.
        for j in range(ndm2.shape[0]):
            # sys.Cn_b = -ndm2[j]*sys.CY_b
            # sys.Cl_b =  ldm[i]*sys.CY_b
            dx = -(lmpn_0 - ndm2[j]*rzzb)
            dz = -(hmpl_0 -  ldm[i]*rxxb)
            # if i == 0:
            #     print(lmpn_0/rzzb,ndm2[j],dx)
            # if j == 0:
            #     print(hmpl_0/rxxb,ldm[i],dz,hmpl_0 + dz)
            sys._buckingham_matrices_and_approximations(cg_shift=[dx,0.,dz])
            sys._lateral_dimensional_properties(sys.b,is_alternate=True)
            # sys._buckingham_matrices_and_approximations()
            # sys._lateral_dimensional_properties(sys.b,is_alternate=True)
            t_sl_exact[j] = sys.b["lat"]["double"][sys.b["lat"]["sl"]]*1.
            z_dr_exact[j] = sys.b["lat"]["zt"][sys.b["lat"]["dr"][0]]*1.
            s_dr_exact[j] = sys.b["lat"]["Sg"][sys.b["lat"]["dr"][0]]*1.
            wn_dr_exact[j] = sys.b["lat"]["wn"][sys.b["lat"]["dr"][0]]*1.
            CAP_dr_exact[j] = -wn_dr_exact[j]**2./sys.CY_b*sys.CW
        #     ###   ###   ###   ###   ###
        #     evals[0:2,j] = sys.b["lat"]["evals"][sys.b["lat"]["rb"]]
        #     evals[2:3,j] = sys.b["lat"]["evals"][sys.b["lat"]["ro"]]
        #     evals[3:4,j] = sys.b["lat"]["evals"][sys.b["lat"]["sl"]]
        #     evals[4:6,j] = sys.b["lat"]["evals"][sys.b["lat"]["dr"]]
        #     ###   ###   ###   ###   ###
        # cs = ["k","k","b","r","g","m"]
        # for l in range(6):
        #     if cs[l] not in []:#"m","k","g"]:
        #         for k in range(len(evals[l])):
        #             alg.plot(np.real(evals[l,k]),np.imag(evals[l,k]),".",c=cs[l],\
        #                 ms=float(k/dm2.shape[0])*5. + 1.0)#,c="k")
        #         alg.plot(np.real(evals[l,:]),np.imag(evals[l,:]),linewidth=0.2,c=cs[l])
        #         alg.set_title(sys.aircraft_name)
        # folder = "hand_qual_plots/"
        # aircraft_name = sys.aircraft_name.lower().replace(" ","_").replace("-","_")
        # name = "eigs_" + aircraft_name
        # savedictr = dict(transparent=True,format="png",dpi=300.0)
        # flg.savefig(folder + "lat_eigvals/" + name + "_{:>02d}".format(i) + ".png",**savedictr)
        # alg.cla()
        # return to normal
        sys.Cl_b = Cl_b_0*1.
        sys.Cn_b = Cn_b_0*1.
        sys._buckingham_matrices_and_approximations()
        sys._lateral_dimensional_properties(sys.b,is_alternate=True)

        # plot spiral
        asl.plot(ndm[:ipo],t_sl_other[:ipo],ls="--",c=str(i/num*0.95))
        asl.plot(ndm[ips:],t_sl[ips:],c=str(i/num*0.95))
        asl.plot(ndm2[2:],t_sl_exact[2:],c=str(i/num*0.95),marker="o",\
            mfc="none",ms=4.,ls="none")

        # plot Dutch-roll damping ratio
        adz.plot(ndm,z_dr_other,ls="--",c=str(i/num*0.95))
        adz.plot(ndm,z_dr,c=str(i/num*0.95))
        adz.plot(ndm2[2:],z_dr_exact[2:],c=str(i/num*0.95),marker="o",\
            mfc="none",ms=4.,ls="none")

        # plot Dutch-roll natural frequency
        adw.plot(ndm,wn_dr_other,ls="--",c=str(i/num*0.95))
        adw.plot(ndm,wn_dr,c=str(i/num*0.95))
        adw.plot(ndm2[2:],wn_dr_exact[2:],c=str(i/num*0.95),marker="o",\
            mfc="none",ms=4.,ls="none")

        # plot Dutch-roll damping rate
        ads.plot(ndm,s_dr_other,ls="--",c=str(i/num*0.95))
        ads.plot(ndm,s_dr,c=str(i/num*0.95))
        ads.plot(ndm2[2:],s_dr_exact[2:],c=str(i/num*0.95),marker="o",\
            mfc="none",ms=4.,ls="none")
        
        # plot Dutch-roll CAP

        adC.loglog(z_dr,CAP_dr,c=str(i/num*0.95))
        for q in range(0,len(CAP_dr)-1,int((num_n-1)/6)): # ind_list: #
            if adC_z_lim[0] <= z_dr[q] <= adC_z_lim[1] and \
                adC_C_lim[0] <= CAP_dr[q] <= adC_C_lim[1]:
                adC.annotate('', xytext=(z_dr[q],CAP_dr[q]),
                xy=(z_dr[q+1],CAP_dr[q+1]), arrowprops=dict(arrowstyle="->",
                color=str(i/num*0.95)), size=10.)
        # if i == 0:
        #     adC.plot(ndm,CAP_dr,c=str(i/num*0.95))
        #     adC.plot(ndm2[2:],CAP_dr_exact[2:],c=str(i/num*0.95),marker="o",\
        #         mfc="none",ms=4.,ls="none")
    
    # shade levels spiral
    sl_l3 = 4.
    alfa = 0.8
    if sys.aircraft_name in ["F-94A"]:
        sl_l1, sl_l2 = 12., 12.
    else:
        if sys.aircraft_name in ["Navion"]:
            sl_l1, sl_l2 = 20., 12.
        elif sys.aircraft_name in ["Lockheed Jetstar","Boeing 747"]:
            sl_l1, sl_l2 = 20., 12.
        else:
            sl_l1, sl_l2 = 20., 12.
        asl.text(0.0,sl_l1,"Level 2",va="bottom",ha="left",c="w",
            bbox=dict(facecolor="0.5",linewidth=0,alpha=alfa,
            boxstyle="Square, pad=0.0"))
    asl.fill_between([0.0,ndm_final],2*[sl_l3],2*[ 0.00],color="0.5",alpha=0.4)
    asl.text(0.2,100.,"Level 1",va="center",ha="right",
        bbox=dict(facecolor="w",linewidth=0,alpha=alfa,
        boxstyle="Square, pad=0.0"))
    asl.text(1.325,sl_l3,"Level 3",va="bottom",ha="right",
        bbox=dict(facecolor="w",linewidth=0,alpha=alfa,
        boxstyle="Square, pad=0.0"))
    asl.text(ndm_final,0.0,"Level 4",va="bottom",ha="right",c="w",
        bbox=dict(facecolor="0.5",linewidth=0,alpha=alfa,
        boxstyle="Square, pad=0.0"))
    # asl.fill_between([0.0,ndm_final],2*[sl_l2],2*[sl_l3],color="0.5",alpha=0.4)
    asl.fill_between([0.0,ndm_final],2*[sl_l1],2*[sl_l2],color="0.5",alpha=0.4)
    # asl.text(dm_final/2.,0.0,"Level 2",va="bottom",ha="center",
    #     bbox=dict(facecolor="w",linewidth=0,alpha=alfa,
    #     boxstyle="Square, pad=0.0"))
    # asl.text(dm_final/2.,(0.04+ph_limit)/2.,"Level 1",va="bottom",ha="center",
    #     bbox=dict(facecolor="w",linewidth=0,alpha=0.6,
    #     boxstyle="Square, pad=0.0"))

    # legend elements
    f_lo = 0         /num
    f_md = int(num/2)/num
    f_hi = (num - 1) /num
    lbl_name = r"$h_{mp_\ell}/r_{xx_b} =$ "
    lgnd_elms_sl = [
        Line2D([0], [0], c='k', ls='none',lw=1,marker="o",\
            mfc="none",ms=4., label='exact'),
        Line2D([0], [0], c='k', ls='-',lw=1, label='Eq. (88)'),
        Line2D([0], [0], c='k', ls='--',lw=1, label='Eq. (89)'),
        Line2D([0], [0], c=str(f_lo*0.95), ls='-',lw=1, label=lbl_name + \
            "{:>3.1f}".format(ldm[int(num*f_lo)])),
        Line2D([0], [0], c=str(f_md*0.95), ls='-',lw=1, label=lbl_name + \
            "{:>3.1f}".format(ldm[int(num*f_md)])),
        Line2D([0], [0], c=str(f_hi*0.95), ls='-',lw=1, label=lbl_name + \
            "{:>3.1f}".format(ldm[int(num*f_hi)]))
    ]
    lgnd_elms_dr = [ # zeta
        Line2D([0], [0], c='k', ls='none',lw=1,marker="o",\
            mfc="none",ms=4., label='exact'),
        Line2D([0], [0], c='k', ls='-',lw=1, label='Eq. (93)'),
        Line2D([0], [0], c='k', ls='--',lw=1, label='Eq. (96)'),
        Line2D([0], [0], c=str(f_lo*0.95), ls='-',lw=1, label=lbl_name + \
            "{:>3.1f}".format(ldm[int(num*f_lo)])),
        Line2D([0], [0], c=str(f_md*0.95), ls='-',lw=1, label=lbl_name + \
            "{:>3.1f}".format(ldm[int(num*f_md)])),
        Line2D([0], [0], c=str(f_hi*0.95), ls='-',lw=1, label=lbl_name + \
            "{:>3.1f}".format(ldm[int(num*f_hi)]))
    ]

    t_2m = 240.
    asl.set_xlim((0.,ndm_final)) # lmpn_sm_0,ndm_final)) # 
    asl.set_ylim((0.,t_2m))
    asl.set_xlabel(r"Yaw dynamic margin, $l_{mp_n}/r_{zz_b}$")
    asl.set_ylabel(r"Spiral time to double, $\tau_{sl}$ [s]")
    asl.legend(handles=lgnd_elms_sl)
    name = "t_sl_" + aircraft_name
    fsl.savefig(folder + "spiral/" + name + ".png",**savedict)

    if sys.aircraft_name == "Navion":
        z_limit = 0.75
        s_limit = 0.6
        w_limit = 3.5
    elif sys.aircraft_name == "F-94A":
        z_limit = 0.4
        s_limit = 0.25
        w_limit = 3.
    elif sys.aircraft_name == "Lockheed Jetstar":
        z_limit = 0.6
        s_limit = 0.25
        w_limit = 2.5
    elif sys.aircraft_name == "Boeing 747":
        z_limit = 0.6
        s_limit = 0.2
        w_limit = 1.5
    else:
        z_limit = 1.
        s_limit = 1.
        w_limit = 5.
    adz.set_xlim((0.,ndm_final))
    adz.set_ylim((0.,z_limit))
    adz.set_xlabel(r"Yaw dynamic margin, $l_{mp_n}/r_{zz_b}$")
    adz.set_ylabel(r"Dutch-roll damping ratio, $\zeta_{dr}$")
    adz.legend(handles=lgnd_elms_dr)
    name = "z_dr_" + aircraft_name
    fdz.savefig(folder + "dutch_roll_z/" + name + ".png",**savedict)

    lgnd_elms_dr[1] = Line2D([0], [0], c='k', ls='-',lw=1, label='Eq. (92)')
    lgnd_elms_dr[2] = Line2D([0], [0], c='k', ls='--',lw=1, label='Eq. (95)')
    adw.set_xlim((0.,ndm_final))
    adw.set_ylim((0.,w_limit))
    adw.set_xlabel(r"Yaw dynamic margin, $l_{mp_n}/r_{zz_b}$")
    adw.set_ylabel(r"Dutch-roll natural frequency, $\omega_{n_{dr}}$ [rad/s]")
    adw.legend(handles=lgnd_elms_dr)
    name = "wn_dr_" + aircraft_name
    fdw.savefig(folder + "dutch_roll_wn/" + name + ".png",**savedict)

    lgnd_elms_dr[1] = Line2D([0], [0], c='k', ls='-',lw=1, label='Eq. (90)')
    lgnd_elms_dr[2] = Line2D([0], [0], c='k', ls='--',lw=1, label='Eq. (94)')
    ads.set_xlim((0.,ndm_final))
    ads.set_ylim((0.,s_limit))
    ads.set_xlabel(r"Yaw dynamic margin, $l_{mp_n}/r_{zz_b}$")
    ads.set_ylabel(r"Dutch-roll damping rate, $\sigma_{dr}$ [s$^{-1}$]")
    ads.legend(handles=lgnd_elms_dr)
    name = "s_dr_" + aircraft_name
    fds.savefig(folder + "dutch_roll_s/" + name + ".png",**savedict)

    # draw limits
    wn_2_CAP = lambda wn : -wn**2./sys.CY_b*sys.CW
    CAP_2_wn = lambda CAP : (-CAP*sys.CY_b/sys.CW)**0.5
    if sys.aircraft_name in ["F-94A"]:
        # maneuvering lines -- A flight phases
        wn_1 = 1.0
        z_1  = 0.4
        s_1  = 0.4
    elif sys.aircraft_name in ["Navion"]:
        # cruise lines -- B flight phases
        wn_1 = 0.4
        z_1  = 0.08
        s_1  = 0.15
    elif sys.aircraft_name in ["Lockheed Jetstar","Boeing 747"]:
        # landing lines -- C flight phases
        wn_1 = 0.4
        z_1  = 0.08
        s_1  = 0.1
    else:
        # plot cruise lines
        wn_1 = 0.4
        z_1  = 0.08
        s_1  = 0.15
    CAP_1 = wn_2_CAP(wn_1)
    wn_3 = 0.4
    CAP_3 = wn_2_CAP(wn_3)
    adC.loglog([5.,0.01],[CAP_3,CAP_3],c="k")
    s_2 = 0.05
    wn_2 = 0.4
    z_2  = 0.02
    CAP_2 = wn_2_CAP(wn_2)
    z_s_2 = s_2/wn_2
    CAP_s_2 = wn_2_CAP(s_2/z_2)
    adC.loglog([5.,z_s_2,z_2,z_2],[CAP_2,CAP_2,CAP_s_2,20.],c="k")
    z_s_1 = s_1/wn_1
    CAP_s_1 = wn_2_CAP(s_1/z_1)
    adC.loglog([5.,z_s_1,z_1,z_1],[CAP_1,CAP_1,CAP_s_1,20.],c="k")
    # plot level labels
    alfa = 0.8
    text = adC.text(2.0,1.0,"Level 1",va="center",ha="center",
        bbox=dict(facecolor="w",linewidth=0,alpha=alfa,
        boxstyle="Square, pad=0.0"))
    text = adC.text(10.**(((np.log10(z_2)+np.log10(z_1)))/2.),
        CAP_s_2,"Level 2",va="bottom",ha="center",
        bbox=dict(facecolor="w",linewidth=0,alpha=alfa,
        boxstyle="Square, pad=0.0"))
    text = adC.text(0.02,1.,"Level 3",va="center",ha="center",
        bbox=dict(facecolor="w",linewidth=0,alpha=alfa,
        boxstyle="Square, pad=0.0"))
    text = adC.text(0.2,0.05,"Level 4",va="center",ha="center",
        bbox=dict(facecolor="w",linewidth=0,alpha=alfa,
        boxstyle="Square, pad=0.0"))
    adC.set_xlim(adC_z_lim)
    adC.set_ylim(adC_C_lim)
    adC.xaxis.set_major_formatter(ScalarFormatter())
    adC.yaxis.set_major_formatter(ScalarFormatter())
    adC.legend(handles=lgnd_elms_dr[3:])
    adC.set_xlabel(r"Dutch-roll damping ratio, $\zeta_{dr}$")
    adC.set_ylabel(r"Dutch-roll CAP [s$^{-2}$]")
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
    for i in range(num_craft): # [1,5,9,19]: # range(1,2): # 
        analyze_aircraft(run_files[i])

    # folder = "aircraft_database/"
    # filename = "Cessna_172.json"
    # analyze_aircraft(file_name=folder+filename)