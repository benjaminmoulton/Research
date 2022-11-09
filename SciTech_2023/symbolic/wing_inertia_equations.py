import sympy as sy
import numpy as np
from matplotlib import pyplot as plt
import time


if __name__ == "__main__":
    start = time.time()
    sym = sy.Symbol
    igr = sy.integrate
    simp = sy.simplify
    exp = sy.expand
    piecewise = sy.Piecewise
    diff = sy.diff
    sin = sy.sin
    cos = sy.cos
    tan = sy.tan
    mat = sy.Matrix
    # sy.init_printing(use_unicode=True)

    # declare variables
    print("declaring variables...")
    xp = sym("xp")
    yp = sym("yp")
    z = sym("z")
    # xcg = sym("xcg")
    # ycg = sym("ycg")
    # zcg = sym("zcg")
    # p = sym("p")
    b = sym("b")
    # Lr = sym("Lr") # sym("\u039B")
    # tanL = tan(Lr)
    ct = sym("ct")
    cr = sym("cr")
    tt = sym("tt")
    tr = sym("tr")

    a0 = sym("a_0")
    a1 = sym("a_1")
    a2 = sym("a_2")
    a3 = sym("a_3")
    a4 = sym("a_4")
    a5 = sym("a_5")

    # create bounds
    print("creating bounds...")
    c_eq = (ct-cr) * yp + cr
    tmax_eq = (tt-tr) * yp + tr
    t_eq = a0 * sy.sqrt(xp) + a1*(xp) + a2*(xp)**2 + a3*(xp)**3 +a4*(xp)**4
    t_eq = a0 * sy.sqrt(xp)
    x = c_eq * ( sy.Rational(1,4) - xp )
    y = b * yp
    # z = t_eq * tmax_eq * c_eq * ( zp - sy.Rational(1,2) )
    dx = c_eq
    dy = b
    dz = 1 # t_eq * tmax_eq * c_eq
    xp_up = 1; xp_lo = 0
    z_up = sy.Rational(1,2) * t_eq * tmax_eq * c_eq; z_lo = - z_up
    yp_up = 1; yp_lo = 0
    xp_bnd = (xp,xp_lo,xp_up)
    yp_bnd = (yp,yp_lo,yp_up)
    zp_bnd = (z,z_lo,z_up)
    # c_eq = (ct-cr) * yp + cr
    # tmax_eq = (tt-tr) * yp + tr
    # t_eq = 1#a0 * sy.sqrt(xp) + a1*(xp) + a2*(xp)**2 + a3*(xp)**3 +a4*(xp)**4
    # x = c_eq * ( sy.Rational(1,4) - xp )
    # y = b * yp
    # z = t_eq * tmax_eq * c_eq * ( sy.Rational(1,2) - zp )
    # dx = c_eq
    # dy = b
    # dz = t_eq * tmax_eq * c_eq
    # xp_up = 1; xp_lo = 0
    # # z_up = sy.Rational(1,2) * t_eq * tmax_eq * c_eq; z_lo = - z_up
    # zp_up = 1; zp_lo = 0
    # yp_up = 1; yp_lo = 0
    # xp_bnd = (xp,xp_lo,xp_up)
    # yp_bnd = (yp,yp_lo,yp_up)
    # zp_bnd = (zp,zp_lo,zp_up)

    # define S integration constants
    S000 = x**0*y**0*z**0 # 0
    S100 = x**1*y**0*z**0 # 1
    S010 = x**0*y**1*z**0 # 2
    S001 = x**0*y**0*z**1 # 3
    S110 = x**1*y**1*z**0 # 4
    S101 = x**1*y**0*z**1 # 5
    S011 = x**0*y**1*z**1 # 6
    S200 = x**2*y**0*z**0 # 7
    S020 = x**0*y**2*z**0 # 8
    S002 = x**0*y**0*z**2 # 9

    R = (40*a0 + 30*a1 + 20*a2 + 15*a3 + 12*a4) / 60 

    S_input = [S000,S100,S010,S001,S110,S101,S011,S200,S020,S002]
    S_input_names = ["S000","S100","S010","S001","S110","S101",
                    "S011","S200","S020","S002"]
    S = []
    print("Calculating...")
    for i in range(len(S_input)):
        print("\t{}...".format(S_input_names[i]))
        z_int = simp( igr(S_input[i],zp_bnd)*dz )
        # print("\t    z_int =",z_int)
        x_int = simp( igr(     z_int,xp_bnd)*dx )
        # print("\t    x_int =",x_int)
        y_int = simp( igr(     x_int,yp_bnd)*dy )
        # print("\t    y_int =",y_int)
        S.append(y_int)

        # if i != 9:
        #     # print(S[-1])
        #     print(simp(S[-1] / R))
    
    # report duration
    duration = time.time() - start
    print()
    print("calculation duration =", duration)
    print()

    
    
    for i in range(len(S)):
        # if i == 9:
        # print(S_input_names[i],"=",S[i])
        Sint = str(S[i]).replace("cr","c_r").replace("ct","c_t")
        Sint = Sint.replace("tr","\\tau_r").replace("tt","\\tau_t")
        Sint = Sint.replace("**","^").replace("*"," ")
        print(S_input_names[i],"=",Sint)
        print()

    subsvals = [
        [cr,1.0],
        [ct,0.5],
        [tr,0.12],
        [tt,0.08],
        # [Lr,0.0],
        [b,2.0],
        [a0,2.969], # 2.969,-1.26,-3.516,2.843,-1.036
        [a1,-1.260],
        [a2,-3.516],
        [a3,2.843],
        [a4,-1.036], # closed trailing edge
        # [a4,-1.015], # open trailing edge
        [a5,0.0]
    ]

    # substitute
    for i in range(len(S)):
        S[i] = S[i].subs(subsvals)
    
    m_eqs = S[0]
    cg_eqs = np.array([
        S[1] / S[0],
        S[2] / S[0],
        float(S[3] / S[0])
    ])
    I_eqs = np.array([
        S[8] + S[9],
        S[7] + S[9],
        S[7] + S[8],
        -S[4],
        float(S[5]),
        -float(S[6]),
    ])

    # properties
    # for wing 12 gc closed TE
    lbm_slug = 32.17405082
    m_wng = 2.66974726 / lbm_slug
    cg_wng = np.array([ 
        -0.13712589,
        0.72910392,
        0.00000000000001
    ])
    I_wng = np.array([
        2.16737992,
        0.15240118,
        2.31726753,
        0.23551551,
        -0.00000002,
        -0.00000100
    ]) / lbm_slug
    I_wng_abt_cg = np.array([
        0.74816221,
        0.10220056,
        0.84784920,
        -0.03140323,
        -0.00000015,
        -0.00000032
    ]) / lbm_slug

    # my eqs
    cr = 1.0 # ],
    ct = 0.5 # ],
    tr = 0.12 # ],
    tt = 0.08 # ],
    b = 2.0 # ],
    a0 = 2.969 # ], # 2.969,-1.26,-3.516,2.843,-1.036
    a1 = -1.260 # ],
    a2 = -3.516 # ],
    a3 = 2.843 # ],
    a4 = -1.036 # 
    cr2, ct2 = cr**2., ct**2.
    crct = cr * ct
    cr3, ct3 = cr**3., ct**3.
    cr2ct, crct2 = cr2 * ct, cr * ct2
    cr4, ct4 = cr**4., ct**4.
    cr3ct, cr2ct2, crct3 = cr3 * ct, cr2 * ct2, cr * ct3

    # calculate kappa values
    ka = tr * (3.*cr2 + 2.*crct + ct2) + tt *(cr2 + 2.*crct + 3.*ct2)
    one = 4.*cr3 + 3.*cr2ct + 2.*crct2 +    ct3
    two =    cr3 + 2.*cr2ct + 3.*crct2 + 4.*ct3
    kb = tr * one + tt * two
    one = 3.*cr2 + 4.*crct + 3.*ct2
    two =    cr2 + 3.*crct + 6.*ct2
    kc = tr * one + 2. * tt * two
    one = (cr + ct) * (2.*cr2 + crct + 2.*ct2)
    two = cr3 + 3.*cr2ct + 6.*crct2 + 10.*ct3
    kd = tr * one + tt * two
    one = 5.*cr4 + 4.*cr3ct + 3.*cr2ct2 + 2.*crct3 +    ct4
    two =    cr4 + 2.*cr3ct + 3.*cr2ct2 + 4.*crct3 + 5.*ct4
    ke = tr * one + tt * two
    one = cr2 + 2.*crct +  2.*ct2
    two = cr2 + 4.*crct + 10.*ct2
    kf = tr * one + tt * two
    one = 35.*cr4 + 20.*cr3ct + 10.*cr2ct2 +  4.*crct3 +     ct4
    two = 15.*cr4 + 20.*cr3ct + 18.*cr2ct2 + 12.*crct3 +  5.*ct4
    thr =  5.*cr4 + 12.*cr3ct + 18.*cr2ct2 + 20.*crct3 + 15.*ct4
    fou =     cr4 +  4.*cr3ct + 10.*cr2ct2 + 20.*crct3 + 35.*ct4
    kg = tr**3. *one + tr**2.*tt *two + tr*tt**2. *thr + tt**3. *fou
    # calculate upsilon values
    u0 = 1./ 60.*( 40.*a0+ 30.*a1+ 20.*a2+ 15.*a3+ 12.*a4)
    u1 = 1./ 60.*( 56.*a0+ 50.*a1+ 40.*a2+ 33.*a3+ 28.*a4)
    u2 = 1./980.*(856.*a0+770.*a1+644.*a2+553.*a3+484.*a4)
    one1 = 2./5.*a0**3. + 1*a0**2.*a1 + 3./4.*a0**2.*a2 + 3./5.*a0**2.*a3
    one2 = 1./2.*a0**2.*a4 + 6./7.*a0*a1**2. + 4./3.*a0*a1*a2
    one3 = 12./11.*a0*a1*a3 + 12./13.*a0*a1*a4
    two1 = 6./11*a0*a2**2. + 12./13.*a0*a2*a3 + 4./5.*a0*a2*a4 +2./5.*a0*a3**2.
    two2 = 12./17.*a0*a3*a4 + 6./19.*a0*a4**2. + 1./4.*a1**3. + 3./5.*a1**2.*a2
    thr1 = 1./2.*a1**2.*a3 + 3./7.*a1**2.*a4 + 1./2.*a1*a2**2. + 6./7.*a1*a2*a3
    thr2 = 3./4.*a1*a2*a4 + 3./8.*a1*a3**2. + 2./3.*a1*a3*a4 + 3./10.*a1*a4**2.
    fou1 = 1./7.*a2**3. + 3./8.*a2**2.*a3 + 1./3.*a2**2.*a4 + 1./3.*a2*a3**2.
    fou2 = 3./5.*a2*a3*a4 + 3./11.*a2*a4**2. + 1./10.*a3**3. + 3./11.*a3**2.*a4
    fou3 = 1./4.*a3*a4**2. + 1./13.*a4**3.
    u3 = one1 + one2 + one3 + two1 + two2 + thr1 + thr2 + fou1 + fou2 + fou3
    SP_S = [
        b/12. * ka * u0,
        -b/80. * kb * u1,
        b**2. / 60. * kc * u0,
        0.0,
        -b**2./240.*kd*u1,
        0.0,
        0.0,
        7.*b/1440.*ke*u2,
        b**3./60.*kf*u0,
        b/3360.*kg*u3
    ]
    # for i in range(len(SP_S)):
    #     print(S[i])
    #     print(SP_S[i])
    #     print()
    
    m_sps = SP_S[0]
    cg_sps = np.array([
        SP_S[1] / SP_S[0], # wrong
        SP_S[2] / SP_S[0],
        float(SP_S[3] / SP_S[0])
    ])
    I_sps = np.array([
        SP_S[8] + SP_S[9],
        SP_S[7] + SP_S[9],
        SP_S[7] + SP_S[8],
        -SP_S[4],
        float(SP_S[5]),
        -float(SP_S[6]),
    ])

    eq_props = np.concatenate(([m_eqs],cg_eqs,I_eqs))
    wn_props = np.concatenate(([m_wng],cg_wng,I_wng))
    sp_props = np.concatenate(([m_sps],cg_sps,I_sps))

    pd_props = (wn_props - eq_props) / wn_props
    pe_props = (wn_props - sp_props) / wn_props
    names = ["mass","cg-x","cg-y","cg-z","Ixx","Iyy","Izz","Ixy","Ixz","Iyz"]

    for i in range(eq_props.shape[0]):
        print(names[i])
        print("eq       = {:< 16.12}".format(eq_props[i]))
        print("SCPP     = {:< 16.12}".format(sp_props[i]))
        print("SW       = {:< 16.12}".format(wn_props[i]))
        print("PD SW eq = {:< 16.12}".format(pd_props[i]))
        print("PD SW SC = {:< 16.12}".format(pe_props[i]))
        print()









    # # for wing: 12 gc closed TE, 18 gc open TE, 6 gc open TE
    # Ixxo_wing = 2.16737992 / 32.17405082 # 2.17819785 / 32.17405082 # 2.19446398 / 32.17405082
    # Iyyo_wing = 0.15240118 / 32.17405082 # 0.15639797 / 32.17405082 # 0.15697218 / 32.17405082
    # Izzo_wing = 2.31726753 / 32.17405082 # 2.33207533 / 32.17405082 # 2.34888376 / 32.17405082
    

    # 6 gc, open TE
    # 0.067633982213,  0.004858720779,  0.072414516344
    # 0.098749335714,  0.012259752381,  0.110843750000
    # 0.098749335714,  0.012259752381,  0.110843750000

    # 0.068206020817,  0.004878844162,  0.073005533967
    # 0.098749371901,  0.012259755920,  0.110843788367

    # 0.008386922411,  0.004124620905,  0.008095518117
    # 0.000000366455,  0.000000288644,  0.000000346135
    # 18 gc, open TE
    # 0.067633982213,  0.004858720779,  0.072414516344
    # 0.098749335714,  0.012259752381,  0.110843750000
    # 0.098749335714,  0.012259752381,  0.110843750000

    # 0.067700454077,  0.004860997171,  0.072483112029
    # 0.098749371901,  0.012259755920,  0.110843788367

    # 0.000981852559,  0.000468297224,  0.000946367815
    # 0.000000366455,  0.000000288644,  0.000000346135
    # 12 gc, closed TE
    # 0.067219395146,  0.004732239963,  0.071873822594
    # 0.098749335714,  0.012259752381,  0.110843750000
    # 0.098749335714,  0.012259752381,  0.110843750000

    # 0.067364222557,  0.004736773148,  0.072022871567
    # 0.098749371901,  0.012259755920,  0.110843788367

    # 0.002149915862,  0.000957019756,  0.002069467234
    # 0.000000366455,  0.000000288644,  0.000000346135
