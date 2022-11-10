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

    r = sym("r",positive=True)

    # create bounds
    print("creating bounds...")
    c_eq = (ct-cr) * yp + cr
    tmax_eq = (tt-tr) * yp + tr
    t_eq = sy.Piecewise((xp/r, xp < r), ((1-xp)/(1-r), xp >= r), (0, True))
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

    # subsvals = [
    #     [cr,1.0],
    #     [ct,0.5],
    #     [tr,0.12],
    #     [tt,0.08],
    #     # [Lr,0.0],
    #     [b,2.0],
    #     [a0,2.969], # 2.969,-1.26,-3.516,2.843,-1.036
    #     [a1,-1.260],
    #     [a2,-3.516],
    #     [a3,2.843],
    #     [a4,-1.036], # closed trailing edge
    #     # [a4,-1.015], # open trailing edge
    #     [a5,0.0]
    # ]

    # # substitute
    # for i in range(len(S)):
    #     S[i] = S[i].subs(subsvals)
    
    # m_eqs = S[0]
    # cg_eqs = np.array([
    #     S[1] / S[0],
    #     S[2] / S[0],
    #     float(S[3] / S[0])
    # ])
    # I_eqs = np.array([
    #     S[8] + S[9],
    #     S[7] + S[9],
    #     S[7] + S[8],
    #     -S[4],
    #     float(S[5]),
    #     -float(S[6]),
    # ])

    # # properties
    # # for wing 12 gc closed TE
    # lbm_slug = 32.17405082
    # m_wng = 2.66974726 / lbm_slug
    # cg_wng = np.array([ 
    #     -0.13712589,
    #     0.72910392,
    #     0.00000000000001
    # ])
    # I_wng = np.array([
    #     2.16737992,
    #     0.15240118,
    #     2.31726753,
    #     -0.23551551,
    #     0.00000002,
    #     0.00000100
    # ]) / lbm_slug

    # eq_props = np.concatenate(([m_eqs],cg_eqs,I_eqs))
    # wn_props = np.concatenate(([m_wng],cg_wng,I_wng))

    # pd_props = (wn_props - eq_props) / wn_props
    # names = ["mass","cg-x","cg-y","cg-z","Ixx","Iyy","Izz","Ixy","Ixz","Iyz"]

    # for i in range(eq_props.shape[0]):
    #     print(names[i])
    #     print("eq = {:< 16.12}".format(eq_props[i]))
    #     print("SW = {:< 16.12}".format(wn_props[i]))
    #     print("PD = {:< 16.12}".format(pd_props[i]))
    #     print()









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
