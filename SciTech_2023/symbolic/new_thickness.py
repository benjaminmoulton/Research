import sympy as sy
import numpy as np
from matplotlib import pyplot as plt


if __name__ == "__main__":
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
    x = sym("x")
    y = sym("y")
    z = sym("z")
    # xcg = sym("xcg")
    # ycg = sym("ycg")
    # zcg = sym("zcg")
    # p = sym("p")
    b = sym("b")
    Lr = sym("Lr") # sym("\u039B")
    tanL = tan(Lr)
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
    x_up =   sy.Rational(1,4) * ( (ct-cr) / b * y + cr )
    x_lo = - sy.Rational(3,4) * ( (ct-cr) / b * y + cr )
    x_c = sy.Rational(1,4) - x / ( (ct-cr) / b * y + cr )
    t_eq = a0 * sy.sqrt(x_c) + a1*(x_c) + a2*(x_c)**2 + a3*(x_c)**3 +a4*(x_c)**4
    z_up =   sy.Rational(1,2) * ((ct-cr) / b * y + cr) * ((tt-tr) / b * y + tr) #* t_eq
    # z_up =   sy.Rational(1,2) * ((ct*tt-cr*tr) / b * y + cr*tr) * t_eq
    z_lo = - z_up
    y_up = b
    y_lo = 0
    x_bnd = (x,x_lo,x_up)
    y_bnd = (y,y_lo,y_up)
    z_bnd = (z,z_lo,z_up)

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

    S_input = [S000,S100,S010,S001,S110,S101,S011,S200,S020,S002]
    S_input_names = ["S000","S100","S010","S001","S110","S101",
                    "S011","S200","S020","S002"]
    S = []
    Sdy = []
    print("Calculating...")
    for i in range(len(S_input)):
        print("\t{}...".format(S_input_names[i]))
        # print(simp(igr(S_input[i],z_bnd)))
        Sdy.append(simp(igr(S_input[i],z_bnd,x_bnd)))
        # print(Sdy[-1])
        S.append(simp(igr(Sdy[-1],y_bnd)))
        # print(S[-1])
    
    for i in range(len(S)):
        # if i == 9:
        print(S_input_names[i],"=",S[i])
        # Sint = str(S[i]).replace("cr","c_r").replace("ct","c_t")
        # Sint = Sint.replace("tr","\\tau_r").replace("tt","\\tau_t")
        # Sint = Sint.replace("**","^").replace("*"," ")
        # print(S_input_names[i],"=",Sint)
        print()
    
    # x_up =   sy.Rational(1,4)
    # x_lo = - sy.Rational(3,4)
    x_up = 1
    x_lo = 0
    # x_c = sy.Rational(1,4) - x
    xp = sym("xp")
    t_eq = a0 * sy.sqrt(xp) + a1*(xp) + a2*(xp)**2 + a3*(xp)**3 +a4*(xp)**4 + a5
    # t_eq = t_eq.subs(x_c,sy.Rational(1,4) - x)
    c_eq = (ct-cr) * y + cr
    tmax_eq = (tt-tr) * y + tr
    tmax_in_z_eq = c_eq * tmax_eq
    # tmax_in_z_eq = (tt*ct-tr*cr) / b * y + tr*cr
    z_up = sy.Rational(1,2) * tmax_in_z_eq * t_eq
    z_lo = - z_up
    y_up = 1
    y_lo = 0
    x_bnd = (xp,x_lo,x_up)
    y_bnd = (y,y_lo,y_up)
    z_bnd = (z,z_lo,z_up)
    S200 = igr(igr(igr((sy.Rational(1,4)-xp)**2*c_eq**2,z_bnd),x_bnd)*c_eq,y_bnd) * b
    S020 = igr(igr(igr(y**2*b**2,z_bnd),x_bnd)*c_eq,y_bnd) * b
    df = igr(z**2,z_bnd) # *c_eq**2
    # print(df)
    func = igr(df,x_bnd)*c_eq
    # print(func)
    igral = simp(igr(func,y_bnd)) * b
    print("S002s=",simp(igral.subs([[a0,0],[a1,0],[a2,0],[a3,0],[a4,0],[a5,1]])))
    print()
    print("S002 =",igral)

    subsvals = [
        [cr,1.0],
        [ct,0.5],
        [tr,0.12],
        [tt,0.08],
        [Lr,0.0],
        [b,2.0]
    ]

    subsavals = [ # 2.969,-1.26,-3.516,2.843,-1.036
        [a0,2.969],
        [a1,-1.260],
        [a2,-3.516],
        [a3,2.843],
        [a4,-1.036], # closed trailing edge
        # [a4,-1.015], # open trailing edge
        [a5,0.0]
    ]
    subsa51vals = [
        [a0,0.0],
        [a1,0.0],
        [a2,0.0],
        [a3,0.0],
        [a4,0.0],
        [a5,1.0]
    ]

    S200 = S200.subs(subsvals)
    S020 = S020.subs(subsvals)
    S002 = igral.subs(subsvals)

    # differences
    S200_as = S200.subs(subsavals)
    S020_as = S020.subs(subsavals)
    S002_as = S002.subs(subsavals)
    S200_a51 = S200.subs(subsa51vals)
    S020_a51 = S020.subs(subsa51vals)
    S002_a51 = S002.subs(subsa51vals)

    Iyyo_as  = S200_as  + S002_as
    Iyyo_a51 = S200_a51 + S002_a51
    Ixxo_as  = S020_as  + S002_as
    Ixxo_a51 = S020_a51 + S002_a51
    Izzo_as  = S200_as  + S020_as
    Izzo_a51 = S200_a51 + S020_a51

    print("{:< 15.12f}, {:< 15.12f}, {:< 15.12f}".format(Ixxo_as,Iyyo_as,Izzo_as))
    print("{:< 15.12f}, {:< 15.12f}, {:< 15.12f}".format(Ixxo_a51,Iyyo_a51,Izzo_a51))
    S7 = S[7].subs(subsvals); S8 = S[8].subs(subsvals); S9 = S[9].subs(subsvals)
    print("{:< 15.12f}, {:< 15.12f}, {:< 15.12f}".format(S8+S9,S7+S9,S7+S8))
    print()


    # properties
    # for wing: 12 gc closed TE, 18 gc open TE, 6 gc open TE
    Ixxo_wing = 2.16737992 / 32.17405082 # 2.17819785 / 32.17405082 # 2.19446398 / 32.17405082
    Iyyo_wing = 0.15240118 / 32.17405082 # 0.15639797 / 32.17405082 # 0.15697218 / 32.17405082
    Izzo_wing = 2.31726753 / 32.17405082 # 2.33207533 / 32.17405082 # 2.34888376 / 32.17405082
    Ixxo_pprism = 3.17716731 / 32.17405082
    Iyyo_pprism = 0.39444601 / 32.17405082
    Izzo_pprism = 3.56629368 / 32.17405082
    print("{:< 15.12f}, {:< 15.12f}, {:< 15.12f}".format(Ixxo_wing,Iyyo_wing,Izzo_wing))
    print("{:< 15.12f}, {:< 15.12f}, {:< 15.12f}".format(Ixxo_pprism,Iyyo_pprism,Izzo_pprism))
    print()
    pdxx_as  = (Ixxo_wing - Ixxo_as) / Ixxo_wing
    pdyy_as  = (Iyyo_wing - Iyyo_as) / Iyyo_wing
    pdzz_as  = (Izzo_wing - Izzo_as) / Izzo_wing
    pdxx_a51 = (Ixxo_pprism - Ixxo_a51) / Ixxo_pprism
    pdyy_a51 = (Iyyo_pprism - Iyyo_a51) / Iyyo_pprism
    pdzz_a51 = (Izzo_pprism - Izzo_a51) / Izzo_pprism
    print("{:< 15.12f}, {:< 15.12f}, {:< 15.12f}".format(pdxx_as,pdyy_as,pdzz_as))
    print("{:< 15.12f}, {:< 15.12f}, {:< 15.12f}".format(pdxx_a51,pdyy_a51,pdzz_a51))
    # 12:40 - 12:45

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


    # a = sym("a")
    # b = sym("b")
    # c = sym("c")
    # d = sym("d")
    # f = sym("f")
    # eq = Cw*Cw * (a*y**4+b*y**3+c*y**2+d*y+f)
    # print(simp(eq.subs(y**4,0).subs(y**3,0)))
    # integ = igr(eq,y_bnd)
    # print(simp(integ))

    # print("Determining untwisted properties...")
    # # mass
    # print("\tmass...")
    # m = S[0]
    # # moments
    # print("\tmoments...")
    # M = simp(mat([[S[1] - tanL * S[2]],
    #         [S[2]],
    #         [S[3]]]))
    # # inertias # Ixx, Iyy, Izz, Ixy, Ixz, Iyz
    # print("\tinertias...")
    # i_main = mat([  [S[8] + S[9]],
    #                 [S[7] + S[9]],
    #                 [S[7] + S[8]],
    #                 [S[4]], ######################################### SIGNS ON prod terms may not be correct
    #                 [S[5]],
    #                 [S[6]]])
    # i_sweep = mat([  [0],
    #                 [tanL**2*S[8] - 2*tanL*S[4]],
    #                 [tanL**2*S[8] - 2*tanL*S[4]],
    #                 [- tanL*S[8]],
    #                 [- tanL*S[6]],
    #                 [0]])
    # I = simp(i_main + i_sweep)
    # # print(I)