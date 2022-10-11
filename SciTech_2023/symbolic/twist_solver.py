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
    xcg = sym("xcg")
    ycg = sym("ycg")
    zcg = sym("zcg")
    p = sym("p")
    b_2 = sym("b_2")
    Lr = sym("Lr") # sym("\u039B")
    tanL = tan(Lr)
    ct = sym("ct")
    cr = sym("cr")
    tt = sym("tt")
    tr = sym("tr")
    # twist variables
    wr = sym("wr")
    w  = sym("w")
    w  = 0
    Sw = sy.sin( w * y / b_2 + wr)
    Cw = sy.cos( w * y / b_2 + wr)

    # create bounds
    print("creating bounds...")
    x_up =   sy.Rational(1,4) * ( (ct-cr) / b_2 * y + cr )
    x_lo = - sy.Rational(3,4) * ( (ct-cr) / b_2 * y + cr )
    z_up =   sy.Rational(1,2) * ((tt*ct-tr*cr) / b_2 * y + tr*cr)
    z_lo = - sy.Rational(1,2) * ((tt*ct-tr*cr) / b_2 * y + tr*cr)
    y_up = b_2
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
        Sdy.append(simp(igr(S_input[i],x_bnd,z_bnd)))
        S.append(simp(igr(Sdy[-1],y_bnd)))
    
    # for i in range(len(S)):
    #     # print(S[i])
    #     # print(simp(S[i]))
    #     print(Sdy[i])
    #     print(simp(Sdy[i]))
    #     print()

    # a = sym("a")
    # b = sym("b")
    # c = sym("c")
    # d = sym("d")
    # f = sym("f")
    # eq = Cw*Cw * (a*y**4+b*y**3+c*y**2+d*y+f)
    # print(simp(eq.subs(y**4,0).subs(y**3,0)))
    # integ = igr(eq,y_bnd)
    # print(simp(integ))

    print("Determining untwisted properties...")
    # mass
    print("\tmass...")
    m = S[0]
    # moments
    print("\tmoments...")
    M = simp(mat([[S[1] - tanL * S[2]],
            [S[2]],
            [S[3]]]))
    # inertias # Ixx, Iyy, Izz, Ixy, Ixz, Iyz
    print("\tinertias...")
    i_main = mat([  [S[8] + S[9]],
                    [S[7] + S[9]],
                    [S[7] + S[8]],
                    [S[4]], ######################################### SIGNS ON prod terms may not be correct
                    [S[5]],
                    [S[6]]])
    i_sweep = mat([  [0],
                    [tanL**2*S[8] - 2*tanL*S[4]],
                    [tanL**2*S[8] - 2*tanL*S[4]],
                    [- tanL*S[8]],
                    [- tanL*S[6]],
                    [0]])
    I = simp(i_main + i_sweep)
    # print(I)

    # twist
    print("Determining twisted properties...")
    # moments
    print("\tmoments...")
    M_twist = simp(igr(mat([[Cw*Sdy[1] + Sw*Sdy[3]],
                        [Sdy[2]],
                        [-Sw*Sdy[1] + Cw*Sdy[3]]]),y_bnd))
    # inertias
    print("\tinertias...")
    i_m = mat([  [Sdy[8] + Sdy[9]],
                    [Sdy[7] + Sdy[9]],
                    [Sdy[7] + Sdy[8]],
                    [-Sdy[4]],
                    [Sdy[5]],
                    [-Sdy[6]]])
    i_twist = simp(igr(mat([[Cw*Cw*i_m[0] - 2*Cw*Sw*i_m[4] + Sw*Sw*i_m[2]],
                    [i_m[1]],
                    [Sw*Sw*i_m[0] + 2*Cw*Sw*i_m[4] + Cw*Cw*i_m[2]],
                    [Cw*i_m[3] + Sw*i_m[5]], # negated
                    [Cw*Sw*(i_m[0] - i_m[2]) + i_m[4]*(Cw*Cw - Sw*Sw)], # negated
                    [Cw*i_m[5] - Sw*i_m[3]] # negated
    ]),y_bnd))
    i_sweep = mat([  [0],
                    [tanL**2*S[8] - 2*tanL*S[4]],
                    [tanL**2*S[8] - 2*tanL*S[4]],
                    [tanL*S[8]],
                    [- tanL*S[6]],
                    [0]])
    I_twist = simp(i_twist + i_sweep)

    print("Reporting values...")
    I_names = ["Ixx","Iyy","Izz","Ixy","Ixz","Iyz"]
    for i in range(6):
        print("\t{}...".format(I_names[i]))
        print(simp(I[i]))
        print("--")
        print(simp(I_twist[i]))#.subs(wr,0)))
        print()



