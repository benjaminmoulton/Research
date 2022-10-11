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

    # create bounds
    print("creating bounds...")
    x_up =   sy.Rational(1,4) * ( (ct-cr) / b * y + cr )
    x_lo = - sy.Rational(3,4) * ( (ct-cr) / b * y + cr )
    x_c = sy.Rational(1,4) - x / ( (ct-cr) / b * y + cr )
    t_eq = a0 * sy.sqrt(x_c) + a1*(x_c) + a2*(x_c)**2 + a3*(x_c)**3 +a4*(x_c)**4
    z_up =   sy.Rational(1,2) * ((ct-cr) / b * y + cr) * ((tt-tr) / b * y + tr) * t_eq
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
        print(simp(igr(S_input[i],z_bnd)))
        Sdy.append(simp(igr(S_input[i],z_bnd,x_bnd)))
        print(Sdy[-1])
        S.append(simp(igr(Sdy[-1],y_bnd)))
        print(S[-1])
    
    for i in range(len(S)):
        print(S_input_names[i],"=",S[i])
        Sint = str(S[i]).replace("cr","c_r").replace("ct","c_t")
        Sint = Sint.replace("tr","\\tau_r").replace("tt","\\tau_t")
        Sint = Sint.replace("**","^").replace("*"," ")
        print(S_input_names[i],"=",Sint)
        # print(simp(S[i]))
        # print(Sdy[i])
        # print(simp(Sdy[i]))
        print()

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