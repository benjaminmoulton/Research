import sympy as sy


if __name__ == "__main__":
    sym = sy.Symbol
    igr = sy.integrate
    simp = sy.simplify
    exp = sy.expand
    piecewise = sy.Piecewise

    # declare variables
    print("declaring variables...")
    xh = sym("xh")
    xhmt = sym("xhmt",positive=True,finite=True,nonzero=True,real=True)

    # naca 4-digit, 4-digit modified
    a0 = sym("a_0")
    a1 = sym("a_1")
    a2 = sym("a_2")
    a3 = sym("a_3")
    a4 = sym("a_4")
    d0 = sym("d_0")
    d1 = sym("d_1")
    d2 = sym("d_2")
    d3 = sym("d_3")

    # create bounds
    print("creating bounds...")
    types = ["naca4digit","naca4digitmodified","diamond"]
    tau_type = types[1]
    # tau equations divided by tau_max
    # naca 4 digit
    if tau_type == "naca4digit":
        tau_eq = a0 * sy.sqrt(xh) + a1*(xh) + a2*(xh)**2 + a3*(xh)**3 +a4*(xh)**4
    # naca 4 digit modified
    if tau_type == "naca4digitmodified":
        LE = a0 * sy.sqrt(xh) + a1*(  xh) + a2*(  xh)**2 + a3*(  xh)**3
        TE = d0               + d1*(1-xh) + d2*(1-xh)**2 + d3*(1-xh)**3
        tau_eq = sy.Piecewise((LE, xh < xhmt), (TE, xh >= xhmt), (0, True))
    # diamond
    if tau_type == "diamond":
        tau_eq = sy.Piecewise((xh/xhmt, xh < xhmt), ((1-xh)/(1-xhmt), xh >= xhmt), (0, True))
    x = ( sy.Rational(1,4) - xh )
    xh_up = 1; xh_lo = 0
    xh_bnd = (xh,xh_lo,xh_up)

    # integrate
    u0 = simp(                     igr(       tau_eq   , xh_bnd, conds="piecewise") )
    u1 = simp( -4                * igr(x    * tau_eq   , xh_bnd, conds="piecewise") )
    u2 = simp( sy.Rational(48,7) * igr(x**2 * tau_eq   , xh_bnd, conds="piecewise") )
    u3 = simp(                     igr(       tau_eq**3, xh_bnd, conds="piecewise") )

    piecewise_types = ["naca4digitmodified","diamond"]
    if tau_type in piecewise_types: # if piecewise equation
        y = sym("y")
        y_bnd = (y,0,1)
        u0 = simp(igr(u0,y_bnd).args[1][0])
        u1 = simp(igr(u1,y_bnd).args[1][0])
        u2 = simp(igr(u2,y_bnd).args[1][0])
        u3 = simp(igr(u3,y_bnd).args[1][0])

    # print("u0 =", u0)
    # print("u1 =", u1)
    # print("u2 =", u2)
    # print("u3 =", u3)
    us = [u0,u1,u2,u3]
    names = ["u0","u1","u2","u3"]
    for i in range(len(names)):
        ueq = str(us[i]).replace("xhmt","x")#r"\hat{x}_{mt}")
        ueq = ueq.replace("**","^").replace("*"," ")
        print(names[i],"=",ueq)
        print()