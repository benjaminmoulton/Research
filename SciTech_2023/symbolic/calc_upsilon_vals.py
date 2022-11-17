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
    # tau equations divided by tau_max
    # # naca 4 digit
    # tau_eq = a0 * sy.sqrt(xh) + a1*(xh) + a2*(xh)**2 + a3*(xh)**3 +a4*(xh)**4
    # # naca 4 digit modified
    # LE = a0 * sy.sqrt(xh) + a1*(  xh) + a2*(  xh)**2 + a3*(  xh)**3
    # TE = d0               + d1*(1-xh) + d2*(1-xh)**2 + d3*(1-xh)**3
    # tau_eq = sy.Piecewise((LE, xh < xhmt), (TE, xh >= xhmt), (0, True))
    # diamond
    tau_eq = sy.Piecewise((xh/xhmt, xh < xhmt), ((1-xh)/(1-xhmt), xh >= xhmt), (0, True))
    x = ( sy.Rational(1,4) - xh )
    xh_up = 1; xh_lo = 0
    xh_bnd = (xh,xh_lo,xh_up)

    # integrate
    u0 = simp(                     igr(       tau_eq   , xh_bnd, conds="piecewise") )
    u1 = simp( -4                * igr(x    * tau_eq   , xh_bnd, conds="piecewise") )
    u2 = simp( sy.Rational(48,7) * igr(x**2 * tau_eq   , xh_bnd, conds="piecewise") )
    u3 = simp(                     igr(       tau_eq**3, xh_bnd, conds="piecewise") )

    # # solve with bounds on xhmt
    # q = sym("q")
    # u0 = sy.solve([xhmt > 0, xhmt < 1, q - u0],q)
    # u1 = sy.solve([xhmt > 0, xhmt < 1, q - u1],q)
    # u2 = sy.solve([xhmt > 0, xhmt < 1, q - u2],q)
    # u3 = sy.solve([xhmt > 0, xhmt < 1, q - u3],q)
    print(u0)
    u0 = sy.refine(u0, sy.Q.negative(xhmt - 1))
    print(u0)

    print("u0 =", u0)
    print("u1 =", u1)
    print("u2 =", u2)
    print("u3 =", u3)