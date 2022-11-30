import sympy as sy


if __name__ == "__main__":
    sym = sy.Symbol
    igr = sy.integrate
    simp = sy.simplify
    exp = sy.expand
    piecewise = sy.Piecewise

    # declare variables
    print("declaring variables...")
    x = sym("x")
    y = sym("y")
    z = sym("z")
    r = sym("r")
    theta = sym("theta")
    Nb = sym("Nb")
    u0 = sym("u0")
    cr = sym("cr")
    ct = sym("ct")
    tt = sym("tt")
    tr = sym("tr")
    rr = sym("rr")
    rt = sym("rt")

    # define thickness and chord functions
    tau_m = (tt-tr)*(r-rr)/(rt-rr) + tr
    c = (ct-cr)*(r-rr)/(rt-rr) + cr

    # define height function
    h = Nb * tau_m * c**2 * u0 / (2 * sy.pi * r)

    # define bounds
    r_up = rt
    r_lo = rr
    theta_up = 2 * sy.pi
    theta_lo = 0
    x_up = sy.Rational(1,2) * h
    x_lo = - sy.Rational(1,2) * h
    x_bnd = (x,x_lo,x_up)
    r_bnd = (r,r_lo,r_up)
    theta_bnd = (theta,theta_lo,theta_up)

    i = x
    j = r * sy.cos(theta)
    k = r * sy.sin(theta)

    S = [
        [0,0,0],
        [1,0,0],
        [0,1,0],
        [0,0,1],
        [1,1,0],
        [1,0,1],
        [0,1,1],
        [2,0,0],
        [0,2,0],
        [0,0,2]
    ]

    for i in range(len(S)):
        f = i**S[i][0] * j**S[i][1] * k**S[i][2]
        iii = simp( igr(f * r, x_bnd,r_bnd,theta_bnd) )
        print("S{}{}{} = {}".format(S[i][0],S[i][1],S[i][2],iii))


    
    ka = 3*cr**2*tr + cr**2*tt + 2*cr*ct*tr + 2*cr*ct*tt + ct**2*tr + 3*ct**2*tt
    ka = tr*(3*cr**2 + 2*cr*ct + ct**2) + tt*(cr**2 + 2*cr*ct + 3*ct**2)
    S000 =    Nb*u0*(-3*cr**2*rr*tr - cr**2*rr*tt + 3*cr**2*rt*tr + cr**2*rt*tt - 2*cr*ct*rr*tr - 2*cr*ct*rr*tt + 2*cr*ct*rt*tr + 2*cr*ct*rt*tt - ct**2*rr*tr - 3*ct**2*rr*tt + ct**2*rt*tr + 3*ct**2*rt*tt)/12
    S100 =    Nb*u0*(-3*cr**2*rr*tr - cr**2*rr*tt + 3*cr**2*rt*tr + cr**2*rt*tt - 2*cr*ct*rr*tr - 2*cr*ct*rr*tt + 2*cr*ct*rt*tr + 2*cr*ct*rt*tt - ct**2*rr*tr - 3*ct**2*rr*tt + ct**2*rt*tr + 3*ct**2*rt*tt)/12
    S010 = 0
    S001 = 0
    S110 = 0
    S101 = 0
    S011 = 0
    S200 = 49*Nb*u0*(-3*cr**2*rr*tr - cr**2*rr*tt + 3*cr**2*rt*tr + cr**2*rt*tt - 2*cr*ct*rr*tr - 2*cr*ct*rr*tt + 2*cr*ct*rt*tr + 2*cr*ct*rt*tt - ct**2*rr*tr - 3*ct**2*rr*tt + ct**2*rt*tr + 3*ct**2*rt*tt)/12
    S020 = Nb*u0*(-10*cr**2*rr**3*tr - 2*cr**2*rr**3*tt + 6*cr**2*rr**2*rt*tr + 3*cr**2*rr*rt**2*tr + cr**2*rr*rt**2*tt + cr**2*rt**3*tr + cr**2*rt**3*tt - 4*cr*ct*rr**3*tr - 2*cr*ct*rr**3*tt - 2*cr*ct*rr**2*rt*tt + 2*cr*ct*rr*rt**2*tr + 2*cr*ct*rt**3*tr + 4*cr*ct*rt**3*tt - ct**2*rr**3*tr - ct**2*rr**3*tt - ct**2*rr**2*rt*tr - 3*ct**2*rr**2*rt*tt - 6*ct**2*rr*rt**2*tt + 2*ct**2*rt**3*tr + 10*ct**2*rt**3*tt)/120
    S002 = Nb*u0*(-10*cr**2*rr**3*tr - 2*cr**2*rr**3*tt + 6*cr**2*rr**2*rt*tr + 3*cr**2*rr*rt**2*tr + cr**2*rr*rt**2*tt + cr**2*rt**3*tr + cr**2*rt**3*tt - 4*cr*ct*rr**3*tr - 2*cr*ct*rr**3*tt - 2*cr*ct*rr**2*rt*tt + 2*cr*ct*rr*rt**2*tr + 2*cr*ct*rt**3*tr + 4*cr*ct*rt**3*tt - ct**2*rr**3*tr - ct**2*rr**3*tt - ct**2*rr**2*rt*tr - 3*ct**2*rr**2*rt*tt - 6*ct**2*rr*rt**2*tt + 2*ct**2*rt**3*tr + 10*ct**2*rt**3*tt)/120
    
    
    
    a = -tr*(10*cr**2 + 4*cr*ct + ct**2) - tt*(2*cr**2 + 2*cr*ct + ct**2)
    b = + tr*(6*cr**2 - ct**2) - tt*(2*cr*ct + 3*ct**2)
    d = + tr*(3*cr**2 + 2*cr*ct) + tt*(cr**2 - 6*ct**2)
    f = + tr*(cr**2 + 2*cr*ct + 2*ct**2) + tt*(cr**2 + 4*cr*ct + 10*ct**2)
    S002 = Nb*u0*( rr**3 * a + rr**2*rt * b + rr*rt**2 * d + rt**3 * f )/120
    
