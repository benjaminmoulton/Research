import sympy as sy


if __name__ == "__main__":
    sym = sy.Symbol
    igr = sy.integrate
    simp = sy.simplify
    exp = sy.expand
    piecewise = sy.Piecewise
    pi = sy.pi

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
    # h = sym("h") # cylinder

    # define bounds
    r_up = rt
    r_lo = rr
    # r_up = sym("R2") # cylinder
    # r_lo = sym("R1") # cylinder
    theta_up = 2 * pi
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

    ka = 3*cr**2*tr + cr**2*tt + 2*cr*ct*tr + 2*cr*ct*tt + ct**2*tr + 3*ct**2*tt

    S_ints = []

    for q in range(len(S)):
        f = i**S[q][0] * j**S[q][1] * k**S[q][2]
        if q in [7,8,9]:
            iii = 10
        else:
            iii = simp( igr(f * r, x_bnd,r_bnd,theta_bnd) )
        S_ints.append(iii)
        print("S{}{}{} = {}".format(S[q][0],S[q][1],S[q][2],iii))
        # if i == 0:
        #     print("or   = {}".format(simp(exp(Nb*(rt-rr)/12*ka*u0))))

    print()
    print()
    m = sym("m")
    V = simp( S_ints[0] )
    print("V =",V)
    xcg = simp(S_ints[1] / S_ints[0])
    print("xcg =",xcg)
    ycg = simp(S_ints[2] / S_ints[0])
    print("ycg =",ycg)
    zcg = simp(S_ints[3] / S_ints[0])
    print("zcg =",zcg)
    Ixx = simp(m * (S_ints[8] + S_ints[9]) / S_ints[0])
    print("Ixx =",Ixx)
    Iyy = simp(m * (S_ints[7] + S_ints[9]) / S_ints[0])
    print("Iyy =",Iyy)
    Izz = simp(m * (S_ints[7] + S_ints[8]) / S_ints[0])
    print("Izz =",Izz)
    Ixy = simp(m * S_ints[4] / S_ints[0])
    print("Ixy =",Ixy)
    Ixz = simp(m * S_ints[5] / S_ints[0])
    print("Ixz =",Ixz)
    Iyz = simp(m * S_ints[6] / S_ints[0])
    print("Iyz =",Iyz)

    # # # code results
    # S000 = Nb*u0*(-3*cr**2*rr*tr - cr**2*rr*tt + 3*cr**2*rt*tr + cr**2*rt*tt - 2*cr*ct*rr*tr - 2*cr*ct*rr*tt + 2*cr*ct*rt*tr + 2*cr*ct*rt*tt - ct**2*rr*tr - 3*ct**2*rr*tt + ct**2*rt*tr + 3*ct**2*rt*tt)/12
    # S100 = 0
    # S010 = 0
    # S001 = 0
    # S110 = 0
    # S101 = 0
    # S011 = 0
    # S200 = Nb**3*u0**3*(-35*cr**6*rr**9*rt*tr**3 - 15*cr**6*rr**9*rt*tr**2*tt - 5*cr**6*rr**9*rt*tr*tt**2 - cr**6*rr**9*rt*tt**3 + 360*cr**6*rr**8*rt**2*tr**3 + 160*cr**6*rr**8*rt**2*tr**2*tt + 56*cr**6*rr**8*rt**2*tr*tt**2 + 12*cr**6*rr**8*rt**2*tt**3 - 1680*cr**6*rr**7*rt**3*tr**3 - 784*cr**6*rr**7*rt**3*tr**2*tt - 294*cr**6*rr**7*rt**3*tr*tt**2 - 70*cr**6*rr**7*rt**3*tt**3 + 4704*cr**6*rr**6*rt**4*tr**3 + 2352*cr**6*rr**6*rt**4*tr**2*tt + 980*cr**6*rr**6*rt**4*tr*tt**2 + 280*cr**6*rr**6*rt**4*tt**3 - 8820*cr**6*rr**5*rt**5*tr**3 - 4900*cr**6*rr**5*rt**5*tr**2*tt - 2450*cr**6*rr**5*rt**5*tr*tt**2 - 1050*cr**6*rr**5*rt**5*tt**3 + 11760*cr**6*rr**4*rt**6*tr**3 + 7840*cr**6*rr**4*rt**6*tr**2*tt + 5880*cr**6*rr**4*rt**6*tr*tt**2 + 1680*cr**6*rr**4*rt**6*tt**3*log(rr) - 1680*cr**6*rr**4*rt**6*tt**3*log(rt) - 924*cr**6*rr**4*rt**6*tt**3 - 11760*cr**6*rr**3*rt**7*tr**3 - 11760*cr**6*rr**3*rt**7*tr**2*tt - 5880*cr**6*rr**3*rt**7*tr*tt**2*log(rr) + 5880*cr**6*rr**3*rt**7*tr*tt**2*log(rt) + 294*cr**6*rr**3*rt**7*tr*tt**2 + 840*cr**6*rr**3*rt**7*tt**3*log(rr) - 840*cr**6*rr**3*rt**7*tt**3*log(rt) + 1638*cr**6*rr**3*rt**7*tt**3 + 10080*cr**6*rr**2*rt**8*tr**3 + 6720*cr**6*rr**2*rt**8*tr**2*tt*log(rr) - 6720*cr**6*rr**2*rt**8*tr**2*tt*log(rt) + 3984*cr**6*rr**2*rt**8*tr**2*tt - 1680*cr**6*rr**2*rt**8*tr*tt**2*log(rr) + 1680*cr**6*rr**2*rt**8*tr*tt**2*log(rt) - 4356*cr**6*rr**2*rt**8*tr*tt**2 + 120*cr**6*rr**2*rt**8*tt**3 - 2520*cr**6*rr*rt**9*tr**3*log(rr) + 2520*cr**6*rr*rt**9*tr**3*log(rt) - 4329*cr**6*rr*rt**9*tr**3 + 840*cr**6*rr*rt**9*tr**2*tt*log(rr) - 840*cr**6*rr*rt**9*tr**2*tt*log(rt) + 3123*cr**6*rr*rt**9*tr**2*tt - 105*cr**6*rr*rt**9*tr*tt**2 - 5*cr**6*rr*rt**9*tt**3 - 280*cr**6*rt**10*tr**3 - 30*cr**5*ct*rr**9*rt*tr**3 - 30*cr**5*ct*rr**9*rt*tr**2*tt - 18*cr**5*ct*rr**9*rt*tr*tt**2 - 6*cr**5*ct*rr**9*rt*tt**3 + 320*cr**5*ct*rr**8*rt**2*tr**3 + 336*cr**5*ct*rr**8*rt**2*tr**2*tt + 216*cr**5*ct*rr**8*rt**2*tr*tt**2 + 80*cr**5*ct*rr**8*rt**2*tt**3 - 1568*cr**5*ct*rr**7*rt**3*tr**3 - 1764*cr**5*ct*rr**7*rt**3*tr**2*tt - 1260*cr**5*ct*rr**7*rt**3*tr*tt**2 - 560*cr**5*ct*rr**7*rt**3*tt**3 + 4704*cr**5*ct*rr**6*rt**4*tr**3 + 5880*cr**5*ct*rr**6*rt**4*tr**2*tt + 5040*cr**5*ct*rr**6*rt**4*tr*tt**2 + 3360*cr**5*ct*rr**6*rt**4*tt**3 - 9800*cr**5*ct*rr**5*rt**5*tr**3 - 14700*cr**5*ct*rr**5*rt**5*tr**2*tt - 18900*cr**5*ct*rr**5*rt**5*tr*tt**2 - 8400*cr**5*ct*rr**5*rt**5*tt**3*log(rr) + 8400*cr**5*ct*rr**5*rt**5*tt**3*log(rt) + 8400*cr**5*ct*rr**5*rt**5*tt**3 + 15680*cr**5*ct*rr**4*rt**6*tr**3 + 35280*cr**5*ct*rr**4*rt**6*tr**2*tt + 30240*cr**5*ct*rr**4*rt**6*tr*tt**2*log(rr) - 30240*cr**5*ct*rr**4*rt**6*tr*tt**2*log(rt) - 16632*cr**5*ct*rr**4*rt**6*tr*tt**2 - 6720*cr**5*ct*rr**4*rt**6*tt**3*log(rr) + 6720*cr**5*ct*rr**4*rt**6*tt**3*log(rt) - 9744*cr**5*ct*rr**4*rt**6*tt**3 - 23520*cr**5*ct*rr**3*rt**7*tr**3 - 35280*cr**5*ct*rr**3*rt**7*tr**2*tt*log(rr) + 35280*cr**5*ct*rr**3*rt**7*tr**2*tt*log(rt) + 1764*cr**5*ct*rr**3*rt**7*tr**2*tt + 15120*cr**5*ct*rr**3*rt**7*tr*tt**2*log(rr) - 15120*cr**5*ct*rr**3*rt**7*tr*tt**2*log(rt) + 29484*cr**5*ct*rr**3*rt**7*tr*tt**2 - 1680*cr**5*ct*rr**3*rt**7*tt**3 + 13440*cr**5*ct*rr**2*rt**8*tr**3*log(rr) - 13440*cr**5*ct*rr**2*rt**8*tr**3*log(rt) + 7968*cr**5*ct*rr**2*rt**8*tr**3 - 10080*cr**5*ct*rr**2*rt**8*tr**2*tt*log(rr) + 10080*cr**5*ct*rr**2*rt**8*tr**2*tt*log(rt) - 26136*cr**5*ct*rr**2*rt**8*tr**2*tt + 2160*cr**5*ct*rr**2*rt**8*tr*tt**2 + 160*cr**5*ct*rr**2*rt**8*tt**3 + 1680*cr**5*ct*rr*rt**9*tr**3*log(rr) - 1680*cr**5*ct*rr*rt**9*tr**3*log(rt) + 6246*cr**5*ct*rr*rt**9*tr**3 - 630*cr**5*ct*rr*rt**9*tr**2*tt - 90*cr**5*ct*rr*rt**9*tr*tt**2 - 10*cr**5*ct*rr*rt**9*tt**3 - 25*cr**4*ct**2*rr**9*rt*tr**3 - 45*cr**4*ct**2*rr**9*rt*tr**2*tt - 45*cr**4*ct**2*rr**9*rt*tr*tt**2 - 25*cr**4*ct**2*rr**9*rt*tt**3 + 280*cr**4*ct**2*rr**8*rt**2*tr**3 + 540*cr**4*ct**2*rr**8*rt**2*tr**2*tt + 600*cr**4*ct**2*rr**8*rt**2*tr*tt**2 + 400*cr**4*ct**2*rr**8*rt**2*tt**3 - 1470*cr**4*ct**2*rr**7*rt**3*tr**3 - 3150*cr**4*ct**2*rr**7*rt**3*tr**2*tt - 4200*cr**4*ct**2*rr**7*rt**3*tr*tt**2 - 4200*cr**4*ct**2*rr**7*rt**3*tt**3 + 4900*cr**4*ct**2*rr**6*rt**4*tr**3 + 12600*cr**4*ct**2*rr**6*rt**4*tr**2*tt + 25200*cr**4*ct**2*rr**6*rt**4*tr*tt**2 + 16800*cr**4*ct**2*rr**6*rt**4*tt**3*log(rr) - 16800*cr**4*ct**2*rr**6*rt**4*tt**3*log(rt) - 24360*cr**4*ct**2*rr**6*rt**4*tt**3 - 12250*cr**4*ct**2*rr**5*rt**5*tr**3 - 47250*cr**4*ct**2*rr**5*rt**5*tr**2*tt - 63000*cr**4*ct**2*rr**5*rt**5*tr*tt**2*log(rr) + 63000*cr**4*ct**2*rr**5*rt**5*tr*tt**2*log(rt) + 63000*cr**4*ct**2*rr**5*rt**5*tr*tt**2 + 21000*cr**4*ct**2*rr**5*rt**5*tt**3*log(rr) - 21000*cr**4*ct**2*rr**5*rt**5*tt**3*log(rt) + 21000*cr**4*ct**2*rr**5*rt**5*tt**3 + 29400*cr**4*ct**2*rr**4*rt**6*tr**3 + 75600*cr**4*ct**2*rr**4*rt**6*tr**2*tt*log(rr) - 75600*cr**4*ct**2*rr**4*rt**6*tr**2*tt*log(rt) - 41580*cr**4*ct**2*rr**4*rt**6*tr**2*tt - 50400*cr**4*ct**2*rr**4*rt**6*tr*tt**2*log(rr) + 50400*cr**4*ct**2*rr**4*rt**6*tr*tt**2*log(rt) - 73080*cr**4*ct**2*rr**4*rt**6*tr*tt**2 + 8400*cr**4*ct**2*rr**4*rt**6*tt**3 - 29400*cr**4*ct**2*rr**3*rt**7*tr**3*log(rr) + 29400*cr**4*ct**2*rr**3*rt**7*tr**3*log(rt) + 1470*cr**4*ct**2*rr**3*rt**7*tr**3 + 37800*cr**4*ct**2*rr**3*rt**7*tr**2*tt*log(rr) - 37800*cr**4*ct**2*rr**3*rt**7*tr**2*tt*log(rt) + 73710*cr**4*ct**2*rr**3*rt**7*tr**2*tt - 12600*cr**4*ct**2*rr**3*rt**7*tr*tt**2 - 1400*cr**4*ct**2*rr**3*rt**7*tt**3 - 8400*cr**4*ct**2*rr**2*rt**8*tr**3*log(rr) + 8400*cr**4*ct**2*rr**2*rt**8*tr**3*log(rt) - 21780*cr**4*ct**2*rr**2*rt**8*tr**3 + 5400*cr**4*ct**2*rr**2*rt**8*tr**2*tt + 1200*cr**4*ct**2*rr**2*rt**8*tr*tt**2 + 200*cr**4*ct**2*rr**2*rt**8*tt**3 - 525*cr**4*ct**2*rr*rt**9*tr**3 - 225*cr**4*ct**2*rr*rt**9*tr**2*tt - 75*cr**4*ct**2*rr*rt**9*tr*tt**2 - 15*cr**4*ct**2*rr*rt**9*tt**3 - 20*cr**3*ct**3*rr**9*rt*tr**3 - 60*cr**3*ct**3*rr**9*rt*tr**2*tt - 100*cr**3*ct**3*rr**9*rt*tr*tt**2 - 100*cr**3*ct**3*rr**9*rt*tt**3 + 240*cr**3*ct**3*rr**8*rt**2*tr**3 + 800*cr**3*ct**3*rr**8*rt**2*tr**2*tt + 1600*cr**3*ct**3*rr**8*rt**2*tr*tt**2 + 2400*cr**3*ct**3*rr**8*rt**2*tt**3 - 1400*cr**3*ct**3*rr**7*rt**3*tr**3 - 5600*cr**3*ct**3*rr**7*rt**3*tr**2*tt - 16800*cr**3*ct**3*rr**7*rt**3*tr*tt**2 - 16800*cr**3*ct**3*rr**7*rt**3*tt**3*log(rr) + 16800*cr**3*ct**3*rr**7*rt**3*tt**3*log(rt) + 32760*cr**3*ct**3*rr**7*rt**3*tt**3 + 5600*cr**3*ct**3*rr**6*rt**4*tr**3 + 33600*cr**3*ct**3*rr**6*rt**4*tr**2*tt + 67200*cr**3*ct**3*rr**6*rt**4*tr*tt**2*log(rr) - 67200*cr**3*ct**3*rr**6*rt**4*tr*tt**2*log(rt) - 97440*cr**3*ct**3*rr**6*rt**4*tr*tt**2 - 33600*cr**3*ct**3*rr**6*rt**4*tt**3*log(rr) + 33600*cr**3*ct**3*rr**6*rt**4*tt**3*log(rt) - 18480*cr**3*ct**3*rr**6*rt**4*tt**3 - 21000*cr**3*ct**3*rr**5*rt**5*tr**3 - 84000*cr**3*ct**3*rr**5*rt**5*tr**2*tt*log(rr) + 84000*cr**3*ct**3*rr**5*rt**5*tr**2*tt*log(rt) + 84000*cr**3*ct**3*rr**5*rt**5*tr**2*tt + 84000*cr**3*ct**3*rr**5*rt**5*tr*tt**2*log(rr) - 84000*cr**3*ct**3*rr**5*rt**5*tr*tt**2*log(rt) + 84000*cr**3*ct**3*rr**5*rt**5*tr*tt**2 - 21000*cr**3*ct**3*rr**5*rt**5*tt**3 + 33600*cr**3*ct**3*rr**4*rt**6*tr**3*log(rr) - 33600*cr**3*ct**3*rr**4*rt**6*tr**3*log(rt) - 18480*cr**3*ct**3*rr**4*rt**6*tr**3 - 67200*cr**3*ct**3*rr**4*rt**6*tr**2*tt*log(rr) + 67200*cr**3*ct**3*rr**4*rt**6*tr**2*tt*log(rt) - 97440*cr**3*ct**3*rr**4*rt**6*tr**2*tt + 33600*cr**3*ct**3*rr**4*rt**6*tr*tt**2 + 5600*cr**3*ct**3*rr**4*rt**6*tt**3 + 16800*cr**3*ct**3*rr**3*rt**7*tr**3*log(rr) - 16800*cr**3*ct**3*rr**3*rt**7*tr**3*log(rt) + 32760*cr**3*ct**3*rr**3*rt**7*tr**3 - 16800*cr**3*ct**3*rr**3*rt**7*tr**2*tt - 5600*cr**3*ct**3*rr**3*rt**7*tr*tt**2 - 1400*cr**3*ct**3*rr**3*rt**7*tt**3 + 2400*cr**3*ct**3*rr**2*rt**8*tr**3 + 1600*cr**3*ct**3*rr**2*rt**8*tr**2*tt + 800*cr**3*ct**3*rr**2*rt**8*tr*tt**2 + 240*cr**3*ct**3*rr**2*rt**8*tt**3 - 100*cr**3*ct**3*rr*rt**9*tr**3 - 100*cr**3*ct**3*rr*rt**9*tr**2*tt - 60*cr**3*ct**3*rr*rt**9*tr*tt**2 - 20*cr**3*ct**3*rr*rt**9*tt**3 - 15*cr**2*ct**4*rr**9*rt*tr**3 - 75*cr**2*ct**4*rr**9*rt*tr**2*tt - 225*cr**2*ct**4*rr**9*rt*tr*tt**2 - 525*cr**2*ct**4*rr**9*rt*tt**3 + 200*cr**2*ct**4*rr**8*rt**2*tr**3 + 1200*cr**2*ct**4*rr**8*rt**2*tr**2*tt + 5400*cr**2*ct**4*rr**8*rt**2*tr*tt**2 + 8400*cr**2*ct**4*rr**8*rt**2*tt**3*log(rr) - 8400*cr**2*ct**4*rr**8*rt**2*tt**3*log(rt) - 21780*cr**2*ct**4*rr**8*rt**2*tt**3 - 1400*cr**2*ct**4*rr**7*rt**3*tr**3 - 12600*cr**2*ct**4*rr**7*rt**3*tr**2*tt - 37800*cr**2*ct**4*rr**7*rt**3*tr*tt**2*log(rr) + 37800*cr**2*ct**4*rr**7*rt**3*tr*tt**2*log(rt) + 73710*cr**2*ct**4*rr**7*rt**3*tr*tt**2 + 29400*cr**2*ct**4*rr**7*rt**3*tt**3*log(rr) - 29400*cr**2*ct**4*rr**7*rt**3*tt**3*log(rt) + 1470*cr**2*ct**4*rr**7*rt**3*tt**3 + 8400*cr**2*ct**4*rr**6*rt**4*tr**3 + 50400*cr**2*ct**4*rr**6*rt**4*tr**2*tt*log(rr) - 50400*cr**2*ct**4*rr**6*rt**4*tr**2*tt*log(rt) - 73080*cr**2*ct**4*rr**6*rt**4*tr**2*tt - 75600*cr**2*ct**4*rr**6*rt**4*tr*tt**2*log(rr) + 75600*cr**2*ct**4*rr**6*rt**4*tr*tt**2*log(rt) - 41580*cr**2*ct**4*rr**6*rt**4*tr*tt**2 + 29400*cr**2*ct**4*rr**6*rt**4*tt**3 - 21000*cr**2*ct**4*rr**5*rt**5*tr**3*log(rr) + 21000*cr**2*ct**4*rr**5*rt**5*tr**3*log(rt) + 21000*cr**2*ct**4*rr**5*rt**5*tr**3 + 63000*cr**2*ct**4*rr**5*rt**5*tr**2*tt*log(rr) - 63000*cr**2*ct**4*rr**5*rt**5*tr**2*tt*log(rt) + 63000*cr**2*ct**4*rr**5*rt**5*tr**2*tt - 47250*cr**2*ct**4*rr**5*rt**5*tr*tt**2 - 12250*cr**2*ct**4*rr**5*rt**5*tt**3 - 16800*cr**2*ct**4*rr**4*rt**6*tr**3*log(rr) + 16800*cr**2*ct**4*rr**4*rt**6*tr**3*log(rt) - 24360*cr**2*ct**4*rr**4*rt**6*tr**3 + 25200*cr**2*ct**4*rr**4*rt**6*tr**2*tt + 12600*cr**2*ct**4*rr**4*rt**6*tr*tt**2 + 4900*cr**2*ct**4*rr**4*rt**6*tt**3 - 4200*cr**2*ct**4*rr**3*rt**7*tr**3 - 4200*cr**2*ct**4*rr**3*rt**7*tr**2*tt - 3150*cr**2*ct**4*rr**3*rt**7*tr*tt**2 - 1470*cr**2*ct**4*rr**3*rt**7*tt**3 + 400*cr**2*ct**4*rr**2*rt**8*tr**3 + 600*cr**2*ct**4*rr**2*rt**8*tr**2*tt + 540*cr**2*ct**4*rr**2*rt**8*tr*tt**2 + 280*cr**2*ct**4*rr**2*rt**8*tt**3 - 25*cr**2*ct**4*rr*rt**9*tr**3 - 45*cr**2*ct**4*rr*rt**9*tr**2*tt - 45*cr**2*ct**4*rr*rt**9*tr*tt**2 - 25*cr**2*ct**4*rr*rt**9*tt**3 - 10*cr*ct**5*rr**9*rt*tr**3 - 90*cr*ct**5*rr**9*rt*tr**2*tt - 630*cr*ct**5*rr**9*rt*tr*tt**2 - 1680*cr*ct**5*rr**9*rt*tt**3*log(rr) + 1680*cr*ct**5*rr**9*rt*tt**3*log(rt) + 6246*cr*ct**5*rr**9*rt*tt**3 + 160*cr*ct**5*rr**8*rt**2*tr**3 + 2160*cr*ct**5*rr**8*rt**2*tr**2*tt + 10080*cr*ct**5*rr**8*rt**2*tr*tt**2*log(rr) - 10080*cr*ct**5*rr**8*rt**2*tr*tt**2*log(rt) - 26136*cr*ct**5*rr**8*rt**2*tr*tt**2 - 13440*cr*ct**5*rr**8*rt**2*tt**3*log(rr) + 13440*cr*ct**5*rr**8*rt**2*tt**3*log(rt) + 7968*cr*ct**5*rr**8*rt**2*tt**3 - 1680*cr*ct**5*rr**7*rt**3*tr**3 - 15120*cr*ct**5*rr**7*rt**3*tr**2*tt*log(rr) + 15120*cr*ct**5*rr**7*rt**3*tr**2*tt*log(rt) + 29484*cr*ct**5*rr**7*rt**3*tr**2*tt + 35280*cr*ct**5*rr**7*rt**3*tr*tt**2*log(rr) - 35280*cr*ct**5*rr**7*rt**3*tr*tt**2*log(rt) + 1764*cr*ct**5*rr**7*rt**3*tr*tt**2 - 23520*cr*ct**5*rr**7*rt**3*tt**3 + 6720*cr*ct**5*rr**6*rt**4*tr**3*log(rr) - 6720*cr*ct**5*rr**6*rt**4*tr**3*log(rt) - 9744*cr*ct**5*rr**6*rt**4*tr**3 - 30240*cr*ct**5*rr**6*rt**4*tr**2*tt*log(rr) + 30240*cr*ct**5*rr**6*rt**4*tr**2*tt*log(rt) - 16632*cr*ct**5*rr**6*rt**4*tr**2*tt + 35280*cr*ct**5*rr**6*rt**4*tr*tt**2 + 15680*cr*ct**5*rr**6*rt**4*tt**3 + 8400*cr*ct**5*rr**5*rt**5*tr**3*log(rr) - 8400*cr*ct**5*rr**5*rt**5*tr**3*log(rt) + 8400*cr*ct**5*rr**5*rt**5*tr**3 - 18900*cr*ct**5*rr**5*rt**5*tr**2*tt - 14700*cr*ct**5*rr**5*rt**5*tr*tt**2 - 9800*cr*ct**5*rr**5*rt**5*tt**3 + 3360*cr*ct**5*rr**4*rt**6*tr**3 + 5040*cr*ct**5*rr**4*rt**6*tr**2*tt + 5880*cr*ct**5*rr**4*rt**6*tr*tt**2 + 4704*cr*ct**5*rr**4*rt**6*tt**3 - 560*cr*ct**5*rr**3*rt**7*tr**3 - 1260*cr*ct**5*rr**3*rt**7*tr**2*tt - 1764*cr*ct**5*rr**3*rt**7*tr*tt**2 - 1568*cr*ct**5*rr**3*rt**7*tt**3 + 80*cr*ct**5*rr**2*rt**8*tr**3 + 216*cr*ct**5*rr**2*rt**8*tr**2*tt + 336*cr*ct**5*rr**2*rt**8*tr*tt**2 + 320*cr*ct**5*rr**2*rt**8*tt**3 - 6*cr*ct**5*rr*rt**9*tr**3 - 18*cr*ct**5*rr*rt**9*tr**2*tt - 30*cr*ct**5*rr*rt**9*tr*tt**2 - 30*cr*ct**5*rr*rt**9*tt**3 - 280*ct**6*rr**10*tt**3 - 5*ct**6*rr**9*rt*tr**3 - 105*ct**6*rr**9*rt*tr**2*tt - 840*ct**6*rr**9*rt*tr*tt**2*log(rr) + 840*ct**6*rr**9*rt*tr*tt**2*log(rt) + 3123*ct**6*rr**9*rt*tr*tt**2 + 2520*ct**6*rr**9*rt*tt**3*log(rr) - 2520*ct**6*rr**9*rt*tt**3*log(rt) - 4329*ct**6*rr**9*rt*tt**3 + 120*ct**6*rr**8*rt**2*tr**3 + 1680*ct**6*rr**8*rt**2*tr**2*tt*log(rr) - 1680*ct**6*rr**8*rt**2*tr**2*tt*log(rt) - 4356*ct**6*rr**8*rt**2*tr**2*tt - 6720*ct**6*rr**8*rt**2*tr*tt**2*log(rr) + 6720*ct**6*rr**8*rt**2*tr*tt**2*log(rt) + 3984*ct**6*rr**8*rt**2*tr*tt**2 + 10080*ct**6*rr**8*rt**2*tt**3 - 840*ct**6*rr**7*rt**3*tr**3*log(rr) + 840*ct**6*rr**7*rt**3*tr**3*log(rt) + 1638*ct**6*rr**7*rt**3*tr**3 + 5880*ct**6*rr**7*rt**3*tr**2*tt*log(rr) - 5880*ct**6*rr**7*rt**3*tr**2*tt*log(rt) + 294*ct**6*rr**7*rt**3*tr**2*tt - 11760*ct**6*rr**7*rt**3*tr*tt**2 - 11760*ct**6*rr**7*rt**3*tt**3 - 1680*ct**6*rr**6*rt**4*tr**3*log(rr) + 1680*ct**6*rr**6*rt**4*tr**3*log(rt) - 924*ct**6*rr**6*rt**4*tr**3 + 5880*ct**6*rr**6*rt**4*tr**2*tt + 7840*ct**6*rr**6*rt**4*tr*tt**2 + 11760*ct**6*rr**6*rt**4*tt**3 - 1050*ct**6*rr**5*rt**5*tr**3 - 2450*ct**6*rr**5*rt**5*tr**2*tt - 4900*ct**6*rr**5*rt**5*tr*tt**2 - 8820*ct**6*rr**5*rt**5*tt**3 + 280*ct**6*rr**4*rt**6*tr**3 + 980*ct**6*rr**4*rt**6*tr**2*tt + 2352*ct**6*rr**4*rt**6*tr*tt**2 + 4704*ct**6*rr**4*rt**6*tt**3 - 70*ct**6*rr**3*rt**7*tr**3 - 294*ct**6*rr**3*rt**7*tr**2*tt - 784*ct**6*rr**3*rt**7*tr*tt**2 - 1680*ct**6*rr**3*rt**7*tt**3 + 12*ct**6*rr**2*rt**8*tr**3 + 56*ct**6*rr**2*rt**8*tr**2*tt + 160*ct**6*rr**2*rt**8*tr*tt**2 + 360*ct**6*rr**2*rt**8*tt**3 - ct**6*rr*rt**9*tr**3 - 5*ct**6*rr*rt**9*tr**2*tt - 15*ct**6*rr*rt**9*tr*tt**2 - 35*ct**6*rr*rt**9*tt**3)/(13440*pi**2*rr*rt*(rr**9 - 9*rr**8*rt + 36*rr**7*rt**2 - 84*rr**6*rt**3 + 126*rr**5*rt**4 - 126*rr**4*rt**5 + 84*rr**3*rt**6 - 36*rr**2*rt**7 + 9*rr*rt**8 - rt**9))
    # S020 = Nb*u0*(-10*cr**2*rr**3*tr - 2*cr**2*rr**3*tt + 6*cr**2*rr**2*rt*tr + 3*cr**2*rr*rt**2*tr + cr**2*rr*rt**2*tt + cr**2*rt**3*tr + cr**2*rt**3*tt - 4*cr*ct*rr**3*tr - 2*cr*ct*rr**3*tt - 2*cr*ct*rr**2*rt*tt + 2*cr*ct*rr*rt**2*tr + 2*cr*ct*rt**3*tr + 4*cr*ct*rt**3*tt - ct**2*rr**3*tr - ct**2*rr**3*tt - ct**2*rr**2*rt*tr - 3*ct**2*rr**2*rt*tt - 6*ct**2*rr*rt**2*tt + 2*ct**2*rt**3*tr + 10*ct**2*rt**3*tt)/120
    # S002 = Nb*u0*(-10*cr**2*rr**3*tr - 2*cr**2*rr**3*tt + 6*cr**2*rr**2*rt*tr + 3*cr**2*rr*rt**2*tr + cr**2*rr*rt**2*tt + cr**2*rt**3*tr + cr**2*rt**3*tt - 4*cr*ct*rr**3*tr - 2*cr*ct*rr**3*tt - 2*cr*ct*rr**2*rt*tt + 2*cr*ct*rr*rt**2*tr + 2*cr*ct*rt**3*tr + 4*cr*ct*rt**3*tt - ct**2*rr**3*tr - ct**2*rr**3*tt - ct**2*rr**2*rt*tr - 3*ct**2*rr**2*rt*tt - 6*ct**2*rr*rt**2*tt + 2*ct**2*rt**3*tr + 10*ct**2*rt**3*tt)/120
    
    ka = 3*cr**2*tr + cr**2*tt + 2*cr*ct*tr + 2*cr*ct*tt + ct**2*tr + 3*ct**2*tt
    ka = tr*(3*cr**2 + 2*cr*ct + ct**2) + tt*(cr**2 + 2*cr*ct + 3*ct**2)
    # S000 = Nb*u0*ka*(rt-rr)/12
    # S100 = 0
    # S010 = 0
    # S001 = 0
    # S110 = 0
    # S101 = 0
    # S011 = 0
    ln = sy.ln

    # >>> a = np.array([16800,67200,33600,50400,75600,6720,30240,1680])
    # >>> print(np.gcd.reduce(a),a / np.gcd.reduce(a))
    
    a = -1*(
         + cr**6*(35*tr**3 + 15*tr**2*tt + 5*tr*tt**2 + tt**3)
         + 6*cr**5*ct*(5*tr**3 + 5*tr**2*tt + 3*tr*tt**2 + tt**3)
         + 5*cr**4*ct**2*(5*tr**3 + 9*tr**2*tt + 9*tr*tt**2 + 5*tt**3)
         + 20*cr**3*ct**3*(tr**3 + 3*tr**2*tt + 5*tr*tt**2 + 5*tt**3)
         + 15*cr**2*ct**4*(tr**3 + 5*tr**2*tt + 15*tr*tt**2 + 35*tt**3)
         + 2*cr*ct**5*(5*tr**3 + 45*tr**2*tt + 315*tr*tt**2 - 3123*tt**3)
         + ct**6*(5*tr**3 + 105*tr**2*tt - 3123*tr*tt**2 + 4329*tt**3)
    )
    al= ln(rr/rt)*(
         - 1680*cr*ct**5*tt**3
         - 840*ct**6*(tr*tt**2 - 3*tt**3)
    )
    b = 4*(
         +    cr**6*(90*tr**3 + 40*tr**2*tt + 14*tr*tt**2 + 3*tt**3)
         +  2*cr**5*ct   *(40*tr**3 +  42*tr**2*tt +   27*tr*tt**2 +   10*tt**3)
         +  5*cr**4*ct**2*(14*tr**3 +  27*tr**2*tt +   30*tr*tt**2 +   20*tt**3)
         + 20*cr**3*ct**3*( 3*tr**3 +  10*tr**2*tt +   20*tr*tt**2 +   30*tt**3)
         +  5*cr**2*ct**4*(10*tr**3 +  60*tr**2*tt +  270*tr*tt**2 - 1089*tt**3)
         +  2*cr   *ct**5*(20*tr**3 + 270*tr**2*tt - 3267*tr*tt**2 +  996*tt**3)
         +  3      *ct**6*(10*tr**3 - 363*tr**2*tt +  332*tr*tt**2 +  840*tt**3)
    )
    bl= 1680*ln(rr/rt)*(
         + 5*cr**2*ct**4* tt**3
         + 2*cr   *ct**5*(3*tr*tt**2 - 4*tt**3)
         +         ct**6*(tr**2*tt - 4*tr*tt**2)
    )
    c = -14*(
         + cr**6*(120*tr**3 + 56*tr**2*tt + 21*tr*tt**2 + 5*tt**3)
         + 2*cr**5*ct*(56*tr**3 + 63*tr**2*tt + 45*tr*tt**2 + 20*tt**3)
         + 15*cr**4*ct**2*(7*tr**3 + 15*tr**2*tt + 20*tr*tt**2 + 20*tt**3)
         + 20*cr**3*ct**3*(5*tr**3 + 20*tr**2*tt + 60*tr*tt**2 - 117*tt**3)
         + 5*cr**2*ct**4*(20*tr**3 + 180*tr**2*tt - 1053*tr*tt**2 - 21*tt**3)
         + 6*cr*ct**5*(20*tr**3 - 351*tr**2*tt - 21*tr*tt**2 + 280*tt**3)
         + 3*ct**6(-39*tr**3 - 7*tr**2*tt + 280*tr*tt**2 + 280*tt**3)
    )
    cl= -840*ln(rr/rt)*(
         + 20*cr**3*ct**3*tt**3
         + 5*cr**2*ct**4*(9*tr*tt**2 - 7*tt**3)
         + 6*cr*ct**5*(3*tr**2*tt - 7*tr*tt**2)
         +   ct**6*(tr**3 - 7*tr**2*tt)
    )
    d = 28*(
         +    cr**6      *( 168*tr**3 +  84*tr**2*tt +  35*tr*tt**2 +  10*tt**3)
         +  6*cr**5*ct   *(  28*tr**3 +  35*tr**2*tt +  30*tr*tt**2 +  20*tt**3)
         +  5*cr**4*ct**2*(  35*tr**3 +  90*tr**2*tt + 180*tr*tt**2 - 174*tt**3)
         + 20*cr**3*ct**3*(  10*tr**3 +  60*tr**2*tt - 174*tr*tt**2 -  33*tt**3)
         + 15*cr**2*ct**4*(  20*tr**3 - 174*tr**2*tt -  99*tr*tt**2 +  70*tt**3)
         +  2*cr   *ct**5*(-174*tr**3 - 297*tr**2*tt + 630*tr*tt**2 + 280*tt**3)
         +          ct**6*(- 33*tr**3 + 210*tr**2*tt + 280*tr*tt**2 + 420*tt**3)
    )
    dl= 1680*ln(rr/rt)*(
         + 10*cr**4*ct**2*tt**3
         + 40*cr**3*ct**3*tr*tt**2
         - 20*cr**3*ct**3*tt**3
         + 30*cr**2*ct**4*tr**2*tt
         - 45*cr**2*ct**4*tr*tt**2
         +  4*cr*ct**5*tr**3
         - 18*cr*ct**5*tr**2*tt
         -    ct**6*tr**3
    )
    e = -70*(
         + cr**6*(126*tr**3 + 70*tr**2*tt + 35*tr*tt**2 + 15*tt**3)
         + 10*cr**5*ct*(14*tr**3 + 21*tr**2*tt + 27*tr*tt**2 - 12*tt**3)
         + 25*cr**4*ct**2*(7*tr**3 + 27*tr**2*tt - 36*tr*tt**2 - 12*tt**3)
         + 300*cr**3*ct**3*(tr**3 - 4*tr**2*tt - 4*tr*tt**2 + tt**3)
         - 25*cr**2*ct**4*(12*tr**3 + 36*tr**2*tt - 27*tr*tt**2 - 7*tt**3)
         - 10*cr*ct**5*(12*tr**3 - 27*tr**2*tt - 21*tr*tt**2 - 14*tt**3)
         + ct**6*(15*tr**3 + 35*tr**2*tt + 70*tr*tt**2 + 126*tt**3)
    )
    el= 4200*ln(rr/rt)* (
         - 2*cr**5*ct*tt**3
         - 15*cr**4*ct**2*tr*tt**2
         + 5*cr**4*ct**2*tt**3
         - 20*cr**3*ct**3*tr**2*tt
         + 20*cr**3*ct**3*tr*tt**2
         - 5*cr**2*ct**4*tr**3
         + 15*cr**2*ct**4*tr**2*tt
         + 2*cr*ct**5*tr**3
    )
    f = (+ 11760*cr**6*tr**3 + 7840*cr**6*tr**2*tt + 5880*cr**6*tr*tt**2 - 924*cr**6*tt**3
         + 15680*cr**5*ct*tr**3 + 35280*cr**5*ct*tr**2*tt - 16632*cr**5*ct*tr*tt**2 - 9744*cr**5*ct*tt**3
         + 29400*cr**4*ct**2*tr**3 - 41580*cr**4*ct**2*tr**2*tt - 73080*cr**4*ct**2*tr*tt**2 + 8400*cr**4*ct**2*tt**3
         - 18480*cr**3*ct**3*tr**3 - 97440*cr**3*ct**3*tr**2*tt + 33600*cr**3*ct**3*tr*tt**2 + 5600*cr**3*ct**3*tt**3
         - 24360*cr**2*ct**4*tr**3 + 25200*cr**2*ct**4*tr**2*tt + 12600*cr**2*ct**4*tr*tt**2 + 4900*cr**2*ct**4*tt**3
         + 3360*cr*ct**5*tr**3 + 5040*cr*ct**5*tr**2*tt + 5880*cr*ct**5*tr*tt**2 + 4704*cr*ct**5*tt**3
         + 280*ct**6*tr**3 + 980*ct**6*tr**2*tt + 2352*ct**6*tr*tt**2 + 4704*ct**6*tt**3
    )
    fl= ln(rr/rt)*(
         + 1680*cr**6*tt**3
         + 30240*cr**5*ct*tr*tt**2
         - 6720*cr**5*ct*tt**3
         + 75600*cr**4*ct**2*tr**2*tt
         - 50400*cr**4*ct**2*tr*tt**2
         + 33600*cr**3*ct**3*tr**3
         - 67200*cr**3*ct**3*tr**2*tt
         - 16800*cr**2*ct**4*tr**3
    )
    g = (- 11760*cr**6*tr**3 - 11760*cr**6*tr**2*tt + 294*cr**6*tr*tt**2 + 1638*cr**6*tt**3
         - 23520*cr**5*ct*tr**3 + 1764*cr**5*ct*tr**2*tt + 29484*cr**5*ct*tr*tt**2 - 1680*cr**5*ct*tt**3
         + 1470*cr**4*ct**2*tr**3 + 73710*cr**4*ct**2*tr**2*tt - 12600*cr**4*ct**2*tr*tt**2 - 1400*cr**4*ct**2*tt**3
         + 32760*cr**3*ct**3*tr**3 - 16800*cr**3*ct**3*tr**2*tt - 5600*cr**3*ct**3*tr*tt**2 - 1400*cr**3*ct**3*tt**3
         - 4200*cr**2*ct**4*tr**3 - 4200*cr**2*ct**4*tr**2*tt - 3150*cr**2*ct**4*tr*tt**2 - 1470*cr**2*ct**4*tt**3
         - 560*cr*ct**5*tr**3 - 1260*cr*ct**5*tr**2*tt - 1764*cr*ct**5*tr*tt**2 - 1568*cr*ct**5*tt**3
         - 70*ct**6*tr**3 - 294*ct**6*tr**2*tt - 784*ct**6*tr*tt**2 - 1680*ct**6*tt**3
    )
    gl= ln(rr/rt)*(
         - 5880*cr**6*tr*tt**2
         + 840*cr**6*tt**3
         - 35280*cr**5*ct*tr**2*tt
         + 15120*cr**5*ct*tr*tt**2
         - 29400*cr**4*ct**2*tr**3
         + 37800*cr**4*ct**2*tr**2*tt
         + 16800*cr**3*ct**3*tr**3
    )
    h = (+ 10080*cr**6*tr**3 + 3984*cr**6*tr**2*tt - 4356*cr**6*tr*tt**2 + 120*cr**6*tt**3
         + 7968*cr**5*ct*tr**3 - 26136*cr**5*ct*tr**2*tt + 2160*cr**5*ct*tr*tt**2 + 160*cr**5*ct*tt**3
         - 21780*cr**4*ct**2*tr**3 + 5400*cr**4*ct**2*tr**2*tt + 1200*cr**4*ct**2*tr*tt**2 + 200*cr**4*ct**2*tt**3
         + 2400*cr**3*ct**3*tr**3 + 1600*cr**3*ct**3*tr**2*tt + 800*cr**3*ct**3*tr*tt**2 + 240*cr**3*ct**3*tt**3
         + 400*cr**2*ct**4*tr**3 + 600*cr**2*ct**4*tr**2*tt + 540*cr**2*ct**4*tr*tt**2 + 280*cr**2*ct**4*tt**3
         + 80*cr*ct**5*tr**3 + 216*cr*ct**5*tr**2*tt + 336*cr*ct**5*tr*tt**2 + 320*cr*ct**5*tt**3
         + 12*ct**6*tr**3 + 56*ct**6*tr**2*tt + 160*ct**6*tr*tt**2 + 360*ct**6*tt**3
    )
    hl= 1680*ln(rr/rt)*(
         + 4*cr**6*tr**2*tt
         - cr**6*tr*tt**2
         + 8*cr**5*ct*tr**3
         - 6*cr**5*ct*tr**2*tt
         - 5*cr**4*ct**2*tr**3
    )
    i = (- 4329*cr**6*tr**3 + 3123*cr**6*tr**2*tt - 105*cr**6*tr*tt**2 - 5*cr**6*tt**3
         + 6246*cr**5*ct*tr**3 - 630*cr**5*ct*tr**2*tt - 90*cr**5*ct*tr*tt**2 - 10*cr**5*ct*tt**3
         - 525*cr**4*ct**2*tr**3 - 225*cr**4*ct**2*tr**2*tt - 75*cr**4*ct**2*tr*tt**2 - 15*cr**4*ct**2*tt**3
         - 100*cr**3*ct**3*tr**3 - 100*cr**3*ct**3*tr**2*tt - 60*cr**3*ct**3*tr*tt**2 - 20*cr**3*ct**3*tt**3
         - 25*cr**2*ct**4*tr**3 - 45*cr**2*ct**4*tr**2*tt - 45*cr**2*ct**4*tr*tt**2 - 25*cr**2*ct**4*tt**3
         - 6*cr*ct**5*tr**3 - 18*cr*ct**5*tr**2*tt - 30*cr*ct**5*tr*tt**2 - 30*cr*ct**5*tt**3
         - ct**6*tr**3 - 5*ct**6*tr**2*tt - 15*ct**6*tr*tt**2 - 35*ct**6*tt**3
    )
    il= 840*ln(rr/rt)*(
         - 3*cr**6*tr**3
         + cr**6*tr**2*tt
         + 2*cr**5*ct*tr**3
    )
    j = - 280*cr**6*tr**3
    k = - 280*ct**6*tt**3

    S200 = Nb**3*u0**3*( 
         + (a+al)*rr**9*rt    + (b+bl)*rr**8*rt**2
         + (c+cl)*rr**7*rt**3 + (d+dl)*rr**6*rt**4
         + (e+el)*rr**5*rt**5 + (f+fl)*rr**4*rt**6
         + (g+gl)*rr**3*rt**7 + (h+hl)*rr**2*rt**8
         + (i+il)*rr   *rt**9 + j*rt**10 + k*rr**10
        )/(13440*pi**2*rr*rt*(rr-rt)**9)
    
    # S020 = Nb*u0*(-10*cr**2*rr**3*tr - 2*cr**2*rr**3*tt + 6*cr**2*rr**2*rt*tr + 3*cr**2*rr*rt**2*tr + cr**2*rr*rt**2*tt + cr**2*rt**3*tr + cr**2*rt**3*tt - 4*cr*ct*rr**3*tr - 2*cr*ct*rr**3*tt - 2*cr*ct*rr**2*rt*tt + 2*cr*ct*rr*rt**2*tr + 2*cr*ct*rt**3*tr + 4*cr*ct*rt**3*tt - ct**2*rr**3*tr - ct**2*rr**3*tt - ct**2*rr**2*rt*tr - 3*ct**2*rr**2*rt*tt - 6*ct**2*rr*rt**2*tt + 2*ct**2*rt**3*tr + 10*ct**2*rt**3*tt)/120
    # S002 = Nb*u0*(-10*cr**2*rr**3*tr - 2*cr**2*rr**3*tt + 6*cr**2*rr**2*rt*tr + 3*cr**2*rr*rt**2*tr + cr**2*rr*rt**2*tt + cr**2*rt**3*tr + cr**2*rt**3*tt - 4*cr*ct*rr**3*tr - 2*cr*ct*rr**3*tt - 2*cr*ct*rr**2*rt*tt + 2*cr*ct*rr*rt**2*tr + 2*cr*ct*rt**3*tr + 4*cr*ct*rt**3*tt - ct**2*rr**3*tr - ct**2*rr**3*tt - ct**2*rr**2*rt*tr - 3*ct**2*rr**2*rt*tt - 6*ct**2*rr*rt**2*tt + 2*ct**2*rt**3*tr + 10*ct**2*rt**3*tt)/120
    
    

    a = -tr*(10*cr**2 + 4*cr*ct + ct**2) - tt*(2*cr**2 + 2*cr*ct + ct**2)
    b = + tr*(6*cr**2 - ct**2) - tt*(2*cr*ct + 3*ct**2)
    d = + tr*(3*cr**2 + 2*cr*ct) + tt*(cr**2 - 6*ct**2)
    f = + tr*(cr**2 + 2*cr*ct + 2*ct**2) + tt*(cr**2 + 4*cr*ct + 10*ct**2)
    S002 = Nb*u0*( rr**3 * a + rr**2*rt * b + rr*rt**2 * d + rt**3 * f )/120