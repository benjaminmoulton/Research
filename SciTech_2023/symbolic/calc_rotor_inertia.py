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

    T = [
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

    T_ints = []

    for q in range(len(T)):
        f = i**T[q][0] * j**T[q][1] * k**T[q][2]
        if q in [100]:#[7,8,9]:
            iii = 10
        else:
            iii = simp( igr(f * r, x_bnd,r_bnd,theta_bnd) )
        T_ints.append(iii)
        print("T{}{}{} = {}".format(T[q][0],T[q][1],T[q][2],iii))
        # if i == 0:
        #     print("or   = {}".format(simp(exp(Nb*(rt-rr)/12*ka*u0))))

    print()
    print()
#     m = sym("m")
#     V = simp( T_ints[0] )
#     print("V =",V)
#     xcg = simp(T_ints[1] / T_ints[0])
#     print("xcg =",xcg)
#     ycg = simp(T_ints[2] / T_ints[0])
#     print("ycg =",ycg)
#     zcg = simp(T_ints[3] / T_ints[0])
#     print("zcg =",zcg)
#     Ixx = simp(m * (T_ints[8] + T_ints[9]) / T_ints[0])
#     print("Ixx =",Ixx)
#     Iyy = simp(m * (T_ints[7] + T_ints[9]) / T_ints[0])
#     print("Iyy =",Iyy)
#     Izz = simp(m * (T_ints[7] + T_ints[8]) / T_ints[0])
#     print("Izz =",Izz)
#     Ixy = simp(m * T_ints[4] / T_ints[0])
#     print("Ixy =",Ixy)
#     Ixz = simp(m * T_ints[5] / T_ints[0])
#     print("Ixz =",Ixz)
#     Iyz = simp(m * T_ints[6] / T_ints[0])
#     print("Iyz =",Iyz)

    
    
    ka = 3*cr**2*tr + cr**2*tt + 2*cr*ct*tr + 2*cr*ct*tt + ct**2*tr + 3*ct**2*tt
    ka = tr*(3*cr**2 + 2*cr*ct + ct**2) + tt*(cr**2 + 2*cr*ct + 3*ct**2)
    # T000 = Nb*u0*ka*(rt-rr)/12
    # T100 = 0
    # T010 = 0
    # T001 = 0
    # T110 = 0
    # T101 = 0
    # T011 = 0
    ln = sy.ln

    # >>> a = np.array([16800,67200,33600,50400,75600,6720,30240,1680])
    # >>> print(np.gcd.reduce(a),a / np.gcd.reduce(a))
    
    a = -1*(
         +    cr**6      *(35*tr**3 +  15*tr**2*tt +    5*tr*tt**2 +      tt**3)
         +  6*cr**5*ct   *( 5*tr**3 +   5*tr**2*tt +    3*tr*tt**2 +      tt**3)
         +  5*cr**4*ct**2*( 5*tr**3 +   9*tr**2*tt +    9*tr*tt**2 +    5*tt**3)
         + 20*cr**3*ct**3*(   tr**3 +   3*tr**2*tt +    5*tr*tt**2 +    5*tt**3)
         + 15*cr**2*ct**4*(   tr**3 +   5*tr**2*tt +   15*tr*tt**2 +   35*tt**3)
         +  2*cr   *ct**5*( 5*tr**3 +  45*tr**2*tt +  315*tr*tt**2 - 3123*tt**3)
         +          ct**6*( 5*tr**3 + 105*tr**2*tt - 3123*tr*tt**2 + 4329*tt**3)
    )
    al= (ln(rr)-ln(rt))*(
         - 1680*cr*ct**5* tt**3
         -  840   *ct**6*(tr*tt**2 - 3*tt**3)
    )
    b = 4*(
         +    cr**6      *(90*tr**3 +  40*tr**2*tt +   14*tr*tt**2 +    3*tt**3)
         +  2*cr**5*ct   *(40*tr**3 +  42*tr**2*tt +   27*tr*tt**2 +   10*tt**3)
         +  5*cr**4*ct**2*(14*tr**3 +  27*tr**2*tt +   30*tr*tt**2 +   20*tt**3)
         + 20*cr**3*ct**3*( 3*tr**3 +  10*tr**2*tt +   20*tr*tt**2 +   30*tt**3)
         +  5*cr**2*ct**4*(10*tr**3 +  60*tr**2*tt +  270*tr*tt**2 - 1089*tt**3)
         +  2*cr   *ct**5*(20*tr**3 + 270*tr**2*tt - 3267*tr*tt**2 +  996*tt**3)
         +  3      *ct**6*(10*tr**3 - 363*tr**2*tt +  332*tr*tt**2 +  840*tt**3)
    )
    bl= 1680*(ln(rr)-ln(rt))*(
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
         + 3*ct**6*(-39*tr**3 - 7*tr**2*tt + 280*tr*tt**2 + 280*tt**3)
    )
    cl= -840*(ln(rr)-ln(rt))*(
         + 20*cr**3*ct**3*tt**3
         +  5*cr**2*ct**4*(9*tr*tt**2 - 7*tt**3)
         +  6*cr   *ct**5*(3*tr**2*tt - 7*tr*tt**2)
         +          ct**6*(tr**3 - 7*tr**2*tt)
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
    dl= 1680*(ln(rr)-ln(rt))*(
         + 10*cr**4*ct**2*tt**3
         + 20*cr**3*ct**3*(2*tr*tt**2 - tt**3)
         + 15*cr**2*ct**4*(2*tr**2*tt - 3*tr*tt**2)
         +  2*cr   *ct**5*(2*tr**3 - 9*tr**2*tt)
         -          ct**6*tr**3
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
    el= 4200*(ln(rr)-ln(rt))* (
         -  2*cr**5*ct   *tt**3
         -  5*cr**4*ct**2*(3*tr*tt**2 - tt**3)
         - 20*cr**3*ct**3*(tr**2*tt - tr*tt**2)
         -  5*cr**2*ct**4*(tr**3 - 3*tr**2*tt)
         +  2*cr   *ct**5*tr**3
    )
    f = 28*(
          + cr**6*(420*tr**3 + 280*tr**2*tt + 210*tr*tt**2 - 33*tt**3)
          + 2*cr**5*ct*(280*tr**3 + 630*tr**2*tt - 297*tr*tt**2 - 174*tt**3)
          + 15*cr**4*ct**2*(70*tr**3 - 99*tr**2*tt - 174*tr*tt**2 + 20*tt**3)
          + 20*cr**3*ct**3*(-33*tr**3 - 174*tr**2*tt + 60*tr*tt**2 + 10*tt**3)
          + 5*cr**2*ct**4*(-174*tr**3 + 180*tr**2*tt + 90*tr*tt**2 + 35*tt**3)
          + 6*cr*ct**5*(20*tr**3 + 30*tr**2*tt + 35*tr*tt**2 + 28*tt**3)
          + ct**6*(10*tr**3 + 35*tr**2*tt + 84*tr*tt**2 + 168*tt**3)
    )
    fl= 1680*(ln(rr)-ln(rt))*(
          + cr**6*tt**3
          + 2*cr**5*ct*(9*tr*tt**2 - 2*tt**3)
          + 15*cr**4*ct**2*(3*tr**2*tt - 2*tr*tt**2)
          + 20*cr**3*ct**3*(tr**3 - 2*tr**2*tt)
          - 10*cr**2*ct**4*tr**3
    )
    g = 14*( # 5*tr**3 is wrong
         + 3*cr**6*(-280*tr**3 - 280*tr**2*tt + 7*tr*tt**2 + 39*tt**3)
         + 6*cr**5*ct*(-280*tr**3 + 21*tr**2*tt + 351*tr*tt**2 - 20*tt**3)
         + 5*cr**4*ct**2*(21*tr**3 + 1053*tr**2*tt - 180*tr*tt**2 - 20*tt**3)
         + 20*cr**3*ct**3*(117*tr**3 - 60*tr**2*tt - 20*tr*tt**2 - 5*tt**3)
         + 15*cr**2*ct**4*(-20*tr**3 - 20*tr**2*tt - 15*tr*tt**2 - 7*tt**3)
         + 2*cr*ct**5*(-20*tr**3 - 45*tr**2*tt - 63*tr*tt**2 - 56*tt**3)
         + ct**6*(-5*tr**3 - 21*tr**2*tt - 56*tr*tt**2 - 120*tt**3)
    ) # through here checking
    gl= -840*(ln(rr)-ln(rt))*(
         + cr**6*(7*tr*tt**2 - tt**3)
         + 6*cr**5*ct*(7*tr**2*tt - 3*tr*tt**2)
         + 5*cr**4*ct**2*(7*tr**3 - 9*tr**2*tt)
         - 20*cr**3*ct**3*tr**3
    )
    h = (+ 12*cr**6*(840*tr**3 + 332*tr**2*tt - 363*tr*tt**2 + 10*tt**3)
         + 8*cr**5*ct*(996*tr**3 - 3267*tr**2*tt + 270*tr*tt**2 + 20*tt**3)
         + 20*cr**4*ct**2*(-1089*tr**3 + 270*tr**2*tt + 60*tr*tt**2 + 10*tt**3)
         + 80*cr**3*ct**3*(30*tr**3 + 20*tr**2*tt + 10*tr*tt**2 + 3*tt**3)
         + 20*cr**2*ct**4*(20*tr**3 + 30*tr**2*tt + 27*tr*tt**2 + 14*tt**3)
         + 8*cr*ct**5*(10*tr**3 + 27*tr**2*tt + 42*tr*tt**2 + 40*tt**3)
         + 4*ct**6*(3*tr**3 + 14*tr**2*tt + 40*tr*tt**2 + 90*tt**3)
    )
    hl= 1680*(ln(rr)-ln(rt))*(
         + 4*cr**6*tr**2*tt
         - cr**6*tr*tt**2
         + 8*cr**5*ct*tr**3
         - 6*cr**5*ct*tr**2*tt
         - 5*cr**4*ct**2*tr**3
    )
    i = (
         - 4329*cr**6*tr**3 + 3123*cr**6*tr**2*tt - 105*cr**6*tr*tt**2 - 5*cr**6*tt**3
         + 6246*cr**5*ct*tr**3 - 630*cr**5*ct*tr**2*tt - 90*cr**5*ct*tr*tt**2 - 10*cr**5*ct*tt**3
         - 525*cr**4*ct**2*tr**3 - 225*cr**4*ct**2*tr**2*tt - 75*cr**4*ct**2*tr*tt**2 - 15*cr**4*ct**2*tt**3
         - 100*cr**3*ct**3*tr**3 - 100*cr**3*ct**3*tr**2*tt - 60*cr**3*ct**3*tr*tt**2 - 20*cr**3*ct**3*tt**3
         - 25*cr**2*ct**4*tr**3 - 45*cr**2*ct**4*tr**2*tt - 45*cr**2*ct**4*tr*tt**2 - 25*cr**2*ct**4*tt**3
         - 6*cr*ct**5*tr**3 - 18*cr*ct**5*tr**2*tt - 30*cr*ct**5*tr*tt**2 - 30*cr*ct**5*tt**3
         - ct**6*tr**3 - 5*ct**6*tr**2*tt - 15*ct**6*tr*tt**2 - 35*ct**6*tt**3
    )
    il= 840*(ln(rr)-ln(rt))*(
         - 3*cr**6*tr**3
         + cr**6*tr**2*tt
         + 2*cr**5*ct*tr**3
    )
    j = - 280*cr**6*tr**3
    k = - 280*ct**6*tt**3

    den = (13440*pi**2*rr*rt*(rr-rt)**9)
    factor = Nb**3*u0**3
    T200 = factor*( 
         + (a+al)*rr**9*rt    + (b+bl)*rr**8*rt**2
         + (c+cl)*rr**7*rt**3 + (d+dl)*rr**6*rt**4
         + (e+el)*rr**5*rt**5 + (f+fl)*rr**4*rt**6
         + (g+gl)*rr**3*rt**7 + (h+hl)*rr**2*rt**8
         + (i+il)*rr   *rt**9 + j*rt**10 + k*rr**10
     ) / den
     
    print("should be zero:",simp((T_ints[7] - sy.expand(T200))*den/factor))
    
    # T020 = Nb*u0*(-10*cr**2*rr**3*tr - 2*cr**2*rr**3*tt + 6*cr**2*rr**2*rt*tr + 3*cr**2*rr*rt**2*tr + cr**2*rr*rt**2*tt + cr**2*rt**3*tr + cr**2*rt**3*tt - 4*cr*ct*rr**3*tr - 2*cr*ct*rr**3*tt - 2*cr*ct*rr**2*rt*tt + 2*cr*ct*rr*rt**2*tr + 2*cr*ct*rt**3*tr + 4*cr*ct*rt**3*tt - ct**2*rr**3*tr - ct**2*rr**3*tt - ct**2*rr**2*rt*tr - 3*ct**2*rr**2*rt*tt - 6*ct**2*rr*rt**2*tt + 2*ct**2*rt**3*tr + 10*ct**2*rt**3*tt)/120
    # T002 = Nb*u0*(-10*cr**2*rr**3*tr - 2*cr**2*rr**3*tt + 6*cr**2*rr**2*rt*tr + 3*cr**2*rr*rt**2*tr + cr**2*rr*rt**2*tt + cr**2*rt**3*tr + cr**2*rt**3*tt - 4*cr*ct*rr**3*tr - 2*cr*ct*rr**3*tt - 2*cr*ct*rr**2*rt*tt + 2*cr*ct*rr*rt**2*tr + 2*cr*ct*rt**3*tr + 4*cr*ct*rt**3*tt - ct**2*rr**3*tr - ct**2*rr**3*tt - ct**2*rr**2*rt*tr - 3*ct**2*rr**2*rt*tt - 6*ct**2*rr*rt**2*tt + 2*ct**2*rt**3*tr + 10*ct**2*rt**3*tt)/120
    
    

    a = -tr*(10*cr**2 + 4*cr*ct + ct**2) - tt*(2*cr**2 + 2*cr*ct + ct**2)
    b = + tr*(6*cr**2 - ct**2) - tt*(2*cr*ct + 3*ct**2)
    d = + tr*(3*cr**2 + 2*cr*ct) + tt*(cr**2 - 6*ct**2)
    f = + tr*(cr**2 + 2*cr*ct + 2*ct**2) + tt*(cr**2 + 4*cr*ct + 10*ct**2)
    T002 = Nb*u0*( rr**3 * a + rr**2*rt * b + rr*rt**2 * d + rt**3 * f )/120