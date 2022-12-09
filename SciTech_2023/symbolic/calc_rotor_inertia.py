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
    
    a = - 280*ct**6*tt**3
    b = -1*(
         +    cr**6      *(35*tr**3 +  15*tr**2*tt +    5*tr*tt**2 +      tt**3)
         +  6*cr**5*ct   *( 5*tr**3 +   5*tr**2*tt +    3*tr*tt**2 +      tt**3)
         +  5*cr**4*ct**2*( 5*tr**3 +   9*tr**2*tt +    9*tr*tt**2 +    5*tt**3)
         + 20*cr**3*ct**3*(   tr**3 +   3*tr**2*tt +    5*tr*tt**2 +    5*tt**3)
         + 15*cr**2*ct**4*(   tr**3 +   5*tr**2*tt +   15*tr*tt**2 +   35*tt**3)
         +  2*cr   *ct**5*( 5*tr**3 +  45*tr**2*tt +  315*tr*tt**2 - 3123*tt**3)
         +          ct**6*( 5*tr**3 + 105*tr**2*tt - 3123*tr*tt**2 + 4329*tt**3)
    )
    bl= (ln(rr)-ln(rt))*(
         - 1680*cr*ct**5* tt**3
         -  840   *ct**6*(tr*tt**2 - 3*tt**3)
    )
    c = 4*(
         +    cr**6      *(90*tr**3 +  40*tr**2*tt +   14*tr*tt**2 +    3*tt**3)
         +  2*cr**5*ct   *(40*tr**3 +  42*tr**2*tt +   27*tr*tt**2 +   10*tt**3)
         +  5*cr**4*ct**2*(14*tr**3 +  27*tr**2*tt +   30*tr*tt**2 +   20*tt**3)
         + 20*cr**3*ct**3*( 3*tr**3 +  10*tr**2*tt +   20*tr*tt**2 +   30*tt**3)
         +  5*cr**2*ct**4*(10*tr**3 +  60*tr**2*tt +  270*tr*tt**2 - 1089*tt**3)
         +  2*cr   *ct**5*(20*tr**3 + 270*tr**2*tt - 3267*tr*tt**2 +  996*tt**3)
         +  3      *ct**6*(10*tr**3 - 363*tr**2*tt +  332*tr*tt**2 +  840*tt**3)
    )
    cl= 1680*(ln(rr)-ln(rt))*(
         + 5*cr**2*ct**4* tt**3
         + 2*cr   *ct**5*(3*tr*tt**2 - 4*tt**3)
         +         ct**6*(tr**2*tt - 4*tr*tt**2)
    )
    d = -14*(
         + cr**6*(120*tr**3 + 56*tr**2*tt + 21*tr*tt**2 + 5*tt**3)
         + 2*cr**5*ct*(56*tr**3 + 63*tr**2*tt + 45*tr*tt**2 + 20*tt**3)
         + 15*cr**4*ct**2*(7*tr**3 + 15*tr**2*tt + 20*tr*tt**2 + 20*tt**3)
         + 20*cr**3*ct**3*(5*tr**3 + 20*tr**2*tt + 60*tr*tt**2 - 117*tt**3)
         + 5*cr**2*ct**4*(20*tr**3 + 180*tr**2*tt - 1053*tr*tt**2 - 21*tt**3)
         + 6*cr*ct**5*(20*tr**3 - 351*tr**2*tt - 21*tr*tt**2 + 280*tt**3)
         + 3*ct**6*(-39*tr**3 - 7*tr**2*tt + 280*tr*tt**2 + 280*tt**3)
    )
    dl= -840*(ln(rr)-ln(rt))*(
         + 20*cr**3*ct**3*tt**3
         +  5*cr**2*ct**4*(9*tr*tt**2 - 7*tt**3)
         +  6*cr   *ct**5*(3*tr**2*tt - 7*tr*tt**2)
         +          ct**6*(tr**3 - 7*tr**2*tt)
    )
    e = 28*(
         +    cr**6      *( 168*tr**3 +  84*tr**2*tt +  35*tr*tt**2 +  10*tt**3)
         +  6*cr**5*ct   *(  28*tr**3 +  35*tr**2*tt +  30*tr*tt**2 +  20*tt**3)
         +  5*cr**4*ct**2*(  35*tr**3 +  90*tr**2*tt + 180*tr*tt**2 - 174*tt**3)
         + 20*cr**3*ct**3*(  10*tr**3 +  60*tr**2*tt - 174*tr*tt**2 -  33*tt**3)
         + 15*cr**2*ct**4*(  20*tr**3 - 174*tr**2*tt -  99*tr*tt**2 +  70*tt**3)
         +  2*cr   *ct**5*(-174*tr**3 - 297*tr**2*tt + 630*tr*tt**2 + 280*tt**3)
         +          ct**6*(- 33*tr**3 + 210*tr**2*tt + 280*tr*tt**2 + 420*tt**3)
    )
    el= 1680*(ln(rr)-ln(rt))*(
         + 10*cr**4*ct**2*tt**3
         + 20*cr**3*ct**3*(2*tr*tt**2 - tt**3)
         + 15*cr**2*ct**4*(2*tr**2*tt - 3*tr*tt**2)
         +  2*cr   *ct**5*(2*tr**3 - 9*tr**2*tt)
         -          ct**6*tr**3
    )
    f = -70*(
         + cr**6*(126*tr**3 + 70*tr**2*tt + 35*tr*tt**2 + 15*tt**3)
         + 10*cr**5*ct*(14*tr**3 + 21*tr**2*tt + 27*tr*tt**2 - 12*tt**3)
         + 25*cr**4*ct**2*(7*tr**3 + 27*tr**2*tt - 36*tr*tt**2 - 12*tt**3)
         + 300*cr**3*ct**3*(tr**3 - 4*tr**2*tt - 4*tr*tt**2 + tt**3)
         - 25*cr**2*ct**4*(12*tr**3 + 36*tr**2*tt - 27*tr*tt**2 - 7*tt**3)
         - 10*cr*ct**5*(12*tr**3 - 27*tr**2*tt - 21*tr*tt**2 - 14*tt**3)
         + ct**6*(15*tr**3 + 35*tr**2*tt + 70*tr*tt**2 + 126*tt**3)
    )
    fl= 4200*(ln(rr)-ln(rt))* (
         -  2*cr**5*ct   *tt**3
         -  5*cr**4*ct**2*(3*tr*tt**2 - tt**3)
         - 20*cr**3*ct**3*(tr**2*tt - tr*tt**2)
         -  5*cr**2*ct**4*(tr**3 - 3*tr**2*tt)
         +  2*cr   *ct**5*tr**3
    )
    g = 28*(
          + cr**6*(420*tr**3 + 280*tr**2*tt + 210*tr*tt**2 - 33*tt**3)
          + 2*cr**5*ct*(280*tr**3 + 630*tr**2*tt - 297*tr*tt**2 - 174*tt**3)
          + 15*cr**4*ct**2*(70*tr**3 - 99*tr**2*tt - 174*tr*tt**2 + 20*tt**3)
          + 20*cr**3*ct**3*(-33*tr**3 - 174*tr**2*tt + 60*tr*tt**2 + 10*tt**3)
          + 5*cr**2*ct**4*(-174*tr**3 + 180*tr**2*tt + 90*tr*tt**2 + 35*tt**3)
          + 6*cr*ct**5*(20*tr**3 + 30*tr**2*tt + 35*tr*tt**2 + 28*tt**3)
          + ct**6*(10*tr**3 + 35*tr**2*tt + 84*tr*tt**2 + 168*tt**3)
    )
    gl= 1680*(ln(rr)-ln(rt))*(
          + cr**6*tt**3
          + 2*cr**5*ct*(9*tr*tt**2 - 2*tt**3)
          + 15*cr**4*ct**2*(3*tr**2*tt - 2*tr*tt**2)
          + 20*cr**3*ct**3*(tr**3 - 2*tr**2*tt)
          - 10*cr**2*ct**4*tr**3
    )
    h = 14*( # 5*tr**3 is wrong
         + 3*cr**6*(-280*tr**3 - 280*tr**2*tt + 7*tr*tt**2 + 39*tt**3)
         + 6*cr**5*ct*(-280*tr**3 + 21*tr**2*tt + 351*tr*tt**2 - 20*tt**3)
         + 5*cr**4*ct**2*(21*tr**3 + 1053*tr**2*tt - 180*tr*tt**2 - 20*tt**3)
         + 20*cr**3*ct**3*(117*tr**3 - 60*tr**2*tt - 20*tr*tt**2 - 5*tt**3)
         + 15*cr**2*ct**4*(-20*tr**3 - 20*tr**2*tt - 15*tr*tt**2 - 7*tt**3)
         + 2*cr*ct**5*(-20*tr**3 - 45*tr**2*tt - 63*tr*tt**2 - 56*tt**3)
         + ct**6*(-5*tr**3 - 21*tr**2*tt - 56*tr*tt**2 - 120*tt**3)
    ) # through here checking
    hl= -840*(ln(rr)-ln(rt))*(
         + cr**6*(7*tr*tt**2 - tt**3)
         + 6*cr**5*ct*(7*tr**2*tt - 3*tr*tt**2)
         + 5*cr**4*ct**2*(7*tr**3 - 9*tr**2*tt)
         - 20*cr**3*ct**3*tr**3
    )
    i = 4*(
         + 3*cr**6*(840*tr**3 + 332*tr**2*tt - 363*tr*tt**2 + 10*tt**3)
         + 2*cr**5*ct*(996*tr**3 - 3267*tr**2*tt + 270*tr*tt**2 + 20*tt**3)
         + 5*cr**4*ct**2*(-1089*tr**3 + 270*tr**2*tt + 60*tr*tt**2 + 10*tt**3)
         + 20*cr**3*ct**3*(30*tr**3 + 20*tr**2*tt + 10*tr*tt**2 + 3*tt**3)
         + 5*cr**2*ct**4*(20*tr**3 + 30*tr**2*tt + 27*tr*tt**2 + 14*tt**3)
         + 2*cr*ct**5*(10*tr**3 + 27*tr**2*tt + 42*tr*tt**2 + 40*tt**3)
         + ct**6*(3*tr**3 + 14*tr**2*tt + 40*tr*tt**2 + 90*tt**3)
    )
    il= 1680*(ln(rr)-ln(rt))*(
         + 4*cr**6*tr**2*tt
         - cr**6*tr*tt**2
         + 8*cr**5*ct*tr**3
         - 6*cr**5*ct*tr**2*tt
         - 5*cr**4*ct**2*tr**3
    )
    j = -1*(
         + cr**6*(4329*tr**3 - 3123*tr**2*tt + 105*tr*tt**2 + 5*tt**3)
         + 2*cr**5*ct*(-3123*tr**3 + 315*tr**2*tt + 45*tr*tt**2 + 5*tt**3)
         + 15*cr**4*ct**2*(35*tr**3 + 15*tr**2*tt + 5*tr*tt**2 + tt**3)
         + 20*cr**3*ct**3*(5*tr**3 + 5*tr**2*tt + 3*tr*tt**2 + tt**3)
         + 5*cr**2*ct**4*(5*tr**3 + 9*tr**2*tt + 9*tr*tt**2 + 5*tt**3)
         + 6*cr*ct**5*(tr**3 + 3*tr**2*tt + 5*tr*tt**2 + 5*tt**3)
         + ct**6*(tr**3 + 5*tr**2*tt + 15*tr*tt**2 + 35*tt**3)
    )
    jl= 840*(ln(rr)-ln(rt))*(
         - 3*cr**6*tr**3
         + cr**6*tr**2*tt
         + 2*cr**5*ct*tr**3
    )
    k = - 280*cr**6*tr**3

    den = (13440*pi**2*rr*rt*(rr-rt)**9)
    factor = Nb**3*u0**3
    T200 = factor*(
         + a*rr**10
         + (b+bl)*rr**9*rt    + (c+cl)*rr**8*rt**2
         + (d+dl)*rr**7*rt**3 + (e+el)*rr**6*rt**4
         + (f+fl)*rr**5*rt**5 + (g+gl)*rr**4*rt**6
         + (h+hl)*rr**3*rt**7 + (i+il)*rr**2*rt**8
         + (j+jl)*rr   *rt**9 + k*rt**10
     ) / den
     
    print("should be zero:",simp((T_ints[7] - sy.expand(T200))*den/factor))
    

    a = -tr*(10*cr**2 + 4*cr*ct + ct**2) - tt*(2*cr**2 + 2*cr*ct + ct**2)
    b = + tr*(6*cr**2 - ct**2) - tt*(2*cr*ct + 3*ct**2)
    d = + tr*(3*cr**2 + 2*cr*ct) + tt*(cr**2 - 6*ct**2)
    f = + tr*(cr**2 + 2*cr*ct + 2*ct**2) + tt*(cr**2 + 4*cr*ct + 10*ct**2)
    T002 = Nb*u0*( rr**3 * a + rr**2*rt * b + rr*rt**2 * d + rt**3 * f )/120



#     a = - 280 c_t^6 t_t^3
#     b = -1 (
#          +    c_r^6       (35 t_r^3 +  15 t_r^2 t_t +    5 t_r t_t^2 +      t_t^3)
#          +  6 c_r^5 c_t    ( 5 t_r^3 +   5 t_r^2 t_t +    3 t_r t_t^2 +      t_t^3)
#          +  5 c_r^4 c_t^2 ( 5 t_r^3 +   9 t_r^2 t_t +    9 t_r t_t^2 +    5 t_t^3)
#          + 20 c_r^3 c_t^3 (   t_r^3 +   3 t_r^2 t_t +    5 t_r t_t^2 +    5 t_t^3)
#          + 15 c_r^2 c_t^4 (   t_r^3 +   5 t_r^2 t_t +   15 t_r t_t^2 +   35 t_t^3)
#          +  2 c_r    c_t^5 ( 5 t_r^3 +  45 t_r^2 t_t +  315 t_r t_t^2 - 3123 t_t^3)
#          +          c_t^6 ( 5 t_r^3 + 105 t_r^2 t_t - 3123 t_r t_t^2 + 4329 t_t^3)
#     )
#     bl= \ln \left( \frac{r_r}{r_t} \right) (
#          - 1680 c_r c_t^5  t_t^3
#          -  840    c_t^6 (t_r t_t^2 - 3 t_t^3)
#     )
#     c = 4 \Big(
#          +    c_r^6       (90 t_r^3 +  40 t_r^2 t_t +   14 t_r t_t^2 +    3 t_t^3)
#          +  2 c_r^5 c_t    (40 t_r^3 +  42 t_r^2 t_t +   27 t_r t_t^2 +   10 t_t^3)
#          +  5 c_r^4 c_t^2 (14 t_r^3 +  27 t_r^2 t_t +   30 t_r t_t^2 +   20 t_t^3)
#          + 20 c_r^3 c_t^3 ( 3 t_r^3 +  10 t_r^2 t_t +   20 t_r t_t^2 +   30 t_t^3)
#          +  5 c_r^2 c_t^4 (10 t_r^3 +  60 t_r^2 t_t +  270 t_r t_t^2 - 1089 t_t^3)
#          +  2 c_r    c_t^5 (20 t_r^3 + 270 t_r^2 t_t - 3267 t_r t_t^2 +  996 t_t^3)
#          +  3       c_t^6 (10 t_r^3 - 363 t_r^2 t_t +  332 t_r t_t^2 +  840 t_t^3)
#     \Big)
#     cl= 1680 \ln \left( \frac{r_r}{r_t} \right) \Big(
#          + 5 c_r^2 c_t^4  t_t^3
#          + 2 c_r    c_t^5 (3 t_r t_t^2 - 4 t_t^3)
#          +         c_t^6 (t_r^2 t_t - 4 t_r t_t^2)
#     \Big)
#     d = -14 \Big(
#          + c_r^6 (120 t_r^3 + 56 t_r^2 t_t + 21 t_r t_t^2 + 5 t_t^3)
#          + 2 c_r^5 c_t (56 t_r^3 + 63 t_r^2 t_t + 45 t_r t_t^2 + 20 t_t^3)
#          + 15 c_r^4 c_t^2 (7 t_r^3 + 15 t_r^2 t_t + 20 t_r t_t^2 + 20 t_t^3)
#          + 20 c_r^3 c_t^3 (5 t_r^3 + 20 t_r^2 t_t + 60 t_r t_t^2 - 117 t_t^3)
#          + 5 c_r^2 c_t^4 (20 t_r^3 + 180 t_r^2 t_t - 1053 t_r t_t^2 - 21 t_t^3)
#          + 6 c_r c_t^5 (20 t_r^3 - 351 t_r^2 t_t - 21 t_r t_t^2 + 280 t_t^3)
#          + 3 c_t^6 (-39 t_r^3 - 7 t_r^2 t_t + 280 t_r t_t^2 + 280 t_t^3)
#     \Big)
#     dl= -840 \ln \left( \frac{r_r}{r_t} \right) \Big(
#          + 20 c_r^3 c_t^3 t_t^3
#          +  5 c_r^2 c_t^4 (9 t_r t_t^2 - 7 t_t^3)
#          +  6 c_r    c_t^5 (3 t_r^2 t_t - 7 t_r t_t^2)
#          +          c_t^6 (t_r^3 - 7 t_r^2 t_t)
#     \Big)
#     e = 28 \Big(
#          +    c_r^6       ( 168 t_r^3 +  84 t_r^2 t_t +  35 t_r t_t^2 +  10 t_t^3)
#          +  6 c_r^5 c_t    (  28 t_r^3 +  35 t_r^2 t_t +  30 t_r t_t^2 +  20 t_t^3)
#          +  5 c_r^4 c_t^2 (  35 t_r^3 +  90 t_r^2 t_t + 180 t_r t_t^2 - 174 t_t^3)
#          + 20 c_r^3 c_t^3 (  10 t_r^3 +  60 t_r^2 t_t - 174 t_r t_t^2 -  33 t_t^3)
#          + 15 c_r^2 c_t^4 (  20 t_r^3 - 174 t_r^2 t_t -  99 t_r t_t^2 +  70 t_t^3)
#          +  2 c_r    c_t^5 (-174 t_r^3 - 297 t_r^2 t_t + 630 t_r t_t^2 + 280 t_t^3)
#          +          c_t^6 (- 33 t_r^3 + 210 t_r^2 t_t + 280 t_r t_t^2 + 420 t_t^3)
#     \Big)
#     el= 1680 \ln \left( \frac{r_r}{r_t} \right) \Big(
#          + 10 c_r^4 c_t^2 t_t^3
#          + 20 c_r^3 c_t^3 (2 t_r t_t^2 - t_t^3)
#          + 15 c_r^2 c_t^4 (2 t_r^2 t_t - 3 t_r t_t^2)
#          +  2 c_r    c_t^5 (2 t_r^3 - 9 t_r^2 t_t)
#          -          c_t^6 t_r^3
#     \Big)
#     f = -70 \Big(
#          + c_r^6 (126 t_r^3 + 70 t_r^2 t_t + 35 t_r t_t^2 + 15 t_t^3)
#          + 10 c_r^5 c_t (14 t_r^3 + 21 t_r^2 t_t + 27 t_r t_t^2 - 12 t_t^3)
#          + 25 c_r^4 c_t^2 (7 t_r^3 + 27 t_r^2 t_t - 36 t_r t_t^2 - 12 t_t^3)
#          + 300 c_r^3 c_t^3 (t_r^3 - 4 t_r^2 t_t - 4 t_r t_t^2 + t_t^3)
#          - 25 c_r^2 c_t^4 (12 t_r^3 + 36 t_r^2 t_t - 27 t_r t_t^2 - 7 t_t^3)
#          - 10 c_r c_t^5 (12 t_r^3 - 27 t_r^2 t_t - 21 t_r t_t^2 - 14 t_t^3)
#          + c_t^6 (15 t_r^3 + 35 t_r^2 t_t + 70 t_r t_t^2 + 126 t_t^3)
#     \Big)
#     fl= 4200 \ln \left( \frac{r_r}{r_t} \right)  \Big(
#          -  2 c_r^5 c_t    t_t^3
#          -  5 c_r^4 c_t^2 (3 t_r t_t^2 - t_t^3)
#          - 20 c_r^3 c_t^3 (t_r^2 t_t - t_r t_t^2)
#          -  5 c_r^2 c_t^4 (t_r^3 - 3 t_r^2 t_t)
#          +  2 c_r    c_t^5 t_r^3
#     \Big)
#     g = 28 \Big(
#           + c_r^6 (420 t_r^3 + 280 t_r^2 t_t + 210 t_r t_t^2 - 33 t_t^3)
#           + 2 c_r^5 c_t (280 t_r^3 + 630 t_r^2 t_t - 297 t_r t_t^2 - 174 t_t^3)
#           + 15 c_r^4 c_t^2 (70 t_r^3 - 99 t_r^2 t_t - 174 t_r t_t^2 + 20 t_t^3)
#           + 20 c_r^3 c_t^3 (-33 t_r^3 - 174 t_r^2 t_t + 60 t_r t_t^2 + 10 t_t^3)
#           + 5 c_r^2 c_t^4 (-174 t_r^3 + 180 t_r^2 t_t + 90 t_r t_t^2 + 35 t_t^3)
#           + 6 c_r c_t^5 (20 t_r^3 + 30 t_r^2 t_t + 35 t_r t_t^2 + 28 t_t^3)
#           + c_t^6 (10 t_r^3 + 35 t_r^2 t_t + 84 t_r t_t^2 + 168 t_t^3)
#     \Big)
#     gl= 1680 \ln \left( \frac{r_r}{r_t} \right) \Big(
#           + c_r^6 t_t^3
#           + 2 c_r^5 c_t (9 t_r t_t^2 - 2 t_t^3)
#           + 15 c_r^4 c_t^2 (3 t_r^2 t_t - 2 t_r t_t^2)
#           + 20 c_r^3 c_t^3 (t_r^3 - 2 t_r^2 t_t)
#           - 10 c_r^2 c_t^4 t_r^3
#     \Big)
#     h = 14 \Big(
#          + 3 c_r^6 (-280 t_r^3 - 280 t_r^2 t_t + 7 t_r t_t^2 + 39 t_t^3)
#          + 6 c_r^5 c_t (-280 t_r^3 + 21 t_r^2 t_t + 351 t_r t_t^2 - 20 t_t^3)
#          + 5 c_r^4 c_t^2 (21 t_r^3 + 1053 t_r^2 t_t - 180 t_r t_t^2 - 20 t_t^3)
#          + 20 c_r^3 c_t^3 (117 t_r^3 - 60 t_r^2 t_t - 20 t_r t_t^2 - 5 t_t^3)
#          + 15 c_r^2 c_t^4 (-20 t_r^3 - 20 t_r^2 t_t - 15 t_r t_t^2 - 7 t_t^3)
#          + 2 c_r c_t^5 (-20 t_r^3 - 45 t_r^2 t_t - 63 t_r t_t^2 - 56 t_t^3)
#          + c_t^6 (-5 t_r^3 - 21 t_r^2 t_t - 56 t_r t_t^2 - 120 t_t^3)
#     \Big)
#     hl= -840 \ln \left( \frac{r_r}{r_t} \right) \Big(
#          + c_r^6 (7 t_r t_t^2 - t_t^3)
#          + 6 c_r^5 c_t (7 t_r^2 t_t - 3 t_r t_t^2)
#          + 5 c_r^4 c_t^2 (7 t_r^3 - 9 t_r^2 t_t)
#          - 20 c_r^3 c_t^3 t_r^3
#     \Big)
#     i = 4 * \Big(
#          + 3 c_r^6 (840 t_r^3 + 332 t_r^2 t_t - 363 t_r t_t^2 + 10 t_t^3)
#          + 2 c_r^5 c_t (996 t_r^3 - 3267 t_r^2 t_t + 270 t_r t_t^2 + 20 t_t^3)
#          + 5 c_r^4 c_t^2 (-1089 t_r^3 + 270 t_r^2 t_t + 60 t_r t_t^2 + 10 t_t^3)
#          + 20 c_r^3 c_t^3 (30 t_r^3 + 20 t_r^2 t_t + 10 t_r t_t^2 + 3 t_t^3)
#          + 5 c_r^2 c_t^4 (20 t_r^3 + 30 t_r^2 t_t + 27 t_r t_t^2 + 14 t_t^3)
#          + 2 c_r c_t^5 (10 t_r^3 + 27 t_r^2 t_t + 42 t_r t_t^2 + 40 t_t^3)
#          + c_t^6 (3 t_r^3 + 14 t_r^2 t_t + 40 t_r t_t^2 + 90 t_t^3)
#     \Big)
#     il= 1680 \ln \left( \frac{r_r}{r_t} \right) \Big(
#          + 4 c_r^6 t_r^2 t_t
#          - c_r^6 t_r t_t^2
#          + 8 c_r^5 c_t t_r^3
#          - 6 c_r^5 c_t t_r^2 t_t
#          - 5 c_r^4 c_t^2 t_r^3
#     \Big)
#     j = -1 \Big(
#          + c_r^6 (4329 t_r^3 - 3123 t_r^2 t_t + 105 t_r t_t^2 + 5 t_t^3)
#          + 2 c_r^5 c_t (-3123 t_r^3 + 315 t_r^2 t_t + 45 t_r t_t^2 + 5 t_t^3)
#          + 15 c_r^4 c_t^2 (35 t_r^3 + 15 t_r^2 t_t + 5 t_r t_t^2 + t_t^3)
#          + 20 c_r^3 c_t^3 (5 t_r^3 + 5 t_r^2 t_t + 3 t_r t_t^2 + t_t^3)
#          + 5 c_r^2 c_t^4 (5 t_r^3 + 9 t_r^2 t_t + 9 t_r t_t^2 + 5 t_t^3)
#          + 6 c_r c_t^5 (t_r^3 + 3 t_r^2 t_t + 5 t_r t_t^2 + 5 t_t^3)
#          + c_t^6 (t_r^3 + 5 t_r^2 t_t + 15 t_r t_t^2 + 35 t_t^3)
#     \Big)
#     jl= 840 \ln \left( \frac{r_r}{r_t} \right) \Big(
#          - 3 c_r^6 t_r^3
#          + c_r^6 t_r^2 t_t
#          + 2 c_r^5 c_t t_r^3
#     \Big)
#     k = - 280 c_r^6 t_r^3