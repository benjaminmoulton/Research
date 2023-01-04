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
    Rt = tt / tr
    Rc = ct / cr
    Rr = rt / rr

    # >>> a = np.array([16800,67200,33600,50400,75600,6720,30240,1680])
    # >>> print(np.gcd.reduce(a),a / np.gcd.reduce(a))
    
    a = Rc**6*Rt**3
    b = -1*(
         +          (35 +  15*Rt +    5*Rt**2 +      Rt**3)
         +  6*Rc   *( 5 +   5*Rt +    3*Rt**2 +      Rt**3)
         +  5*Rc**2*( 5 +   9*Rt +    9*Rt**2 +    5*Rt**3)
         + 20*Rc**3*( 1 +   3*Rt +    5*Rt**2 +    5*Rt**3)
         + 15*Rc**4*( 1 +   5*Rt +   15*Rt**2 +   35*Rt**3)
         +  2*Rc**5*( 5 +  45*Rt +  315*Rt**2 - 3123*Rt**3)
         +    Rc**6*( 5 + 105*Rt - 3123*Rt**2 + 4329*Rt**3)
    )
    j = -1*(
         +          ( 4329 - 3123*Rt + 105*Rt**2 +  5*Rt**3)
         +  2*Rc   *(-3123 +  315*Rt +  45*Rt**2 +  5*Rt**3)
         + 15*Rc**2*(   35 +   15*Rt +   5*Rt**2 +    Rt**3)
         + 20*Rc**3*(    5 +    5*Rt +   3*Rt**2 +    Rt**3)
         +  5*Rc**4*(    5 +    9*Rt +   9*Rt**2 +  5*Rt**3)
         +  6*Rc**5*(    1 +    3*Rt +   5*Rt**2 +  5*Rt**3)
         +    Rc**6*(    1 +    5*Rt +  15*Rt**2 + 35*Rt**3)
    )
    bl= -840*ln(Rr)*(
         - 2*Rc**5* Rt**3
         -   Rc**6*(Rt**2 - 3*Rt**3)
    )
    jl= -840*ln(Rr)*(
         + (Rt - 3)
         + 2*Rc
    )
    c = 4*(
         +          (90 +  40*Rt +   14*Rt**2 +    3*Rt**3)
         +  2*Rc   *(40 +  42*Rt +   27*Rt**2 +   10*Rt**3)
         +  5*Rc**2*(14 +  27*Rt +   30*Rt**2 +   20*Rt**3)
         + 20*Rc**3*( 3 +  10*Rt +   20*Rt**2 +   30*Rt**3)
         +  5*Rc**4*(10 +  60*Rt +  270*Rt**2 - 1089*Rt**3)
         +  2*Rc**5*(20 + 270*Rt - 3267*Rt**2 +  996*Rt**3)
         +  3*Rc**6*(10 - 363*Rt +  332*Rt**2 +  840*Rt**3)
    )
    i = 4*(
         +  3*      (  840 +  332*Rt - 363*Rt**2 + 10*Rt**3)
         +  2*Rc   *(  996 - 3267*Rt + 270*Rt**2 + 20*Rt**3)
         +  5*Rc**2*(-1089 +  270*Rt +  60*Rt**2 + 10*Rt**3)
         + 20*Rc**3*(   30 +   20*Rt +  10*Rt**2 +  3*Rt**3)
         +  5*Rc**4*(   20 +   30*Rt +  27*Rt**2 + 14*Rt**3)
         +  2*Rc**5*(   10 +   27*Rt +  42*Rt**2 + 40*Rt**3)
         +    Rc**6*(    3 +   14*Rt +  40*Rt**2 + 90*Rt**3)
    )
    cl= -1680*ln(Rr)*(
         + 5*Rc**4* Rt**3
         + 2*Rc**5*(3*Rt**2 - 4*Rt**3)
         +   Rc**6*(Rt - 4*Rt**2)
    )
    il= -1680*ln(Rr)*(
         - 5*Rc**2
         - 2*Rc   *(3*Rt - 4)
         -         (Rt**2 - 4*Rt)
    )
    d = -14*(
         +          (120 +  56*Rt +   21*Rt**2 +   5*Rt**3)
         +  2*Rc   *( 56 +  63*Rt +   45*Rt**2 +  20*Rt**3)
         + 15*Rc**2*(  7 +  15*Rt +   20*Rt**2 +  20*Rt**3)
         + 20*Rc**3*(  5 +  20*Rt +   60*Rt**2 - 117*Rt**3)
         +  5*Rc**4*( 20 + 180*Rt - 1053*Rt**2 -  21*Rt**3)
         +  6*Rc**5*( 20 - 351*Rt -   21*Rt**2 + 280*Rt**3)
         +  3*Rc**6*(-39 -   7*Rt +  280*Rt**2 + 280*Rt**3)
    )
    h = -14*(
         +  3      *( 280 +  280*Rt -   7*Rt**2 -  39*Rt**3)
         +  6*Rc   *( 280 -   21*Rt - 351*Rt**2 +  20*Rt**3)
         +  5*Rc**2*( -21 - 1053*Rt + 180*Rt**2 +  20*Rt**3)
         + 20*Rc**3*(-117 +   60*Rt +  20*Rt**2 +   5*Rt**3)
         + 15*Rc**4*(  20 +   20*Rt +  15*Rt**2 +   7*Rt**3)
         +  2*Rc**5*(  20 +   45*Rt +  63*Rt**2 +  56*Rt**3)
         +    Rc**6*(   5 +   21*Rt +  56*Rt**2 + 120*Rt**3)
    )
    dl= 840*ln(Rr)*(
         + 20*Rc**3*Rt**3
         +  5*Rc**4*(9*Rt**2 - 7*Rt**3)
         +  6*Rc**5*(3*Rt - 7*Rt**2)
         +    Rc**6*(1 - 7*Rt)
    )
    hl= 840*ln(Rr)*(
         +          (7*Rt**2 - Rt**3)
         +  6*Rc   *(7*Rt - 3*Rt**2)
         +  5*Rc**2*(7 - 9*Rt)
         - 20*Rc**3
    )
    e = 28*(
         +          ( 168 +  84*Rt +  35*Rt**2 +  10*Rt**3)
         +  6*Rc   *(  28 +  35*Rt +  30*Rt**2 +  20*Rt**3)
         +  5*Rc**2*(  35 +  90*Rt + 180*Rt**2 - 174*Rt**3)
         + 20*Rc**3*(  10 +  60*Rt - 174*Rt**2 -  33*Rt**3)
         + 15*Rc**4*(  20 - 174*Rt -  99*Rt**2 +  70*Rt**3)
         +  2*Rc**5*(-174 - 297*Rt + 630*Rt**2 + 280*Rt**3)
         +    Rc**6*(- 33 + 210*Rt + 280*Rt**2 + 420*Rt**3)
    )
    g = 28*(
          +          ( 420 + 280*Rt + 210*Rt**2 -  33*Rt**3)
          +  2*Rc   *( 280 + 630*Rt - 297*Rt**2 - 174*Rt**3)
          + 15*Rc**2*(  70 -  99*Rt - 174*Rt**2 +  20*Rt**3)
          + 20*Rc**3*( -33 - 174*Rt +  60*Rt**2 +  10*Rt**3)
          +  5*Rc**4*(-174 + 180*Rt +  90*Rt**2 +  35*Rt**3)
          +  6*Rc**5*(  20 +  30*Rt +  35*Rt**2 +  28*Rt**3)
          +    Rc**6*(  10 +  35*Rt +  84*Rt**2 + 168*Rt**3)
    )
    el= -1680*ln(Rr)*(
         + 10*Rc**2*Rt**3
         + 20*Rc**3*(2*Rt**2 - Rt**3)
         + 15*Rc**4*(2*Rt - 3*Rt**2)
         +  2*Rc**5*(2 - 9*Rt)
         -    Rc**6
    )
    gl= -1680*ln(Rr)*(
          +           Rt**3
          -  2*Rc   *(2*Rt**3 - 9*Rt**2)
          - 15*Rc**2*(2*Rt**2 - 3*Rt)
          - 20*Rc**3*(2*Rt - 1)
          - 10*Rc**4
    )
    f = -70*(
         +           (126 + 70*Rt + 35*Rt**2 +  15*Rt**3)
         +  10*Rc   *( 14 + 21*Rt + 27*Rt**2 -  12*Rt**3)
         +  25*Rc**2*(  7 + 27*Rt - 36*Rt**2 -  12*Rt**3)
         + 300*Rc**3*(  1 -  4*Rt -  4*Rt**2 +     Rt**3)
         -  25*Rc**4*( 12 + 36*Rt - 27*Rt**2 -   7*Rt**3)
         -  10*Rc**5*( 12 - 27*Rt - 21*Rt**2 -  14*Rt**3)
         +     Rc**6*( 15 + 35*Rt + 70*Rt**2 + 126*Rt**3)
    )
    fl= -4200*ln(Rr)* (
         -  2*Rc   * Rt**3
         -  5*Rc**2*(3*Rt**2 - Rt**3)
         - 20*Rc**3*(Rt - Rt**2)
         -  5*Rc**4*(1 - 3*Rt)
         +  2*Rc**5
    )

    den = (13440*pi**2*rr*Rr*(1-Rr)**9)
    factor = Nb**3*u0**3
    T200 = factor  * cr**6*tr**3*(
         - 280*(a + Rr**10)
         + (b+bl)*Rr + (j+jl)*Rr**9 
         + (c+cl)*Rr**2 + (i+il)*Rr**8 
         + (d+dl)*Rr**3 + (h+hl)*Rr**7 
         + (e+el)*Rr**4 + (g+gl)*Rr**6 
         + (f+fl)*Rr**5 
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