import sympy as sy
import numpy as np
from matplotlib import pyplot as plt
import time


if __name__ == "__main__":
    start = time.time()
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
    pi = sy.pi
    frac = sy.Rational

    # declare variables
    print("declaring variables...")
    x = sym("x")
    r0 = sym("r0")
    r1 = sym("r1")
    t = sym("t")
    t0 = sym("t0")
    t1 = sym("t1")
    h = sym("h")
    rp = sym("r")
    p = sym("p")
    mp = sym("m")
    ay0 = sym("b0")
    ay1 = sym("b1")
    az0 = sym("c0")
    az1 = sym("c1")

    # create bounds
    print("creating bounds...")
    r = (r1-r0) * x/h + r0
    ay = ( (ay1-ay0) * x/h + ay0 )
    az = ( (az1-az0) * x/h + az0 )
    x_up = h; x_lo = 0
    theta_up = t1 # 2 * pi # 
    theta_lo = t0 # 0 # 
    # r_eq = ( (ay*cos(t))**2 + (az*sin(t))**2 )**frac(1,2)
    r_up = 1; r_lo = 0
    x_bnd = (x,x_lo,x_up)
    r_bnd = (rp,r_lo,r_up)
    theta_bnd = (t,theta_lo,theta_up)

    i = x
    j = ay * rp * cos(t)
    k = az * rp * sin(t)

    V = [
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
    V_ints = []

    print("tapered cylinder eqs...")

    for q in range(len(V)):
        f = p * i**V[q][0] * j**V[q][1] * k**V[q][2]
        iii = simp( igr(ay * az * f * rp, r_bnd,x_bnd,theta_bnd) )
        V_ints.append(iii)
        print("V{}{}{} = {}".format(V[q][0],V[q][1],V[q][2],iii))
        
    print()
    print()
    quit()

    print("tapered elliptic cylinder eqs...")
    m = simp(V_ints[0])
    cglist = ["x","y","z"]
    prodlist = ["xy","xz","yz"]
    prodcglist = [(0,1),(0,2),(1,2)]
    momlist = ["xx","yy","zz"]
    momlisti = [(8,9),(7,9),(7,8)]
    momcglist = [(1,2),(0,2),(0,1)]
    cg = []

    for q in range(len(V)):
        iii = simp(V_ints[q])
        if q == 0:
            print("V = {}".format(simp(m/p)))
            print("m = {}".format(m))
        elif q < 4:
            cg.append(simp(iii/m))
            print("{}-cg = {}".format(cglist[q-1],cg[-1]))
        elif q < 7:
            cgprod = cg[prodcglist[q-4][0]] * cg[prodcglist[q-4][1]]
            I = simp(mp*(iii/m - cgprod))
            print("I{}  = {}".format(prodlist[q-4],I))
            I = simp(mp*(iii/m))
            print("I{}o = {}".format(prodlist[q-4],I))
        else:
            ii0 = simp(V_ints[momlisti[q-7][0]])
            ii1 = simp(V_ints[momlisti[q-7][1]])
            cg0 = cg[momcglist[q-7][0]]
            cg1 = cg[momcglist[q-7][1]]
            I = simp(mp*( (ii0+ii1)/m - cg0**2 - cg1**2) )
            print("I{}  = {}".format(momlist[q-7],I))
            I = simp(mp*( (ii0+ii1)/m) )
            print("I{}o = {}".format(momlist[q-7],I))
        
    print()
    print()
    quit()

    # print("elliptic cylinder eqs...")
    # m = simp( V_ints[0].subs({ay0:ay1,az0:az1}) )
    # cglist = ["x","y","z"]
    # prodlist = ["xy","xz","yz"]
    # prodcglist = [(0,1),(0,2),(1,2)]
    # momlist = ["xx","yy","zz"]
    # momlisti = [(8,9),(7,9),(7,8)]
    # momcglist = [(1,2),(0,2),(0,1)]
    # cg = []

    # for q in range(len(V)):
    #     iii = simp( V_ints[q].subs({ay0:ay1,az0:az1}) )
    #     if q == 0:
    #         print("V = {}".format(simp(m/p)))
    #         print("m = {}".format(m))
    #     elif q < 4:
    #         cg.append(simp(iii/m))
    #         print("{}-cg = {}".format(cglist[q-1],cg[-1]))
    #     elif q < 7:
    #         cgprod = cg[prodcglist[q-4][0]] * cg[prodcglist[q-4][1]]
    #         I = simp(mp*(iii/m - cgprod))
    #         print("I{} = {}".format(prodlist[q-4],I))
    #     else:
    #         ii0 = simp( V_ints[momlisti[q-7][0]].subs({ay0:ay1,az0:az1}) )
    #         ii1 = simp( V_ints[momlisti[q-7][1]].subs({ay0:ay1,az0:az1}) )
    #         cg0 = cg[momcglist[q-7][0]]
    #         cg1 = cg[momcglist[q-7][1]]
    #         I = simp(mp*( (ii0+ii1)/m - cg0**2 - cg1**2) )
    #         print("I{} = {}".format(momlist[q-7],I))
        
    # print()
    # print()

    # print("tapered cylinder eqs...")
    # m = simp( V_ints[0].subs({ay0:r0,az0:r0,ay1:r1,az1:r1}) )
    # cglist = ["x","y","z"]
    # prodlist = ["xy","xz","yz"]
    # prodcglist = [(0,1),(0,2),(1,2)]
    # momlist = ["xx","yy","zz"]
    # momlisti = [(8,9),(7,9),(7,8)]
    # momcglist = [(1,2),(0,2),(0,1)]
    # cg = []

    # for q in range(len(V)):
    #     iii = simp( V_ints[q].subs({ay0:r0,az0:r0,ay1:r1,az1:r1}) )
    #     if q == 0:
    #         print("V = {}".format(simp(m/p)))
    #         print("m = {}".format(m))
    #     elif q < 4:
    #         cg.append(simp(iii/m))
    #         print("{}-cg = {}".format(cglist[q-1],cg[-1]))
    #     elif q < 7:
    #         cgprod = cg[prodcglist[q-4][0]] * cg[prodcglist[q-4][1]]
    #         I = simp(mp*(iii/m - cgprod))
    #         print("I{} = {}".format(prodlist[q-4],I))
    #     else:
    #         ii0 = simp( V_ints[momlisti[q-7][0]].subs({ay0:r0,az0:r0,ay1:r1,az1:r1}) )
    #         ii1 = simp( V_ints[momlisti[q-7][1]].subs({ay0:r0,az0:r0,ay1:r1,az1:r1}) )
    #         cg0 = cg[momcglist[q-7][0]]
    #         cg1 = cg[momcglist[q-7][1]]
    #         I = simp(mp*( (ii0+ii1)/m - cg0**2 - cg1**2) )
    #         print("I{} = {}".format(momlist[q-7],I))
        
    # print()
    # print()

    # print("tapered cylinder about origin eqs...")
    # m = simp( V_ints[0].subs({ay0:r0,az0:r0,ay1:r1,az1:r1}) )
    # cglist = ["x","y","z"]
    # prodlist = ["xy","xz","yz"]
    # prodcglist = [(0,1),(0,2),(1,2)]
    # momlist = ["xx","yy","zz"]
    # momlisti = [(8,9),(7,9),(7,8)]
    # momcglist = [(1,2),(0,2),(0,1)]
    # cg = []

    # for q in range(len(V)):
    #     iii = simp( V_ints[q].subs({ay0:r0,az0:r0,ay1:r1,az1:r1}) )
    #     if q == 0:
    #         print("V = {}".format(simp(m/p)))
    #         print("m = {}".format(m))
    #     elif q < 4:
    #         cg.append(simp(iii/m))
    #         print("{}-cg = {}".format(cglist[q-1],cg[-1]))
    #     elif q < 7:
    #         cgprod = cg[prodcglist[q-4][0]] * cg[prodcglist[q-4][1]]
    #         I = simp(mp*(iii/m ))#- cgprod))
    #         print("I{} = {}".format(prodlist[q-4],I))
    #     else:
    #         ii0 = simp( V_ints[momlisti[q-7][0]].subs({ay0:r0,az0:r0,ay1:r1,az1:r1}) )
    #         ii1 = simp( V_ints[momlisti[q-7][1]].subs({ay0:r0,az0:r0,ay1:r1,az1:r1}) )
    #         cg0 = cg[momcglist[q-7][0]]
    #         cg1 = cg[momcglist[q-7][1]]
    #         I = simp(mp*( (ii0+ii1)/m ))#- cg0**2 - cg1**2) )
    #         print("I{} = {}".format(momlist[q-7],I))
        
    # print()
    # print()

    # print("cylinder eqs...")
    # m = simp( V_ints[0].subs({ay0:rp,ay1:rp,az0:rp,az1:rp}) )
    # cglist = ["x","y","z"]
    # prodlist = ["xy","xz","yz"]
    # prodcglist = [(0,1),(0,2),(1,2)]
    # momlist = ["xx","yy","zz"]
    # momlisti = [(8,9),(7,9),(7,8)]
    # momcglist = [(1,2),(0,2),(0,1)]
    # cg = []

    # for q in range(len(V)):
    #     iii = simp( V_ints[q].subs({ay0:rp,ay1:rp,az0:rp,az1:rp}) )
    #     if q == 0:
    #         print("V = {}".format(simp(m/p)))
    #         print("m = {}".format(m))
    #     elif q < 4:
    #         cg.append(simp(iii/m))
    #         print("{}-cg = {}".format(cglist[q-1],cg[-1]))
    #     elif q < 7:
    #         cgprod = cg[prodcglist[q-4][0]] * cg[prodcglist[q-4][1]]
    #         I = simp(mp*(iii/m - cgprod))
    #         print("I{} = {}".format(prodlist[q-4],I))
    #     else:
    #         ii0 = simp( V_ints[momlisti[q-7][0]].subs({ay0:rp,ay1:rp,az0:rp,az1:rp}) )
    #         ii1 = simp( V_ints[momlisti[q-7][1]].subs({ay0:rp,ay1:rp,az0:rp,az1:rp}) )
    #         cg0 = cg[momcglist[q-7][0]]
    #         cg1 = cg[momcglist[q-7][1]]
    #         I = simp(mp*( (ii0+ii1)/m - cg0**2 - cg1**2) )
    #         print("I{} = {}".format(momlist[q-7],I))
        
    # print()
    # print()

    # print("cone eqs, inertia about apex at origin...")
    # m = simp( V_ints[0].subs({ay0:0,az0:0,ay1:rp,az1:rp}) )
    # cglist = ["x","y","z"]
    # prodlist = ["xy","xz","yz"]
    # prodcglist = [(0,1),(0,2),(1,2)]
    # momlist = ["xx","yy","zz"]
    # momlisti = [(8,9),(7,9),(7,8)]
    # momcglist = [(1,2),(0,2),(0,1)]
    # cg = []

    # for q in range(len(V)):
    #     iii = simp( V_ints[q].subs({ay0:0,az0:0,ay1:rp,az1:rp}) )
    #     if q == 0:
    #         print("V = {}".format(simp(m/p)))
    #         print("m = {}".format(m))
    #     elif q < 4:
    #         cg.append(simp(iii/m))
    #         print("{}-cg = {}".format(cglist[q-1],cg[-1]))
    #     elif q < 7:
    #         cgprod = cg[prodcglist[q-4][0]] * cg[prodcglist[q-4][1]]
    #         I = simp(mp*(iii/m ))#- cgprod))
    #         print("I{} = {}".format(prodlist[q-4],I))
    #     else:
    #         ii0 = simp( V_ints[momlisti[q-7][0]].subs({ay0:0,az0:0,ay1:rp,az1:rp}) )
    #         ii1 = simp( V_ints[momlisti[q-7][1]].subs({ay0:0,az0:0,ay1:rp,az1:rp}) )
    #         cg0 = cg[momcglist[q-7][0]]
    #         cg1 = cg[momcglist[q-7][1]]
    #         I = simp(mp*( (ii0+ii1)/m ))#- cg0**2 - cg1**2) )
    #         print("I{} = {}".format(momlist[q-7],I))
        
    # print()
    # print()

    # print("cone eqs, inertia about cg location...")
    # m = simp( V_ints[0].subs({ay0:0,az0:0,ay1:rp,az1:rp}) )
    # cglist = ["x","y","z"]
    # prodlist = ["xy","xz","yz"]
    # prodcglist = [(0,1),(0,2),(1,2)]
    # momlist = ["xx","yy","zz"]
    # momlisti = [(8,9),(7,9),(7,8)]
    # momcglist = [(1,2),(0,2),(0,1)]
    # cg = []

    # for q in range(len(V)):
    #     iii = simp( V_ints[q].subs({ay0:0,az0:0,ay1:rp,az1:rp}) )
    #     if q == 0:
    #         print("V = {}".format(simp(m/p)))
    #         print("m = {}".format(m))
    #     elif q < 4:
    #         cg.append(simp(iii/m))
    #         print("{}-cg = {}".format(cglist[q-1],cg[-1]))
    #     elif q < 7:
    #         cgprod = cg[prodcglist[q-4][0]] * cg[prodcglist[q-4][1]]
    #         I = simp(mp*(iii/m - cgprod))
    #         print("I{} = {}".format(prodlist[q-4],I))
    #     else:
    #         ii0 = simp( V_ints[momlisti[q-7][0]].subs({ay0:0,az0:0,ay1:rp,az1:rp}) )
    #         ii1 = simp( V_ints[momlisti[q-7][1]].subs({ay0:0,az0:0,ay1:rp,az1:rp}) )
    #         cg0 = cg[momcglist[q-7][0]]
    #         cg1 = cg[momcglist[q-7][1]]
    #         I = simp(mp*( (ii0+ii1)/m - cg0**2 - cg1**2) )
    #         print("I{} = {}".format(momlist[q-7],I))
    
    # print()
    # print()

    print("tapered half-cylinder (+y) eqs...")
    theta_up =   pi / 2
    theta_lo = - pi / 2
    theta_bnd = (t,theta_lo,theta_up)
    H = V
    H_ints = []
    for q in range(len(H)):
        f = p * i**H[q][0] * j**H[q][1] * k**H[q][2]
        iii = simp( igr(ay * az * f * rp, r_bnd,x_bnd,theta_bnd) )
        H_ints.append(iii)
        print("H{}{}{} = {}".format(H[q][0],H[q][1],H[q][2],iii))
        
    print()
    print()

    print("tapered elliptic half-cylinder (+y) eqs...")
    m = simp(H_ints[0])
    cglist = ["x","y","z"]
    prodlist = ["xy","xz","yz"]
    prodcglist = [(0,1),(0,2),(1,2)]
    momlist = ["xx","yy","zz"]
    momlisti = [(8,9),(7,9),(7,8)]
    momcglist = [(1,2),(0,2),(0,1)]
    cg = []

    for q in range(len(H)):
        iii = simp(H_ints[q])
        if q == 0:
            print("V = {}".format(simp(m/p)))
            print("m = {}".format(m))
        elif q < 4:
            cg.append(simp(iii/m))
            print("{}-cg = {}".format(cglist[q-1],cg[-1]))
        elif q < 7:
            cgprod = cg[prodcglist[q-4][0]] * cg[prodcglist[q-4][1]]
            I = simp(mp*(iii/m - cgprod))
            print("I{}  = {}".format(prodlist[q-4],I))
            I = simp(mp*(iii/m))
            print("I{}o = {}".format(prodlist[q-4],I))
        else:
            ii0 = simp(H_ints[momlisti[q-7][0]])
            ii1 = simp(H_ints[momlisti[q-7][1]])
            cg0 = cg[momcglist[q-7][0]]
            cg1 = cg[momcglist[q-7][1]]
            I = simp(mp*( (ii0+ii1)/m - cg0**2 - cg1**2) )
            print("I{}  = {}".format(momlist[q-7],I))
            I = simp(mp*( (ii0+ii1)/m) )
            print("I{}o = {}".format(momlist[q-7],I))
        
    print()
    print()

    print("tapered quarter-cylinder (+y+z) eqs...")
    theta_up =   pi / 2
    theta_lo = 0
    theta_bnd = (t,theta_lo,theta_up)
    Q = V
    Q_ints = []
    for q in range(len(Q)):
        f = p * i**Q[q][0] * j**Q[q][1] * k**Q[q][2]
        iii = simp( igr(ay * az * f * rp, r_bnd,x_bnd,theta_bnd) )
        Q_ints.append(iii)
        print("Q{}{}{} = {}".format(Q[q][0],Q[q][1],Q[q][2],iii))
        
    print()
    print()

    print("tapered elliptic quarter-cylinder (+y+z) eqs...")
    m = simp(Q_ints[0])
    cglist = ["x","y","z"]
    prodlist = ["xy","xz","yz"]
    prodcglist = [(0,1),(0,2),(1,2)]
    momlist = ["xx","yy","zz"]
    momlisti = [(8,9),(7,9),(7,8)]
    momcglist = [(1,2),(0,2),(0,1)]
    cg = []

    for q in range(len(Q)):
        iii = simp(Q_ints[q])
        if q == 0:
            print("V = {}".format(simp(m/p)))
            print("m = {}".format(m))
        elif q < 4:
            cg.append(simp(iii/m))
            print("{}-cg = {}".format(cglist[q-1],cg[-1]))
        elif q < 7:
            cgprod = cg[prodcglist[q-4][0]] * cg[prodcglist[q-4][1]]
            I = simp(mp*(iii/m - cgprod))
            print("I{}  = {}".format(prodlist[q-4],I))
            I = simp(mp*(iii/m))
            print("I{}o = {}".format(prodlist[q-4],I))
        else:
            ii0 = simp(Q_ints[momlisti[q-7][0]])
            ii1 = simp(Q_ints[momlisti[q-7][1]])
            cg0 = cg[momcglist[q-7][0]]
            cg1 = cg[momcglist[q-7][1]]
            I = simp(mp*( (ii0+ii1)/m - cg0**2 - cg1**2) )
            print("I{}  = {}".format(momlist[q-7],I))
            I = simp(mp*( (ii0+ii1)/m) )
            print("I{}o = {}".format(momlist[q-7],I))
        
    print()
    print()
    