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
    sqrt = sy.sqrt
    sin = sy.sin
    cos = sy.cos
    tan = sy.tan
    mat = sy.Matrix
    pi = sy.pi

    # declare variables
    print("declaring variables...")
    ax = sym("ax")
    ay = sym("ay")
    az = sym("az")
    t = sym("t")
    f = sym("f")
    rp = sym("r")
    p = sym("p")
    mp = sym("m")

    # create bounds
    print("creating bounds...")
    st = sin(t); ct = cos(t)
    sf = sin(f); cf = cos(f)
    r = sqrt( 1 - (ct*cf/ax)**2 - (ct*sf/ay)**2 - (st/az)**2 )
    f_up = 2 * pi; f_lo = 0
    theta_up = pi / 2
    theta_lo = - pi / 2
    r_up = r; r_lo = 0
    f_bnd = (f,f_lo,f_up)
    r_bnd = (rp,r_lo,r_up)
    theta_bnd = (t,theta_lo,theta_up)

    i = ax * cos(t) * cos(f)
    j = ay * cos(t) * sin(f)
    k = az * sin(t)

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

    print("ellipsoid eqs...")

    for q in range(len(V)):
        f = p * i**V[q][0] * j**V[q][1] * k**V[q][2]
        iii = simp( igr(f * rp**2 * sin(t), r_bnd) )
        print("i   = ", iii)
        iii = simp( igr(iii,theta_bnd) )
        print("ii  = ", iii)
        iii = simp( igr(iii,f_bnd) )
        print("iii = ", iii)
        V_ints.append(iii)
        print("V{}{}{} = {}".format(V[q][0],V[q][1],V[q][2],iii))
        
    print()
    print()

    print("ellipsoid inertia eqs...")
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
            print("I{} = {}".format(prodlist[q-4],I))
        else:
            ii0 = simp(V_ints[momlisti[q-7][0]])
            ii1 = simp(V_ints[momlisti[q-7][1]])
            cg0 = cg[momcglist[q-7][0]]
            cg1 = cg[momcglist[q-7][1]]
            I = simp(mp*( (ii0+ii1)/m - cg0**2 - cg1**2) )
            print("I{} = {}".format(momlist[q-7],I))
        
    print()
    print()

    print("sphere eqs...")
    m = simp(V_ints[0].replace(ax,rp).replace(ay,rp).replace(az,rp))
    cglist = ["x","y","z"]
    prodlist = ["xy","xz","yz"]
    prodcglist = [(0,1),(0,2),(1,2)]
    momlist = ["xx","yy","zz"]
    momlisti = [(8,9),(7,9),(7,8)]
    momcglist = [(1,2),(0,2),(0,1)]
    cg = []

    for q in range(len(V)):
        iii = simp(V_ints[q].replace(ax,rp).replace(ay,rp).replace(az,rp))
        if q == 0:
            print("V = {}".format(simp(m/p)))
            print("m = {}".format(m))
        elif q < 4:
            cg.append(simp(iii/m))
            print("{}-cg = {}".format(cglist[q-1],cg[-1]))
        elif q < 7:
            cgprod = cg[prodcglist[q-4][0]] * cg[prodcglist[q-4][1]]
            I = simp(mp*(iii/m - cgprod))
            print("I{} = {}".format(prodlist[q-4],I))
        else:
            ii0 = simp(V_ints[momlisti[q-7][0]].replace(ax,rp).replace(ay,rp).replace(az,rp))
            ii1 = simp(V_ints[momlisti[q-7][1]].replace(ax,rp).replace(ay,rp).replace(az,rp))
            cg0 = cg[momcglist[q-7][0]]
            cg1 = cg[momcglist[q-7][1]]
            I = simp(mp*( (ii0+ii1)/m - cg0**2 - cg1**2) )
            print("I{} = {}".format(momlist[q-7],I))
        
    print()
    print()