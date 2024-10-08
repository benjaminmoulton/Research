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

    # declare variables
    print("declaring variables...")
    x = sym("x")
    h = sym("h")
    yp = sym("y")
    zp = sym("z")
    p = sym("p")
    mp = sym("m")
    y0 = sym("y0")
    z0 = sym("z0")
    y1 = sym("y1")
    z1 = sym("z1")
    lx = sym("lx")
    ly = sym("ly")
    lz = sym("lz")    

    # create bounds
    print("creating bounds...")
    y = (y1-y0) * x/h + y0
    z = (z1-z0) * x/h + z0
    x_up = h; x_lo = 0
    y_up = y/2; y_lo = -y/2
    z_up = z/2; z_lo = -z/2
    x_bnd = ( x,x_lo,x_up)
    y_bnd = (yp,y_lo,y_up)
    z_bnd = (zp,z_lo,z_up)

    i = x
    j = yp
    k = zp

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

    print("tapered prismoid eqs...")

    for q in range(len(V)):
        f = p * i**V[q][0] * j**V[q][1] * k**V[q][2]
        iii = simp( igr(f, z_bnd,y_bnd,x_bnd) )
        V_ints.append(iii)
        print("V{}{}{} = {}".format(V[q][0],V[q][1],V[q][2],iii))
        
    print()
    print()

    print("tapered prismoid inertia eqs...")
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

    # print("prismoid eqs...")
    # m = simp( V_ints[0].subs({z1:z0,y1:y0}).subs({z0:lz,y0:ly,h:lx}) )
    # cglist = ["x","y","z"]
    # prodlist = ["xy","xz","yz"]
    # prodcglist = [(0,1),(0,2),(1,2)]
    # momlist = ["xx","yy","zz"]
    # momlisti = [(8,9),(7,9),(7,8)]
    # momcglist = [(1,2),(0,2),(0,1)]
    # cg = []

    # for q in range(len(V)):
    #     iii = simp( V_ints[q].subs({z1:z0,y1:y0}).subs({z0:lz,y0:ly,h:lx}) )
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
    #         ii0 = simp( V_ints[momlisti[q-7][0]].subs({z1:z0,y1:y0}).subs({z0:lz,y0:ly,h:lx}) )
    #         ii1 = simp( V_ints[momlisti[q-7][1]].subs({z1:z0,y1:y0}).subs({z0:lz,y0:ly,h:lx}) )
    #         cg0 = cg[momcglist[q-7][0]]
    #         cg1 = cg[momcglist[q-7][1]]
    #         I = simp(mp*( (ii0+ii1)/m - cg0**2 - cg1**2) )
    #         print("I{} = {}".format(momlist[q-7],I))
        
    # print()
    # print()

    # print("rectangular pyramid eqs...")
    # m = simp( V_ints[0].subs({z1:0,y1:0}).subs({z0:lz,y0:ly,h:lx}) )
    # cglist = ["x","y","z"]
    # prodlist = ["xy","xz","yz"]
    # prodcglist = [(0,1),(0,2),(1,2)]
    # momlist = ["xx","yy","zz"]
    # momlisti = [(8,9),(7,9),(7,8)]
    # momcglist = [(1,2),(0,2),(0,1)]
    # cg = []

    # for q in range(len(V)):
    #     iii = simp( V_ints[q].subs({z1:0,y1:0}).subs({z0:lz,y0:ly,h:lx}) )
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
    #         ii0 = simp( V_ints[momlisti[q-7][0]].subs({z1:0,y1:0}).subs({z0:lz,y0:ly,h:lx}) )
    #         ii1 = simp( V_ints[momlisti[q-7][1]].subs({z1:0,y1:0}).subs({z0:lz,y0:ly,h:lx}) )
    #         cg0 = cg[momcglist[q-7][0]]
    #         cg1 = cg[momcglist[q-7][1]]
    #         I = simp(mp*( (ii0+ii1)/m - cg0**2 - cg1**2) )
    #         print("I{} = {}".format(momlist[q-7],I))
        
    # print()
    # print()

    print("tapered half-prismoid (+y) eqs...")
    y_up = y/2; y_lo = 0
    # z_up = z/2; z_lo = -z/2
    y_bnd = (yp,y_lo,y_up)
    # z_bnd = (zp,z_lo,z_up)
    H = V
    H_ints = []
    for q in range(len(H)):
        f = p * i**H[q][0] * j**H[q][1] * k**H[q][2]
        iii = simp( igr(f, z_bnd,y_bnd,x_bnd) )
        H_ints.append(iii)
        print("H{}{}{} = {}".format(H[q][0],H[q][1],H[q][2],iii))
        
    print()
    print()

    print("tapered half-prismoid inertia (+y) eqs...")
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

    print("tapered quarter-prismoid (+y+z) eqs...")
    y_up = y/2; y_lo = 0
    z_up = z/2; z_lo = 0
    y_bnd = (yp,y_lo,y_up)
    z_bnd = (zp,z_lo,z_up)
    Q = V
    Q_ints = []
    for q in range(len(Q)):
        f = p * i**Q[q][0] * j**Q[q][1] * k**Q[q][2]
        iii = simp( igr(f, z_bnd,y_bnd,x_bnd) )
        Q_ints.append(iii)
        print("Q{}{}{} = {}".format(Q[q][0],Q[q][1],Q[q][2],iii))
        
    print()
    print()

    print("tapered quarter-prismoid inertia (+y+z) eqs...")
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
