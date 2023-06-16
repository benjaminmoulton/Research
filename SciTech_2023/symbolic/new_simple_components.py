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
    r0 = sym("r0")
    r1 = sym("r1")
    t = sym("t")
    h = sym("h")
    rp = sym("rp")

    # create bounds
    print("creating bounds...")
    r = (r1-r0) * x + r0
    x_up = h; x_lo = 0
    theta_up = 2 * pi
    theta_lo = 0
    r_up = r; r_lo = 0
    x_bnd = (x,x_lo,x_up)
    r_bnd = (rp,r_lo,r_up)
    theta_bnd = (t,theta_lo,theta_up)

    i = x
    j = r * cos(t)
    k = r * sin(t)

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

    for q in range(len(V)):
        f = i**V[q][0] * j**V[q][1] * k**V[q][2]
        iii = simp( igr(f * r, r_bnd,x_bnd,theta_bnd) )
        V_ints.append(iii)
        print("T{}{}{} = {}".format(V[q][0],V[q][1],V[q][2],iii))
        
    print()
    print()