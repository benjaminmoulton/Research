import sympy as sy
import numpy as np
from matplotlib import pyplot as plt
import time as t
from organize import organize_equation

# current time
start_time = t.time()

sym = sy.Symbol
igr = sy.integrate
simp = sy.simplify
exp = sy.expand
piecewise = sy.Piecewise
diff = sy.diff
tan = sy.tan
# sy.init_printing(use_unicode=True)

# declare variables
print("declaring variables...")
x = sym("x")
y = sym("y")
z = sym("z")
xcg = sym("xcg")
ycg = sym("ycg")
zcg = sym("zcg")
p = sym("p")
b_2 = sym("b_2")
Lr = sym("Lr") # sym("\u039B")
tanL = tan(Lr)
ct = sym("ct")
cr = sym("cr")
tt = sym("tt")
tr = sym("tr")
ts = sym("ts")
twb = sym("twb")
Xle = sym("Xle")
Xte = sym("Xte")


# create bounds
print("creating bounds...")
x_up =   (sy.Rational(1,4) - Xle)*( (ct-cr)/b_2 * y + cr ) - ts - twb
x_lo = - (sy.Rational(3,4) - Xte)*( (ct-cr)/b_2 * y + cr ) + ts + twb
z_up =   sy.Rational(1,2) * ((tt*ct-tr*cr)/b_2 * y + tr*cr) - ts - twb
z_lo = - sy.Rational(1,2) * ((tt*ct-tr*cr)/b_2 * y + tr*cr) + ts + twb
y_up = b_2
y_lo = 0
x_bnd = (x,x_lo,x_up)
y_bnd = (y,y_lo,y_up)
z_bnd = (z,z_lo,z_up)

# c = sym("c"); t = sym("t")
# A = ((x_up - x_lo) * (z_up - z_lo)).subs([[y,0],[cr,c],[tr,t]])
# A_sout = A.subs([[ts,0],[twb,0],[Xle,0],[Xte,0]])
# A_sinn = A.subs([[twb,0],[Xle,0],[Xte,0]])
# A_wout = A.subs([[twb,0]])
# A_winn = A

# eq_A = simp(A_sout - A_sinn + A_wout - A_winn)
# print(eq_A)
# eq_B = -2*c*twb*(Xle+Xte) + 2*c*(ts+twb)*(1+t) -4*(ts+twb)**2
# print(sy.expand(eq_B))
# print(simp(eq_A-eq_B))
# quit()

# create simplification constants
ka = tt*ct*(2*ct+cr) + tr*cr*(ct+2*cr)
kb = tt*ct*(3*ct+cr) + tr*cr*(ct+cr)
kc = tt*ct*(3*ct**2+2*ct*cr+cr**2) + tr*cr*(ct**2+2*ct*cr+3*cr**2)
kd = 4 *  b_2 * ( tt*ct*(12*ct+3*cr) + tr*cr*(3*ct+2*cr) )
ke = 2*tt*ct*(6*ct**2+3*ct*cr+cr**2) + tr*cr*(3*ct**2+4*ct*cr+3*cr**2)
kf = 7*tt*ct*(4*ct**3+3*ct**2*cr+2*ct*cr**2+cr**3) + \
    7*tr*cr*(ct**3+2*ct**2*cr+3*ct*cr**2+4*cr**3)
kg = tt**3*ct**3*(4*ct+cr) + tt**2*ct**2*tr*cr*(3*ct+2*cr) + \
    tt*ct*tr**2*cr**2*(2*ct+3*cr) + tr**3*cr**3*(ct+4*cr)

# expand and multiply for convention
ka = exp(ka)
kb = exp(kb)
kc = exp(kc)
kd = exp(kd)
ke = exp(ke)
kf = exp(kf)
kg = exp(kg)
subexps = [(ka,"k_a"),(kb,"k_b"),(kc,"k_c"),(kd,"k_d"),(ke,"k_e"),(kf,"k_f"),
(kg,"k_g")]

def ksubs(expr):
    for i in range(len(subexps)):
        expr = simp(simp(expr).subs(subexps[0][0],subexps[0][1]))
        # expr = simp(exp(expr).subs(subexps[0][0],subexps[0][1]))
    
    return expr

# perform integrations
if False:
    print("integrating...")
    # mass
    print("\t mass...")
    m = igr(  p,  x_bnd,z_bnd,y_bnd)
    # moments
    print("\t moments...")
    myz = p*x
    mxz = p*y
    mxy = p*z
    Lyz = p*-tanL*y
    Lxz = p*0
    Lxy = p*0

    print("\t\t Myz...")
    Myz = igr(myz, x_bnd,z_bnd,y_bnd) + igr(Lyz, x_bnd,z_bnd,y_bnd)
    print("\t\t Mxz...")
    Mxz = igr(mxz, x_bnd,z_bnd,y_bnd) + igr(Lxz, x_bnd,z_bnd,y_bnd)
    print("\t\t Mxy...")
    Mxy = igr(mxy, x_bnd,z_bnd,y_bnd) + igr(Lxy, x_bnd,z_bnd,y_bnd)

    # sub inertias
    print("\t sub inertias...")
    ixx = igr( p*(y**2+z**2), x_bnd,z_bnd)
    iyy = igr( p*(x**2+z**2), x_bnd,z_bnd)
    izz = igr( p*(x**2+y**2), x_bnd,z_bnd)
    ixy = igr( p*(x*y), x_bnd,z_bnd)
    ixz = igr( p*(x*z), x_bnd,z_bnd)
    iyz = igr( p*(y*z), x_bnd,z_bnd)

    # new inertias
    print("\t full inertias...")
    print("\t\t Ixx...")
    Ixx_new = igr(ixx,y_bnd) - m * (ycg**2 + zcg**2)
    print("\t\t Iyy...")
    Iyy_new = igr(iyy,y_bnd) \
        - p*igr( -y**2*tanL**2 + 2*y*x*tanL, x_bnd,z_bnd,y_bnd) \
        - m * (xcg**2 + zcg**2)
    print("\t\t Izz...")
    Izz_new = igr(izz,y_bnd) \
        - p*igr( -y**2*tanL**2 + 2*y*x*tanL, x_bnd,z_bnd,y_bnd) \
        - m * (xcg**2 + ycg**2)
    print("\t\t Ixy...")
    Ixy_new = igr(ixy,y_bnd) - p*igr( y**2*tanL, x_bnd,z_bnd,y_bnd) \
        - m * (xcg * ycg)
    print("\t\t Ixz...")
    Ixz_new = igr(ixz,y_bnd)  - p*igr( y*z*tanL, x_bnd,z_bnd,y_bnd) \
        - m * (xcg * zcg)
    print("\t\t Iyz...")
    Iyz_new = (igr(iyz,y_bnd) / m - (ycg*zcg)) * m

    # simplify expressions
    print("skipping simplifying...")
    # print("simplifying...")
    # print("\t mass...")
    # m = simp(m)#ksubs(m))
    # print("\t moments...")
    # print("\t\t Myz...")
    # Myz = simp(Myz)
    # print("\t\t Mxz...")
    # Mxz = simp(Mxz)
    # print("\t\t Mxy...")
    # Mxy = simp(Mxy)
    # print("\t moments of inertia...")
    # print("\t\t Ixx...")
    # Ixx_new = simp(Ixx_new)
    # print("\t\t Iyy...")
    # Iyy_new = simp(Iyy_new)
    # print("\t\t Izz...")
    # Izz_new = simp(Izz_new)
    # print("\t products of inertia...")
    # print("\t\t Ixy...")
    # Ixy_new = simp(Ixy_new)
    # print("\t\t Ixz...")
    # Ixz_new = simp(Ixz_new)
    # print("\t\t Iyz...")
    # Iyz_new = simp(Iyz_new)

    # assign results
    # print("skipping reporting...")
    print("reporting...\n\n")
    print("m =", m, "\n")
    print("Myz =", Myz, "\n")
    print("Mxz =", Mxz, "\n")
    print("Mxy =", Mxy, "\n")

    print("Ixx_new =", Ixx_new, "\n")
    print("Iyy_new =", Iyy_new, "\n")
    print("Izz_new =", Izz_new, "\n")
    print("Ixy_new =", Ixy_new, "\n")
    print("Ixz_new =", Ixz_new, "\n")
    print("Iyz_new =", Iyz_new, "\n")

    print("Time elapse:",t.time()-start_time,"seconds\n")
else:
    print("skipped integrating...")

# p_val    = 0.5
# b_2_val  = 2.0
# Lr_val   = np.deg2rad(10.0)
# ct_val   = 1.0
# cr_val   = 1.5
# tt_val   = 0.12
# tr_val   = 0.12
# ts_val   = 0.125 / 12.0
# twb_val  = 0.25 / 12.0
# Xle_val = 0.3
# Xte_val = 0.4
# p_val    = 2.565855969
# b_2_val  = 1.285349055
# Lr_val   = np.deg2rad(11.81445049)
# ct_val   = 1.07501694
# cr_val   = 1.145833333
# tt_val   = 0.12
# tr_val   = 0.12
# ts_val   = 0.8 / 10.0 / 2.54 / 12.0
# twb_val  = 3.2 / 10.0 / 2.54 / 12.0
# Xle_val = 0.39 # from 0.39 to 0.49 x/c
# Xte_val = 0.51

# replace_vals = [
#     [b_2,b_2_val],
#     [cr,cr_val],
#     [ct,ct_val],
#     [tr,tr_val],
#     [tt,tt_val],
#     [Lr,Lr_val],
#     [p,p_val]
# ]

# ckanging_vals = [
# [   [ts,0.0],
#     [twb,0.0],
#     [Xle,0.0],
#     [Xte,0.0]],
# [   [ts,ts_val],
#     [twb,0.0],
#     [Xle,0.0],
#     [Xte,0.0]],
# [   [ts,ts_val],
#     [twb,0.0],
#     [Xle,Xle_val],
#     [Xte,Xte_val]],
# [   [ts,ts_val],
#     [twb,twb_val],
#     [Xle,Xle_val],
#     [Xte,Xte_val]]
# ]

# # initialize arrays
# m_arr   = np.zeros(4)
# Myz_arr = np.zeros(4)
# Mxz_arr = np.zeros(4)
# Mxy_arr = np.zeros(4)
# Ixx_arr = np.zeros(4)
# Iyy_arr = np.zeros(4)
# Izz_arr = np.zeros(4)
# Ixy_arr = np.zeros(4)
# Ixz_arr = np.zeros(4)
# Iyz_arr = np.zeros(4)

# for i in range(4):

#     # replace mass and moments
#     m_arr[i]   =   m.subs(replace_vals).subs(ckanging_vals[i])
#     Myz_arr[i] = Myz.subs(replace_vals).subs(ckanging_vals[i])
#     Mxz_arr[i] = Mxz.subs(replace_vals).subs(ckanging_vals[i])
#     Mxy_arr[i] = Mxy.subs(replace_vals).subs(ckanging_vals[i])

# # negate some
# index = [1,3]
# m_arr[index] *= -1.0
# Myz_arr[index] *= -1.0
# Mxz_arr[index] *= -1.0
# Mxy_arr[index] *= -1.0

# m       = np.sum(  m_arr)
# xcg_val = np.sum(Myz_arr) / m
# ycg_val = np.sum(Mxz_arr) / m
# zcg_val = np.sum(Mxy_arr) / m

# # affix cg location
# replace_vals2 = [
#     [xcg,xcg_val],
#     [ycg,ycg_val],
#     [zcg,zcg_val]
# ]

# # replace inertias
# for i in range(4):
#     Ixx_arr[i] = Ixx_new.subs(replace_vals).subs(ckanging_vals[i]).subs(replace_vals2)
#     Iyy_arr[i] = Iyy_new.subs(replace_vals).subs(ckanging_vals[i]).subs(replace_vals2)
#     Izz_arr[i] = Izz_new.subs(replace_vals).subs(ckanging_vals[i]).subs(replace_vals2)
#     Ixy_arr[i] = Ixy_new.subs(replace_vals).subs(ckanging_vals[i]).subs(replace_vals2)
#     Ixz_arr[i] = Ixz_new.subs(replace_vals).subs(ckanging_vals[i]).subs(replace_vals2)
#     Iyz_arr[i] = Iyz_new.subs(replace_vals).subs(ckanging_vals[i]).subs(replace_vals2)

# # negate some
# Ixy_arr *= -1.0
# Iyz_arr *= -1.0
# Ixx_arr[index] *= -1.0
# Iyy_arr[index] *= -1.0
# Izz_arr[index] *= -1.0
# Ixy_arr[index] *= -1.0
# Ixz_arr[index] *= -1.0
# Iyz_arr[index] *= -1.0

# Ixx_new = np.sum(Ixx_arr)
# Iyy_new = np.sum(Iyy_arr)
# Izz_new = np.sum(Izz_arr)
# Ixy_new = np.sum(Ixy_arr)
# Ixz_new = np.sum(Ixz_arr)
# Iyz_new = np.sum(Iyz_arr)

# # print(m_arr)
# # print(Myz_arr / m_arr)
# # print(Mxz_arr / m_arr)
# # print(Mxy_arr / m_arr)

# print("reporting with values...")
# print("m =", m)
# print("Myz =", xcg_val / m)
# print("Mxz =", ycg_val / m)
# print("Mxy =", zcg_val / m)
# print("xcg =", xcg_val)
# print("ycg =", ycg_val)
# print("zcg =", zcg_val)
# print("Ixx_new =", Ixx_new)
# print("Iyy_new =", Iyy_new)
# print("Izz_new =", Izz_new)
# print("Ixy_new =", Ixy_new)
# print("Ixz_new =", Ixz_new)
# print("Iyz_new =", Iyz_new)
# print()
# for inn in [Ixx_new,Ixy_new,Ixz_new]:
#     print("\t{:< 19.16f}".format(inn),end="")
# print()
# for inn in [Ixy_new,Iyy_new,Iyz_new]:
#     print("\t{:< 19.16f}".format(inn),end="")
# print()
# for inn in [Ixz_new,Iyz_new,Izz_new]:
#     print("\t{:< 19.16f}".format(inn),end="")
# print()

# Mass properties of shellwing
#      Configuration: Default
#      Coordinate system: -- default --

# Density = 16.08699420 pounds per cubic foot

# Mass = 1.51653435 pounds
# Mass = 0.0471353997 slugs

# Volume = 0.09427083 cubic feet

# Surface area = 14.68557315 square feet

# Center of mass: ( feet )
# 	X = -0.45582142
# 	Y = 0.92780847
# 	Z = 0.00000000

# Principal axes of inertia and principal moments of inertia: ( pounds * square feet )
# Taken at the center of mass.
# 	 Ix = (-0.17069455,  0.98532399,  0.00000000)   	Px = 0.16112601
# 	 Iy = (-0.98532399, -0.17069455,  0.00000000)   	Py = 0.51386477
# 	 Iz = ( 0.00000000,  0.00000000,  1.00000000)   	Pz = 0.66303218

# Moments of inertia: ( pounds * square feet )
# Taken at the center of mass and aligned with the output coordinate system. (Using positive tensor notation.)
	# Lxx = 0.50358715	Lxy = -0.05932693	Lxz = 0.00000000
	# Lyx = -0.05932693	Lyy = 0.17140362	Lyz = 0.00000000
	# Lzx = 0.00000000	Lzy = 0.00000000	Lzz = 0.66303218
# slugs * ft^2
#     Lxx = 0.015651990740	Lxy = 0.001843940138	Lxz = 0.00000000000
#     Lyx = 0.001843940138	Lyy = 0.005327395412	Lyz = 0.00000000000
#     Lzx = 0.000000000000	Lzy = 0.000000000000	Lzz = 0.02060770125

# Moments of inertia: ( pounds * square feet )
# Taken at the output coordinate system. (Using positive tensor notation.)
# 	Ixx = 1.80906323	Ixy = -0.70069202	Ixz = 0.00000000
# 	Iyx = -0.70069202	Iyy = 0.48649877	Iyz = 0.00000000
# 	Izx = 0.00000000	Izy = 0.00000000	Izz = 2.28360341
# slugs * ft^2
# Ixx = 0.05622748897	Ixy = 0.02177820663	Ixz = 0.00000000000
# Iyx = 0.02177820663	Iyy = 0.01512086685	Iyz = 0.00000000000
# Izx = 0.00000000000	Izy = 0.00000000000	Izz = 0.07097667091




p_a = tr*cr*(2*cr+ct) + tt*ct*(cr+2*ct)
X_0 = 1-Xle-Xte
X_1 = 1+2*Xle-2*Xte

h_ao = p_a*X_0 - 6*(cr+ct)*(ts + twb)*X_0 \
    - 6*(tr*cr+tt*ct)*(ts + twb) + 24*(ts + twb)**2
h_a = (tr*cr+tt*ct)*(ts + twb) + (cr+ct)*( ts + twb*X_0 ) - 4*(ts + twb)**2

m = p/6*b_2*( h_ao ) 
m_all = p*b_2*( h_a )


p_b  = tr*cr*(cr+ct) + tt*ct*(cr+3*ct)
p_c  = tr*cr*(3*cr**2+2*cr*ct+ct**2) + tt*ct*(cr**2+2*cr*ct+3*ct**2)
p_h  = cr**2 + cr*ct + ct**2
h_bo = p_b*X_0 + 4*(ts + twb)*( -X_0*(cr + 2*ct) - (tr*cr+2*tt*ct) \
    + 6*(ts + twb))
h_b  = (ts + twb)*(tr*cr+2*tt*ct) - 6*(ts + twb)**2 + ts*(cr + 2*ct) \
    + twb*X_0*(cr + 2*ct)
h_co = X_1*( X_0*( p_c - 8*p_h*(ts + twb) ) \
    + 4*(ts + twb)*(6*(cr + ct)*(ts + twb) - p_a) )
h_c  = ts*(2*p_h + p_a - 6*ts*(cr + ct)) - twb*X_1*( 6*(cr + ct)*(2*ts + twb) \
    - p_a - 2*p_h*X_0 )

Myz = -p/48*b_2*( 4*b_2*tan(Lr)*h_bo + h_co )
Myz_all = -p/12*b_2*( 4*b_2*tan(Lr)*h_b + h_c )


# - delta in front for left / righ_t wings (delta is - for righ_t wings)
Mxz = p/12*b_2**2*h_bo
Mxz_all = p/3*b_2**2*( h_b )


Mxy = 0
Mxy_all = 0


p_g  = tr**3*cr**3*(4*cr+ct) + tr**2*cr**2*tt*ct*(3*cr+2*ct) \
    + tr*cr*tt**2*ct**2*(2*cr+3*ct) + tt**3*ct**3*(cr+4*ct)
p_d  = tr*cr*(2*cr+3*ct) + tt*ct*(3*cr+12*ct)
p_i  = tr**2*cr**2*(3*cr+ct) + tr*cr*tt*ct*(2*cr+2*ct) + tt**2*ct**2*(cr+3*ct)
p_j  = tr**3*cr**3 + tr**2*cr**2*tt*ct + tr*cr*tt**2*ct**2 + tt**3*ct**3
p_k  = tr**2*cr**2 + tr*cr*tt*ct + tt**2*ct**2
h_do = p_d*X_0 + 10*(ts+twb)*(8*(ts+twb) - (cr+3*ct)*X_0 - (tr*cr+3*tt*ct))
h_d  = (cr+3*ct)*( ts + twb*X_0 ) + (ts+twb)*((tr*cr+3*tt*ct) - 8*(ts+twb) )
h_go = p_g*X_0 - 10*(ts+twb)*( p_i*X_0 + p_j ) \
     + 40*(ts+twb)**2*( 2*p_k + p_a*X_0 ) \
     - 80*(ts+twb)**3*( 3*(tr*cr+tt*ct) + (cr+ct)*X_0 ) + 320*(ts+twb)**4
h_g  =  + ts*( p_i + p_j) + twb*( p_i*X_0 + p_j) \
     - 4*p_a*ts**2*(Xle+Xte) - 4*(ts + twb)**2*( p_a*X_0 + 2*p_k) \
     + 8*(cr+ct)*ts**3*(Xle+Xte) \
     + 8*(ts + twb)**3*( (cr+ct)*X_0 + 3*(tr*cr+tt*ct)) \
     - 32*(ts + twb)**4

Ixx_new = p/240*b_2*( 4*b_2**2*h_do + h_go ) - (ycg**2+zcg**2)*p/6*b_2*h_ao

Ixx_all = p/24*b_2*(  4*b_2**2*h_d  + h_g  ) - (ycg**2+zcg**2)*p  *b_2*h_a


X_2 = 8*Xle**2 - 4*Xle + 8*Xte**2 - 12*Xte + 5
X_3 = X_0 * (16*Xle**2 - 16*Xle*Xte + 4*Xle + 16*Xte**2 - 20*Xte + 7)

p_f  = tr*cr*(4*cr**3+3*cr**2*ct+2*cr*ct**2+ct**3) \
    + tt*ct*(cr**3+2*cr**2*ct+3*cr*ct**2+4*ct**3)
p_l  = cr**3 + cr**2*ct + cr*ct**2 + ct**3
p_e  = tr*cr*(3*cr**2+4*cr*ct+3*ct**2) + 2*tt*ct*(cr**2+3*cr*ct+6*ct**2)
p_m  = cr**2 + 2*cr*ct + 3*ct**2
h_eo = X_1*( + p_e*X_0 - 10*p_m*X_0*(ts+twb) - 10*p_b*(ts+twb) \
    + 40*(cr+2*ct)*(ts+twb)**2 )
h_e = X_1*( + (p_m*X_0+p_b)*twb \
        - 4*(cr+2*ct)*(2*ts*twb + twb**2) ) \
        + (p_m+p_b)*ts - 4*(cr+2*ct)*ts**2
h_fo =  + p_f*X_3 - 10*(ts+twb)*( + p_l*X_3 + p_c*X_2 ) \
     + 40*(ts+twb)**2*( 2*p_h*X_2 + 4*p_a*X_0 ) \
     - 320*(ts+twb)**3*( + 3*(cr+ct)*X_0 + (tr*cr+tt*ct) ) \
     + 1280*(ts+twb)**4
h_f = + p_l*(7*ts+twb*X_3) \
     + p_c*(5*ts+twb*X_2) \
     + 8*ts**2*(+ p_h*(X_2-5) + 2*p_a*(X_0-1) ) \
     - 8*(ts + twb)**2*( + p_h*X_2 + 2*p_a*X_0) \
     + 96*ts**3*(cr+ct)*(1-X_0) \
     + 32*(ts + twb)**3*( + 3*(cr+ct)*X_0 + (tr*cr+tt*ct)) \
     - 128*(ts + twb)**4


Iyy_new = p/960*b_2*( + 4*b_2*( 4*b_2*tan(Lr)**2*h_do + 2*tan(Lr)*h_eo ) \
     + h_fo + 4*h_go ) - (xcg**2+zcg**2)*p/6*b_2*h_ao

Iyy_all = p/96 *b_2*( + 4*b_2*( 4*b_2*tan(Lr)**2*h_d  + 2*tan(Lr)*h_e  ) \
     + h_f  + 4*h_g  ) - (xcg**2+zcg**2)*p  *b_2*h_a


Izz_new = p/960*b_2*( + 4*b_2*( + 4*b_2*(tan(Lr)**2 + 1)*h_do + 2*tan(Lr)*h_eo ) \
     + h_fo )  - (xcg**2+ycg**2)*p/6*b_2*( h_ao )
# something wrong with the way organize_equation handles Izz (one coeff off by factor)

Izz_all = p/96 *b_2*( + 4*b_2*( + 4*b_2*(tan(Lr)**2 + 1)*h_d  + 2*tan(Lr)*h_e  ) \
     + h_f  )  - (xcg**2+ycg**2)*p  *b_2*h_a



# print(simp(sy.expand(h_do)-sy.expand(h_do2)))

Ixy_new = -p/240*b_2*(  \
     + b_2*( + 4*b_2*tan(Lr)*h_do + h_eo ) \
     ) - xcg*ycg*p/6*b_2*( h_ao )

Ixy_all = -p/24*b_2*( \
     + b_2*( + 4*b_2*tan(Lr)*h_d + h_e ) \
     ) - xcg*ycg*p*b_2*( h_a )


# Ixz = -ycg*zcg*m
Ixz_new = -p/6*b_2*xcg*zcg*( h_ao ) 
Ixz_all = -p*b_2*xcg*zcg*( h_a )


# Iyz = -ycg*zcg*m
Iyz_new = -p/6*b_2*ycg*zcg*( h_ao ) 
Iyz_all = -p*b_2*ycg*zcg*( h_a )

# kvals = [ka,kb,kc,kd1,ke,kf1,kg]
# hvals = [ka,ha,kb,kc,he,kg,kd,hh,hi,hj,kf,hl,ke,hn]
# for i in range(len(kvals)):
#     for j in range(len(hvals)):
#         if sy.expand(kvals[i]) == sy.expand(hvals[j]):
#             print("h",chr(ord("a") + j),"==","k",chr(ord("a") + i))


# m = -b_2*p*(2*cr**2*Xle*tr + 2*cr**2*Xte*tr - 2*cr**2*tr + cr*ct*Xle*tr + cr*ct*Xle*tt + cr*ct*Xte*tr + cr*ct*Xte*tt - cr*ct*tr - cr*ct*tt - 6*cr*ts*Xle - 6*cr*ts*Xte + 6*cr*ts*tr + 6*cr*ts - 6*cr*Xle*twb - 6*cr*Xte*twb + 6*cr*twb*tr + 6*cr*twb + 2*ct**2*Xle*tt + 2*ct**2*Xte*tt - 2*ct**2*tt - 6*ct*ts*Xle - 6*ct*ts*Xte + 6*ct*ts*tt + 6*ct*ts - 6*ct*Xle*twb - 6*ct*Xte*twb + 6*ct*twb*tt + 6*ct*twb - 24*ts**2 - 48*ts*twb - 24*twb**2)/6 

# Myz = b_2*p*(4*b_2*cr**2*Xle*tr*tan(Lr) + 4*b_2*cr**2*Xte*tr*tan(Lr) - 4*b_2*cr**2*tr*tan(Lr) + 4*b_2*cr*ct*Xle*tr*tan(Lr) + 4*b_2*cr*ct*Xle*tt*tan(Lr) + 4*b_2*cr*ct*Xte*tr*tan(Lr) + 4*b_2*cr*ct*Xte*tt*tan(Lr) - 4*b_2*cr*ct*tr*tan(Lr) - 4*b_2*cr*ct*tt*tan(Lr) - 16*b_2*cr*ts*Xle*tan(Lr) - 16*b_2*cr*ts*Xte*tan(Lr) + 16*b_2*cr*ts*tr*tan(Lr) + 16*b_2*cr*ts*tan(Lr) - 16*b_2*cr*Xle*twb*tan(Lr) - 16*b_2*cr*Xte*twb*tan(Lr) + 16*b_2*cr*twb*tr*tan(Lr) + 16*b_2*cr*twb*tan(Lr) + 12*b_2*ct**2*Xle*tt*tan(Lr) + 12*b_2*ct**2*Xte*tt*tan(Lr) - 12*b_2*ct**2*tt*tan(Lr) - 32*b_2*ct*ts*Xle*tan(Lr) - 32*b_2*ct*ts*Xte*tan(Lr) + 32*b_2*ct*ts*tt*tan(Lr) + 32*b_2*ct*ts*tan(Lr) - 32*b_2*ct*Xle*twb*tan(Lr) - 32*b_2*ct*Xte*twb*tan(Lr) + 32*b_2*ct*twb*tt*tan(Lr) + 32*b_2*ct*twb*tan(Lr) - 96*b_2*ts**2*tan(Lr) - 192*b_2*ts*twb*tan(Lr) - 96*b_2*twb**2*tan(Lr) + 6*cr**3*Xle**2*tr - 3*cr**3*Xle*tr - 6*cr**3*Xte**2*tr + 9*cr**3*Xte*tr - 3*cr**3*tr + 4*cr**2*ct*Xle**2*tr + 2*cr**2*ct*Xle**2*tt - 2*cr**2*ct*Xle*tr - cr**2*ct*Xle*tt - 4*cr**2*ct*Xte**2*tr - 2*cr**2*ct*Xte**2*tt + 6*cr**2*ct*Xte*tr + 3*cr**2*ct*Xte*tt - 2*cr**2*ct*tr - cr**2*ct*tt - 16*cr**2*ts*Xle**2 + 16*cr**2*ts*Xle*tr + 8*cr**2*ts*Xle + 16*cr**2*ts*Xte**2 - 16*cr**2*ts*Xte*tr - 24*cr**2*ts*Xte + 8*cr**2*ts*tr + 8*cr**2*ts - 16*cr**2*Xle**2*twb + 16*cr**2*Xle*twb*tr + 8*cr**2*Xle*twb + 16*cr**2*Xte**2*twb - 16*cr**2*Xte*twb*tr - 24*cr**2*Xte*twb + 8*cr**2*twb*tr + 8*cr**2*twb + 2*cr*ct**2*Xle**2*tr + 4*cr*ct**2*Xle**2*tt - cr*ct**2*Xle*tr - 2*cr*ct**2*Xle*tt - 2*cr*ct**2*Xte**2*tr - 4*cr*ct**2*Xte**2*tt + 3*cr*ct**2*Xte*tr + 6*cr*ct**2*Xte*tt - cr*ct**2*tr - 2*cr*ct**2*tt - 16*cr*ct*ts*Xle**2 + 8*cr*ct*ts*Xle*tr + 8*cr*ct*ts*Xle*tt + 8*cr*ct*ts*Xle + 16*cr*ct*ts*Xte**2 - 8*cr*ct*ts*Xte*tr - 8*cr*ct*ts*Xte*tt - 24*cr*ct*ts*Xte + 4*cr*ct*ts*tr + 4*cr*ct*ts*tt + 8*cr*ct*ts - 16*cr*ct*Xle**2*twb + 8*cr*ct*Xle*twb*tr + 8*cr*ct*Xle*twb*tt + 8*cr*ct*Xle*twb + 16*cr*ct*Xte**2*twb - 8*cr*ct*Xte*twb*tr - 8*cr*ct*Xte*twb*tt - 24*cr*ct*Xte*twb + 4*cr*ct*twb*tr + 4*cr*ct*twb*tt + 8*cr*ct*twb - 48*cr*ts**2*Xle + 48*cr*ts**2*Xte - 24*cr*ts**2 - 96*cr*ts*Xle*twb + 96*cr*ts*Xte*twb - 48*cr*ts*twb - 48*cr*Xle*twb**2 + 48*cr*Xte*twb**2 - 24*cr*twb**2 + 6*ct**3*Xle**2*tt - 3*ct**3*Xle*tt - 6*ct**3*Xte**2*tt + 9*ct**3*Xte*tt - 3*ct**3*tt - 16*ct**2*ts*Xle**2 + 16*ct**2*ts*Xle*tt + 8*ct**2*ts*Xle + 16*ct**2*ts*Xte**2 - 16*ct**2*ts*Xte*tt - 24*ct**2*ts*Xte + 8*ct**2*ts*tt + 8*ct**2*ts - 16*ct**2*Xle**2*twb + 16*ct**2*Xle*twb*tt + 8*ct**2*Xle*twb + 16*ct**2*Xte**2*twb - 16*ct**2*Xte*twb*tt - 24*ct**2*Xte*twb + 8*ct**2*twb*tt + 8*ct**2*twb - 48*ct*ts**2*Xle + 48*ct*ts**2*Xte - 24*ct*ts**2 - 96*ct*ts*Xle*twb + 96*ct*ts*Xte*twb - 48*ct*ts*twb - 48*ct*Xle*twb**2 + 48*ct*Xte*twb**2 - 24*ct*twb**2)/48 

# Mxz = -b_2**2*p*(cr**2*Xle*tr + cr**2*Xte*tr - cr**2*tr + cr*ct*Xle*tr + cr*ct*Xle*tt + cr*ct*Xte*tr + cr*ct*Xte*tt - cr*ct*tr - cr*ct*tt - 4*cr*ts*Xle - 4*cr*ts*Xte + 4*cr*ts*tr + 4*cr*ts - 4*cr*Xle*twb - 4*cr*Xte*twb + 4*cr*twb*tr + 4*cr*twb + 3*ct**2*Xle*tt + 3*ct**2*Xte*tt - 3*ct**2*tt - 8*ct*ts*Xle - 8*ct*ts*Xte + 8*ct*ts*tt + 8*ct*ts - 8*ct*Xle*twb - 8*ct*Xte*twb + 8*ct*twb*tt + 8*ct*twb - 24*ts**2 - 48*ts*twb - 24*twb**2)/12 

# Mxy = 0 

# Ixx_new = -b_2*p*(8*b_2**2*cr**2*Xle*tr + 8*b_2**2*cr**2*Xte*tr - 8*b_2**2*cr**2*tr + 12*b_2**2*cr*ct*Xle*tr + 12*b_2**2*cr*ct*Xle*tt + 12*b_2**2*cr*ct*Xte*tr + 12*b_2**2*cr*ct*Xte*tt - 12*b_2**2*cr*ct*tr - 12*b_2**2*cr*ct*tt - 40*b_2**2*cr*ts*Xle - 40*b_2**2*cr*ts*Xte + 40*b_2**2*cr*ts*tr + 40*b_2**2*cr*ts - 40*b_2**2*cr*Xle*twb - 40*b_2**2*cr*Xte*twb + 40*b_2**2*cr*twb*tr + 40*b_2**2*cr*twb + 48*b_2**2*ct**2*Xle*tt + 48*b_2**2*ct**2*Xte*tt - 48*b_2**2*ct**2*tt - 120*b_2**2*ct*ts*Xle - 120*b_2**2*ct*ts*Xte + 120*b_2**2*ct*ts*tt + 120*b_2**2*ct*ts - 120*b_2**2*ct*Xle*twb - 120*b_2**2*ct*Xte*twb + 120*b_2**2*ct*twb*tt + 120*b_2**2*ct*twb - 320*b_2**2*ts**2 - 640*b_2**2*ts*twb - 320*b_2**2*twb**2 + 4*cr**4*Xle*tr**3 + 4*cr**4*Xte*tr**3 - 4*cr**4*tr**3 + cr**3*ct*Xle*tr**3 + 3*cr**3*ct*Xle*tr**2*tt + cr**3*ct*Xte*tr**3 + 3*cr**3*ct*Xte*tr**2*tt - cr**3*ct*tr**3 - 3*cr**3*ct*tr**2*tt - 30*cr**3*ts*Xle*tr**2 - 30*cr**3*ts*Xte*tr**2 + 10*cr**3*ts*tr**3 + 30*cr**3*ts*tr**2 - 30*cr**3*Xle*twb*tr**2 - 30*cr**3*Xte*twb*tr**2 + 10*cr**3*twb*tr**3 + 30*cr**3*twb*tr**2 + 2*cr**2*ct**2*Xle*tr**2*tt + 2*cr**2*ct**2*Xle*tr*tt**2 + 2*cr**2*ct**2*Xte*tr**2*tt + 2*cr**2*ct**2*Xte*tr*tt**2 - 2*cr**2*ct**2*tr**2*tt - 2*cr**2*ct**2*tr*tt**2 - 10*cr**2*ct*ts*Xle*tr**2 - 20*cr**2*ct*ts*Xle*tr*tt - 10*cr**2*ct*ts*Xte*tr**2 - 20*cr**2*ct*ts*Xte*tr*tt + 10*cr**2*ct*ts*tr**2*tt + 10*cr**2*ct*ts*tr**2 + 20*cr**2*ct*ts*tr*tt - 10*cr**2*ct*Xle*twb*tr**2 - 20*cr**2*ct*Xle*twb*tr*tt - 10*cr**2*ct*Xte*twb*tr**2 - 20*cr**2*ct*Xte*twb*tr*tt + 10*cr**2*ct*twb*tr**2*tt + 10*cr**2*ct*twb*tr**2 + 20*cr**2*ct*twb*tr*tt + 80*cr**2*ts**2*Xle*tr + 80*cr**2*ts**2*Xte*tr - 80*cr**2*ts**2*tr**2 - 80*cr**2*ts**2*tr + 160*cr**2*ts*Xle*twb*tr + 160*cr**2*ts*Xte*twb*tr - 160*cr**2*ts*twb*tr**2 - 160*cr**2*ts*twb*tr + 80*cr**2*Xle*twb**2*tr - 80*cr**2*Xle*tr*ycg**2 - 80*cr**2*Xle*tr*zcg**2 + 80*cr**2*Xte*twb**2*tr - 80*cr**2*Xte*tr*ycg**2 - 80*cr**2*Xte*tr*zcg**2 - 80*cr**2*twb**2*tr**2 - 80*cr**2*twb**2*tr + 80*cr**2*tr*ycg**2 + 80*cr**2*tr*zcg**2 + 3*cr*ct**3*Xle*tr*tt**2 + cr*ct**3*Xle*tt**3 + 3*cr*ct**3*Xte*tr*tt**2 + cr*ct**3*Xte*tt**3 - 3*cr*ct**3*tr*tt**2 - cr*ct**3*tt**3 - 20*cr*ct**2*ts*Xle*tr*tt - 10*cr*ct**2*ts*Xle*tt**2 - 20*cr*ct**2*ts*Xte*tr*tt - 10*cr*ct**2*ts*Xte*tt**2 + 10*cr*ct**2*ts*tr*tt**2 + 20*cr*ct**2*ts*tr*tt + 10*cr*ct**2*ts*tt**2 - 20*cr*ct**2*Xle*twb*tr*tt - 10*cr*ct**2*Xle*twb*tt**2 - 20*cr*ct**2*Xte*twb*tr*tt - 10*cr*ct**2*Xte*twb*tt**2 + 10*cr*ct**2*twb*tr*tt**2 + 20*cr*ct**2*twb*tr*tt + 10*cr*ct**2*twb*tt**2 + 40*cr*ct*ts**2*Xle*tr + 40*cr*ct*ts**2*Xle*tt + 40*cr*ct*ts**2*Xte*tr + 40*cr*ct*ts**2*Xte*tt - 80*cr*ct*ts**2*tr*tt - 40*cr*ct*ts**2*tr - 40*cr*ct*ts**2*tt + 80*cr*ct*ts*Xle*twb*tr + 80*cr*ct*ts*Xle*twb*tt + 80*cr*ct*ts*Xte*twb*tr + 80*cr*ct*ts*Xte*twb*tt - 160*cr*ct*ts*twb*tr*tt - 80*cr*ct*ts*twb*tr - 80*cr*ct*ts*twb*tt + 40*cr*ct*Xle*twb**2*tr + 40*cr*ct*Xle*twb**2*tt - 40*cr*ct*Xle*tr*ycg**2 - 40*cr*ct*Xle*tr*zcg**2 - 40*cr*ct*Xle*tt*ycg**2 - 40*cr*ct*Xle*tt*zcg**2 + 40*cr*ct*Xte*twb**2*tr + 40*cr*ct*Xte*twb**2*tt - 40*cr*ct*Xte*tr*ycg**2 - 40*cr*ct*Xte*tr*zcg**2 - 40*cr*ct*Xte*tt*ycg**2 - 40*cr*ct*Xte*tt*zcg**2 - 80*cr*ct*twb**2*tr*tt - 40*cr*ct*twb**2*tr - 40*cr*ct*twb**2*tt + 40*cr*ct*tr*ycg**2 + 40*cr*ct*tr*zcg**2 + 40*cr*ct*tt*ycg**2 + 40*cr*ct*tt*zcg**2 - 80*cr*ts**3*Xle - 80*cr*ts**3*Xte + 240*cr*ts**3*tr + 80*cr*ts**3 - 240*cr*ts**2*Xle*twb - 240*cr*ts**2*Xte*twb + 720*cr*ts**2*twb*tr + 240*cr*ts**2*twb - 240*cr*ts*Xle*twb**2 + 240*cr*ts*Xle*ycg**2 + 240*cr*ts*Xle*zcg**2 - 240*cr*ts*Xte*twb**2 + 240*cr*ts*Xte*ycg**2 + 240*cr*ts*Xte*zcg**2 + 720*cr*ts*twb**2*tr + 240*cr*ts*twb**2 - 240*cr*ts*tr*ycg**2 - 240*cr*ts*tr*zcg**2 - 240*cr*ts*ycg**2 - 240*cr*ts*zcg**2 - 80*cr*Xle*twb**3 + 240*cr*Xle*twb*ycg**2 + 240*cr*Xle*twb*zcg**2 - 80*cr*Xte*twb**3 + 240*cr*Xte*twb*ycg**2 + 240*cr*Xte*twb*zcg**2 + 240*cr*twb**3*tr + 80*cr*twb**3 - 240*cr*twb*tr*ycg**2 - 240*cr*twb*tr*zcg**2 - 240*cr*twb*ycg**2 - 240*cr*twb*zcg**2 + 4*ct**4*Xle*tt**3 + 4*ct**4*Xte*tt**3 - 4*ct**4*tt**3 - 30*ct**3*ts*Xle*tt**2 - 30*ct**3*ts*Xte*tt**2 + 10*ct**3*ts*tt**3 + 30*ct**3*ts*tt**2 - 30*ct**3*Xle*twb*tt**2 - 30*ct**3*Xte*twb*tt**2 + 10*ct**3*twb*tt**3 + 30*ct**3*twb*tt**2 + 80*ct**2*ts**2*Xle*tt + 80*ct**2*ts**2*Xte*tt - 80*ct**2*ts**2*tt**2 - 80*ct**2*ts**2*tt + 160*ct**2*ts*Xle*twb*tt + 160*ct**2*ts*Xte*twb*tt - 160*ct**2*ts*twb*tt**2 - 160*ct**2*ts*twb*tt + 80*ct**2*Xle*twb**2*tt - 80*ct**2*Xle*tt*ycg**2 - 80*ct**2*Xle*tt*zcg**2 + 80*ct**2*Xte*twb**2*tt - 80*ct**2*Xte*tt*ycg**2 - 80*ct**2*Xte*tt*zcg**2 - 80*ct**2*twb**2*tt**2 - 80*ct**2*twb**2*tt + 80*ct**2*tt*ycg**2 + 80*ct**2*tt*zcg**2 - 80*ct*ts**3*Xle - 80*ct*ts**3*Xte + 240*ct*ts**3*tt + 80*ct*ts**3 - 240*ct*ts**2*Xle*twb - 240*ct*ts**2*Xte*twb + 720*ct*ts**2*twb*tt + 240*ct*ts**2*twb - 240*ct*ts*Xle*twb**2 + 240*ct*ts*Xle*ycg**2 + 240*ct*ts*Xle*zcg**2 - 240*ct*ts*Xte*twb**2 + 240*ct*ts*Xte*ycg**2 + 240*ct*ts*Xte*zcg**2 + 720*ct*ts*twb**2*tt + 240*ct*ts*twb**2 - 240*ct*ts*tt*ycg**2 - 240*ct*ts*tt*zcg**2 - 240*ct*ts*ycg**2 - 240*ct*ts*zcg**2 - 80*ct*Xle*twb**3 + 240*ct*Xle*twb*ycg**2 + 240*ct*Xle*twb*zcg**2 - 80*ct*Xte*twb**3 + 240*ct*Xte*twb*ycg**2 + 240*ct*Xte*twb*zcg**2 + 240*ct*twb**3*tt + 80*ct*twb**3 - 240*ct*twb*tt*ycg**2 - 240*ct*twb*tt*zcg**2 - 240*ct*twb*ycg**2 - 240*ct*twb*zcg**2 - 320*ts**4 - 1280*ts**3*twb - 1920*ts**2*twb**2 + 960*ts**2*ycg**2 + 960*ts**2*zcg**2 - 1280*ts*twb**3 + 1920*ts*twb*ycg**2 + 1920*ts*twb*zcg**2 - 320*twb**4 + 960*twb**2*ycg**2 + 960*twb**2*zcg**2)/240 

# Iyy_new = -b_2*p*(32*b_2**2*cr**2*Xle*tr*tan(Lr)**2 + 32*b_2**2*cr**2*Xte*tr*tan(Lr)**2 - 32*b_2**2*cr**2*tr*tan(Lr)**2 + 48*b_2**2*cr*ct*Xle*tr*tan(Lr)**2 + 48*b_2**2*cr*ct*Xle*tt*tan(Lr)**2 + 48*b_2**2*cr*ct*Xte*tr*tan(Lr)**2 + 48*b_2**2*cr*ct*Xte*tt*tan(Lr)**2 - 48*b_2**2*cr*ct*tr*tan(Lr)**2 - 48*b_2**2*cr*ct*tt*tan(Lr)**2 - 160*b_2**2*cr*ts*Xle*tan(Lr)**2 - 160*b_2**2*cr*ts*Xte*tan(Lr)**2 + 160*b_2**2*cr*ts*tr*tan(Lr)**2 + 160*b_2**2*cr*ts*tan(Lr)**2 - 160*b_2**2*cr*Xle*twb*tan(Lr)**2 - 160*b_2**2*cr*Xte*twb*tan(Lr)**2 + 160*b_2**2*cr*twb*tr*tan(Lr)**2 + 160*b_2**2*cr*twb*tan(Lr)**2 + 192*b_2**2*ct**2*Xle*tt*tan(Lr)**2 + 192*b_2**2*ct**2*Xte*tt*tan(Lr)**2 - 192*b_2**2*ct**2*tt*tan(Lr)**2 - 480*b_2**2*ct*ts*Xle*tan(Lr)**2 - 480*b_2**2*ct*ts*Xte*tan(Lr)**2 + 480*b_2**2*ct*ts*tt*tan(Lr)**2 + 480*b_2**2*ct*ts*tan(Lr)**2 - 480*b_2**2*ct*Xle*twb*tan(Lr)**2 - 480*b_2**2*ct*Xte*twb*tan(Lr)**2 + 480*b_2**2*ct*twb*tt*tan(Lr)**2 + 480*b_2**2*ct*twb*tan(Lr)**2 - 1280*b_2**2*ts**2*tan(Lr)**2 - 2560*b_2**2*ts*twb*tan(Lr)**2 - 1280*b_2**2*twb**2*tan(Lr)**2 + 48*b_2*cr**3*Xle**2*tr*tan(Lr) - 24*b_2*cr**3*Xle*tr*tan(Lr) - 48*b_2*cr**3*Xte**2*tr*tan(Lr) + 72*b_2*cr**3*Xte*tr*tan(Lr) - 24*b_2*cr**3*tr*tan(Lr) + 64*b_2*cr**2*ct*Xle**2*tr*tan(Lr) + 32*b_2*cr**2*ct*Xle**2*tt*tan(Lr) - 32*b_2*cr**2*ct*Xle*tr*tan(Lr) - 16*b_2*cr**2*ct*Xle*tt*tan(Lr) - 64*b_2*cr**2*ct*Xte**2*tr*tan(Lr) - 32*b_2*cr**2*ct*Xte**2*tt*tan(Lr) + 96*b_2*cr**2*ct*Xte*tr*tan(Lr) + 48*b_2*cr**2*ct*Xte*tt*tan(Lr) - 32*b_2*cr**2*ct*tr*tan(Lr) - 16*b_2*cr**2*ct*tt*tan(Lr) - 160*b_2*cr**2*ts*Xle**2*tan(Lr) + 160*b_2*cr**2*ts*Xle*tr*tan(Lr) + 80*b_2*cr**2*ts*Xle*tan(Lr) + 160*b_2*cr**2*ts*Xte**2*tan(Lr) - 160*b_2*cr**2*ts*Xte*tr*tan(Lr) - 240*b_2*cr**2*ts*Xte*tan(Lr) + 80*b_2*cr**2*ts*tr*tan(Lr) + 80*b_2*cr**2*ts*tan(Lr) - 160*b_2*cr**2*Xle**2*twb*tan(Lr) + 160*b_2*cr**2*Xle*twb*tr*tan(Lr) + 80*b_2*cr**2*Xle*twb*tan(Lr) + 160*b_2*cr**2*Xte**2*twb*tan(Lr) - 160*b_2*cr**2*Xte*twb*tr*tan(Lr) - 240*b_2*cr**2*Xte*twb*tan(Lr) + 80*b_2*cr**2*twb*tr*tan(Lr) + 80*b_2*cr**2*twb*tan(Lr) + 48*b_2*cr*ct**2*Xle**2*tr*tan(Lr) + 96*b_2*cr*ct**2*Xle**2*tt*tan(Lr) - 24*b_2*cr*ct**2*Xle*tr*tan(Lr) - 48*b_2*cr*ct**2*Xle*tt*tan(Lr) - 48*b_2*cr*ct**2*Xte**2*tr*tan(Lr) - 96*b_2*cr*ct**2*Xte**2*tt*tan(Lr) + 72*b_2*cr*ct**2*Xte*tr*tan(Lr) + 144*b_2*cr*ct**2*Xte*tt*tan(Lr) - 24*b_2*cr*ct**2*tr*tan(Lr) - 48*b_2*cr*ct**2*tt*tan(Lr) - 320*b_2*cr*ct*ts*Xle**2*tan(Lr) + 160*b_2*cr*ct*ts*Xle*tr*tan(Lr) + 160*b_2*cr*ct*ts*Xle*tt*tan(Lr) + 160*b_2*cr*ct*ts*Xle*tan(Lr) + 320*b_2*cr*ct*ts*Xte**2*tan(Lr) - 160*b_2*cr*ct*ts*Xte*tr*tan(Lr) - 160*b_2*cr*ct*ts*Xte*tt*tan(Lr) - 480*b_2*cr*ct*ts*Xte*tan(Lr) + 80*b_2*cr*ct*ts*tr*tan(Lr) + 80*b_2*cr*ct*ts*tt*tan(Lr) + 160*b_2*cr*ct*ts*tan(Lr) - 320*b_2*cr*ct*Xle**2*twb*tan(Lr) + 160*b_2*cr*ct*Xle*twb*tr*tan(Lr) + 160*b_2*cr*ct*Xle*twb*tt*tan(Lr) + 160*b_2*cr*ct*Xle*twb*tan(Lr) + 320*b_2*cr*ct*Xte**2*twb*tan(Lr) - 160*b_2*cr*ct*Xte*twb*tr*tan(Lr) - 160*b_2*cr*ct*Xte*twb*tt*tan(Lr) - 480*b_2*cr*ct*Xte*twb*tan(Lr) + 80*b_2*cr*ct*twb*tr*tan(Lr) + 80*b_2*cr*ct*twb*tt*tan(Lr) + 160*b_2*cr*ct*twb*tan(Lr) - 640*b_2*cr*ts**2*Xle*tan(Lr) + 640*b_2*cr*ts**2*Xte*tan(Lr) - 320*b_2*cr*ts**2*tan(Lr) - 1280*b_2*cr*ts*Xle*twb*tan(Lr) + 1280*b_2*cr*ts*Xte*twb*tan(Lr) - 640*b_2*cr*ts*twb*tan(Lr) - 640*b_2*cr*Xle*twb**2*tan(Lr) + 640*b_2*cr*Xte*twb**2*tan(Lr) - 320*b_2*cr*twb**2*tan(Lr) + 192*b_2*ct**3*Xle**2*tt*tan(Lr) - 96*b_2*ct**3*Xle*tt*tan(Lr) - 192*b_2*ct**3*Xte**2*tt*tan(Lr) + 288*b_2*ct**3*Xte*tt*tan(Lr) - 96*b_2*ct**3*tt*tan(Lr) - 480*b_2*ct**2*ts*Xle**2*tan(Lr) + 480*b_2*ct**2*ts*Xle*tt*tan(Lr) + 240*b_2*ct**2*ts*Xle*tan(Lr) + 480*b_2*ct**2*ts*Xte**2*tan(Lr) - 480*b_2*ct**2*ts*Xte*tt*tan(Lr) - 720*b_2*ct**2*ts*Xte*tan(Lr) + 240*b_2*ct**2*ts*tt*tan(Lr) + 240*b_2*ct**2*ts*tan(Lr) - 480*b_2*ct**2*Xle**2*twb*tan(Lr) + 480*b_2*ct**2*Xle*twb*tt*tan(Lr) + 240*b_2*ct**2*Xle*twb*tan(Lr) + 480*b_2*ct**2*Xte**2*twb*tan(Lr) - 480*b_2*ct**2*Xte*twb*tt*tan(Lr) - 720*b_2*ct**2*Xte*twb*tan(Lr) + 240*b_2*ct**2*twb*tt*tan(Lr) + 240*b_2*ct**2*twb*tan(Lr) - 1280*b_2*ct*ts**2*Xle*tan(Lr) + 1280*b_2*ct*ts**2*Xte*tan(Lr) - 640*b_2*ct*ts**2*tan(Lr) - 2560*b_2*ct*ts*Xle*twb*tan(Lr) + 2560*b_2*ct*ts*Xte*twb*tan(Lr) - 1280*b_2*ct*ts*twb*tan(Lr) - 1280*b_2*ct*Xle*twb**2*tan(Lr) + 1280*b_2*ct*Xte*twb**2*tan(Lr) - 640*b_2*ct*twb**2*tan(Lr) + 64*cr**4*Xle**3*tr - 48*cr**4*Xle**2*tr + 16*cr**4*Xle*tr**3 + 12*cr**4*Xle*tr + 64*cr**4*Xte**3*tr - 144*cr**4*Xte**2*tr + 16*cr**4*Xte*tr**3 + 108*cr**4*Xte*tr - 16*cr**4*tr**3 - 28*cr**4*tr + 48*cr**3*ct*Xle**3*tr + 16*cr**3*ct*Xle**3*tt - 36*cr**3*ct*Xle**2*tr - 12*cr**3*ct*Xle**2*tt + 4*cr**3*ct*Xle*tr**3 + 12*cr**3*ct*Xle*tr**2*tt + 9*cr**3*ct*Xle*tr + 3*cr**3*ct*Xle*tt + 48*cr**3*ct*Xte**3*tr + 16*cr**3*ct*Xte**3*tt - 108*cr**3*ct*Xte**2*tr - 36*cr**3*ct*Xte**2*tt + 4*cr**3*ct*Xte*tr**3 + 12*cr**3*ct*Xte*tr**2*tt + 81*cr**3*ct*Xte*tr + 27*cr**3*ct*Xte*tt - 4*cr**3*ct*tr**3 - 12*cr**3*ct*tr**2*tt - 21*cr**3*ct*tr - 7*cr**3*ct*tt - 160*cr**3*ts*Xle**3 + 240*cr**3*ts*Xle**2*tr + 120*cr**3*ts*Xle**2 - 120*cr**3*ts*Xle*tr**2 - 120*cr**3*ts*Xle*tr - 30*cr**3*ts*Xle - 160*cr**3*ts*Xte**3 + 240*cr**3*ts*Xte**2*tr + 360*cr**3*ts*Xte**2 - 120*cr**3*ts*Xte*tr**2 - 360*cr**3*ts*Xte*tr - 270*cr**3*ts*Xte + 40*cr**3*ts*tr**3 + 120*cr**3*ts*tr**2 + 150*cr**3*ts*tr + 70*cr**3*ts - 160*cr**3*Xle**3*twb + 240*cr**3*Xle**2*twb*tr + 120*cr**3*Xle**2*twb - 120*cr**3*Xle*twb*tr**2 - 120*cr**3*Xle*twb*tr - 30*cr**3*Xle*twb - 160*cr**3*Xte**3*twb + 240*cr**3*Xte**2*twb*tr + 360*cr**3*Xte**2*twb - 120*cr**3*Xte*twb*tr**2 - 360*cr**3*Xte*twb*tr - 270*cr**3*Xte*twb + 40*cr**3*twb*tr**3 + 120*cr**3*twb*tr**2 + 150*cr**3*twb*tr + 70*cr**3*twb + 32*cr**2*ct**2*Xle**3*tr + 32*cr**2*ct**2*Xle**3*tt - 24*cr**2*ct**2*Xle**2*tr - 24*cr**2*ct**2*Xle**2*tt + 8*cr**2*ct**2*Xle*tr**2*tt + 8*cr**2*ct**2*Xle*tr*tt**2 + 6*cr**2*ct**2*Xle*tr + 6*cr**2*ct**2*Xle*tt + 32*cr**2*ct**2*Xte**3*tr + 32*cr**2*ct**2*Xte**3*tt - 72*cr**2*ct**2*Xte**2*tr - 72*cr**2*ct**2*Xte**2*tt + 8*cr**2*ct**2*Xte*tr**2*tt + 8*cr**2*ct**2*Xte*tr*tt**2 + 54*cr**2*ct**2*Xte*tr + 54*cr**2*ct**2*Xte*tt - 8*cr**2*ct**2*tr**2*tt - 8*cr**2*ct**2*tr*tt**2 - 14*cr**2*ct**2*tr - 14*cr**2*ct**2*tt - 160*cr**2*ct*ts*Xle**3 + 160*cr**2*ct*ts*Xle**2*tr + 80*cr**2*ct*ts*Xle**2*tt + 120*cr**2*ct*ts*Xle**2 - 40*cr**2*ct*ts*Xle*tr**2 - 80*cr**2*ct*ts*Xle*tr*tt - 80*cr**2*ct*ts*Xle*tr - 40*cr**2*ct*ts*Xle*tt - 30*cr**2*ct*ts*Xle - 160*cr**2*ct*ts*Xte**3 + 160*cr**2*ct*ts*Xte**2*tr + 80*cr**2*ct*ts*Xte**2*tt + 360*cr**2*ct*ts*Xte**2 - 40*cr**2*ct*ts*Xte*tr**2 - 80*cr**2*ct*ts*Xte*tr*tt - 240*cr**2*ct*ts*Xte*tr - 120*cr**2*ct*ts*Xte*tt - 270*cr**2*ct*ts*Xte + 40*cr**2*ct*ts*tr**2*tt + 40*cr**2*ct*ts*tr**2 + 80*cr**2*ct*ts*tr*tt + 100*cr**2*ct*ts*tr + 50*cr**2*ct*ts*tt + 70*cr**2*ct*ts - 160*cr**2*ct*Xle**3*twb + 160*cr**2*ct*Xle**2*twb*tr + 80*cr**2*ct*Xle**2*twb*tt + 120*cr**2*ct*Xle**2*twb - 40*cr**2*ct*Xle*twb*tr**2 - 80*cr**2*ct*Xle*twb*tr*tt - 80*cr**2*ct*Xle*twb*tr - 40*cr**2*ct*Xle*twb*tt - 30*cr**2*ct*Xle*twb - 160*cr**2*ct*Xte**3*twb + 160*cr**2*ct*Xte**2*twb*tr + 80*cr**2*ct*Xte**2*twb*tt + 360*cr**2*ct*Xte**2*twb - 40*cr**2*ct*Xte*twb*tr**2 - 80*cr**2*ct*Xte*twb*tr*tt - 240*cr**2*ct*Xte*twb*tr - 120*cr**2*ct*Xte*twb*tt - 270*cr**2*ct*Xte*twb + 40*cr**2*ct*twb*tr**2*tt + 40*cr**2*ct*twb*tr**2 + 80*cr**2*ct*twb*tr*tt + 100*cr**2*ct*twb*tr + 50*cr**2*ct*twb*tt + 70*cr**2*ct*twb - 640*cr**2*ts**2*Xle**2 + 640*cr**2*ts**2*Xle*tr + 320*cr**2*ts**2*Xle - 640*cr**2*ts**2*Xte**2 + 640*cr**2*ts**2*Xte*tr + 960*cr**2*ts**2*Xte - 320*cr**2*ts**2*tr**2 - 640*cr**2*ts**2*tr - 400*cr**2*ts**2 - 1280*cr**2*ts*Xle**2*twb + 1280*cr**2*ts*Xle*twb*tr + 640*cr**2*ts*Xle*twb - 1280*cr**2*ts*Xte**2*twb + 1280*cr**2*ts*Xte*twb*tr + 1920*cr**2*ts*Xte*twb - 640*cr**2*ts*twb*tr**2 - 1280*cr**2*ts*twb*tr - 800*cr**2*ts*twb - 640*cr**2*Xle**2*twb**2 + 640*cr**2*Xle*twb**2*tr + 320*cr**2*Xle*twb**2 - 320*cr**2*Xle*tr*xcg**2 - 320*cr**2*Xle*tr*zcg**2 - 640*cr**2*Xte**2*twb**2 + 640*cr**2*Xte*twb**2*tr + 960*cr**2*Xte*twb**2 - 320*cr**2*Xte*tr*xcg**2 - 320*cr**2*Xte*tr*zcg**2 - 320*cr**2*twb**2*tr**2 - 640*cr**2*twb**2*tr - 400*cr**2*twb**2 + 320*cr**2*tr*xcg**2 + 320*cr**2*tr*zcg**2 + 16*cr*ct**3*Xle**3*tr + 48*cr*ct**3*Xle**3*tt - 12*cr*ct**3*Xle**2*tr - 36*cr*ct**3*Xle**2*tt + 12*cr*ct**3*Xle*tr*tt**2 + 3*cr*ct**3*Xle*tr + 4*cr*ct**3*Xle*tt**3 + 9*cr*ct**3*Xle*tt + 16*cr*ct**3*Xte**3*tr + 48*cr*ct**3*Xte**3*tt - 36*cr*ct**3*Xte**2*tr - 108*cr*ct**3*Xte**2*tt + 12*cr*ct**3*Xte*tr*tt**2 + 27*cr*ct**3*Xte*tr + 4*cr*ct**3*Xte*tt**3 + 81*cr*ct**3*Xte*tt - 12*cr*ct**3*tr*tt**2 - 7*cr*ct**3*tr - 4*cr*ct**3*tt**3 - 21*cr*ct**3*tt - 160*cr*ct**2*ts*Xle**3 + 80*cr*ct**2*ts*Xle**2*tr + 160*cr*ct**2*ts*Xle**2*tt + 120*cr*ct**2*ts*Xle**2 - 80*cr*ct**2*ts*Xle*tr*tt - 40*cr*ct**2*ts*Xle*tr - 40*cr*ct**2*ts*Xle*tt**2 - 80*cr*ct**2*ts*Xle*tt - 30*cr*ct**2*ts*Xle \
#   - 160*cr*ct**2*ts*Xte**3 + 80*cr*ct**2*ts*Xte**2*tr + 160*cr*ct**2*ts*Xte**2*tt + 360*cr*ct**2*ts*Xte**2 - 80*cr*ct**2*ts*Xte*tr*tt - 120*cr*ct**2*ts*Xte*tr - 40*cr*ct**2*ts*Xte*tt**2 - 240*cr*ct**2*ts*Xte*tt - 270*cr*ct**2*ts*Xte + 40*cr*ct**2*ts*tr*tt**2 + 80*cr*ct**2*ts*tr*tt + 50*cr*ct**2*ts*tr + 40*cr*ct**2*ts*tt**2 + 100*cr*ct**2*ts*tt + 70*cr*ct**2*ts - 160*cr*ct**2*Xle**3*twb + 80*cr*ct**2*Xle**2*twb*tr + 160*cr*ct**2*Xle**2*twb*tt + 120*cr*ct**2*Xle**2*twb - 80*cr*ct**2*Xle*twb*tr*tt - 40*cr*ct**2*Xle*twb*tr - 40*cr*ct**2*Xle*twb*tt**2 - 80*cr*ct**2*Xle*twb*tt - 30*cr*ct**2*Xle*twb - 160*cr*ct**2*Xte**3*twb + 80*cr*ct**2*Xte**2*twb*tr + 160*cr*ct**2*Xte**2*twb*tt + 360*cr*ct**2*Xte**2*twb - 80*cr*ct**2*Xte*twb*tr*tt - 120*cr*ct**2*Xte*twb*tr - 40*cr*ct**2*Xte*twb*tt**2 - 240*cr*ct**2*Xte*twb*tt - 270*cr*ct**2*Xte*twb + 40*cr*ct**2*twb*tr*tt**2 + 80*cr*ct**2*twb*tr*tt + 50*cr*ct**2*twb*tr + 40*cr*ct**2*twb*tt**2 + 100*cr*ct**2*twb*tt + 70*cr*ct**2*twb - 640*cr*ct*ts**2*Xle**2 + 320*cr*ct*ts**2*Xle*tr + 320*cr*ct*ts**2*Xle*tt + 320*cr*ct*ts**2*Xle - 640*cr*ct*ts**2*Xte**2 + 320*cr*ct*ts**2*Xte*tr + 320*cr*ct*ts**2*Xte*tt + 960*cr*ct*ts**2*Xte - 320*cr*ct*ts**2*tr*tt - 320*cr*ct*ts**2*tr - 320*cr*ct*ts**2*tt - 400*cr*ct*ts**2 - 1280*cr*ct*ts*Xle**2*twb + 640*cr*ct*ts*Xle*twb*tr + 640*cr*ct*ts*Xle*twb*tt + 640*cr*ct*ts*Xle*twb - 1280*cr*ct*ts*Xte**2*twb + 640*cr*ct*ts*Xte*twb*tr + 640*cr*ct*ts*Xte*twb*tt + 1920*cr*ct*ts*Xte*twb - 640*cr*ct*ts*twb*tr*tt - 640*cr*ct*ts*twb*tr - 640*cr*ct*ts*twb*tt - 800*cr*ct*ts*twb - 640*cr*ct*Xle**2*twb**2 + 320*cr*ct*Xle*twb**2*tr + 320*cr*ct*Xle*twb**2*tt + 320*cr*ct*Xle*twb**2 - 160*cr*ct*Xle*tr*xcg**2 - 160*cr*ct*Xle*tr*zcg**2 - 160*cr*ct*Xle*tt*xcg**2 - 160*cr*ct*Xle*tt*zcg**2 - 640*cr*ct*Xte**2*twb**2 + 320*cr*ct*Xte*twb**2*tr + 320*cr*ct*Xte*twb**2*tt + 960*cr*ct*Xte*twb**2 - 160*cr*ct*Xte*tr*xcg**2 - 160*cr*ct*Xte*tr*zcg**2 - 160*cr*ct*Xte*tt*xcg**2 - 160*cr*ct*Xte*tt*zcg**2 - 320*cr*ct*twb**2*tr*tt - 320*cr*ct*twb**2*tr - 320*cr*ct*twb**2*tt - 400*cr*ct*twb**2 + 160*cr*ct*tr*xcg**2 + 160*cr*ct*tr*zcg**2 + 160*cr*ct*tt*xcg**2 + 160*cr*ct*tt*zcg**2 - 1280*cr*ts**3*Xle - 1280*cr*ts**3*Xte + 1280*cr*ts**3*tr + 1280*cr*ts**3 - 3840*cr*ts**2*Xle*twb - 3840*cr*ts**2*Xte*twb + 3840*cr*ts**2*twb*tr + 3840*cr*ts**2*twb - 3840*cr*ts*Xle*twb**2 + 960*cr*ts*Xle*xcg**2 + 960*cr*ts*Xle*zcg**2 - 3840*cr*ts*Xte*twb**2 + 960*cr*ts*Xte*xcg**2 + 960*cr*ts*Xte*zcg**2 + 3840*cr*ts*twb**2*tr + 3840*cr*ts*twb**2 - 960*cr*ts*tr*xcg**2 - 960*cr*ts*tr*zcg**2 - 960*cr*ts*xcg**2 - 960*cr*ts*zcg**2 - 1280*cr*Xle*twb**3 + 960*cr*Xle*twb*xcg**2 + 960*cr*Xle*twb*zcg**2 - 1280*cr*Xte*twb**3 + 960*cr*Xte*twb*xcg**2 + 960*cr*Xte*twb*zcg**2 + 1280*cr*twb**3*tr + 1280*cr*twb**3 - 960*cr*twb*tr*xcg**2 - 960*cr*twb*tr*zcg**2 - 960*cr*twb*xcg**2 - 960*cr*twb*zcg**2 + 64*ct**4*Xle**3*tt - 48*ct**4*Xle**2*tt + 16*ct**4*Xle*tt**3 + 12*ct**4*Xle*tt + 64*ct**4*Xte**3*tt - 144*ct**4*Xte**2*tt + 16*ct**4*Xte*tt**3 + 108*ct**4*Xte*tt - 16*ct**4*tt**3 - 28*ct**4*tt - 160*ct**3*ts*Xle**3 + 240*ct**3*ts*Xle**2*tt + 120*ct**3*ts*Xle**2 - 120*ct**3*ts*Xle*tt**2 - 120*ct**3*ts*Xle*tt - 30*ct**3*ts*Xle - 160*ct**3*ts*Xte**3 + 240*ct**3*ts*Xte**2*tt + 360*ct**3*ts*Xte**2 - 120*ct**3*ts*Xte*tt**2 - 360*ct**3*ts*Xte*tt - 270*ct**3*ts*Xte + 40*ct**3*ts*tt**3 + 120*ct**3*ts*tt**2 + 150*ct**3*ts*tt + 70*ct**3*ts - 160*ct**3*Xle**3*twb + 240*ct**3*Xle**2*twb*tt + 120*ct**3*Xle**2*twb - 120*ct**3*Xle*twb*tt**2 - 120*ct**3*Xle*twb*tt - 30*ct**3*Xle*twb - 160*ct**3*Xte**3*twb + 240*ct**3*Xte**2*twb*tt + 360*ct**3*Xte**2*twb - 120*ct**3*Xte*twb*tt**2 - 360*ct**3*Xte*twb*tt - 270*ct**3*Xte*twb + 40*ct**3*twb*tt**3 + 120*ct**3*twb*tt**2 + 150*ct**3*twb*tt + 70*ct**3*twb - 640*ct**2*ts**2*Xle**2 + 640*ct**2*ts**2*Xle*tt + 320*ct**2*ts**2*Xle - 640*ct**2*ts**2*Xte**2 + 640*ct**2*ts**2*Xte*tt + 960*ct**2*ts**2*Xte - 320*ct**2*ts**2*tt**2 - 640*ct**2*ts**2*tt - 400*ct**2*ts**2 - 1280*ct**2*ts*Xle**2*twb + 1280*ct**2*ts*Xle*twb*tt + 640*ct**2*ts*Xle*twb - 1280*ct**2*ts*Xte**2*twb + 1280*ct**2*ts*Xte*twb*tt + 1920*ct**2*ts*Xte*twb - 640*ct**2*ts*twb*tt**2 - 1280*ct**2*ts*twb*tt - 800*ct**2*ts*twb - 640*ct**2*Xle**2*twb**2 + 640*ct**2*Xle*twb**2*tt + 320*ct**2*Xle*twb**2 - 320*ct**2*Xle*tt*xcg**2 - 320*ct**2*Xle*tt*zcg**2 - 640*ct**2*Xte**2*twb**2 + 640*ct**2*Xte*twb**2*tt + 960*ct**2*Xte*twb**2 - 320*ct**2*Xte*tt*xcg**2 - 320*ct**2*Xte*tt*zcg**2 - 320*ct**2*twb**2*tt**2 - 640*ct**2*twb**2*tt - 400*ct**2*twb**2 + 320*ct**2*tt*xcg**2 + 320*ct**2*tt*zcg**2 - 1280*ct*ts**3*Xle - 1280*ct*ts**3*Xte + 1280*ct*ts**3*tt + 1280*ct*ts**3 - 3840*ct*ts**2*Xle*twb - 3840*ct*ts**2*Xte*twb + 3840*ct*ts**2*twb*tt + 3840*ct*ts**2*twb - 3840*ct*ts*Xle*twb**2 + 960*ct*ts*Xle*xcg**2 + 960*ct*ts*Xle*zcg**2 - 3840*ct*ts*Xte*twb**2 + 960*ct*ts*Xte*xcg**2 + 960*ct*ts*Xte*zcg**2 + 3840*ct*ts*twb**2*tt + 3840*ct*ts*twb**2 - 960*ct*ts*tt*xcg**2 - 960*ct*ts*tt*zcg**2 - 960*ct*ts*xcg**2 - 960*ct*ts*zcg**2 - 1280*ct*Xle*twb**3 + 960*ct*Xle*twb*xcg**2 + 960*ct*Xle*twb*zcg**2 - 1280*ct*Xte*twb**3 + 960*ct*Xte*twb*xcg**2 + 960*ct*Xte*twb*zcg**2 + 1280*ct*twb**3*tt + 1280*ct*twb**3 - 960*ct*twb*tt*xcg**2 - 960*ct*twb*tt*zcg**2 - 960*ct*twb*xcg**2 - 960*ct*twb*zcg**2 - 2560*ts**4 - 10240*ts**3*twb - 15360*ts**2*twb**2 + 3840*ts**2*xcg**2 + 3840*ts**2*zcg**2 - 10240*ts*twb**3 + 7680*ts*twb*xcg**2 + 7680*ts*twb*zcg**2 - 2560*twb**4 + 3840*twb**2*xcg**2 + 3840*twb**2*zcg**2)/960 

# Izz_new = -b_2*p*(32*b_2**2*cr**2*Xle*tr*tan(Lr)**2 + 32*b_2**2*cr**2*Xle*tr + 32*b_2**2*cr**2*Xte*tr*tan(Lr)**2 + 32*b_2**2*cr**2*Xte*tr - 32*b_2**2*cr**2*tr*tan(Lr)**2 - 32*b_2**2*cr**2*tr + 48*b_2**2*cr*ct*Xle*tr*tan(Lr)**2 + 48*b_2**2*cr*ct*Xle*tr + 48*b_2**2*cr*ct*Xle*tt*tan(Lr)**2 + 48*b_2**2*cr*ct*Xle*tt + 48*b_2**2*cr*ct*Xte*tr*tan(Lr)**2 + 48*b_2**2*cr*ct*Xte*tr + 48*b_2**2*cr*ct*Xte*tt*tan(Lr)**2 + 48*b_2**2*cr*ct*Xte*tt - 48*b_2**2*cr*ct*tr*tan(Lr)**2 - 48*b_2**2*cr*ct*tr - 48*b_2**2*cr*ct*tt*tan(Lr)**2 - 48*b_2**2*cr*ct*tt - 160*b_2**2*cr*ts*Xle*tan(Lr)**2 - 160*b_2**2*cr*ts*Xle - 160*b_2**2*cr*ts*Xte*tan(Lr)**2 - 160*b_2**2*cr*ts*Xte + 160*b_2**2*cr*ts*tr*tan(Lr)**2 + 160*b_2**2*cr*ts*tr + 160*b_2**2*cr*ts*tan(Lr)**2 + 160*b_2**2*cr*ts - 160*b_2**2*cr*Xle*twb*tan(Lr)**2 - 160*b_2**2*cr*Xle*twb - 160*b_2**2*cr*Xte*twb*tan(Lr)**2 - 160*b_2**2*cr*Xte*twb + 160*b_2**2*cr*twb*tr*tan(Lr)**2 + 160*b_2**2*cr*twb*tr + 160*b_2**2*cr*twb*tan(Lr)**2 + 160*b_2**2*cr*twb + 192*b_2**2*ct**2*Xle*tt*tan(Lr)**2 + 192*b_2**2*ct**2*Xle*tt + 192*b_2**2*ct**2*Xte*tt*tan(Lr)**2 + 192*b_2**2*ct**2*Xte*tt - 192*b_2**2*ct**2*tt*tan(Lr)**2 - 192*b_2**2*ct**2*tt - 480*b_2**2*ct*ts*Xle*tan(Lr)**2 - 480*b_2**2*ct*ts*Xle - 480*b_2**2*ct*ts*Xte*tan(Lr)**2 - 480*b_2**2*ct*ts*Xte + 480*b_2**2*ct*ts*tt*tan(Lr)**2 + 480*b_2**2*ct*ts*tt + 480*b_2**2*ct*ts*tan(Lr)**2 + 480*b_2**2*ct*ts - 480*b_2**2*ct*Xle*twb*tan(Lr)**2 - 480*b_2**2*ct*Xle*twb - 480*b_2**2*ct*Xte*twb*tan(Lr)**2 - 480*b_2**2*ct*Xte*twb + 480*b_2**2*ct*twb*tt*tan(Lr)**2 + 480*b_2**2*ct*twb*tt + 480*b_2**2*ct*twb*tan(Lr)**2 + 480*b_2**2*ct*twb - 1280*b_2**2*ts**2*tan(Lr)**2 - 1280*b_2**2*ts**2 - 2560*b_2**2*ts*twb*tan(Lr)**2 - 2560*b_2**2*ts*twb - 1280*b_2**2*twb**2*tan(Lr)**2 - 1280*b_2**2*twb**2 + 48*b_2*cr**3*Xle**2*tr*tan(Lr) - 24*b_2*cr**3*Xle*tr*tan(Lr) - 48*b_2*cr**3*Xte**2*tr*tan(Lr) + 72*b_2*cr**3*Xte*tr*tan(Lr) - 24*b_2*cr**3*tr*tan(Lr) + 64*b_2*cr**2*ct*Xle**2*tr*tan(Lr) + 32*b_2*cr**2*ct*Xle**2*tt*tan(Lr) - 32*b_2*cr**2*ct*Xle*tr*tan(Lr) - 16*b_2*cr**2*ct*Xle*tt*tan(Lr) - 64*b_2*cr**2*ct*Xte**2*tr*tan(Lr) - 32*b_2*cr**2*ct*Xte**2*tt*tan(Lr) + 96*b_2*cr**2*ct*Xte*tr*tan(Lr) + 48*b_2*cr**2*ct*Xte*tt*tan(Lr) - 32*b_2*cr**2*ct*tr*tan(Lr) - 16*b_2*cr**2*ct*tt*tan(Lr) - 160*b_2*cr**2*ts*Xle**2*tan(Lr) + 160*b_2*cr**2*ts*Xle*tr*tan(Lr) + 80*b_2*cr**2*ts*Xle*tan(Lr) + 160*b_2*cr**2*ts*Xte**2*tan(Lr) - 160*b_2*cr**2*ts*Xte*tr*tan(Lr) - 240*b_2*cr**2*ts*Xte*tan(Lr) + 80*b_2*cr**2*ts*tr*tan(Lr) + 80*b_2*cr**2*ts*tan(Lr) - 160*b_2*cr**2*Xle**2*twb*tan(Lr) + 160*b_2*cr**2*Xle*twb*tr*tan(Lr) + 80*b_2*cr**2*Xle*twb*tan(Lr) + 160*b_2*cr**2*Xte**2*twb*tan(Lr) - 160*b_2*cr**2*Xte*twb*tr*tan(Lr) - 240*b_2*cr**2*Xte*twb*tan(Lr) + 80*b_2*cr**2*twb*tr*tan(Lr) + 80*b_2*cr**2*twb*tan(Lr) + 48*b_2*cr*ct**2*Xle**2*tr*tan(Lr) + 96*b_2*cr*ct**2*Xle**2*tt*tan(Lr) - 24*b_2*cr*ct**2*Xle*tr*tan(Lr) - 48*b_2*cr*ct**2*Xle*tt*tan(Lr) - 48*b_2*cr*ct**2*Xte**2*tr*tan(Lr) - 96*b_2*cr*ct**2*Xte**2*tt*tan(Lr) + 72*b_2*cr*ct**2*Xte*tr*tan(Lr) + 144*b_2*cr*ct**2*Xte*tt*tan(Lr) - 24*b_2*cr*ct**2*tr*tan(Lr) - 48*b_2*cr*ct**2*tt*tan(Lr) - 320*b_2*cr*ct*ts*Xle**2*tan(Lr) + 160*b_2*cr*ct*ts*Xle*tr*tan(Lr) + 160*b_2*cr*ct*ts*Xle*tt*tan(Lr) + 160*b_2*cr*ct*ts*Xle*tan(Lr) + 320*b_2*cr*ct*ts*Xte**2*tan(Lr) - 160*b_2*cr*ct*ts*Xte*tr*tan(Lr) - 160*b_2*cr*ct*ts*Xte*tt*tan(Lr) - 480*b_2*cr*ct*ts*Xte*tan(Lr) + 80*b_2*cr*ct*ts*tr*tan(Lr) + 80*b_2*cr*ct*ts*tt*tan(Lr) + 160*b_2*cr*ct*ts*tan(Lr) - 320*b_2*cr*ct*Xle**2*twb*tan(Lr) + 160*b_2*cr*ct*Xle*twb*tr*tan(Lr) + 160*b_2*cr*ct*Xle*twb*tt*tan(Lr) + 160*b_2*cr*ct*Xle*twb*tan(Lr) + 320*b_2*cr*ct*Xte**2*twb*tan(Lr) - 160*b_2*cr*ct*Xte*twb*tr*tan(Lr) - 160*b_2*cr*ct*Xte*twb*tt*tan(Lr) - 480*b_2*cr*ct*Xte*twb*tan(Lr) + 80*b_2*cr*ct*twb*tr*tan(Lr) + 80*b_2*cr*ct*twb*tt*tan(Lr) + 160*b_2*cr*ct*twb*tan(Lr) - 640*b_2*cr*ts**2*Xle*tan(Lr) + 640*b_2*cr*ts**2*Xte*tan(Lr) - 320*b_2*cr*ts**2*tan(Lr) - 1280*b_2*cr*ts*Xle*twb*tan(Lr) + 1280*b_2*cr*ts*Xte*twb*tan(Lr) - 640*b_2*cr*ts*twb*tan(Lr) - 640*b_2*cr*Xle*twb**2*tan(Lr) + 640*b_2*cr*Xte*twb**2*tan(Lr) - 320*b_2*cr*twb**2*tan(Lr) + 192*b_2*ct**3*Xle**2*tt*tan(Lr) - 96*b_2*ct**3*Xle*tt*tan(Lr) - 192*b_2*ct**3*Xte**2*tt*tan(Lr) + 288*b_2*ct**3*Xte*tt*tan(Lr) - 96*b_2*ct**3*tt*tan(Lr) - 480*b_2*ct**2*ts*Xle**2*tan(Lr) + 480*b_2*ct**2*ts*Xle*tt*tan(Lr) + 240*b_2*ct**2*ts*Xle*tan(Lr) + 480*b_2*ct**2*ts*Xte**2*tan(Lr) - 480*b_2*ct**2*ts*Xte*tt*tan(Lr) - 720*b_2*ct**2*ts*Xte*tan(Lr) + 240*b_2*ct**2*ts*tt*tan(Lr) + 240*b_2*ct**2*ts*tan(Lr) - 480*b_2*ct**2*Xle**2*twb*tan(Lr) + 480*b_2*ct**2*Xle*twb*tt*tan(Lr) + 240*b_2*ct**2*Xle*twb*tan(Lr) + 480*b_2*ct**2*Xte**2*twb*tan(Lr) - 480*b_2*ct**2*Xte*twb*tt*tan(Lr) - 720*b_2*ct**2*Xte*twb*tan(Lr) + 240*b_2*ct**2*twb*tt*tan(Lr) + 240*b_2*ct**2*twb*tan(Lr) - 1280*b_2*ct*ts**2*Xle*tan(Lr) + 1280*b_2*ct*ts**2*Xte*tan(Lr) - 640*b_2*ct*ts**2*tan(Lr) - 2560*b_2*ct*ts*Xle*twb*tan(Lr) + 2560*b_2*ct*ts*Xte*twb*tan(Lr) - 1280*b_2*ct*ts*twb*tan(Lr) - 1280*b_2*ct*Xle*twb**2*tan(Lr) + 1280*b_2*ct*Xte*twb**2*tan(Lr) - 640*b_2*ct*twb**2*tan(Lr) + 64*cr**4*Xle**3*tr - 48*cr**4*Xle**2*tr + 12*cr**4*Xle*tr + 64*cr**4*Xte**3*tr - 144*cr**4*Xte**2*tr + 108*cr**4*Xte*tr - 28*cr**4*tr + 48*cr**3*ct*Xle**3*tr + 16*cr**3*ct*Xle**3*tt - 36*cr**3*ct*Xle**2*tr - 12*cr**3*ct*Xle**2*tt + 9*cr**3*ct*Xle*tr + 3*cr**3*ct*Xle*tt + 48*cr**3*ct*Xte**3*tr + 16*cr**3*ct*Xte**3*tt - 108*cr**3*ct*Xte**2*tr - 36*cr**3*ct*Xte**2*tt + 81*cr**3*ct*Xte*tr + 27*cr**3*ct*Xte*tt - 21*cr**3*ct*tr - 7*cr**3*ct*tt - 160*cr**3*ts*Xle**3 + 240*cr**3*ts*Xle**2*tr + 120*cr**3*ts*Xle**2 - 120*cr**3*ts*Xle*tr - 30*cr**3*ts*Xle - 160*cr**3*ts*Xte**3 + 240*cr**3*ts*Xte**2*tr + 360*cr**3*ts*Xte**2 - 360*cr**3*ts*Xte*tr - 270*cr**3*ts*Xte + 150*cr**3*ts*tr + 70*cr**3*ts - 160*cr**3*Xle**3*twb + 240*cr**3*Xle**2*twb*tr + 120*cr**3*Xle**2*twb - 120*cr**3*Xle*twb*tr - 30*cr**3*Xle*twb - 160*cr**3*Xte**3*twb + 240*cr**3*Xte**2*twb*tr + 360*cr**3*Xte**2*twb - 360*cr**3*Xte*twb*tr - 270*cr**3*Xte*twb + 150*cr**3*twb*tr + 70*cr**3*twb + 32*cr**2*ct**2*Xle**3*tr + 32*cr**2*ct**2*Xle**3*tt - 24*cr**2*ct**2*Xle**2*tr - 24*cr**2*ct**2*Xle**2*tt + 6*cr**2*ct**2*Xle*tr + 6*cr**2*ct**2*Xle*tt + 32*cr**2*ct**2*Xte**3*tr + 32*cr**2*ct**2*Xte**3*tt - 72*cr**2*ct**2*Xte**2*tr - 72*cr**2*ct**2*Xte**2*tt + 54*cr**2*ct**2*Xte*tr + 54*cr**2*ct**2*Xte*tt - 14*cr**2*ct**2*tr - 14*cr**2*ct**2*tt - 160*cr**2*ct*ts*Xle**3 + 160*cr**2*ct*ts*Xle**2*tr + 80*cr**2*ct*ts*Xle**2*tt + 120*cr**2*ct*ts*Xle**2 - 80*cr**2*ct*ts*Xle*tr - 40*cr**2*ct*ts*Xle*tt - 30*cr**2*ct*ts*Xle - 160*cr**2*ct*ts*Xte**3 + 160*cr**2*ct*ts*Xte**2*tr + 80*cr**2*ct*ts*Xte**2*tt + 360*cr**2*ct*ts*Xte**2 - 240*cr**2*ct*ts*Xte*tr - 120*cr**2*ct*ts*Xte*tt - 270*cr**2*ct*ts*Xte + 100*cr**2*ct*ts*tr + 50*cr**2*ct*ts*tt + 70*cr**2*ct*ts - 160*cr**2*ct*Xle**3*twb + 160*cr**2*ct*Xle**2*twb*tr + 80*cr**2*ct*Xle**2*twb*tt + 120*cr**2*ct*Xle**2*twb - 80*cr**2*ct*Xle*twb*tr - 40*cr**2*ct*Xle*twb*tt - 30*cr**2*ct*Xle*twb - 160*cr**2*ct*Xte**3*twb + 160*cr**2*ct*Xte**2*twb*tr + 80*cr**2*ct*Xte**2*twb*tt + 360*cr**2*ct*Xte**2*twb - 240*cr**2*ct*Xte*twb*tr - 120*cr**2*ct*Xte*twb*tt - 270*cr**2*ct*Xte*twb + 100*cr**2*ct*twb*tr + 50*cr**2*ct*twb*tt + 70*cr**2*ct*twb - 640*cr**2*ts**2*Xle**2 + 320*cr**2*ts**2*Xle*tr + 320*cr**2*ts**2*Xle - 640*cr**2*ts**2*Xte**2 + 320*cr**2*ts**2*Xte*tr + 960*cr**2*ts**2*Xte - 320*cr**2*ts**2*tr - 400*cr**2*ts**2 - 1280*cr**2*ts*Xle**2*twb + 640*cr**2*ts*Xle*twb*tr + 640*cr**2*ts*Xle*twb - 1280*cr**2*ts*Xte**2*twb + 640*cr**2*ts*Xte*twb*tr + 1920*cr**2*ts*Xte*twb - 640*cr**2*ts*twb*tr - 800*cr**2*ts*twb - 640*cr**2*Xle**2*twb**2 + 320*cr**2*Xle*twb**2*tr + 320*cr**2*Xle*twb**2 - 320*cr**2*Xle*tr*xcg**2 - 320*cr**2*Xle*tr*ycg**2 - 640*cr**2*Xte**2*twb**2 + 320*cr**2*Xte*twb**2*tr + 960*cr**2*Xte*twb**2 - 320*cr**2*Xte*tr*xcg**2 - 320*cr**2*Xte*tr*ycg**2 - 320*cr**2*twb**2*tr - 400*cr**2*twb**2 + 320*cr**2*tr*xcg**2 + 320*cr**2*tr*ycg**2 + 16*cr*ct**3*Xle**3*tr + 48*cr*ct**3*Xle**3*tt - 12*cr*ct**3*Xle**2*tr - 36*cr*ct**3*Xle**2*tt + 3*cr*ct**3*Xle*tr + 9*cr*ct**3*Xle*tt + 16*cr*ct**3*Xte**3*tr + 48*cr*ct**3*Xte**3*tt - 36*cr*ct**3*Xte**2*tr - 108*cr*ct**3*Xte**2*tt + 27*cr*ct**3*Xte*tr + 81*cr*ct**3*Xte*tt - 7*cr*ct**3*tr - 21*cr*ct**3*tt - 160*cr*ct**2*ts*Xle**3 + 80*cr*ct**2*ts*Xle**2*tr + 160*cr*ct**2*ts*Xle**2*tt + 120*cr*ct**2*ts*Xle**2 - 40*cr*ct**2*ts*Xle*tr - 80*cr*ct**2*ts*Xle*tt - 30*cr*ct**2*ts*Xle - 160*cr*ct**2*ts*Xte**3 + 80*cr*ct**2*ts*Xte**2*tr + 160*cr*ct**2*ts*Xte**2*tt + 360*cr*ct**2*ts*Xte**2 - 120*cr*ct**2*ts*Xte*tr - 240*cr*ct**2*ts*Xte*tt - 270*cr*ct**2*ts*Xte + 50*cr*ct**2*ts*tr + 100*cr*ct**2*ts*tt + 70*cr*ct**2*ts - 160*cr*ct**2*Xle**3*twb + 80*cr*ct**2*Xle**2*twb*tr + 160*cr*ct**2*Xle**2*twb*tt + 120*cr*ct**2*Xle**2*twb - 40*cr*ct**2*Xle*twb*tr - 80*cr*ct**2*Xle*twb*tt - 30*cr*ct**2*Xle*twb - 160*cr*ct**2*Xte**3*twb + 80*cr*ct**2*Xte**2*twb*tr + 160*cr*ct**2*Xte**2*twb*tt + 360*cr*ct**2*Xte**2*twb \
#   - 120*cr*ct**2*Xte*twb*tr - 240*cr*ct**2*Xte*twb*tt - 270*cr*ct**2*Xte*twb + 50*cr*ct**2*twb*tr + 100*cr*ct**2*twb*tt + 70*cr*ct**2*twb - 640*cr*ct*ts**2*Xle**2 + 160*cr*ct*ts**2*Xle*tr + 160*cr*ct*ts**2*Xle*tt + 320*cr*ct*ts**2*Xle - 640*cr*ct*ts**2*Xte**2 + 160*cr*ct*ts**2*Xte*tr + 160*cr*ct*ts**2*Xte*tt + 960*cr*ct*ts**2*Xte - 160*cr*ct*ts**2*tr - 160*cr*ct*ts**2*tt - 400*cr*ct*ts**2 - 1280*cr*ct*ts*Xle**2*twb + 320*cr*ct*ts*Xle*twb*tr + 320*cr*ct*ts*Xle*twb*tt + 640*cr*ct*ts*Xle*twb - 1280*cr*ct*ts*Xte**2*twb + 320*cr*ct*ts*Xte*twb*tr + 320*cr*ct*ts*Xte*twb*tt + 1920*cr*ct*ts*Xte*twb - 320*cr*ct*ts*twb*tr - 320*cr*ct*ts*twb*tt - 800*cr*ct*ts*twb - 640*cr*ct*Xle**2*twb**2 + 160*cr*ct*Xle*twb**2*tr + 160*cr*ct*Xle*twb**2*tt + 320*cr*ct*Xle*twb**2 - 160*cr*ct*Xle*tr*xcg**2 - 160*cr*ct*Xle*tr*ycg**2 - 160*cr*ct*Xle*tt*xcg**2 - 160*cr*ct*Xle*tt*ycg**2 - 640*cr*ct*Xte**2*twb**2 + 160*cr*ct*Xte*twb**2*tr + 160*cr*ct*Xte*twb**2*tt + 960*cr*ct*Xte*twb**2 - 160*cr*ct*Xte*tr*xcg**2 - 160*cr*ct*Xte*tr*ycg**2 - 160*cr*ct*Xte*tt*xcg**2 - 160*cr*ct*Xte*tt*ycg**2 - 160*cr*ct*twb**2*tr - 160*cr*ct*twb**2*tt - 400*cr*ct*twb**2 + 160*cr*ct*tr*xcg**2 + 160*cr*ct*tr*ycg**2 + 160*cr*ct*tt*xcg**2 + 160*cr*ct*tt*ycg**2 - 960*cr*ts**3*Xle - 960*cr*ts**3*Xte + 320*cr*ts**3*tr + 960*cr*ts**3 - 2880*cr*ts**2*Xle*twb - 2880*cr*ts**2*Xte*twb + 960*cr*ts**2*twb*tr + 2880*cr*ts**2*twb - 2880*cr*ts*Xle*twb**2 + 960*cr*ts*Xle*xcg**2 + 960*cr*ts*Xle*ycg**2 - 2880*cr*ts*Xte*twb**2 + 960*cr*ts*Xte*xcg**2 + 960*cr*ts*Xte*ycg**2 + 960*cr*ts*twb**2*tr + 2880*cr*ts*twb**2 - 960*cr*ts*tr*xcg**2 - 960*cr*ts*tr*ycg**2 - 960*cr*ts*xcg**2 - 960*cr*ts*ycg**2 - 960*cr*Xle*twb**3 + 960*cr*Xle*twb*xcg**2 + 960*cr*Xle*twb*ycg**2 - 960*cr*Xte*twb**3 + 960*cr*Xte*twb*xcg**2 + 960*cr*Xte*twb*ycg**2 + 320*cr*twb**3*tr + 960*cr*twb**3 - 960*cr*twb*tr*xcg**2 - 960*cr*twb*tr*ycg**2 - 960*cr*twb*xcg**2 - 960*cr*twb*ycg**2 + 64*ct**4*Xle**3*tt - 48*ct**4*Xle**2*tt + 12*ct**4*Xle*tt + 64*ct**4*Xte**3*tt - 144*ct**4*Xte**2*tt + 108*ct**4*Xte*tt - 28*ct**4*tt - 160*ct**3*ts*Xle**3 + 240*ct**3*ts*Xle**2*tt + 120*ct**3*ts*Xle**2 - 120*ct**3*ts*Xle*tt - 30*ct**3*ts*Xle - 160*ct**3*ts*Xte**3 + 240*ct**3*ts*Xte**2*tt + 360*ct**3*ts*Xte**2 - 360*ct**3*ts*Xte*tt - 270*ct**3*ts*Xte + 150*ct**3*ts*tt + 70*ct**3*ts - 160*ct**3*Xle**3*twb + 240*ct**3*Xle**2*twb*tt + 120*ct**3*Xle**2*twb - 120*ct**3*Xle*twb*tt - 30*ct**3*Xle*twb - 160*ct**3*Xte**3*twb + 240*ct**3*Xte**2*twb*tt + 360*ct**3*Xte**2*twb - 360*ct**3*Xte*twb*tt - 270*ct**3*Xte*twb + 150*ct**3*twb*tt + 70*ct**3*twb - 640*ct**2*ts**2*Xle**2 + 320*ct**2*ts**2*Xle*tt + 320*ct**2*ts**2*Xle - 640*ct**2*ts**2*Xte**2 + 320*ct**2*ts**2*Xte*tt + 960*ct**2*ts**2*Xte - 320*ct**2*ts**2*tt - 400*ct**2*ts**2 - 1280*ct**2*ts*Xle**2*twb + 640*ct**2*ts*Xle*twb*tt + 640*ct**2*ts*Xle*twb - 1280*ct**2*ts*Xte**2*twb + 640*ct**2*ts*Xte*twb*tt + 1920*ct**2*ts*Xte*twb - 640*ct**2*ts*twb*tt - 800*ct**2*ts*twb - 640*ct**2*Xle**2*twb**2 + 320*ct**2*Xle*twb**2*tt + 320*ct**2*Xle*twb**2 - 320*ct**2*Xle*tt*xcg**2 - 320*ct**2*Xle*tt*ycg**2 - 640*ct**2*Xte**2*twb**2 + 320*ct**2*Xte*twb**2*tt + 960*ct**2*Xte*twb**2 - 320*ct**2*Xte*tt*xcg**2 - 320*ct**2*Xte*tt*ycg**2 - 320*ct**2*twb**2*tt - 400*ct**2*twb**2 + 320*ct**2*tt*xcg**2 + 320*ct**2*tt*ycg**2 - 960*ct*ts**3*Xle - 960*ct*ts**3*Xte + 320*ct*ts**3*tt + 960*ct*ts**3 - 2880*ct*ts**2*Xle*twb - 2880*ct*ts**2*Xte*twb + 960*ct*ts**2*twb*tt + 2880*ct*ts**2*twb - 2880*ct*ts*Xle*twb**2 + 960*ct*ts*Xle*xcg**2 + 960*ct*ts*Xle*ycg**2 - 2880*ct*ts*Xte*twb**2 + 960*ct*ts*Xte*xcg**2 + 960*ct*ts*Xte*ycg**2 + 960*ct*ts*twb**2*tt + 2880*ct*ts*twb**2 - 960*ct*ts*tt*xcg**2 - 960*ct*ts*tt*ycg**2 - 960*ct*ts*xcg**2 - 960*ct*ts*ycg**2 - 960*ct*Xle*twb**3 + 960*ct*Xle*twb*xcg**2 + 960*ct*Xle*twb*ycg**2 - 960*ct*Xte*twb**3 + 960*ct*Xte*twb*xcg**2 + 960*ct*Xte*twb*ycg**2 + 320*ct*twb**3*tt + 960*ct*twb**3 - 960*ct*twb*tt*xcg**2 - 960*ct*twb*tt*ycg**2 - 960*ct*twb*xcg**2 - 960*ct*twb*ycg**2 - 1280*ts**4 - 5120*ts**3*twb - 7680*ts**2*twb**2 + 3840*ts**2*xcg**2 + 3840*ts**2*ycg**2 - 5120*ts*twb**3 + 7680*ts*twb*xcg**2 + 7680*ts*twb*ycg**2 - 1280*twb**4 + 3840*twb**2*xcg**2 + 3840*twb**2*ycg**2)/960 

# Ixy_new = b_2*p*(8*b_2**2*cr**2*Xle*tr*tan(Lr) + 8*b_2**2*cr**2*Xte*tr*tan(Lr) - 8*b_2**2*cr**2*tr*tan(Lr) + 12*b_2**2*cr*ct*Xle*tr*tan(Lr) + 12*b_2**2*cr*ct*Xle*tt*tan(Lr) + 12*b_2**2*cr*ct*Xte*tr*tan(Lr) + 12*b_2**2*cr*ct*Xte*tt*tan(Lr) - 12*b_2**2*cr*ct*tr*tan(Lr) - 12*b_2**2*cr*ct*tt*tan(Lr) - 40*b_2**2*cr*ts*Xle*tan(Lr) - 40*b_2**2*cr*ts*Xte*tan(Lr) + 40*b_2**2*cr*ts*tr*tan(Lr) + 40*b_2**2*cr*ts*tan(Lr) - 40*b_2**2*cr*Xle*twb*tan(Lr) - 40*b_2**2*cr*Xte*twb*tan(Lr) + 40*b_2**2*cr*twb*tr*tan(Lr) + 40*b_2**2*cr*twb*tan(Lr) + 48*b_2**2*ct**2*Xle*tt*tan(Lr) + 48*b_2**2*ct**2*Xte*tt*tan(Lr) - 48*b_2**2*ct**2*tt*tan(Lr) - 120*b_2**2*ct*ts*Xle*tan(Lr) - 120*b_2**2*ct*ts*Xte*tan(Lr) + 120*b_2**2*ct*ts*tt*tan(Lr) + 120*b_2**2*ct*ts*tan(Lr) - 120*b_2**2*ct*Xle*twb*tan(Lr) - 120*b_2**2*ct*Xte*twb*tan(Lr) + 120*b_2**2*ct*twb*tt*tan(Lr) + 120*b_2**2*ct*twb*tan(Lr) - 320*b_2**2*ts**2*tan(Lr) - 640*b_2**2*ts*twb*tan(Lr) - 320*b_2**2*twb**2*tan(Lr) + 6*b_2*cr**3*Xle**2*tr - 3*b_2*cr**3*Xle*tr - 6*b_2*cr**3*Xte**2*tr + 9*b_2*cr**3*Xte*tr - 3*b_2*cr**3*tr + 8*b_2*cr**2*ct*Xle**2*tr + 4*b_2*cr**2*ct*Xle**2*tt - 4*b_2*cr**2*ct*Xle*tr - 2*b_2*cr**2*ct*Xle*tt - 8*b_2*cr**2*ct*Xte**2*tr - 4*b_2*cr**2*ct*Xte**2*tt + 12*b_2*cr**2*ct*Xte*tr + 6*b_2*cr**2*ct*Xte*tt - 4*b_2*cr**2*ct*tr - 2*b_2*cr**2*ct*tt - 20*b_2*cr**2*ts*Xle**2 + 20*b_2*cr**2*ts*Xle*tr + 10*b_2*cr**2*ts*Xle + 20*b_2*cr**2*ts*Xte**2 - 20*b_2*cr**2*ts*Xte*tr - 30*b_2*cr**2*ts*Xte + 10*b_2*cr**2*ts*tr + 10*b_2*cr**2*ts - 20*b_2*cr**2*Xle**2*twb + 20*b_2*cr**2*Xle*twb*tr + 10*b_2*cr**2*Xle*twb + 20*b_2*cr**2*Xte**2*twb - 20*b_2*cr**2*Xte*twb*tr - 30*b_2*cr**2*Xte*twb + 10*b_2*cr**2*twb*tr + 10*b_2*cr**2*twb + 6*b_2*cr*ct**2*Xle**2*tr + 12*b_2*cr*ct**2*Xle**2*tt - 3*b_2*cr*ct**2*Xle*tr - 6*b_2*cr*ct**2*Xle*tt - 6*b_2*cr*ct**2*Xte**2*tr - 12*b_2*cr*ct**2*Xte**2*tt + 9*b_2*cr*ct**2*Xte*tr + 18*b_2*cr*ct**2*Xte*tt - 3*b_2*cr*ct**2*tr - 6*b_2*cr*ct**2*tt - 40*b_2*cr*ct*ts*Xle**2 + 20*b_2*cr*ct*ts*Xle*tr + 20*b_2*cr*ct*ts*Xle*tt + 20*b_2*cr*ct*ts*Xle + 40*b_2*cr*ct*ts*Xte**2 - 20*b_2*cr*ct*ts*Xte*tr - 20*b_2*cr*ct*ts*Xte*tt - 60*b_2*cr*ct*ts*Xte + 10*b_2*cr*ct*ts*tr + 10*b_2*cr*ct*ts*tt + 20*b_2*cr*ct*ts - 40*b_2*cr*ct*Xle**2*twb + 20*b_2*cr*ct*Xle*twb*tr + 20*b_2*cr*ct*Xle*twb*tt + 20*b_2*cr*ct*Xle*twb + 40*b_2*cr*ct*Xte**2*twb - 20*b_2*cr*ct*Xte*twb*tr - 20*b_2*cr*ct*Xte*twb*tt - 60*b_2*cr*ct*Xte*twb + 10*b_2*cr*ct*twb*tr + 10*b_2*cr*ct*twb*tt + 20*b_2*cr*ct*twb - 80*b_2*cr*ts**2*Xle + 80*b_2*cr*ts**2*Xte - 40*b_2*cr*ts**2 - 160*b_2*cr*ts*Xle*twb + 160*b_2*cr*ts*Xte*twb - 80*b_2*cr*ts*twb - 80*b_2*cr*Xle*twb**2 + 80*b_2*cr*Xte*twb**2 - 40*b_2*cr*twb**2 + 24*b_2*ct**3*Xle**2*tt - 12*b_2*ct**3*Xle*tt - 24*b_2*ct**3*Xte**2*tt + 36*b_2*ct**3*Xte*tt - 12*b_2*ct**3*tt - 60*b_2*ct**2*ts*Xle**2 + 60*b_2*ct**2*ts*Xle*tt + 30*b_2*ct**2*ts*Xle + 60*b_2*ct**2*ts*Xte**2 - 60*b_2*ct**2*ts*Xte*tt - 90*b_2*ct**2*ts*Xte + 30*b_2*ct**2*ts*tt + 30*b_2*ct**2*ts - 60*b_2*ct**2*Xle**2*twb + 60*b_2*ct**2*Xle*twb*tt + 30*b_2*ct**2*Xle*twb + 60*b_2*ct**2*Xte**2*twb - 60*b_2*ct**2*Xte*twb*tt - 90*b_2*ct**2*Xte*twb + 30*b_2*ct**2*twb*tt + 30*b_2*ct**2*twb - 160*b_2*ct*ts**2*Xle + 160*b_2*ct*ts**2*Xte - 80*b_2*ct*ts**2 - 320*b_2*ct*ts*Xle*twb + 320*b_2*ct*ts*Xte*twb - 160*b_2*ct*ts*twb - 160*b_2*ct*Xle*twb**2 + 160*b_2*ct*Xte*twb**2 - 80*b_2*ct*twb**2 + 80*cr**2*Xle*tr*xcg*ycg + 80*cr**2*Xte*tr*xcg*ycg - 80*cr**2*tr*xcg*ycg + 40*cr*ct*Xle*tr*xcg*ycg + 40*cr*ct*Xle*tt*xcg*ycg + 40*cr*ct*Xte*tr*xcg*ycg + 40*cr*ct*Xte*tt*xcg*ycg - 40*cr*ct*tr*xcg*ycg - 40*cr*ct*tt*xcg*ycg - 240*cr*ts*Xle*xcg*ycg - 240*cr*ts*Xte*xcg*ycg + 240*cr*ts*tr*xcg*ycg + 240*cr*ts*xcg*ycg - 240*cr*Xle*twb*xcg*ycg - 240*cr*Xte*twb*xcg*ycg + 240*cr*twb*tr*xcg*ycg + 240*cr*twb*xcg*ycg + 80*ct**2*Xle*tt*xcg*ycg + 80*ct**2*Xte*tt*xcg*ycg - 80*ct**2*tt*xcg*ycg - 240*ct*ts*Xle*xcg*ycg - 240*ct*ts*Xte*xcg*ycg + 240*ct*ts*tt*xcg*ycg + 240*ct*ts*xcg*ycg - 240*ct*Xle*twb*xcg*ycg - 240*ct*Xte*twb*xcg*ycg + 240*ct*twb*tt*xcg*ycg + 240*ct*twb*xcg*ycg - 960*ts**2*xcg*ycg - 1920*ts*twb*xcg*ycg - 960*twb**2*xcg*ycg)/240 

# Ixz_new = b_2*p*xcg*zcg*(2*cr**2*Xle*tr + 2*cr**2*Xte*tr - 2*cr**2*tr + cr*ct*Xle*tr + cr*ct*Xle*tt + cr*ct*Xte*tr + cr*ct*Xte*tt - cr*ct*tr - cr*ct*tt - 6*cr*ts*Xle - 6*cr*ts*Xte + 6*cr*ts*tr + 6*cr*ts - 6*cr*Xle*twb - 6*cr*Xte*twb + 6*cr*twb*tr + 6*cr*twb + 2*ct**2*Xle*tt + 2*ct**2*Xte*tt - 2*ct**2*tt - 6*ct*ts*Xle - 6*ct*ts*Xte + 6*ct*ts*tt + 6*ct*ts - 6*ct*Xle*twb - 6*ct*Xte*twb + 6*ct*twb*tt + 6*ct*twb - 24*ts**2 - 48*ts*twb - 24*twb**2)/6 

# Iyz_new = b_2*p*ycg*zcg*(2*cr**2*Xle*tr + 2*cr**2*Xte*tr - 2*cr**2*tr + cr*ct*Xle*tr + cr*ct*Xle*tt + cr*ct*Xte*tr + cr*ct*Xte*tt - cr*ct*tr - cr*ct*tt - 6*cr*ts*Xle - 6*cr*ts*Xte + 6*cr*ts*tr + 6*cr*ts - 6*cr*Xle*twb - 6*cr*Xte*twb + 6*cr*twb*tr + 6*cr*twb + 2*ct**2*Xle*tt + 2*ct**2*Xte*tt - 2*ct**2*tt - 6*ct*ts*Xle - 6*ct*ts*Xte + 6*ct*ts*tt + 6*ct*ts - 6*ct*Xle*twb - 6*ct*Xte*twb + 6*ct*twb*tt + 6*ct*twb - 24*ts**2 - 48*ts*twb - 24*twb**2)/6 




p_val    = 0.5
b_2_val  = 2.0
Lr_val   = np.deg2rad(10.0)
ct_val   = 1.0
cr_val   = 1.5
tt_val   = 0.12
tr_val   = 0.12
ts_val   = 0.125 / 12.0
twb_val  = 0.25 / 12.0
Xle_val = 0.3
Xte_val = 0.4
# p_val    = 2.67851396616085
# b_2_val  = 1.28534905532406
# Lr_val   = np.deg2rad(11.8144504850503)
# ct_val   = 1.07501694006085
# cr_val   = 1.14583333333333
# tt_val   = 0.12
# tr_val   = 0.12
# ts_val   = 0.8 / 10.0 / 2.54 / 12.0
# twb_val  = 3.2 / 10.0 / 2.54 / 12.0
# Xle_val = 0.39 # from 0.39 to 0.49 x/c
# Xte_val = 0.51

replace_vals = [
    [p,p_val],
    [b_2,b_2_val],
    [cr,cr_val],
    [ct,ct_val],
    [tr,tr_val],
    [tt,tt_val],
    [Lr,Lr_val],
    [p,p_val],
    [ts,ts_val],
    [twb,twb_val],
    [Xle,Xle_val],
    [Xte,Xte_val]
]

array = [p_a,p_b,p_c,p_d,p_e,p_f,p_g,p_h,p_i,p_j,p_k,p_l,p_m]
for i in range(len(array)):
    print("p",chr(ord("a") + i),"=",array[i].subs(replace_vals))
array = [X_0,X_1,X_2,X_3]
for i in range(len(array)):
    print("X",i,"=",array[i].subs(replace_vals))
array = [h_a,h_b,h_c,h_d,h_e,h_f,h_g]
for i in range(len(array)):
    print("h",chr(ord("a") + i),"=",array[i].subs(replace_vals))



# m_sout = m.subs([[ts,0],[twb,0],[Xle,0],[Xte,0]])
# m_sinn = m.subs([[twb,0],[Xle,0],[Xte,0]])
# m_wout = m.subs([[twb,0]])
# m_winn = m

# eq_m = simp(m_sout - m_sinn + m_wout - m_winn)
m_val = m_all.subs(replace_vals)
m_tru = 0.0471353997
m_pdf = (m_tru-m_val) / m_tru * 100.0
m_co_val = m.subs(replace_vals)
m_co_tru = 0.02809375
m_co_pdf = float((m_co_tru-m_co_val) / m_co_tru * 100.0)

# print(eq_m)
# print(eq_m.subs([[ts,0],[twb,0],[Xle,0],[Xte,0]]))
# print(eq_m.subs(replace_vals))

# Myz_sout = Myz.subs([[ts,0],[twb,0],[Xle,0],[Xte,0]])
# Myz_sinn = Myz.subs([[twb,0],[Xle,0],[Xte,0]])
# Myz_wout = Myz.subs([[twb,0]])
# Myz_winn = Myz

# eq_Myz = simp(Myz_sout - Myz_sinn + Myz_wout - Myz_winn)
xcg_val = Myz_all.subs(replace_vals) / m_val
xcg_tru = -0.45582142
xcg_pdf = (xcg_tru-xcg_val) / xcg_tru * 100.0
xcg_co_val = Myz.subs(replace_vals) / m_co_val
xcg_co_tru = -0.4024386751239950
xcg_co_pdf = (xcg_co_tru-xcg_co_val) / xcg_co_tru * 100.0

# print(Myz)
# print()
# print(eq_Myz)
# print(eq_Myz.subs([[ts,0],[twb,0],[Xle,0],[Xte,0]]))
# print(eq_Myz.subs(replace_vals))
# print(Myz_all.subs(replace_vals))

# Mxz_sout = Mxz.subs([[ts,0],[twb,0],[Xle,0],[Xte,0]])
# Mxz_sinn = Mxz.subs([[twb,0],[Xle,0],[Xte,0]])
# Mxz_wout = Mxz.subs([[twb,0]])
# Mxz_winn = Mxz

# eq_Mxz = simp(Mxz_sout - Mxz_sinn + Mxz_wout - Mxz_winn)
ycg_val = Mxz_all.subs(replace_vals) / m_val
ycg_tru = 0.92780847
ycg_pdf = (ycg_tru-ycg_val) / ycg_tru * 100.0
ycg_co_val = Mxz.subs(replace_vals) / m_co_val
ycg_co_tru = 0.8109010011123470
ycg_co_pdf = (ycg_co_tru-ycg_co_val) / ycg_co_tru * 100.0

# print(eq_Mxz)
# print(eq_Mxz.subs([[ts,0],[twb,0],[Xle,0],[Xte,0]]))
# print(eq_Mxz.subs(replace_vals))
# print(Mxz_all.subs(replace_vals))

try:
    zcg_val = float( float(Mxy_all.subs(replace_vals)) / m_val )
except:
    zcg_val = float( float(Mxy_all) / m_val )
zcg_tru = 0.0
zcg_pdf = float(zcg_tru-zcg_val)
try:
    zcg_co_val = float( float(Mxy.subs(replace_vals)) / m_co_val )
except:
    zcg_co_val = float( float(Mxy) / m_co_val )
zcg_co_tru = 0.0
zcg_co_pdf = float(zcg_co_tru-zcg_co_val)

# settle replacement vals for cg location
cgloc = [[xcg,float(xcg_val)],[ycg,float(ycg_val)],[zcg,float(zcg_val)]]
cgloc_co = [[xcg,float(xcg_co_val)],[ycg,float(ycg_co_val)],[zcg,float(zcg_co_val)]]
# print(replace_vals + cgloc)
# replace_vals.append([xcg,float(xcg_val)])
# replace_vals.append([ycg,float(ycg_val)])
# replace_vals.append([zcg,float(zcg_val)])


# Ixx_sout = Ixx_new.subs([[ts,0],[twb,0],[Xle,0],[Xte,0]])
# Ixx_sinn = Ixx_new.subs([[twb,0],[Xle,0],[Xte,0]])
# Ixx_wout = Ixx_new.subs([[twb,0]])
# Ixx_winn = Ixx_new

# eq_Ixx = simp(Ixx_sout - Ixx_sinn + Ixx_wout - Ixx_winn)
Ixx_val = Ixx_all.subs(replace_vals + cgloc)
Ixx_tru = 0.015651990740
Ixx_pdf = (Ixx_tru-Ixx_val) / Ixx_tru * 100.0
Ixx_co_val = Ixx_new.subs(replace_vals + cgloc_co)
Ixx_co_tru = 0.0085809608620292
Ixx_co_pdf = (Ixx_co_tru-Ixx_co_val) / Ixx_co_tru * 100.0


# Iyy_sout = Iyy_new.subs([[ts,0],[twb,0],[Xle,0],[Xte,0]])
# Iyy_sinn = Iyy_new.subs([[twb,0],[Xle,0],[Xte,0]])
# Iyy_wout = Iyy_new.subs([[twb,0]])
# Iyy_winn = Iyy_new

# eq_Iyy = simp(Iyy_sout - Iyy_sinn + Iyy_wout - Iyy_winn)
Iyy_val = Iyy_all.subs(replace_vals + cgloc)
Iyy_tru = 0.005327395412 #### -0.0008757176284707 ####
Iyy_pdf = (Iyy_tru-Iyy_val) / Iyy_tru * 100.0
Iyy_co_val = Iyy_new.subs(replace_vals + cgloc_co)
Iyy_co_tru = 0.0004114335651049
Iyy_co_pdf = (Iyy_co_tru-Iyy_co_val) / Iyy_co_tru * 100.0


# Izz_sout = Izz_new.subs([[ts,0],[twb,0],[Xle,0],[Xte,0]])
# Izz_sinn = Izz_new.subs([[twb,0],[Xle,0],[Xte,0]])
# Izz_wout = Izz_new.subs([[twb,0]])
# Izz_winn = Izz_new

# eq_Izz = simp(Izz_sout - Izz_sinn + Izz_wout - Izz_winn)
Izz_val = Izz_all.subs(replace_vals + cgloc)
Izz_tru = 0.02060770125 #### 0.0019527473954077 ####
Izz_pdf = (Izz_tru-Izz_val) / Izz_tru * 100.0
Izz_co_val = Izz_new.subs(replace_vals + cgloc_co)
Izz_co_tru = 0.0089504625065612
Izz_co_pdf = (Izz_co_tru-Izz_co_val) / Izz_co_tru * 100.0


# Ixy_sout = Ixy_new.subs([[ts,0],[twb,0],[Xle,0],[Xte,0]])
# Ixy_sinn = Ixy_new.subs([[twb,0],[Xle,0],[Xte,0]])
# Ixy_wout = Ixy_new.subs([[twb,0]])
# Ixy_winn = Ixy_new

# eq_Ixy = simp(Ixy_sout - Ixy_sinn + Ixy_wout - Ixy_winn)
Ixy_val = Ixy_all.subs(replace_vals + cgloc)
Ixy_tru = -0.001843940138
Ixy_pdf = (Ixy_tru-Ixy_val) / Ixy_tru * 100.0
Ixy_co_val = Ixy_new.subs(replace_vals + cgloc_co)
Ixy_co_tru = -0.0010813583108170
Ixy_co_pdf = (Ixy_co_tru-Ixy_co_val) / Ixy_co_tru * 100.0

# Ixz_sout = Ixz_new.subs([[ts,0],[twb,0],[Xle,0],[Xte,0]])
# Ixz_sinn = Ixz_new.subs([[twb,0],[Xle,0],[Xte,0]])
# Ixz_wout = Ixz_new.subs([[twb,0]])
# Ixz_winn = Ixz_new

# eq_Ixz = simp(Ixz_sout - Ixz_sinn + Ixz_wout - Ixz_winn)
Ixz_val = float(Ixz_all.subs(replace_vals + cgloc))
Ixz_tru = 0.0
Ixz_pdf = Ixz_tru-Ixz_val
Ixz_co_val = float(Ixz_new.subs(replace_vals + cgloc_co))
Ixz_co_tru = 0.0
Ixz_co_pdf = Ixz_co_tru-Ixz_co_val


# Iyz_sout = Iyz_new.subs([[ts,0],[twb,0],[Xle,0],[Xte,0]])
# Iyz_sinn = Iyz_new.subs([[twb,0],[Xle,0],[Xte,0]])
# Iyz_wout = Iyz_new.subs([[twb,0]])
# Iyz_winn = Iyz_new

# eq_Iyz = simp(Iyz_sout - Iyz_sinn + Iyz_wout - Iyz_winn)
Iyz_val = float(Iyz_all.subs(replace_vals + cgloc))
Iyz_tru = 0.0
Iyz_pdf = Iyz_tru-Iyz_val
Iyz_co_val = float(Iyz_new.subs(replace_vals + cgloc_co))
Iyz_co_tru = 0.0
Iyz_co_pdf = Iyz_co_tru-Iyz_co_val

# print(eq_Iyz)
# print(eq_Iyz.subs([[ts,0],[twb,0],[Xle,0],[Xte,0]]))
# print(eq_Iyz.subs(replace_vals))
# print(Iyz_all.subs(replace_vals))


print("reporting...")
print("  m_all= {:< 19.16f} ?? {:< 11.8f}, %= {:> 13.8f}".format(\
    m_val,m_tru,m_pdf))
print("     co  {:< 19.16f} ?? {:< 11.8f}, %= {:> 13.8f}\n".format(\
    m_co_val,m_co_tru,m_co_pdf))

print("xcg_all= {:< 19.16f} ?? {:< 11.8f}, %= {:> 13.8f}".format(\
    xcg_val,xcg_tru,xcg_pdf))
print("     co  {:< 19.16f} ?? {:< 11.8f}, %= {:> 13.8f}".format(\
    xcg_co_val,xcg_co_tru,xcg_co_pdf))
print("ycg_all= {:< 19.16f} ?? {:< 11.8f}, %= {:> 13.8f}".format(\
    ycg_val,ycg_tru,ycg_pdf))
print("     co  {:< 19.16f} ?? {:< 11.8f}, %= {:> 13.8f}".format(\
    ycg_co_val,ycg_co_tru,ycg_co_pdf))
print("zcg_all= {:< 19.16f} ?? {:< 11.8f}, %= {:> 13.8f}".format(\
    zcg_val,zcg_tru,zcg_pdf))
print("     co  {:< 19.16f} ?? {:< 11.8f}, %= {:> 13.8f}\n".format(\
    zcg_co_val,zcg_co_tru,zcg_co_pdf))

print("Ixx_all= {:< 19.16f} ?? {:< 11.8f}, %= {:> 13.8f}".format(\
    Ixx_val,Ixx_tru,Ixx_pdf))
print("     co  {:< 19.16f} ?? {:< 11.8f}, %= {:> 13.8f}".format(\
    Ixx_co_val,Ixx_co_tru,Ixx_co_pdf))
print("Iyy_all= {:< 19.16f} ?? {:< 11.8f}, %= {:> 13.8f}".format(\
    Iyy_val,Iyy_tru,Iyy_pdf))
print("     co  {:< 19.16f} ?? {:< 11.8f}, %= {:> 13.8f}".format(\
    Iyy_co_val,Iyy_co_tru,Iyy_co_pdf))
print("Izz_all= {:< 19.16f} ?? {:< 11.8f}, %= {:> 13.8f}".format(\
    Izz_val,Izz_tru,Izz_pdf))
print("     co  {:< 19.16f} ?? {:< 11.8f}, %= {:> 13.8f}\n".format(\
    Izz_co_val,Izz_co_tru,Izz_co_pdf))

print("Ixy_all= {:< 19.16f} ?? {:< 11.8f}, %= {:> 13.8f}".format(\
    Ixy_val,Ixy_tru,Ixy_pdf))
print("     co  {:< 19.16f} ?? {:< 11.8f}, %= {:> 13.8f}".format(\
    Ixy_co_val,Ixy_co_tru,Ixy_co_pdf))
print("Ixz_all= {:< 19.16f} ?? {:< 11.8f}, %= {:> 13.8f}".format(\
    Ixz_val,Ixz_tru,Ixz_pdf))
print("     co  {:< 19.16f} ?? {:< 11.8f}, %= {:> 13.8f}".format(\
    Ixz_co_val,Ixz_co_tru,Ixz_co_pdf))
print("Iyz_all= {:< 19.16f} ?? {:< 11.8f}, %= {:> 13.8f}".format(\
    Iyz_val,Iyz_tru,Iyz_pdf))
print("     co  {:< 19.16f} ?? {:< 11.8f}, %= {:> 13.8f}".format(\
    Iyz_co_val,Iyz_co_tru,Iyz_co_pdf))