import sympy as sy
import numpy as np
from matplotlib import pyplot as plt
import time as t

# current time
start_time = t.time()

sym = sy.Symbol
igr = sy.integrate
simp = sy.simplify
exp = sy.expand
piecewise = sy.Piecewise
diff = sy.diff
# sy.init_printing(use_unicode=True)

# declare variables
print("declaring variables...")
x = sym("x")
y = sym("y")
z = sym("z")
xcg = sym("x_cg")
ycg = sym("y_cg")
zcg = sym("z_cg")
p = sym("p")
b_2 = sym("b/2")
L = sym("\u039B")
Lr = sym("\u039B_r")
tanL = sy.tan(Lr)
# tanL = sy.tan(L * y / b_2 + Lr)
ct = sym("c_t")
cr = sym("c_r")
tt = sym("t_t")
tr = sym("t_r")
a0 = sym("a_0")
a1 = sym("a_1")
a2 = sym("a_2")
a3 = sym("a_3")
a4 = sym("a_4")
w = sym("\u03C9")
wr = sym("\u03C9_r")
S = sy.sin( w * y / b_2 + wr)
C = sy.cos( w * y / b_2 + wr)
sk = sym("sk")
spk = sym("spk")
spLE = sym("spLE")
spTE = sym("spTE")
# sk = 0


# make only mounting angle
S = S.subs(w,0); C = C.subs(w,0)
# make no twist
S = S.subs(wr,0); C = C.subs(wr,0)

# create bounds
print("creating bounds...")
# x_up_wL = (-tanL + sy.Rational(1,4)*(ct-cr)/b_2) * y + sy.Rational(1,4)*cr - sk
# x_lo_wL = (-tanL - sy.Rational(3,4)*(ct-cr)/b_2) * y - sy.Rational(3,4)*cr + sk
x_up =   (sy.Rational(1,4) - spLE)*( (ct-cr)/b_2 * y + cr ) - sk - spk
x_lo = - (sy.Rational(3,4) - spTE)*( (ct-cr)/b_2 * y + cr ) + sk + spk
ch = simp(x_up - x_lo)
t_eq = a0 * (x/ch)**0.5 + a1*(x/ch) + a2*(x/ch)**2 + a3*(x/ch)**3 +a4*(x/ch)**4
z_up =   sy.Rational(1,2) * ((tt*ct-tr*cr)/b_2 * y + tr*cr) - sk - spk # * t_eq
z_lo = - sy.Rational(1,2) * ((tt*ct-tr*cr)/b_2 * y + tr*cr) + sk + spk # * t_eq
y_up = b_2
y_lo = 0
x_bnd = (x,x_lo,x_up)
# x_bnwL = (x,x_lo_wL,x_up_wL)
y_bnd = (y,y_lo,y_up)
z_bnd = (z,z_lo,z_up)

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
Sm = S * 1; Cm = C * 1
# Sm = Sm.subs(w,0); Cm = Cm.subs(w,0)
# myz_twist = myz # # # Cm*myz + Sm*mxy
# mxz_twist = mxz # # # mxz
# mxy_twist = mxy # # # -Sm*myz + Cm*mxy
myz_twist = Cm*myz + Sm*mxy
mxz_twist = mxz
mxy_twist = -Sm*myz + Cm*mxy
Lyz_twist = Lyz # # # Cm*Lyz + Sm*Lxy
Lxz_twist = Lxz # # # Lxz
Lxy_twist = Lxy # # # -Sm*Lyz + Cm*Lxy
# Lyz_twist = Cm*Lyz + Sm*Lxy
# Lxz_twist = Lxz
# Lxy_twist = -Sm*Lyz + Cm*Lxy
Myz = igr(myz_twist, x_bnd,z_bnd,y_bnd) + igr(Lyz_twist, x_bnd,z_bnd,y_bnd)
Mxz = igr(mxz_twist, x_bnd,z_bnd,y_bnd) + igr(Lxz_twist, x_bnd,z_bnd,y_bnd)
Mxy = igr(mxy_twist, x_bnd,z_bnd,y_bnd) + igr(Lxy_twist, x_bnd,z_bnd,y_bnd)
Myz_no_twist = igr(myz, x_bnd,z_bnd,y_bnd) + igr(Lyz, x_bnd,z_bnd,y_bnd)
Mxz_no_twist = igr(mxz, x_bnd,z_bnd,y_bnd) + igr(Lxz, x_bnd,z_bnd,y_bnd)
Mxy_no_twist = igr(mxy, x_bnd,z_bnd,y_bnd) + igr(Lxy, x_bnd,z_bnd,y_bnd)

m = sy.factor(m)
Myz = sy.factor(Myz)
Mxz = sy.factor(Mxz)
Mxy = sy.factor(Mxy)
Myz_no_twist = sy.factor(Myz_no_twist)
Mxz_no_twist = sy.factor(Mxz_no_twist)
Mxy_no_twist = sy.factor(Mxy_no_twist)
# print("Myz =", Myz, "\n")
# print("Mxz =", Mxz, "\n")
# print("Mxy =", Mxy, "\n")
# print("Myzn =", Myz_no_twist, "\n")
# print("Mxzn =", Mxz_no_twist, "\n")
# print("Mxyn =", Mxy_no_twist, "\n")


# # replacement variables
# print("replacing...")
# b_2_val = 8.0
# cr_val = 2.0
# ct_val = 1.50
# tr_val = tt_val = 0.12
# Lr_val = np.deg2rad(15.0)
# p_val = 0.0175
# wr_val = np.deg2rad(10.0)
# w_val = np.deg2rad(10.0)
# sk_val = 0.0

# # organize for plugging in values and printing
# names = ["all", "only mount", "no twist","no sweep","noSwonlyM"]
# types = ["me","SW","PD"]
# val_single = [
#     [b_2,b_2_val],
#     [cr,cr_val],
#     [ct,ct_val],
#     [tr,tr_val],
#     [tt,tt_val],
#     [Lr,Lr_val],
#     [p,p_val],
#     [w,w_val],
#     [wr,wr_val],
#     [sk,sk_val]
# ]
# vals = [val_single] * len(names)
# vals = np.array(vals)
# vals[1,8,1] = vals[4,8,1] = np.deg2rad(20.0)
# vals[1,7,1] =vals[4,7,1] =vals[2,7,1] =vals[2,8,1] =vals[3,5,1]=vals[4,5,1]=0.0
# xcgs = np.zeros((len(names),len(types)))
# ycgs = np.zeros((len(names),len(types)))
# zcgs = np.zeros((len(names),len(types)))
# xcgs[0,1] = -1.40572654; ycgs[0,1] = 3.62162115; zcgs[0,1] = 0.07463028
# xcgs[1,1] = -1.38708173; ycgs[1,1] = 3.62162162; zcgs[1,1] = 0.15165589
# xcgs[2,1] = -1.41382275; ycgs[2,1] = 3.62162162; zcgs[2,1] = 0.0
# xcgs[3,1] = -0.43531313; ycgs[3,1] = 3.62162081; zcgs[3,1] = 0.07465744

# for i in range(len(names)):
#     m_i = m.subs(vals[i])
#     Myz_i = Myz.subs(vals[i])
#     Mxz_i = Mxz.subs(vals[i])
#     Mxy_i = Mxy.subs(vals[i])

#     xcgs[i,0] = Myz_i / m_i
#     ycgs[i,0] = Mxz_i / m_i
#     zcgs[i,0] = Mxy_i / m_i
    
#     # percent difference
#     if xcgs[i,1] != 0.0:
#         xcgs[i,2] = (xcgs[i,1]-xcgs[i,0]) / xcgs[i,1] * 100.0
#     if ycgs[i,1] != 0.0:
#         ycgs[i,2] = (ycgs[i,1]-ycgs[i,0]) / ycgs[i,1] * 100.0
#     if zcgs[i,1] != 0.0:
#         zcgs[i,2] = (zcgs[i,1]-zcgs[i,0]) / zcgs[i,1] * 100.0

def print_table(name,namelist,typelist,values):
    print("{:=^5s}".format(name),end="")

    widths = 17

    # print column headers
    for i in range(len(namelist)):
        print(" {:-^{}s}".format(namelist[i],widths),end="")
    print()

    # print values
    for i in range(len(typelist)):
        print("{:^5s}".format(typelist[i]),end="")
        for j in range(len(namelist)):
            print(" {:> {}.12f}".format(values[j][i],widths),end="")
        print()

    print("{: ^5s}".format(""),end="")
    for j in range(len(namelist)):
        if values[j,2] <= -1.0E-4:
            err = "-" * 5
        elif values[j,2] >=  1.0E-4:
            err = "+" * 5
        else:
            err = ""
        print(" {: ^{}s}".format(err,widths),end="")
    print()

# print("m   =",   m, "\n")
# print("Myz =", Myz, "\n")
# print("Mxz =", Mxz, "\n")
# print("Mxy =", Mxy, "\n")

# print_table("xcg",names,types,xcgs)
# print_table("ycg",names,types,ycgs)
# print_table("zcg",names,types,zcgs)

# print()


# sub inertias
print("\t sub inertias...")
ixx = igr( p*(y**2+z**2), x_bnd,z_bnd)
iyy = igr( p*(x**2+z**2), x_bnd,z_bnd)
izz = igr( p*(x**2+y**2), x_bnd,z_bnd)
ixy = igr( p*(x*y), x_bnd,z_bnd)
ixz = igr( p*(x*z), x_bnd,z_bnd)
iyz = igr( p*(y*z), x_bnd,z_bnd)
# # twisted sub inertias
# ixx_twist = ixx
# iyy_twist = iyy
# izz_twist = izz
# ixy_twist = ixy
# ixz_twist = ixz
# iyz_twist = iyz
# twist equations
ixx_twist = C*C*ixx - 2*C*S*ixz + S*S*izz
iyy_twist = iyy
izz_twist = S*S*ixx + 2*C*S*ixz + C*C*izz
ixy_twist = C*ixy + S*iyz # negated
ixz_twist = C*S*(ixx - izz) + ixz*(C*C - S*S) # negated
iyz_twist = C*iyz - S*ixy # negated

# print("i__ twist before y integration")
# print(ixx_twist)
# print(iyy_twist)
# print(izz_twist)
# print(ixy_twist)
# print(ixz_twist)
# print(iyz_twist)
# print()

# # old inertias
# Ixx = (igr( p*(y**2+z**2), x_bnwL,z_bnd,y_bnd) / m - ycg**2 - zcg**2) * m
# Iyy = (igr( p*(x**2+z**2), x_bnwL,z_bnd,y_bnd) / m - xcg**2 - zcg**2) * m
# Izz = (igr( p*(x**2+y**2), x_bnwL,z_bnd,y_bnd) / m - xcg**2 - ycg**2) * m
# Ixy = (igr( p*(x*y), x_bnwL,z_bnd,y_bnd) / m - (xcg*ycg)) * m
# Ixz = (igr( p*(x*z), x_bnwL,z_bnd,y_bnd) / m - (xcg*zcg)) * m
# Iyz = (igr( p*(y*z), x_bnwL,z_bnd,y_bnd) / m - (ycg*zcg)) * m

# new inertias
print("\t full inertias...")
print("\t\t Ixx...")
Ixx_new = igr(ixx_twist,y_bnd) - m * (ycg**2 + zcg**2)
print("\t\t Iyy...")
Iyy_new = igr(iyy_twist,y_bnd) \
    - p*igr( -y**2*tanL**2 + 2*y*x*tanL, x_bnd,z_bnd,y_bnd) \
    - m * (xcg**2 + zcg**2)
print("\t\t Izz...")
Izz_new = igr(izz_twist,y_bnd) \
    - p*igr( -y**2*tanL**2 + 2*y*x*tanL, x_bnd,z_bnd,y_bnd) \
    - m * (xcg**2 + ycg**2)
print("\t\t Ixy...")
Ixy_new = igr(ixy_twist,y_bnd) - p*igr( y**2*tanL, x_bnd,z_bnd,y_bnd) \
    - m * (xcg * ycg)
print("\t\t Ixz...")
Ixz_new = igr(ixz_twist,y_bnd)  - p*igr( y*z*tanL, x_bnd,z_bnd,y_bnd) \
    - m * (xcg * zcg)
print("\t\t Iyz...")
Iyz_new = (igr(iyz_twist,y_bnd) / m - (ycg*zcg)) * m
# # differences
# Ixx_diff = simp(Ixx - Ixx_new)
# Iyy_diff = simp(Iyy - Iyy_new)
# Izz_diff = simp(Izz - Izz_new)
# Ixy_diff = simp(Ixy - Ixy_new)
# Ixz_diff = simp(Ixz - Ixz_new)
# Iyz_diff = simp(Iyz - Iyz_new)

# simplify expressions
# print("skipping simplifying...")
print("simplifying...")
print("\t mass and moments...")
m = sy.factor(m)#ksubs(m))
Myz = sy.factor(Myz)
Mxz = sy.factor(Mxz)
Mxy = sy.factor(Mxy)
print("\t moments of inertia...")
Ixx_new = sy.factor(simp(sy.trigsimp(Ixx_new)))
Iyy_new = sy.factor(simp(sy.trigsimp(Iyy_new)))
Izz_new = sy.factor(simp(sy.trigsimp(Izz_new)))
print("\t products of inertia...")
Ixy_new = sy.factor(simp(sy.trigsimp(Ixy_new)))
Ixz_new = sy.factor(simp(sy.trigsimp(Ixz_new)))
Iyz_new = sy.factor(simp(sy.trigsimp(Iyz_new)))

# assign results
print("reporting...\n\n")
print("m =", m, "\n")
print("Myz =", Myz, "\n")
print("Mxz =", Mxz, "\n")
print("Mxy =", Mxy, "\n")
# print("Ixx_old =", Ixx, "\n")
# print("Iyy_old =", Iyy, "\n")
# print("Izz_old =", Izz, "\n")
# print("Ixy_old =", Ixy, "\n")
# print("Ixz_old =", Ixz, "\n")
# print("Iyz_old =", Iyz, "\n")

# print("diff =", Ixx_diff, "\n")
# print("diff =", Iyy_diff, "\n")
# print("diff =", Izz_diff, "\n")
# print("diff =", Ixy_diff, "\n")
# print("diff =", Ixz_diff, "\n")
# print("diff =", Iyz_diff, "\n")

print("Ixx_new =", Ixx_new, "\n")
print("Iyy_new =", Iyy_new, "\n")
print("Izz_new =", Izz_new, "\n")
print("Ixy_new =", Ixy_new, "\n")
print("Ixz_new =", Ixz_new, "\n")
print("Iyz_new =", Iyz_new, "\n")

print("Time elapse:",t.time()-start_time,"seconds")

# p_val    = 0.5
# b_2_val  = 2.0
# Lr_val   = np.deg2rad(10.0)
# ct_val   = 1.0
# cr_val   = 1.5
# tt_val   = 0.12
# tr_val   = 0.12
# sk_val   = 0.125 / 12.0
# spk_val  = 0.25 / 12.0
# spLE_val = 0.3
# spTE_val = 0.4
p_val    = 2.565855969
b_2_val  = 1.285349055
Lr_val   = np.deg2rad(11.81445049)
ct_val   = 1.07501694
cr_val   = 1.145833333
tt_val   = 0.12
tr_val   = 0.12
sk_val   = 0.8 / 10.0 / 2.54 / 12.0
spk_val  = 3.2 / 10.0 / 2.54 / 12.0
spLE_val = 0.39 # from 0.39 to 0.49 x/c
spTE_val = 0.51

replace_vals = [
    [b_2,b_2_val],
    [cr,cr_val],
    [ct,ct_val],
    [tr,tr_val],
    [tt,tt_val],
    [Lr,Lr_val],
    [p,p_val]
]

changing_vals = [
[   [sk,0.0],
    [spk,0.0],
    [spLE,0.0],
    [spTE,0.0]],
[   [sk,sk_val],
    [spk,0.0],
    [spLE,0.0],
    [spTE,0.0]],
[   [sk,sk_val],
    [spk,0.0],
    [spLE,spLE_val],
    [spTE,spTE_val]],
[   [sk,sk_val],
    [spk,spk_val],
    [spLE,spLE_val],
    [spTE,spTE_val]]
]

# initialize arrays
m_arr   = np.zeros(4)
Myz_arr = np.zeros(4)
Mxz_arr = np.zeros(4)
Mxy_arr = np.zeros(4)
Ixx_arr = np.zeros(4)
Iyy_arr = np.zeros(4)
Izz_arr = np.zeros(4)
Ixy_arr = np.zeros(4)
Ixz_arr = np.zeros(4)
Iyz_arr = np.zeros(4)

for i in range(4):

    # replace mass and moments
    m_arr[i]   =   m.subs(replace_vals).subs(changing_vals[i])
    Myz_arr[i] = Myz.subs(replace_vals).subs(changing_vals[i])
    Mxz_arr[i] = Mxz.subs(replace_vals).subs(changing_vals[i])
    Mxy_arr[i] = Mxy.subs(replace_vals).subs(changing_vals[i])

# negate some
index = [1,3]
m_arr[index] *= -1.0
Myz_arr[index] *= -1.0
Mxz_arr[index] *= -1.0
Mxy_arr[index] *= -1.0

m       = np.sum(  m_arr)
xcg_val = np.sum(Myz_arr) / m
ycg_val = np.sum(Mxz_arr) / m
zcg_val = np.sum(Mxy_arr) / m

# affix cg location
replace_vals2 = [
    [xcg,xcg_val],
    [ycg,ycg_val],
    [zcg,zcg_val]
]

# replace inertias
for i in range(4):
    Ixx_arr[i] = Ixx_new.subs(replace_vals).subs(changing_vals[i]).subs(replace_vals2)
    Iyy_arr[i] = Iyy_new.subs(replace_vals).subs(changing_vals[i]).subs(replace_vals2)
    Izz_arr[i] = Izz_new.subs(replace_vals).subs(changing_vals[i]).subs(replace_vals2)
    Ixy_arr[i] = Ixy_new.subs(replace_vals).subs(changing_vals[i]).subs(replace_vals2)
    Ixz_arr[i] = Ixz_new.subs(replace_vals).subs(changing_vals[i]).subs(replace_vals2)
    Iyz_arr[i] = Iyz_new.subs(replace_vals).subs(changing_vals[i]).subs(replace_vals2)

# negate some
Ixy_arr *= -1.0
Iyz_arr *= -1.0
Ixx_arr[index] *= -1.0
Iyy_arr[index] *= -1.0
Izz_arr[index] *= -1.0
Ixy_arr[index] *= -1.0
Ixz_arr[index] *= -1.0
Iyz_arr[index] *= -1.0

Ixx_new = np.sum(Ixx_arr)
Iyy_new = np.sum(Iyy_arr)
Izz_new = np.sum(Izz_arr)
Ixy_new = np.sum(Ixy_arr)
Ixz_new = np.sum(Ixz_arr)
Iyz_new = np.sum(Iyz_arr)

print(m_arr)
print(Myz_arr / m_arr)
print(Mxz_arr / m_arr)
print(Mxy_arr / m_arr)

print("reporting with values...")
print("m =", m)
print("Myz =", xcg_val / m)
print("Mxz =", ycg_val / m)
print("Mxy =", zcg_val / m)
print("xcg =", xcg_val)
print("ycg =", ycg_val)
print("zcg =", zcg_val)
print("Ixx_new =", Ixx_new)
print("Iyy_new =", Iyy_new)
print("Izz_new =", Izz_new)
print("Ixy_new =", Ixy_new)
print("Ixz_new =", Ixz_new)
print("Iyz_new =", Iyz_new)
print()
for inn in [Ixx_new,Ixy_new,Ixz_new]:
    print("\t{:< 19.16f}".format(inn),end="")
print()
for inn in [Ixy_new,Iyy_new,Iyz_new]:
    print("\t{:< 19.16f}".format(inn),end="")
print()
for inn in [Ixz_new,Iyz_new,Izz_new]:
    print("\t{:< 19.16f}".format(inn),end="")
print()

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
# 	Lxx = 0.50358715	Lxy = -0.05932693	Lxz = 0.00000000
# 	Lyx = -0.05932693	Lyy = 0.17140362	Lyz = 0.00000000
# 	Lzx = 0.00000000	Lzy = 0.00000000	Lzz = 0.66303218
# slugs * ft^2
#     Lxx = 0.015651990740	Lxy = 0.001843940138	Lxz = 0.00000000000
#     Lyx = 0.001843940138	Lyy = 0.005327395412	Lyz = 0.00000000000
#     Lzx = 0.000000000000	Lzy = 0.000000000000	Lzz = 0.02060770125

# Moments of inertia: ( pounds * square feet )
# Taken at the output coordinate system. (Using positive tensor notation.)
# 	Ixx = 1.80906323	Ixy = -0.70069202	Ixz = 0.00000000
# 	Iyx = -0.70069202	Iyy = 0.48649877	Iyz = 0.00000000
# 	Izx = 0.00000000	Izy = 0.00000000	Izz = 2.28360341


# Ba = tr*cr*(2*cr+ct) + tt*ct*(cr+2*ct)
# Bb = cr*(1 + tr) + ct*(1 + tt)

# m = b/2*p*( (1 - spLE - spTE)*Ba + 
#     6*(cr+ct)*(sk + spk)*(spLE + spTE)
#     - 6*Bb*(sk + spk) + 24*(sk + spk)**2 )/6 






# m = -b/2*p*(2*c_r**2*spLE*t_r + 2*c_r**2*spTE*t_r - 2*c_r**2*t_r + c_r*c_t*spLE*t_r + c_r*c_t*spLE*t_t + c_r*c_t*spTE*t_r + c_r*c_t*spTE*t_t - c_r*c_t*t_r - c_r*c_t*t_t - 6*c_r*sk*spLE - 6*c_r*sk*spTE + 6*c_r*sk*t_r + 6*c_r*sk - 6*c_r*spLE*spk - 6*c_r*spTE*spk + 6*c_r*spk*t_r + 6*c_r*spk + 2*c_t**2*spLE*t_t + 2*c_t**2*spTE*t_t - 2*c_t**2*t_t - 6*c_t*sk*spLE - 6*c_t*sk*spTE + 6*c_t*sk*t_t + 6*c_t*sk - 6*c_t*spLE*spk - 6*c_t*spTE*spk + 6*c_t*spk*t_t + 6*c_t*spk - 24*sk**2 - 48*sk*spk - 24*spk**2)/6 

# Myz = b/2*p*(4*b/2*c_r**2*spLE*t_r*tan(Λ_r) + 4*b/2*c_r**2*spTE*t_r*tan(Λ_r) - 4*b/2*c_r**2*t_r*tan(Λ_r) + 4*b/2*c_r*c_t*spLE*t_r*tan(Λ_r) + 4*b/2*c_r*c_t*spLE*t_t*tan(Λ_r) + 4*b/2*c_r*c_t*spTE*t_r*tan(Λ_r) + 4*b/2*c_r*c_t*spTE*t_t*tan(Λ_r) - 4*b/2*c_r*c_t*t_r*tan(Λ_r) - 4*b/2*c_r*c_t*t_t*tan(Λ_r) - 16*b/2*c_r*sk*spLE*tan(Λ_r) - 16*b/2*c_r*sk*spTE*tan(Λ_r) + 16*b/2*c_r*sk*t_r*tan(Λ_r) + 16*b/2*c_r*sk*tan(Λ_r) - 16*b/2*c_r*spLE*spk*tan(Λ_r) - 16*b/2*c_r*spTE*spk*tan(Λ_r) + 16*b/2*c_r*spk*t_r*tan(Λ_r) + 16*b/2*c_r*spk*tan(Λ_r) + 12*b/2*c_t**2*spLE*t_t*tan(Λ_r) + 12*b/2*c_t**2*spTE*t_t*tan(Λ_r) - 12*b/2*c_t**2*t_t*tan(Λ_r) - 32*b/2*c_t*sk*spLE*tan(Λ_r) - 32*b/2*c_t*sk*spTE*tan(Λ_r) + 32*b/2*c_t*sk*t_t*tan(Λ_r) + 32*b/2*c_t*sk*tan(Λ_r) - 32*b/2*c_t*spLE*spk*tan(Λ_r) - 32*b/2*c_t*spTE*spk*tan(Λ_r) + 32*b/2*c_t*spk*t_t*tan(Λ_r) + 32*b/2*c_t*spk*tan(Λ_r) - 96*b/2*sk**2*tan(Λ_r) - 192*b/2*sk*spk*tan(Λ_r) - 96*b/2*spk**2*tan(Λ_r) + 6*c_r**3*spLE**2*t_r - 3*c_r**3*spLE*t_r - 6*c_r**3*spTE**2*t_r + 9*c_r**3*spTE*t_r - 3*c_r**3*t_r + 4*c_r**2*c_t*spLE**2*t_r + 2*c_r**2*c_t*spLE**2*t_t - 2*c_r**2*c_t*spLE*t_r - c_r**2*c_t*spLE*t_t - 4*c_r**2*c_t*spTE**2*t_r - 2*c_r**2*c_t*spTE**2*t_t + 6*c_r**2*c_t*spTE*t_r + 3*c_r**2*c_t*spTE*t_t - 2*c_r**2*c_t*t_r - c_r**2*c_t*t_t - 16*c_r**2*sk*spLE**2 + 16*c_r**2*sk*spLE*t_r + 8*c_r**2*sk*spLE + 16*c_r**2*sk*spTE**2 - 16*c_r**2*sk*spTE*t_r - 24*c_r**2*sk*spTE + 8*c_r**2*sk*t_r + 8*c_r**2*sk - 16*c_r**2*spLE**2*spk + 16*c_r**2*spLE*spk*t_r + 8*c_r**2*spLE*spk + 16*c_r**2*spTE**2*spk - 16*c_r**2*spTE*spk*t_r - 24*c_r**2*spTE*spk + 8*c_r**2*spk*t_r + 8*c_r**2*spk + 2*c_r*c_t**2*spLE**2*t_r + 4*c_r*c_t**2*spLE**2*t_t - c_r*c_t**2*spLE*t_r - 2*c_r*c_t**2*spLE*t_t - 2*c_r*c_t**2*spTE**2*t_r - 4*c_r*c_t**2*spTE**2*t_t + 3*c_r*c_t**2*spTE*t_r + 6*c_r*c_t**2*spTE*t_t - c_r*c_t**2*t_r - 2*c_r*c_t**2*t_t - 16*c_r*c_t*sk*spLE**2 + 8*c_r*c_t*sk*spLE*t_r + 8*c_r*c_t*sk*spLE*t_t + 8*c_r*c_t*sk*spLE + 16*c_r*c_t*sk*spTE**2 - 8*c_r*c_t*sk*spTE*t_r - 8*c_r*c_t*sk*spTE*t_t - 24*c_r*c_t*sk*spTE + 4*c_r*c_t*sk*t_r + 4*c_r*c_t*sk*t_t + 8*c_r*c_t*sk - 16*c_r*c_t*spLE**2*spk + 8*c_r*c_t*spLE*spk*t_r + 8*c_r*c_t*spLE*spk*t_t + 8*c_r*c_t*spLE*spk + 16*c_r*c_t*spTE**2*spk - 8*c_r*c_t*spTE*spk*t_r - 8*c_r*c_t*spTE*spk*t_t - 24*c_r*c_t*spTE*spk + 4*c_r*c_t*spk*t_r + 4*c_r*c_t*spk*t_t + 8*c_r*c_t*spk - 48*c_r*sk**2*spLE + 48*c_r*sk**2*spTE - 24*c_r*sk**2 - 96*c_r*sk*spLE*spk + 96*c_r*sk*spTE*spk - 48*c_r*sk*spk - 48*c_r*spLE*spk**2 + 48*c_r*spTE*spk**2 - 24*c_r*spk**2 + 6*c_t**3*spLE**2*t_t - 3*c_t**3*spLE*t_t - 6*c_t**3*spTE**2*t_t + 9*c_t**3*spTE*t_t - 3*c_t**3*t_t - 16*c_t**2*sk*spLE**2 + 16*c_t**2*sk*spLE*t_t + 8*c_t**2*sk*spLE + 16*c_t**2*sk*spTE**2 - 16*c_t**2*sk*spTE*t_t - 24*c_t**2*sk*spTE + 8*c_t**2*sk*t_t + 8*c_t**2*sk - 16*c_t**2*spLE**2*spk + 16*c_t**2*spLE*spk*t_t + 8*c_t**2*spLE*spk + 16*c_t**2*spTE**2*spk - 16*c_t**2*spTE*spk*t_t - 24*c_t**2*spTE*spk + 8*c_t**2*spk*t_t + 8*c_t**2*spk - 48*c_t*sk**2*spLE + 48*c_t*sk**2*spTE - 24*c_t*sk**2 - 96*c_t*sk*spLE*spk + 96*c_t*sk*spTE*spk - 48*c_t*sk*spk - 48*c_t*spLE*spk**2 + 48*c_t*spTE*spk**2 - 24*c_t*spk**2)/48 

# Mxz = -b/2**2*p*(c_r**2*spLE*t_r + c_r**2*spTE*t_r - c_r**2*t_r + c_r*c_t*spLE*t_r + c_r*c_t*spLE*t_t + c_r*c_t*spTE*t_r + c_r*c_t*spTE*t_t - c_r*c_t*t_r - c_r*c_t*t_t - 4*c_r*sk*spLE - 4*c_r*sk*spTE + 4*c_r*sk*t_r + 4*c_r*sk - 4*c_r*spLE*spk - 4*c_r*spTE*spk + 4*c_r*spk*t_r + 4*c_r*spk + 3*c_t**2*spLE*t_t + 3*c_t**2*spTE*t_t - 3*c_t**2*t_t - 8*c_t*sk*spLE - 8*c_t*sk*spTE + 8*c_t*sk*t_t + 8*c_t*sk - 8*c_t*spLE*spk - 8*c_t*spTE*spk + 8*c_t*spk*t_t + 8*c_t*spk - 24*sk**2 - 48*sk*spk - 24*spk**2)/12 

# Mxy = 0 

# Ixx_new = -b/2*p*(8*b/2**2*c_r**2*spLE*t_r + 8*b/2**2*c_r**2*spTE*t_r - 8*b/2**2*c_r**2*t_r + 12*b/2**2*c_r*c_t*spLE*t_r + 12*b/2**2*c_r*c_t*spLE*t_t + 12*b/2**2*c_r*c_t*spTE*t_r + 12*b/2**2*c_r*c_t*spTE*t_t - 12*b/2**2*c_r*c_t*t_r - 12*b/2**2*c_r*c_t*t_t - 40*b/2**2*c_r*sk*spLE - 40*b/2**2*c_r*sk*spTE + 40*b/2**2*c_r*sk*t_r + 40*b/2**2*c_r*sk - 40*b/2**2*c_r*spLE*spk - 40*b/2**2*c_r*spTE*spk + 40*b/2**2*c_r*spk*t_r + 40*b/2**2*c_r*spk + 48*b/2**2*c_t**2*spLE*t_t + 48*b/2**2*c_t**2*spTE*t_t - 48*b/2**2*c_t**2*t_t - 120*b/2**2*c_t*sk*spLE - 120*b/2**2*c_t*sk*spTE + 120*b/2**2*c_t*sk*t_t + 120*b/2**2*c_t*sk - 120*b/2**2*c_t*spLE*spk - 120*b/2**2*c_t*spTE*spk + 120*b/2**2*c_t*spk*t_t + 120*b/2**2*c_t*spk - 320*b/2**2*sk**2 - 640*b/2**2*sk*spk - 320*b/2**2*spk**2 + 4*c_r**4*spLE*t_r**3 + 4*c_r**4*spTE*t_r**3 - 4*c_r**4*t_r**3 + c_r**3*c_t*spLE*t_r**3 + 3*c_r**3*c_t*spLE*t_r**2*t_t + c_r**3*c_t*spTE*t_r**3 + 3*c_r**3*c_t*spTE*t_r**2*t_t - c_r**3*c_t*t_r**3 - 3*c_r**3*c_t*t_r**2*t_t - 30*c_r**3*sk*spLE*t_r**2 - 30*c_r**3*sk*spTE*t_r**2 + 10*c_r**3*sk*t_r**3 + 30*c_r**3*sk*t_r**2 - 30*c_r**3*spLE*spk*t_r**2 - 30*c_r**3*spTE*spk*t_r**2 + 10*c_r**3*spk*t_r**3 + 30*c_r**3*spk*t_r**2 + 2*c_r**2*c_t**2*spLE*t_r**2*t_t + 2*c_r**2*c_t**2*spLE*t_r*t_t**2 + 2*c_r**2*c_t**2*spTE*t_r**2*t_t + 2*c_r**2*c_t**2*spTE*t_r*t_t**2 - 2*c_r**2*c_t**2*t_r**2*t_t - 2*c_r**2*c_t**2*t_r*t_t**2 - 10*c_r**2*c_t*sk*spLE*t_r**2 - 20*c_r**2*c_t*sk*spLE*t_r*t_t - 10*c_r**2*c_t*sk*spTE*t_r**2 - 20*c_r**2*c_t*sk*spTE*t_r*t_t + 10*c_r**2*c_t*sk*t_r**2*t_t + 10*c_r**2*c_t*sk*t_r**2 + 20*c_r**2*c_t*sk*t_r*t_t - 10*c_r**2*c_t*spLE*spk*t_r**2 - 20*c_r**2*c_t*spLE*spk*t_r*t_t - 10*c_r**2*c_t*spTE*spk*t_r**2 - 20*c_r**2*c_t*spTE*spk*t_r*t_t + 10*c_r**2*c_t*spk*t_r**2*t_t + 10*c_r**2*c_t*spk*t_r**2 + 20*c_r**2*c_t*spk*t_r*t_t + 80*c_r**2*sk**2*spLE*t_r + 80*c_r**2*sk**2*spTE*t_r - 80*c_r**2*sk**2*t_r**2 - 80*c_r**2*sk**2*t_r + 160*c_r**2*sk*spLE*spk*t_r + 160*c_r**2*sk*spTE*spk*t_r - 160*c_r**2*sk*spk*t_r**2 - 160*c_r**2*sk*spk*t_r + 80*c_r**2*spLE*spk**2*t_r - 80*c_r**2*spLE*t_r*y_cg**2 - 80*c_r**2*spLE*t_r*z_cg**2 + 80*c_r**2*spTE*spk**2*t_r - 80*c_r**2*spTE*t_r*y_cg**2 - 80*c_r**2*spTE*t_r*z_cg**2 - 80*c_r**2*spk**2*t_r**2 - 80*c_r**2*spk**2*t_r + 80*c_r**2*t_r*y_cg**2 + 80*c_r**2*t_r*z_cg**2 + 3*c_r*c_t**3*spLE*t_r*t_t**2 + c_r*c_t**3*spLE*t_t**3 + 3*c_r*c_t**3*spTE*t_r*t_t**2 + c_r*c_t**3*spTE*t_t**3 - 3*c_r*c_t**3*t_r*t_t**2 - c_r*c_t**3*t_t**3 - 20*c_r*c_t**2*sk*spLE*t_r*t_t - 10*c_r*c_t**2*sk*spLE*t_t**2 - 20*c_r*c_t**2*sk*spTE*t_r*t_t - 10*c_r*c_t**2*sk*spTE*t_t**2 + 10*c_r*c_t**2*sk*t_r*t_t**2 + 20*c_r*c_t**2*sk*t_r*t_t + 10*c_r*c_t**2*sk*t_t**2 - 20*c_r*c_t**2*spLE*spk*t_r*t_t - 10*c_r*c_t**2*spLE*spk*t_t**2 - 20*c_r*c_t**2*spTE*spk*t_r*t_t - 10*c_r*c_t**2*spTE*spk*t_t**2 + 10*c_r*c_t**2*spk*t_r*t_t**2 + 20*c_r*c_t**2*spk*t_r*t_t + 10*c_r*c_t**2*spk*t_t**2 + 40*c_r*c_t*sk**2*spLE*t_r + 40*c_r*c_t*sk**2*spLE*t_t + 40*c_r*c_t*sk**2*spTE*t_r + 40*c_r*c_t*sk**2*spTE*t_t - 80*c_r*c_t*sk**2*t_r*t_t - 40*c_r*c_t*sk**2*t_r - 40*c_r*c_t*sk**2*t_t + 80*c_r*c_t*sk*spLE*spk*t_r + 80*c_r*c_t*sk*spLE*spk*t_t + 80*c_r*c_t*sk*spTE*spk*t_r + 80*c_r*c_t*sk*spTE*spk*t_t - 160*c_r*c_t*sk*spk*t_r*t_t - 80*c_r*c_t*sk*spk*t_r - 80*c_r*c_t*sk*spk*t_t + 40*c_r*c_t*spLE*spk**2*t_r + 40*c_r*c_t*spLE*spk**2*t_t - 40*c_r*c_t*spLE*t_r*y_cg**2 - 40*c_r*c_t*spLE*t_r*z_cg**2 - 40*c_r*c_t*spLE*t_t*y_cg**2 - 40*c_r*c_t*spLE*t_t*z_cg**2 + 40*c_r*c_t*spTE*spk**2*t_r + 40*c_r*c_t*spTE*spk**2*t_t - 40*c_r*c_t*spTE*t_r*y_cg**2 - 40*c_r*c_t*spTE*t_r*z_cg**2 - 40*c_r*c_t*spTE*t_t*y_cg**2 - 40*c_r*c_t*spTE*t_t*z_cg**2 - 80*c_r*c_t*spk**2*t_r*t_t - 40*c_r*c_t*spk**2*t_r - 40*c_r*c_t*spk**2*t_t + 40*c_r*c_t*t_r*y_cg**2 + 40*c_r*c_t*t_r*z_cg**2 + 40*c_r*c_t*t_t*y_cg**2 + 40*c_r*c_t*t_t*z_cg**2 - 80*c_r*sk**3*spLE - 80*c_r*sk**3*spTE + 240*c_r*sk**3*t_r + 80*c_r*sk**3 - 240*c_r*sk**2*spLE*spk - 240*c_r*sk**2*spTE*spk + 720*c_r*sk**2*spk*t_r + 240*c_r*sk**2*spk - 240*c_r*sk*spLE*spk**2 + 240*c_r*sk*spLE*y_cg**2 + 240*c_r*sk*spLE*z_cg**2 - 240*c_r*sk*spTE*spk**2 + 240*c_r*sk*spTE*y_cg**2 + 240*c_r*sk*spTE*z_cg**2 + 720*c_r*sk*spk**2*t_r + 240*c_r*sk*spk**2 - 240*c_r*sk*t_r*y_cg**2 - 240*c_r*sk*t_r*z_cg**2 - 240*c_r*sk*y_cg**2 - 240*c_r*sk*z_cg**2 - 80*c_r*spLE*spk**3 + 240*c_r*spLE*spk*y_cg**2 + 240*c_r*spLE*spk*z_cg**2 - 80*c_r*spTE*spk**3 + 240*c_r*spTE*spk*y_cg**2 + 240*c_r*spTE*spk*z_cg**2 + 240*c_r*spk**3*t_r + 80*c_r*spk**3 - 240*c_r*spk*t_r*y_cg**2 - 240*c_r*spk*t_r*z_cg**2 - 240*c_r*spk*y_cg**2 - 240*c_r*spk*z_cg**2 + 4*c_t**4*spLE*t_t**3 + 4*c_t**4*spTE*t_t**3 - 4*c_t**4*t_t**3 - 30*c_t**3*sk*spLE*t_t**2 - 30*c_t**3*sk*spTE*t_t**2 + 10*c_t**3*sk*t_t**3 + 30*c_t**3*sk*t_t**2 - 30*c_t**3*spLE*spk*t_t**2 - 30*c_t**3*spTE*spk*t_t**2 + 10*c_t**3*spk*t_t**3 + 30*c_t**3*spk*t_t**2 + 80*c_t**2*sk**2*spLE*t_t + 80*c_t**2*sk**2*spTE*t_t - 80*c_t**2*sk**2*t_t**2 - 80*c_t**2*sk**2*t_t + 160*c_t**2*sk*spLE*spk*t_t + 160*c_t**2*sk*spTE*spk*t_t - 160*c_t**2*sk*spk*t_t**2 - 160*c_t**2*sk*spk*t_t + 80*c_t**2*spLE*spk**2*t_t - 80*c_t**2*spLE*t_t*y_cg**2 - 80*c_t**2*spLE*t_t*z_cg**2 + 80*c_t**2*spTE*spk**2*t_t - 80*c_t**2*spTE*t_t*y_cg**2 - 80*c_t**2*spTE*t_t*z_cg**2 - 80*c_t**2*spk**2*t_t**2 - 80*c_t**2*spk**2*t_t + 80*c_t**2*t_t*y_cg**2 + 80*c_t**2*t_t*z_cg**2 - 80*c_t*sk**3*spLE - 80*c_t*sk**3*spTE + 240*c_t*sk**3*t_t + 80*c_t*sk**3 - 240*c_t*sk**2*spLE*spk - 240*c_t*sk**2*spTE*spk + 720*c_t*sk**2*spk*t_t + 240*c_t*sk**2*spk - 240*c_t*sk*spLE*spk**2 + 240*c_t*sk*spLE*y_cg**2 + 240*c_t*sk*spLE*z_cg**2 - 240*c_t*sk*spTE*spk**2 + 240*c_t*sk*spTE*y_cg**2 + 240*c_t*sk*spTE*z_cg**2 + 720*c_t*sk*spk**2*t_t + 240*c_t*sk*spk**2 - 240*c_t*sk*t_t*y_cg**2 - 240*c_t*sk*t_t*z_cg**2 - 240*c_t*sk*y_cg**2 - 240*c_t*sk*z_cg**2 - 80*c_t*spLE*spk**3 + 240*c_t*spLE*spk*y_cg**2 + 240*c_t*spLE*spk*z_cg**2 - 80*c_t*spTE*spk**3 + 240*c_t*spTE*spk*y_cg**2 + 240*c_t*spTE*spk*z_cg**2 + 240*c_t*spk**3*t_t + 80*c_t*spk**3 - 240*c_t*spk*t_t*y_cg**2 - 240*c_t*spk*t_t*z_cg**2 - 240*c_t*spk*y_cg**2 - 240*c_t*spk*z_cg**2 - 320*sk**4 - 1280*sk**3*spk - 1920*sk**2*spk**2 + 960*sk**2*y_cg**2 + 960*sk**2*z_cg**2 - 1280*sk*spk**3 + 1920*sk*spk*y_cg**2 + 1920*sk*spk*z_cg**2 - 320*spk**4 + 960*spk**2*y_cg**2 + 960*spk**2*z_cg**2)/240 

# Iyy_new = -b/2*p*(32*b/2**2*c_r**2*spLE*t_r*tan(Λ_r)**2 + 32*b/2**2*c_r**2*spTE*t_r*tan(Λ_r)**2 - 32*b/2**2*c_r**2*t_r*tan(Λ_r)**2 + 48*b/2**2*c_r*c_t*spLE*t_r*tan(Λ_r)**2 + 48*b/2**2*c_r*c_t*spLE*t_t*tan(Λ_r)**2 + 48*b/2**2*c_r*c_t*spTE*t_r*tan(Λ_r)**2 + 48*b/2**2*c_r*c_t*spTE*t_t*tan(Λ_r)**2 - 48*b/2**2*c_r*c_t*t_r*tan(Λ_r)**2 - 48*b/2**2*c_r*c_t*t_t*tan(Λ_r)**2 - 160*b/2**2*c_r*sk*spLE*tan(Λ_r)**2 - 160*b/2**2*c_r*sk*spTE*tan(Λ_r)**2 + 160*b/2**2*c_r*sk*t_r*tan(Λ_r)**2 + 160*b/2**2*c_r*sk*tan(Λ_r)**2 - 160*b/2**2*c_r*spLE*spk*tan(Λ_r)**2 - 160*b/2**2*c_r*spTE*spk*tan(Λ_r)**2 + 160*b/2**2*c_r*spk*t_r*tan(Λ_r)**2 + 160*b/2**2*c_r*spk*tan(Λ_r)**2 + 192*b/2**2*c_t**2*spLE*t_t*tan(Λ_r)**2 + 192*b/2**2*c_t**2*spTE*t_t*tan(Λ_r)**2 - 192*b/2**2*c_t**2*t_t*tan(Λ_r)**2 - 480*b/2**2*c_t*sk*spLE*tan(Λ_r)**2 - 480*b/2**2*c_t*sk*spTE*tan(Λ_r)**2 + 480*b/2**2*c_t*sk*t_t*tan(Λ_r)**2 + 480*b/2**2*c_t*sk*tan(Λ_r)**2 - 480*b/2**2*c_t*spLE*spk*tan(Λ_r)**2 - 480*b/2**2*c_t*spTE*spk*tan(Λ_r)**2 + 480*b/2**2*c_t*spk*t_t*tan(Λ_r)**2 + 480*b/2**2*c_t*spk*tan(Λ_r)**2 - 1280*b/2**2*sk**2*tan(Λ_r)**2 - 2560*b/2**2*sk*spk*tan(Λ_r)**2 - 1280*b/2**2*spk**2*tan(Λ_r)**2 + 48*b/2*c_r**3*spLE**2*t_r*tan(Λ_r) - 24*b/2*c_r**3*spLE*t_r*tan(Λ_r) - 48*b/2*c_r**3*spTE**2*t_r*tan(Λ_r) + 72*b/2*c_r**3*spTE*t_r*tan(Λ_r) - 24*b/2*c_r**3*t_r*tan(Λ_r) + 64*b/2*c_r**2*c_t*spLE**2*t_r*tan(Λ_r) + 32*b/2*c_r**2*c_t*spLE**2*t_t*tan(Λ_r) - 32*b/2*c_r**2*c_t*spLE*t_r*tan(Λ_r) - 16*b/2*c_r**2*c_t*spLE*t_t*tan(Λ_r) - 64*b/2*c_r**2*c_t*spTE**2*t_r*tan(Λ_r) - 32*b/2*c_r**2*c_t*spTE**2*t_t*tan(Λ_r) + 96*b/2*c_r**2*c_t*spTE*t_r*tan(Λ_r) + 48*b/2*c_r**2*c_t*spTE*t_t*tan(Λ_r) - 32*b/2*c_r**2*c_t*t_r*tan(Λ_r) - 16*b/2*c_r**2*c_t*t_t*tan(Λ_r) - 160*b/2*c_r**2*sk*spLE**2*tan(Λ_r) + 160*b/2*c_r**2*sk*spLE*t_r*tan(Λ_r) + 80*b/2*c_r**2*sk*spLE*tan(Λ_r) + 160*b/2*c_r**2*sk*spTE**2*tan(Λ_r) - 160*b/2*c_r**2*sk*spTE*t_r*tan(Λ_r) - 240*b/2*c_r**2*sk*spTE*tan(Λ_r) + 80*b/2*c_r**2*sk*t_r*tan(Λ_r) + 80*b/2*c_r**2*sk*tan(Λ_r) - 160*b/2*c_r**2*spLE**2*spk*tan(Λ_r) + 160*b/2*c_r**2*spLE*spk*t_r*tan(Λ_r) + 80*b/2*c_r**2*spLE*spk*tan(Λ_r) + 160*b/2*c_r**2*spTE**2*spk*tan(Λ_r) - 160*b/2*c_r**2*spTE*spk*t_r*tan(Λ_r) - 240*b/2*c_r**2*spTE*spk*tan(Λ_r) + 80*b/2*c_r**2*spk*t_r*tan(Λ_r) + 80*b/2*c_r**2*spk*tan(Λ_r) + 48*b/2*c_r*c_t**2*spLE**2*t_r*tan(Λ_r) + 96*b/2*c_r*c_t**2*spLE**2*t_t*tan(Λ_r) - 24*b/2*c_r*c_t**2*spLE*t_r*tan(Λ_r) - 48*b/2*c_r*c_t**2*spLE*t_t*tan(Λ_r) - 48*b/2*c_r*c_t**2*spTE**2*t_r*tan(Λ_r) - 96*b/2*c_r*c_t**2*spTE**2*t_t*tan(Λ_r) + 72*b/2*c_r*c_t**2*spTE*t_r*tan(Λ_r) + 144*b/2*c_r*c_t**2*spTE*t_t*tan(Λ_r) - 24*b/2*c_r*c_t**2*t_r*tan(Λ_r) - 48*b/2*c_r*c_t**2*t_t*tan(Λ_r) - 320*b/2*c_r*c_t*sk*spLE**2*tan(Λ_r) + 160*b/2*c_r*c_t*sk*spLE*t_r*tan(Λ_r) + 160*b/2*c_r*c_t*sk*spLE*t_t*tan(Λ_r) + 160*b/2*c_r*c_t*sk*spLE*tan(Λ_r) + 320*b/2*c_r*c_t*sk*spTE**2*tan(Λ_r) - 160*b/2*c_r*c_t*sk*spTE*t_r*tan(Λ_r) - 160*b/2*c_r*c_t*sk*spTE*t_t*tan(Λ_r) - 480*b/2*c_r*c_t*sk*spTE*tan(Λ_r) + 80*b/2*c_r*c_t*sk*t_r*tan(Λ_r) + 80*b/2*c_r*c_t*sk*t_t*tan(Λ_r) + 160*b/2*c_r*c_t*sk*tan(Λ_r) - 320*b/2*c_r*c_t*spLE**2*spk*tan(Λ_r) + 160*b/2*c_r*c_t*spLE*spk*t_r*tan(Λ_r) + 160*b/2*c_r*c_t*spLE*spk*t_t*tan(Λ_r) + 160*b/2*c_r*c_t*spLE*spk*tan(Λ_r) + 320*b/2*c_r*c_t*spTE**2*spk*tan(Λ_r) - 160*b/2*c_r*c_t*spTE*spk*t_r*tan(Λ_r) - 160*b/2*c_r*c_t*spTE*spk*t_t*tan(Λ_r) - 480*b/2*c_r*c_t*spTE*spk*tan(Λ_r) + 80*b/2*c_r*c_t*spk*t_r*tan(Λ_r) + 80*b/2*c_r*c_t*spk*t_t*tan(Λ_r) + 160*b/2*c_r*c_t*spk*tan(Λ_r) - 640*b/2*c_r*sk**2*spLE*tan(Λ_r) + 640*b/2*c_r*sk**2*spTE*tan(Λ_r) - 320*b/2*c_r*sk**2*tan(Λ_r) - 1280*b/2*c_r*sk*spLE*spk*tan(Λ_r) + 1280*b/2*c_r*sk*spTE*spk*tan(Λ_r) - 640*b/2*c_r*sk*spk*tan(Λ_r) - 640*b/2*c_r*spLE*spk**2*tan(Λ_r) + 640*b/2*c_r*spTE*spk**2*tan(Λ_r) - 320*b/2*c_r*spk**2*tan(Λ_r) + 192*b/2*c_t**3*spLE**2*t_t*tan(Λ_r) - 96*b/2*c_t**3*spLE*t_t*tan(Λ_r) - 192*b/2*c_t**3*spTE**2*t_t*tan(Λ_r) + 288*b/2*c_t**3*spTE*t_t*tan(Λ_r) - 96*b/2*c_t**3*t_t*tan(Λ_r) - 480*b/2*c_t**2*sk*spLE**2*tan(Λ_r) + 480*b/2*c_t**2*sk*spLE*t_t*tan(Λ_r) + 240*b/2*c_t**2*sk*spLE*tan(Λ_r) + 480*b/2*c_t**2*sk*spTE**2*tan(Λ_r) - 480*b/2*c_t**2*sk*spTE*t_t*tan(Λ_r) - 720*b/2*c_t**2*sk*spTE*tan(Λ_r) + 240*b/2*c_t**2*sk*t_t*tan(Λ_r) + 240*b/2*c_t**2*sk*tan(Λ_r) - 480*b/2*c_t**2*spLE**2*spk*tan(Λ_r) + 480*b/2*c_t**2*spLE*spk*t_t*tan(Λ_r) + 240*b/2*c_t**2*spLE*spk*tan(Λ_r) + 480*b/2*c_t**2*spTE**2*spk*tan(Λ_r) - 480*b/2*c_t**2*spTE*spk*t_t*tan(Λ_r) - 720*b/2*c_t**2*spTE*spk*tan(Λ_r) + 240*b/2*c_t**2*spk*t_t*tan(Λ_r) + 240*b/2*c_t**2*spk*tan(Λ_r) - 1280*b/2*c_t*sk**2*spLE*tan(Λ_r) + 1280*b/2*c_t*sk**2*spTE*tan(Λ_r) - 640*b/2*c_t*sk**2*tan(Λ_r) - 2560*b/2*c_t*sk*spLE*spk*tan(Λ_r) + 2560*b/2*c_t*sk*spTE*spk*tan(Λ_r) - 1280*b/2*c_t*sk*spk*tan(Λ_r) - 1280*b/2*c_t*spLE*spk**2*tan(Λ_r) + 1280*b/2*c_t*spTE*spk**2*tan(Λ_r) - 640*b/2*c_t*spk**2*tan(Λ_r) + 64*c_r**4*spLE**3*t_r - 48*c_r**4*spLE**2*t_r + 16*c_r**4*spLE*t_r**3 + 12*c_r**4*spLE*t_r + 64*c_r**4*spTE**3*t_r - 144*c_r**4*spTE**2*t_r + 16*c_r**4*spTE*t_r**3 + 108*c_r**4*spTE*t_r - 16*c_r**4*t_r**3 - 28*c_r**4*t_r + 48*c_r**3*c_t*spLE**3*t_r + 16*c_r**3*c_t*spLE**3*t_t - 36*c_r**3*c_t*spLE**2*t_r - 12*c_r**3*c_t*spLE**2*t_t + 4*c_r**3*c_t*spLE*t_r**3 + 12*c_r**3*c_t*spLE*t_r**2*t_t + 9*c_r**3*c_t*spLE*t_r + 3*c_r**3*c_t*spLE*t_t + 48*c_r**3*c_t*spTE**3*t_r + 16*c_r**3*c_t*spTE**3*t_t - 108*c_r**3*c_t*spTE**2*t_r - 36*c_r**3*c_t*spTE**2*t_t + 4*c_r**3*c_t*spTE*t_r**3 + 12*c_r**3*c_t*spTE*t_r**2*t_t + 81*c_r**3*c_t*spTE*t_r + 27*c_r**3*c_t*spTE*t_t - 4*c_r**3*c_t*t_r**3 - 12*c_r**3*c_t*t_r**2*t_t - 21*c_r**3*c_t*t_r - 7*c_r**3*c_t*t_t - 160*c_r**3*sk*spLE**3 + 240*c_r**3*sk*spLE**2*t_r + 120*c_r**3*sk*spLE**2 - 120*c_r**3*sk*spLE*t_r**2 - 120*c_r**3*sk*spLE*t_r - 30*c_r**3*sk*spLE - 160*c_r**3*sk*spTE**3 + 240*c_r**3*sk*spTE**2*t_r + 360*c_r**3*sk*spTE**2 - 120*c_r**3*sk*spTE*t_r**2 - 360*c_r**3*sk*spTE*t_r - 270*c_r**3*sk*spTE + 40*c_r**3*sk*t_r**3 + 120*c_r**3*sk*t_r**2 + 150*c_r**3*sk*t_r + 70*c_r**3*sk - 160*c_r**3*spLE**3*spk + 240*c_r**3*spLE**2*spk*t_r + 120*c_r**3*spLE**2*spk - 120*c_r**3*spLE*spk*t_r**2 - 120*c_r**3*spLE*spk*t_r - 30*c_r**3*spLE*spk - 160*c_r**3*spTE**3*spk + 240*c_r**3*spTE**2*spk*t_r + 360*c_r**3*spTE**2*spk - 120*c_r**3*spTE*spk*t_r**2 - 360*c_r**3*spTE*spk*t_r - 270*c_r**3*spTE*spk + 40*c_r**3*spk*t_r**3 + 120*c_r**3*spk*t_r**2 + 150*c_r**3*spk*t_r + 70*c_r**3*spk + 32*c_r**2*c_t**2*spLE**3*t_r + 32*c_r**2*c_t**2*spLE**3*t_t - 24*c_r**2*c_t**2*spLE**2*t_r - 24*c_r**2*c_t**2*spLE**2*t_t + 8*c_r**2*c_t**2*spLE*t_r**2*t_t + 8*c_r**2*c_t**2*spLE*t_r*t_t**2 + 6*c_r**2*c_t**2*spLE*t_r + 6*c_r**2*c_t**2*spLE*t_t + 32*c_r**2*c_t**2*spTE**3*t_r + 32*c_r**2*c_t**2*spTE**3*t_t - 72*c_r**2*c_t**2*spTE**2*t_r - 72*c_r**2*c_t**2*spTE**2*t_t + 8*c_r**2*c_t**2*spTE*t_r**2*t_t + 8*c_r**2*c_t**2*spTE*t_r*t_t**2 + 54*c_r**2*c_t**2*spTE*t_r + 54*c_r**2*c_t**2*spTE*t_t - 8*c_r**2*c_t**2*t_r**2*t_t - 8*c_r**2*c_t**2*t_r*t_t**2 - 14*c_r**2*c_t**2*t_r - 14*c_r**2*c_t**2*t_t - 160*c_r**2*c_t*sk*spLE**3 + 160*c_r**2*c_t*sk*spLE**2*t_r + 80*c_r**2*c_t*sk*spLE**2*t_t + 120*c_r**2*c_t*sk*spLE**2 - 40*c_r**2*c_t*sk*spLE*t_r**2 - 80*c_r**2*c_t*sk*spLE*t_r*t_t - 80*c_r**2*c_t*sk*spLE*t_r - 40*c_r**2*c_t*sk*spLE*t_t - 30*c_r**2*c_t*sk*spLE - 160*c_r**2*c_t*sk*spTE**3 + 160*c_r**2*c_t*sk*spTE**2*t_r + 80*c_r**2*c_t*sk*spTE**2*t_t + 360*c_r**2*c_t*sk*spTE**2 - 40*c_r**2*c_t*sk*spTE*t_r**2 - 80*c_r**2*c_t*sk*spTE*t_r*t_t - 240*c_r**2*c_t*sk*spTE*t_r - 120*c_r**2*c_t*sk*spTE*t_t - 270*c_r**2*c_t*sk*spTE + 40*c_r**2*c_t*sk*t_r**2*t_t + 40*c_r**2*c_t*sk*t_r**2 + 80*c_r**2*c_t*sk*t_r*t_t + 100*c_r**2*c_t*sk*t_r + 50*c_r**2*c_t*sk*t_t + 70*c_r**2*c_t*sk - 160*c_r**2*c_t*spLE**3*spk + 160*c_r**2*c_t*spLE**2*spk*t_r + 80*c_r**2*c_t*spLE**2*spk*t_t + 120*c_r**2*c_t*spLE**2*spk - 40*c_r**2*c_t*spLE*spk*t_r**2 - 80*c_r**2*c_t*spLE*spk*t_r*t_t - 80*c_r**2*c_t*spLE*spk*t_r - 40*c_r**2*c_t*spLE*spk*t_t - 30*c_r**2*c_t*spLE*spk - 160*c_r**2*c_t*spTE**3*spk + 160*c_r**2*c_t*spTE**2*spk*t_r + 80*c_r**2*c_t*spTE**2*spk*t_t + 360*c_r**2*c_t*spTE**2*spk - 40*c_r**2*c_t*spTE*spk*t_r**2 - 80*c_r**2*c_t*spTE*spk*t_r*t_t - 240*c_r**2*c_t*spTE*spk*t_r - 120*c_r**2*c_t*spTE*spk*t_t - 270*c_r**2*c_t*spTE*spk + 40*c_r**2*c_t*spk*t_r**2*t_t + 40*c_r**2*c_t*spk*t_r**2 + 80*c_r**2*c_t*spk*t_r*t_t + 100*c_r**2*c_t*spk*t_r + 50*c_r**2*c_t*spk*t_t + 70*c_r**2*c_t*spk - 640*c_r**2*sk**2*spLE**2 + 640*c_r**2*sk**2*spLE*t_r + 320*c_r**2*sk**2*spLE - 640*c_r**2*sk**2*spTE**2 + 640*c_r**2*sk**2*spTE*t_r + 960*c_r**2*sk**2*spTE - 320*c_r**2*sk**2*t_r**2 - 640*c_r**2*sk**2*t_r - 400*c_r**2*sk**2 - 1280*c_r**2*sk*spLE**2*spk + 1280*c_r**2*sk*spLE*spk*t_r + 640*c_r**2*sk*spLE*spk - 1280*c_r**2*sk*spTE**2*spk + 1280*c_r**2*sk*spTE*spk*t_r + 1920*c_r**2*sk*spTE*spk - 640*c_r**2*sk*spk*t_r**2 - 1280*c_r**2*sk*spk*t_r - 800*c_r**2*sk*spk - 640*c_r**2*spLE**2*spk**2 + 640*c_r**2*spLE*spk**2*t_r + 320*c_r**2*spLE*spk**2 - 320*c_r**2*spLE*t_r*x_cg**2 - 320*c_r**2*spLE*t_r*z_cg**2 - 640*c_r**2*spTE**2*spk**2 + 640*c_r**2*spTE*spk**2*t_r + 960*c_r**2*spTE*spk**2 - 320*c_r**2*spTE*t_r*x_cg**2 - 320*c_r**2*spTE*t_r*z_cg**2 - 320*c_r**2*spk**2*t_r**2 - 640*c_r**2*spk**2*t_r - 400*c_r**2*spk**2 + 320*c_r**2*t_r*x_cg**2 + 320*c_r**2*t_r*z_cg**2 + 16*c_r*c_t**3*spLE**3*t_r + 48*c_r*c_t**3*spLE**3*t_t - 12*c_r*c_t**3*spLE**2*t_r - 36*c_r*c_t**3*spLE**2*t_t + 12*c_r*c_t**3*spLE*t_r*t_t**2 + 3*c_r*c_t**3*spLE*t_r + 4*c_r*c_t**3*spLE*t_t**3 + 9*c_r*c_t**3*spLE*t_t + 16*c_r*c_t**3*spTE**3*t_r + 48*c_r*c_t**3*spTE**3*t_t - 36*c_r*c_t**3*spTE**2*t_r - 108*c_r*c_t**3*spTE**2*t_t + 12*c_r*c_t**3*spTE*t_r*t_t**2 + 27*c_r*c_t**3*spTE*t_r + 4*c_r*c_t**3*spTE*t_t**3 + 81*c_r*c_t**3*spTE*t_t - 12*c_r*c_t**3*t_r*t_t**2 - 7*c_r*c_t**3*t_r - 4*c_r*c_t**3*t_t**3 - 21*c_r*c_t**3*t_t - 160*c_r*c_t**2*sk*spLE**3 + 80*c_r*c_t**2*sk*spLE**2*t_r + 160*c_r*c_t**2*sk*spLE**2*t_t + 120*c_r*c_t**2*sk*spLE**2 - 80*c_r*c_t**2*sk*spLE*t_r*t_t - 40*c_r*c_t**2*sk*spLE*t_r - 40*c_r*c_t**2*sk*spLE*t_t**2 - 80*c_r*c_t**2*sk*spLE*t_t - 30*c_r*c_t**2*sk*spLE 
#   - 160*c_r*c_t**2*sk*spTE**3 + 80*c_r*c_t**2*sk*spTE**2*t_r + 160*c_r*c_t**2*sk*spTE**2*t_t + 360*c_r*c_t**2*sk*spTE**2 - 80*c_r*c_t**2*sk*spTE*t_r*t_t - 120*c_r*c_t**2*sk*spTE*t_r - 40*c_r*c_t**2*sk*spTE*t_t**2 - 240*c_r*c_t**2*sk*spTE*t_t - 270*c_r*c_t**2*sk*spTE + 40*c_r*c_t**2*sk*t_r*t_t**2 + 80*c_r*c_t**2*sk*t_r*t_t + 50*c_r*c_t**2*sk*t_r + 40*c_r*c_t**2*sk*t_t**2 + 100*c_r*c_t**2*sk*t_t + 70*c_r*c_t**2*sk - 160*c_r*c_t**2*spLE**3*spk + 80*c_r*c_t**2*spLE**2*spk*t_r + 160*c_r*c_t**2*spLE**2*spk*t_t + 120*c_r*c_t**2*spLE**2*spk - 80*c_r*c_t**2*spLE*spk*t_r*t_t - 40*c_r*c_t**2*spLE*spk*t_r - 40*c_r*c_t**2*spLE*spk*t_t**2 - 80*c_r*c_t**2*spLE*spk*t_t - 30*c_r*c_t**2*spLE*spk - 160*c_r*c_t**2*spTE**3*spk + 80*c_r*c_t**2*spTE**2*spk*t_r + 160*c_r*c_t**2*spTE**2*spk*t_t + 360*c_r*c_t**2*spTE**2*spk - 80*c_r*c_t**2*spTE*spk*t_r*t_t - 120*c_r*c_t**2*spTE*spk*t_r - 40*c_r*c_t**2*spTE*spk*t_t**2 - 240*c_r*c_t**2*spTE*spk*t_t - 270*c_r*c_t**2*spTE*spk + 40*c_r*c_t**2*spk*t_r*t_t**2 + 80*c_r*c_t**2*spk*t_r*t_t + 50*c_r*c_t**2*spk*t_r + 40*c_r*c_t**2*spk*t_t**2 + 100*c_r*c_t**2*spk*t_t + 70*c_r*c_t**2*spk - 640*c_r*c_t*sk**2*spLE**2 + 320*c_r*c_t*sk**2*spLE*t_r + 320*c_r*c_t*sk**2*spLE*t_t + 320*c_r*c_t*sk**2*spLE - 640*c_r*c_t*sk**2*spTE**2 + 320*c_r*c_t*sk**2*spTE*t_r + 320*c_r*c_t*sk**2*spTE*t_t + 960*c_r*c_t*sk**2*spTE - 320*c_r*c_t*sk**2*t_r*t_t - 320*c_r*c_t*sk**2*t_r - 320*c_r*c_t*sk**2*t_t - 400*c_r*c_t*sk**2 - 1280*c_r*c_t*sk*spLE**2*spk + 640*c_r*c_t*sk*spLE*spk*t_r + 640*c_r*c_t*sk*spLE*spk*t_t + 640*c_r*c_t*sk*spLE*spk - 1280*c_r*c_t*sk*spTE**2*spk + 640*c_r*c_t*sk*spTE*spk*t_r + 640*c_r*c_t*sk*spTE*spk*t_t + 1920*c_r*c_t*sk*spTE*spk - 640*c_r*c_t*sk*spk*t_r*t_t - 640*c_r*c_t*sk*spk*t_r - 640*c_r*c_t*sk*spk*t_t - 800*c_r*c_t*sk*spk - 640*c_r*c_t*spLE**2*spk**2 + 320*c_r*c_t*spLE*spk**2*t_r + 320*c_r*c_t*spLE*spk**2*t_t + 320*c_r*c_t*spLE*spk**2 - 160*c_r*c_t*spLE*t_r*x_cg**2 - 160*c_r*c_t*spLE*t_r*z_cg**2 - 160*c_r*c_t*spLE*t_t*x_cg**2 - 160*c_r*c_t*spLE*t_t*z_cg**2 - 640*c_r*c_t*spTE**2*spk**2 + 320*c_r*c_t*spTE*spk**2*t_r + 320*c_r*c_t*spTE*spk**2*t_t + 960*c_r*c_t*spTE*spk**2 - 160*c_r*c_t*spTE*t_r*x_cg**2 - 160*c_r*c_t*spTE*t_r*z_cg**2 - 160*c_r*c_t*spTE*t_t*x_cg**2 - 160*c_r*c_t*spTE*t_t*z_cg**2 - 320*c_r*c_t*spk**2*t_r*t_t - 320*c_r*c_t*spk**2*t_r - 320*c_r*c_t*spk**2*t_t - 400*c_r*c_t*spk**2 + 160*c_r*c_t*t_r*x_cg**2 + 160*c_r*c_t*t_r*z_cg**2 + 160*c_r*c_t*t_t*x_cg**2 + 160*c_r*c_t*t_t*z_cg**2 - 1280*c_r*sk**3*spLE - 1280*c_r*sk**3*spTE + 1280*c_r*sk**3*t_r + 1280*c_r*sk**3 - 3840*c_r*sk**2*spLE*spk - 3840*c_r*sk**2*spTE*spk + 3840*c_r*sk**2*spk*t_r + 3840*c_r*sk**2*spk - 3840*c_r*sk*spLE*spk**2 + 960*c_r*sk*spLE*x_cg**2 + 960*c_r*sk*spLE*z_cg**2 - 3840*c_r*sk*spTE*spk**2 + 960*c_r*sk*spTE*x_cg**2 + 960*c_r*sk*spTE*z_cg**2 + 3840*c_r*sk*spk**2*t_r + 3840*c_r*sk*spk**2 - 960*c_r*sk*t_r*x_cg**2 - 960*c_r*sk*t_r*z_cg**2 - 960*c_r*sk*x_cg**2 - 960*c_r*sk*z_cg**2 - 1280*c_r*spLE*spk**3 + 960*c_r*spLE*spk*x_cg**2 + 960*c_r*spLE*spk*z_cg**2 - 1280*c_r*spTE*spk**3 + 960*c_r*spTE*spk*x_cg**2 + 960*c_r*spTE*spk*z_cg**2 + 1280*c_r*spk**3*t_r + 1280*c_r*spk**3 - 960*c_r*spk*t_r*x_cg**2 - 960*c_r*spk*t_r*z_cg**2 - 960*c_r*spk*x_cg**2 - 960*c_r*spk*z_cg**2 + 64*c_t**4*spLE**3*t_t - 48*c_t**4*spLE**2*t_t + 16*c_t**4*spLE*t_t**3 + 12*c_t**4*spLE*t_t + 64*c_t**4*spTE**3*t_t - 144*c_t**4*spTE**2*t_t + 16*c_t**4*spTE*t_t**3 + 108*c_t**4*spTE*t_t - 16*c_t**4*t_t**3 - 28*c_t**4*t_t - 160*c_t**3*sk*spLE**3 + 240*c_t**3*sk*spLE**2*t_t + 120*c_t**3*sk*spLE**2 - 120*c_t**3*sk*spLE*t_t**2 - 120*c_t**3*sk*spLE*t_t - 30*c_t**3*sk*spLE - 160*c_t**3*sk*spTE**3 + 240*c_t**3*sk*spTE**2*t_t + 360*c_t**3*sk*spTE**2 - 120*c_t**3*sk*spTE*t_t**2 - 360*c_t**3*sk*spTE*t_t - 270*c_t**3*sk*spTE + 40*c_t**3*sk*t_t**3 + 120*c_t**3*sk*t_t**2 + 150*c_t**3*sk*t_t + 70*c_t**3*sk - 160*c_t**3*spLE**3*spk + 240*c_t**3*spLE**2*spk*t_t + 120*c_t**3*spLE**2*spk - 120*c_t**3*spLE*spk*t_t**2 - 120*c_t**3*spLE*spk*t_t - 30*c_t**3*spLE*spk - 160*c_t**3*spTE**3*spk + 240*c_t**3*spTE**2*spk*t_t + 360*c_t**3*spTE**2*spk - 120*c_t**3*spTE*spk*t_t**2 - 360*c_t**3*spTE*spk*t_t - 270*c_t**3*spTE*spk + 40*c_t**3*spk*t_t**3 + 120*c_t**3*spk*t_t**2 + 150*c_t**3*spk*t_t + 70*c_t**3*spk - 640*c_t**2*sk**2*spLE**2 + 640*c_t**2*sk**2*spLE*t_t + 320*c_t**2*sk**2*spLE - 640*c_t**2*sk**2*spTE**2 + 640*c_t**2*sk**2*spTE*t_t + 960*c_t**2*sk**2*spTE - 320*c_t**2*sk**2*t_t**2 - 640*c_t**2*sk**2*t_t - 400*c_t**2*sk**2 - 1280*c_t**2*sk*spLE**2*spk + 1280*c_t**2*sk*spLE*spk*t_t + 640*c_t**2*sk*spLE*spk - 1280*c_t**2*sk*spTE**2*spk + 1280*c_t**2*sk*spTE*spk*t_t + 1920*c_t**2*sk*spTE*spk - 640*c_t**2*sk*spk*t_t**2 - 1280*c_t**2*sk*spk*t_t - 800*c_t**2*sk*spk - 640*c_t**2*spLE**2*spk**2 + 640*c_t**2*spLE*spk**2*t_t + 320*c_t**2*spLE*spk**2 - 320*c_t**2*spLE*t_t*x_cg**2 - 320*c_t**2*spLE*t_t*z_cg**2 - 640*c_t**2*spTE**2*spk**2 + 640*c_t**2*spTE*spk**2*t_t + 960*c_t**2*spTE*spk**2 - 320*c_t**2*spTE*t_t*x_cg**2 - 320*c_t**2*spTE*t_t*z_cg**2 - 320*c_t**2*spk**2*t_t**2 - 640*c_t**2*spk**2*t_t - 400*c_t**2*spk**2 + 320*c_t**2*t_t*x_cg**2 + 320*c_t**2*t_t*z_cg**2 - 1280*c_t*sk**3*spLE - 1280*c_t*sk**3*spTE + 1280*c_t*sk**3*t_t + 1280*c_t*sk**3 - 3840*c_t*sk**2*spLE*spk - 3840*c_t*sk**2*spTE*spk + 3840*c_t*sk**2*spk*t_t + 3840*c_t*sk**2*spk - 3840*c_t*sk*spLE*spk**2 + 960*c_t*sk*spLE*x_cg**2 + 960*c_t*sk*spLE*z_cg**2 - 3840*c_t*sk*spTE*spk**2 + 960*c_t*sk*spTE*x_cg**2 + 960*c_t*sk*spTE*z_cg**2 + 3840*c_t*sk*spk**2*t_t + 3840*c_t*sk*spk**2 - 960*c_t*sk*t_t*x_cg**2 - 960*c_t*sk*t_t*z_cg**2 - 960*c_t*sk*x_cg**2 - 960*c_t*sk*z_cg**2 - 1280*c_t*spLE*spk**3 + 960*c_t*spLE*spk*x_cg**2 + 960*c_t*spLE*spk*z_cg**2 - 1280*c_t*spTE*spk**3 + 960*c_t*spTE*spk*x_cg**2 + 960*c_t*spTE*spk*z_cg**2 + 1280*c_t*spk**3*t_t + 1280*c_t*spk**3 - 960*c_t*spk*t_t*x_cg**2 - 960*c_t*spk*t_t*z_cg**2 - 960*c_t*spk*x_cg**2 - 960*c_t*spk*z_cg**2 - 2560*sk**4 - 10240*sk**3*spk - 15360*sk**2*spk**2 + 3840*sk**2*x_cg**2 + 3840*sk**2*z_cg**2 - 10240*sk*spk**3 + 7680*sk*spk*x_cg**2 + 7680*sk*spk*z_cg**2 - 2560*spk**4 + 3840*spk**2*x_cg**2 + 3840*spk**2*z_cg**2)/960 

# Izz_new = -b/2*p*(32*b/2**2*c_r**2*spLE*t_r*tan(Λ_r)**2 + 32*b/2**2*c_r**2*spLE*t_r + 32*b/2**2*c_r**2*spTE*t_r*tan(Λ_r)**2 + 32*b/2**2*c_r**2*spTE*t_r - 32*b/2**2*c_r**2*t_r*tan(Λ_r)**2 - 32*b/2**2*c_r**2*t_r + 48*b/2**2*c_r*c_t*spLE*t_r*tan(Λ_r)**2 + 48*b/2**2*c_r*c_t*spLE*t_r + 48*b/2**2*c_r*c_t*spLE*t_t*tan(Λ_r)**2 + 48*b/2**2*c_r*c_t*spLE*t_t + 48*b/2**2*c_r*c_t*spTE*t_r*tan(Λ_r)**2 + 48*b/2**2*c_r*c_t*spTE*t_r + 48*b/2**2*c_r*c_t*spTE*t_t*tan(Λ_r)**2 + 48*b/2**2*c_r*c_t*spTE*t_t - 48*b/2**2*c_r*c_t*t_r*tan(Λ_r)**2 - 48*b/2**2*c_r*c_t*t_r - 48*b/2**2*c_r*c_t*t_t*tan(Λ_r)**2 - 48*b/2**2*c_r*c_t*t_t - 160*b/2**2*c_r*sk*spLE*tan(Λ_r)**2 - 160*b/2**2*c_r*sk*spLE - 160*b/2**2*c_r*sk*spTE*tan(Λ_r)**2 - 160*b/2**2*c_r*sk*spTE + 160*b/2**2*c_r*sk*t_r*tan(Λ_r)**2 + 160*b/2**2*c_r*sk*t_r + 160*b/2**2*c_r*sk*tan(Λ_r)**2 + 160*b/2**2*c_r*sk - 160*b/2**2*c_r*spLE*spk*tan(Λ_r)**2 - 160*b/2**2*c_r*spLE*spk - 160*b/2**2*c_r*spTE*spk*tan(Λ_r)**2 - 160*b/2**2*c_r*spTE*spk + 160*b/2**2*c_r*spk*t_r*tan(Λ_r)**2 + 160*b/2**2*c_r*spk*t_r + 160*b/2**2*c_r*spk*tan(Λ_r)**2 + 160*b/2**2*c_r*spk + 192*b/2**2*c_t**2*spLE*t_t*tan(Λ_r)**2 + 192*b/2**2*c_t**2*spLE*t_t + 192*b/2**2*c_t**2*spTE*t_t*tan(Λ_r)**2 + 192*b/2**2*c_t**2*spTE*t_t - 192*b/2**2*c_t**2*t_t*tan(Λ_r)**2 - 192*b/2**2*c_t**2*t_t - 480*b/2**2*c_t*sk*spLE*tan(Λ_r)**2 - 480*b/2**2*c_t*sk*spLE - 480*b/2**2*c_t*sk*spTE*tan(Λ_r)**2 - 480*b/2**2*c_t*sk*spTE + 480*b/2**2*c_t*sk*t_t*tan(Λ_r)**2 + 480*b/2**2*c_t*sk*t_t + 480*b/2**2*c_t*sk*tan(Λ_r)**2 + 480*b/2**2*c_t*sk - 480*b/2**2*c_t*spLE*spk*tan(Λ_r)**2 - 480*b/2**2*c_t*spLE*spk - 480*b/2**2*c_t*spTE*spk*tan(Λ_r)**2 - 480*b/2**2*c_t*spTE*spk + 480*b/2**2*c_t*spk*t_t*tan(Λ_r)**2 + 480*b/2**2*c_t*spk*t_t + 480*b/2**2*c_t*spk*tan(Λ_r)**2 + 480*b/2**2*c_t*spk - 1280*b/2**2*sk**2*tan(Λ_r)**2 - 1280*b/2**2*sk**2 - 2560*b/2**2*sk*spk*tan(Λ_r)**2 - 2560*b/2**2*sk*spk - 1280*b/2**2*spk**2*tan(Λ_r)**2 - 1280*b/2**2*spk**2 + 48*b/2*c_r**3*spLE**2*t_r*tan(Λ_r) - 24*b/2*c_r**3*spLE*t_r*tan(Λ_r) - 48*b/2*c_r**3*spTE**2*t_r*tan(Λ_r) + 72*b/2*c_r**3*spTE*t_r*tan(Λ_r) - 24*b/2*c_r**3*t_r*tan(Λ_r) + 64*b/2*c_r**2*c_t*spLE**2*t_r*tan(Λ_r) + 32*b/2*c_r**2*c_t*spLE**2*t_t*tan(Λ_r) - 32*b/2*c_r**2*c_t*spLE*t_r*tan(Λ_r) - 16*b/2*c_r**2*c_t*spLE*t_t*tan(Λ_r) - 64*b/2*c_r**2*c_t*spTE**2*t_r*tan(Λ_r) - 32*b/2*c_r**2*c_t*spTE**2*t_t*tan(Λ_r) + 96*b/2*c_r**2*c_t*spTE*t_r*tan(Λ_r) + 48*b/2*c_r**2*c_t*spTE*t_t*tan(Λ_r) - 32*b/2*c_r**2*c_t*t_r*tan(Λ_r) - 16*b/2*c_r**2*c_t*t_t*tan(Λ_r) - 160*b/2*c_r**2*sk*spLE**2*tan(Λ_r) + 160*b/2*c_r**2*sk*spLE*t_r*tan(Λ_r) + 80*b/2*c_r**2*sk*spLE*tan(Λ_r) + 160*b/2*c_r**2*sk*spTE**2*tan(Λ_r) - 160*b/2*c_r**2*sk*spTE*t_r*tan(Λ_r) - 240*b/2*c_r**2*sk*spTE*tan(Λ_r) + 80*b/2*c_r**2*sk*t_r*tan(Λ_r) + 80*b/2*c_r**2*sk*tan(Λ_r) - 160*b/2*c_r**2*spLE**2*spk*tan(Λ_r) + 160*b/2*c_r**2*spLE*spk*t_r*tan(Λ_r) + 80*b/2*c_r**2*spLE*spk*tan(Λ_r) + 160*b/2*c_r**2*spTE**2*spk*tan(Λ_r) - 160*b/2*c_r**2*spTE*spk*t_r*tan(Λ_r) - 240*b/2*c_r**2*spTE*spk*tan(Λ_r) + 80*b/2*c_r**2*spk*t_r*tan(Λ_r) + 80*b/2*c_r**2*spk*tan(Λ_r) + 48*b/2*c_r*c_t**2*spLE**2*t_r*tan(Λ_r) + 96*b/2*c_r*c_t**2*spLE**2*t_t*tan(Λ_r) - 24*b/2*c_r*c_t**2*spLE*t_r*tan(Λ_r) - 48*b/2*c_r*c_t**2*spLE*t_t*tan(Λ_r) - 48*b/2*c_r*c_t**2*spTE**2*t_r*tan(Λ_r) - 96*b/2*c_r*c_t**2*spTE**2*t_t*tan(Λ_r) + 72*b/2*c_r*c_t**2*spTE*t_r*tan(Λ_r) + 144*b/2*c_r*c_t**2*spTE*t_t*tan(Λ_r) - 24*b/2*c_r*c_t**2*t_r*tan(Λ_r) - 48*b/2*c_r*c_t**2*t_t*tan(Λ_r) - 320*b/2*c_r*c_t*sk*spLE**2*tan(Λ_r) + 160*b/2*c_r*c_t*sk*spLE*t_r*tan(Λ_r) + 160*b/2*c_r*c_t*sk*spLE*t_t*tan(Λ_r) + 160*b/2*c_r*c_t*sk*spLE*tan(Λ_r) + 320*b/2*c_r*c_t*sk*spTE**2*tan(Λ_r) - 160*b/2*c_r*c_t*sk*spTE*t_r*tan(Λ_r) - 160*b/2*c_r*c_t*sk*spTE*t_t*tan(Λ_r) - 480*b/2*c_r*c_t*sk*spTE*tan(Λ_r) + 80*b/2*c_r*c_t*sk*t_r*tan(Λ_r) + 80*b/2*c_r*c_t*sk*t_t*tan(Λ_r) + 160*b/2*c_r*c_t*sk*tan(Λ_r) - 320*b/2*c_r*c_t*spLE**2*spk*tan(Λ_r) + 160*b/2*c_r*c_t*spLE*spk*t_r*tan(Λ_r) + 160*b/2*c_r*c_t*spLE*spk*t_t*tan(Λ_r) + 160*b/2*c_r*c_t*spLE*spk*tan(Λ_r) + 320*b/2*c_r*c_t*spTE**2*spk*tan(Λ_r) - 160*b/2*c_r*c_t*spTE*spk*t_r*tan(Λ_r) - 160*b/2*c_r*c_t*spTE*spk*t_t*tan(Λ_r) - 480*b/2*c_r*c_t*spTE*spk*tan(Λ_r) + 80*b/2*c_r*c_t*spk*t_r*tan(Λ_r) + 80*b/2*c_r*c_t*spk*t_t*tan(Λ_r) + 160*b/2*c_r*c_t*spk*tan(Λ_r) - 640*b/2*c_r*sk**2*spLE*tan(Λ_r) + 640*b/2*c_r*sk**2*spTE*tan(Λ_r) - 320*b/2*c_r*sk**2*tan(Λ_r) - 1280*b/2*c_r*sk*spLE*spk*tan(Λ_r) + 1280*b/2*c_r*sk*spTE*spk*tan(Λ_r) - 640*b/2*c_r*sk*spk*tan(Λ_r) - 640*b/2*c_r*spLE*spk**2*tan(Λ_r) + 640*b/2*c_r*spTE*spk**2*tan(Λ_r) - 320*b/2*c_r*spk**2*tan(Λ_r) + 192*b/2*c_t**3*spLE**2*t_t*tan(Λ_r) - 96*b/2*c_t**3*spLE*t_t*tan(Λ_r) - 192*b/2*c_t**3*spTE**2*t_t*tan(Λ_r) + 288*b/2*c_t**3*spTE*t_t*tan(Λ_r) - 96*b/2*c_t**3*t_t*tan(Λ_r) - 480*b/2*c_t**2*sk*spLE**2*tan(Λ_r) + 480*b/2*c_t**2*sk*spLE*t_t*tan(Λ_r) + 240*b/2*c_t**2*sk*spLE*tan(Λ_r) + 480*b/2*c_t**2*sk*spTE**2*tan(Λ_r) - 480*b/2*c_t**2*sk*spTE*t_t*tan(Λ_r) - 720*b/2*c_t**2*sk*spTE*tan(Λ_r) + 240*b/2*c_t**2*sk*t_t*tan(Λ_r) + 240*b/2*c_t**2*sk*tan(Λ_r) - 480*b/2*c_t**2*spLE**2*spk*tan(Λ_r) + 480*b/2*c_t**2*spLE*spk*t_t*tan(Λ_r) + 240*b/2*c_t**2*spLE*spk*tan(Λ_r) + 480*b/2*c_t**2*spTE**2*spk*tan(Λ_r) - 480*b/2*c_t**2*spTE*spk*t_t*tan(Λ_r) - 720*b/2*c_t**2*spTE*spk*tan(Λ_r) + 240*b/2*c_t**2*spk*t_t*tan(Λ_r) + 240*b/2*c_t**2*spk*tan(Λ_r) - 1280*b/2*c_t*sk**2*spLE*tan(Λ_r) + 1280*b/2*c_t*sk**2*spTE*tan(Λ_r) - 640*b/2*c_t*sk**2*tan(Λ_r) - 2560*b/2*c_t*sk*spLE*spk*tan(Λ_r) + 2560*b/2*c_t*sk*spTE*spk*tan(Λ_r) - 1280*b/2*c_t*sk*spk*tan(Λ_r) - 1280*b/2*c_t*spLE*spk**2*tan(Λ_r) + 1280*b/2*c_t*spTE*spk**2*tan(Λ_r) - 640*b/2*c_t*spk**2*tan(Λ_r) + 64*c_r**4*spLE**3*t_r - 48*c_r**4*spLE**2*t_r + 12*c_r**4*spLE*t_r + 64*c_r**4*spTE**3*t_r - 144*c_r**4*spTE**2*t_r + 108*c_r**4*spTE*t_r - 28*c_r**4*t_r + 48*c_r**3*c_t*spLE**3*t_r + 16*c_r**3*c_t*spLE**3*t_t - 36*c_r**3*c_t*spLE**2*t_r - 12*c_r**3*c_t*spLE**2*t_t + 9*c_r**3*c_t*spLE*t_r + 3*c_r**3*c_t*spLE*t_t + 48*c_r**3*c_t*spTE**3*t_r + 16*c_r**3*c_t*spTE**3*t_t - 108*c_r**3*c_t*spTE**2*t_r - 36*c_r**3*c_t*spTE**2*t_t + 81*c_r**3*c_t*spTE*t_r + 27*c_r**3*c_t*spTE*t_t - 21*c_r**3*c_t*t_r - 7*c_r**3*c_t*t_t - 160*c_r**3*sk*spLE**3 + 240*c_r**3*sk*spLE**2*t_r + 120*c_r**3*sk*spLE**2 - 120*c_r**3*sk*spLE*t_r - 30*c_r**3*sk*spLE - 160*c_r**3*sk*spTE**3 + 240*c_r**3*sk*spTE**2*t_r + 360*c_r**3*sk*spTE**2 - 360*c_r**3*sk*spTE*t_r - 270*c_r**3*sk*spTE + 150*c_r**3*sk*t_r + 70*c_r**3*sk - 160*c_r**3*spLE**3*spk + 240*c_r**3*spLE**2*spk*t_r + 120*c_r**3*spLE**2*spk - 120*c_r**3*spLE*spk*t_r - 30*c_r**3*spLE*spk - 160*c_r**3*spTE**3*spk + 240*c_r**3*spTE**2*spk*t_r + 360*c_r**3*spTE**2*spk - 360*c_r**3*spTE*spk*t_r - 270*c_r**3*spTE*spk + 150*c_r**3*spk*t_r + 70*c_r**3*spk + 32*c_r**2*c_t**2*spLE**3*t_r + 32*c_r**2*c_t**2*spLE**3*t_t - 24*c_r**2*c_t**2*spLE**2*t_r - 24*c_r**2*c_t**2*spLE**2*t_t + 6*c_r**2*c_t**2*spLE*t_r + 6*c_r**2*c_t**2*spLE*t_t + 32*c_r**2*c_t**2*spTE**3*t_r + 32*c_r**2*c_t**2*spTE**3*t_t - 72*c_r**2*c_t**2*spTE**2*t_r - 72*c_r**2*c_t**2*spTE**2*t_t + 54*c_r**2*c_t**2*spTE*t_r + 54*c_r**2*c_t**2*spTE*t_t - 14*c_r**2*c_t**2*t_r - 14*c_r**2*c_t**2*t_t - 160*c_r**2*c_t*sk*spLE**3 + 160*c_r**2*c_t*sk*spLE**2*t_r + 80*c_r**2*c_t*sk*spLE**2*t_t + 120*c_r**2*c_t*sk*spLE**2 - 80*c_r**2*c_t*sk*spLE*t_r - 40*c_r**2*c_t*sk*spLE*t_t - 30*c_r**2*c_t*sk*spLE - 160*c_r**2*c_t*sk*spTE**3 + 160*c_r**2*c_t*sk*spTE**2*t_r + 80*c_r**2*c_t*sk*spTE**2*t_t + 360*c_r**2*c_t*sk*spTE**2 - 240*c_r**2*c_t*sk*spTE*t_r - 120*c_r**2*c_t*sk*spTE*t_t - 270*c_r**2*c_t*sk*spTE + 100*c_r**2*c_t*sk*t_r + 50*c_r**2*c_t*sk*t_t + 70*c_r**2*c_t*sk - 160*c_r**2*c_t*spLE**3*spk + 160*c_r**2*c_t*spLE**2*spk*t_r + 80*c_r**2*c_t*spLE**2*spk*t_t + 120*c_r**2*c_t*spLE**2*spk - 80*c_r**2*c_t*spLE*spk*t_r - 40*c_r**2*c_t*spLE*spk*t_t - 30*c_r**2*c_t*spLE*spk - 160*c_r**2*c_t*spTE**3*spk + 160*c_r**2*c_t*spTE**2*spk*t_r + 80*c_r**2*c_t*spTE**2*spk*t_t + 360*c_r**2*c_t*spTE**2*spk - 240*c_r**2*c_t*spTE*spk*t_r - 120*c_r**2*c_t*spTE*spk*t_t - 270*c_r**2*c_t*spTE*spk + 100*c_r**2*c_t*spk*t_r + 50*c_r**2*c_t*spk*t_t + 70*c_r**2*c_t*spk - 640*c_r**2*sk**2*spLE**2 + 320*c_r**2*sk**2*spLE*t_r + 320*c_r**2*sk**2*spLE - 640*c_r**2*sk**2*spTE**2 + 320*c_r**2*sk**2*spTE*t_r + 960*c_r**2*sk**2*spTE - 320*c_r**2*sk**2*t_r - 400*c_r**2*sk**2 - 1280*c_r**2*sk*spLE**2*spk + 640*c_r**2*sk*spLE*spk*t_r + 640*c_r**2*sk*spLE*spk - 1280*c_r**2*sk*spTE**2*spk + 640*c_r**2*sk*spTE*spk*t_r + 1920*c_r**2*sk*spTE*spk - 640*c_r**2*sk*spk*t_r - 800*c_r**2*sk*spk - 640*c_r**2*spLE**2*spk**2 + 320*c_r**2*spLE*spk**2*t_r + 320*c_r**2*spLE*spk**2 - 320*c_r**2*spLE*t_r*x_cg**2 - 320*c_r**2*spLE*t_r*y_cg**2 - 640*c_r**2*spTE**2*spk**2 + 320*c_r**2*spTE*spk**2*t_r + 960*c_r**2*spTE*spk**2 - 320*c_r**2*spTE*t_r*x_cg**2 - 320*c_r**2*spTE*t_r*y_cg**2 - 320*c_r**2*spk**2*t_r - 400*c_r**2*spk**2 + 320*c_r**2*t_r*x_cg**2 + 320*c_r**2*t_r*y_cg**2 + 16*c_r*c_t**3*spLE**3*t_r + 48*c_r*c_t**3*spLE**3*t_t - 12*c_r*c_t**3*spLE**2*t_r - 36*c_r*c_t**3*spLE**2*t_t + 3*c_r*c_t**3*spLE*t_r + 9*c_r*c_t**3*spLE*t_t + 16*c_r*c_t**3*spTE**3*t_r + 48*c_r*c_t**3*spTE**3*t_t - 36*c_r*c_t**3*spTE**2*t_r - 108*c_r*c_t**3*spTE**2*t_t + 27*c_r*c_t**3*spTE*t_r + 81*c_r*c_t**3*spTE*t_t - 7*c_r*c_t**3*t_r - 21*c_r*c_t**3*t_t - 160*c_r*c_t**2*sk*spLE**3 + 80*c_r*c_t**2*sk*spLE**2*t_r + 160*c_r*c_t**2*sk*spLE**2*t_t + 120*c_r*c_t**2*sk*spLE**2 - 40*c_r*c_t**2*sk*spLE*t_r - 80*c_r*c_t**2*sk*spLE*t_t - 30*c_r*c_t**2*sk*spLE - 160*c_r*c_t**2*sk*spTE**3 + 80*c_r*c_t**2*sk*spTE**2*t_r + 160*c_r*c_t**2*sk*spTE**2*t_t + 360*c_r*c_t**2*sk*spTE**2 - 120*c_r*c_t**2*sk*spTE*t_r - 240*c_r*c_t**2*sk*spTE*t_t - 270*c_r*c_t**2*sk*spTE + 50*c_r*c_t**2*sk*t_r + 100*c_r*c_t**2*sk*t_t + 70*c_r*c_t**2*sk - 160*c_r*c_t**2*spLE**3*spk + 80*c_r*c_t**2*spLE**2*spk*t_r + 160*c_r*c_t**2*spLE**2*spk*t_t + 120*c_r*c_t**2*spLE**2*spk - 40*c_r*c_t**2*spLE*spk*t_r - 80*c_r*c_t**2*spLE*spk*t_t - 30*c_r*c_t**2*spLE*spk - 160*c_r*c_t**2*spTE**3*spk + 80*c_r*c_t**2*spTE**2*spk*t_r + 160*c_r*c_t**2*spTE**2*spk*t_t + 360*c_r*c_t**2*spTE**2*spk 
#   - 120*c_r*c_t**2*spTE*spk*t_r - 240*c_r*c_t**2*spTE*spk*t_t - 270*c_r*c_t**2*spTE*spk + 50*c_r*c_t**2*spk*t_r + 100*c_r*c_t**2*spk*t_t + 70*c_r*c_t**2*spk - 640*c_r*c_t*sk**2*spLE**2 + 160*c_r*c_t*sk**2*spLE*t_r + 160*c_r*c_t*sk**2*spLE*t_t + 320*c_r*c_t*sk**2*spLE - 640*c_r*c_t*sk**2*spTE**2 + 160*c_r*c_t*sk**2*spTE*t_r + 160*c_r*c_t*sk**2*spTE*t_t + 960*c_r*c_t*sk**2*spTE - 160*c_r*c_t*sk**2*t_r - 160*c_r*c_t*sk**2*t_t - 400*c_r*c_t*sk**2 - 1280*c_r*c_t*sk*spLE**2*spk + 320*c_r*c_t*sk*spLE*spk*t_r + 320*c_r*c_t*sk*spLE*spk*t_t + 640*c_r*c_t*sk*spLE*spk - 1280*c_r*c_t*sk*spTE**2*spk + 320*c_r*c_t*sk*spTE*spk*t_r + 320*c_r*c_t*sk*spTE*spk*t_t + 1920*c_r*c_t*sk*spTE*spk - 320*c_r*c_t*sk*spk*t_r - 320*c_r*c_t*sk*spk*t_t - 800*c_r*c_t*sk*spk - 640*c_r*c_t*spLE**2*spk**2 + 160*c_r*c_t*spLE*spk**2*t_r + 160*c_r*c_t*spLE*spk**2*t_t + 320*c_r*c_t*spLE*spk**2 - 160*c_r*c_t*spLE*t_r*x_cg**2 - 160*c_r*c_t*spLE*t_r*y_cg**2 - 160*c_r*c_t*spLE*t_t*x_cg**2 - 160*c_r*c_t*spLE*t_t*y_cg**2 - 640*c_r*c_t*spTE**2*spk**2 + 160*c_r*c_t*spTE*spk**2*t_r + 160*c_r*c_t*spTE*spk**2*t_t + 960*c_r*c_t*spTE*spk**2 - 160*c_r*c_t*spTE*t_r*x_cg**2 - 160*c_r*c_t*spTE*t_r*y_cg**2 - 160*c_r*c_t*spTE*t_t*x_cg**2 - 160*c_r*c_t*spTE*t_t*y_cg**2 - 160*c_r*c_t*spk**2*t_r - 160*c_r*c_t*spk**2*t_t - 400*c_r*c_t*spk**2 + 160*c_r*c_t*t_r*x_cg**2 + 160*c_r*c_t*t_r*y_cg**2 + 160*c_r*c_t*t_t*x_cg**2 + 160*c_r*c_t*t_t*y_cg**2 - 960*c_r*sk**3*spLE - 960*c_r*sk**3*spTE + 320*c_r*sk**3*t_r + 960*c_r*sk**3 - 2880*c_r*sk**2*spLE*spk - 2880*c_r*sk**2*spTE*spk + 960*c_r*sk**2*spk*t_r + 2880*c_r*sk**2*spk - 2880*c_r*sk*spLE*spk**2 + 960*c_r*sk*spLE*x_cg**2 + 960*c_r*sk*spLE*y_cg**2 - 2880*c_r*sk*spTE*spk**2 + 960*c_r*sk*spTE*x_cg**2 + 960*c_r*sk*spTE*y_cg**2 + 960*c_r*sk*spk**2*t_r + 2880*c_r*sk*spk**2 - 960*c_r*sk*t_r*x_cg**2 - 960*c_r*sk*t_r*y_cg**2 - 960*c_r*sk*x_cg**2 - 960*c_r*sk*y_cg**2 - 960*c_r*spLE*spk**3 + 960*c_r*spLE*spk*x_cg**2 + 960*c_r*spLE*spk*y_cg**2 - 960*c_r*spTE*spk**3 + 960*c_r*spTE*spk*x_cg**2 + 960*c_r*spTE*spk*y_cg**2 + 320*c_r*spk**3*t_r + 960*c_r*spk**3 - 960*c_r*spk*t_r*x_cg**2 - 960*c_r*spk*t_r*y_cg**2 - 960*c_r*spk*x_cg**2 - 960*c_r*spk*y_cg**2 + 64*c_t**4*spLE**3*t_t - 48*c_t**4*spLE**2*t_t + 12*c_t**4*spLE*t_t + 64*c_t**4*spTE**3*t_t - 144*c_t**4*spTE**2*t_t + 108*c_t**4*spTE*t_t - 28*c_t**4*t_t - 160*c_t**3*sk*spLE**3 + 240*c_t**3*sk*spLE**2*t_t + 120*c_t**3*sk*spLE**2 - 120*c_t**3*sk*spLE*t_t - 30*c_t**3*sk*spLE - 160*c_t**3*sk*spTE**3 + 240*c_t**3*sk*spTE**2*t_t + 360*c_t**3*sk*spTE**2 - 360*c_t**3*sk*spTE*t_t - 270*c_t**3*sk*spTE + 150*c_t**3*sk*t_t + 70*c_t**3*sk - 160*c_t**3*spLE**3*spk + 240*c_t**3*spLE**2*spk*t_t + 120*c_t**3*spLE**2*spk - 120*c_t**3*spLE*spk*t_t - 30*c_t**3*spLE*spk - 160*c_t**3*spTE**3*spk + 240*c_t**3*spTE**2*spk*t_t + 360*c_t**3*spTE**2*spk - 360*c_t**3*spTE*spk*t_t - 270*c_t**3*spTE*spk + 150*c_t**3*spk*t_t + 70*c_t**3*spk - 640*c_t**2*sk**2*spLE**2 + 320*c_t**2*sk**2*spLE*t_t + 320*c_t**2*sk**2*spLE - 640*c_t**2*sk**2*spTE**2 + 320*c_t**2*sk**2*spTE*t_t + 960*c_t**2*sk**2*spTE - 320*c_t**2*sk**2*t_t - 400*c_t**2*sk**2 - 1280*c_t**2*sk*spLE**2*spk + 640*c_t**2*sk*spLE*spk*t_t + 640*c_t**2*sk*spLE*spk - 1280*c_t**2*sk*spTE**2*spk + 640*c_t**2*sk*spTE*spk*t_t + 1920*c_t**2*sk*spTE*spk - 640*c_t**2*sk*spk*t_t - 800*c_t**2*sk*spk - 640*c_t**2*spLE**2*spk**2 + 320*c_t**2*spLE*spk**2*t_t + 320*c_t**2*spLE*spk**2 - 320*c_t**2*spLE*t_t*x_cg**2 - 320*c_t**2*spLE*t_t*y_cg**2 - 640*c_t**2*spTE**2*spk**2 + 320*c_t**2*spTE*spk**2*t_t + 960*c_t**2*spTE*spk**2 - 320*c_t**2*spTE*t_t*x_cg**2 - 320*c_t**2*spTE*t_t*y_cg**2 - 320*c_t**2*spk**2*t_t - 400*c_t**2*spk**2 + 320*c_t**2*t_t*x_cg**2 + 320*c_t**2*t_t*y_cg**2 - 960*c_t*sk**3*spLE - 960*c_t*sk**3*spTE + 320*c_t*sk**3*t_t + 960*c_t*sk**3 - 2880*c_t*sk**2*spLE*spk - 2880*c_t*sk**2*spTE*spk + 960*c_t*sk**2*spk*t_t + 2880*c_t*sk**2*spk - 2880*c_t*sk*spLE*spk**2 + 960*c_t*sk*spLE*x_cg**2 + 960*c_t*sk*spLE*y_cg**2 - 2880*c_t*sk*spTE*spk**2 + 960*c_t*sk*spTE*x_cg**2 + 960*c_t*sk*spTE*y_cg**2 + 960*c_t*sk*spk**2*t_t + 2880*c_t*sk*spk**2 - 960*c_t*sk*t_t*x_cg**2 - 960*c_t*sk*t_t*y_cg**2 - 960*c_t*sk*x_cg**2 - 960*c_t*sk*y_cg**2 - 960*c_t*spLE*spk**3 + 960*c_t*spLE*spk*x_cg**2 + 960*c_t*spLE*spk*y_cg**2 - 960*c_t*spTE*spk**3 + 960*c_t*spTE*spk*x_cg**2 + 960*c_t*spTE*spk*y_cg**2 + 320*c_t*spk**3*t_t + 960*c_t*spk**3 - 960*c_t*spk*t_t*x_cg**2 - 960*c_t*spk*t_t*y_cg**2 - 960*c_t*spk*x_cg**2 - 960*c_t*spk*y_cg**2 - 1280*sk**4 - 5120*sk**3*spk - 7680*sk**2*spk**2 + 3840*sk**2*x_cg**2 + 3840*sk**2*y_cg**2 - 5120*sk*spk**3 + 7680*sk*spk*x_cg**2 + 7680*sk*spk*y_cg**2 - 1280*spk**4 + 3840*spk**2*x_cg**2 + 3840*spk**2*y_cg**2)/960 

# Ixy_new = b/2*p*(8*b/2**2*c_r**2*spLE*t_r*tan(Λ_r) + 8*b/2**2*c_r**2*spTE*t_r*tan(Λ_r) - 8*b/2**2*c_r**2*t_r*tan(Λ_r) + 12*b/2**2*c_r*c_t*spLE*t_r*tan(Λ_r) + 12*b/2**2*c_r*c_t*spLE*t_t*tan(Λ_r) + 12*b/2**2*c_r*c_t*spTE*t_r*tan(Λ_r) + 12*b/2**2*c_r*c_t*spTE*t_t*tan(Λ_r) - 12*b/2**2*c_r*c_t*t_r*tan(Λ_r) - 12*b/2**2*c_r*c_t*t_t*tan(Λ_r) - 40*b/2**2*c_r*sk*spLE*tan(Λ_r) - 40*b/2**2*c_r*sk*spTE*tan(Λ_r) + 40*b/2**2*c_r*sk*t_r*tan(Λ_r) + 40*b/2**2*c_r*sk*tan(Λ_r) - 40*b/2**2*c_r*spLE*spk*tan(Λ_r) - 40*b/2**2*c_r*spTE*spk*tan(Λ_r) + 40*b/2**2*c_r*spk*t_r*tan(Λ_r) + 40*b/2**2*c_r*spk*tan(Λ_r) + 48*b/2**2*c_t**2*spLE*t_t*tan(Λ_r) + 48*b/2**2*c_t**2*spTE*t_t*tan(Λ_r) - 48*b/2**2*c_t**2*t_t*tan(Λ_r) - 120*b/2**2*c_t*sk*spLE*tan(Λ_r) - 120*b/2**2*c_t*sk*spTE*tan(Λ_r) + 120*b/2**2*c_t*sk*t_t*tan(Λ_r) + 120*b/2**2*c_t*sk*tan(Λ_r) - 120*b/2**2*c_t*spLE*spk*tan(Λ_r) - 120*b/2**2*c_t*spTE*spk*tan(Λ_r) + 120*b/2**2*c_t*spk*t_t*tan(Λ_r) + 120*b/2**2*c_t*spk*tan(Λ_r) - 320*b/2**2*sk**2*tan(Λ_r) - 640*b/2**2*sk*spk*tan(Λ_r) - 320*b/2**2*spk**2*tan(Λ_r) + 6*b/2*c_r**3*spLE**2*t_r - 3*b/2*c_r**3*spLE*t_r - 6*b/2*c_r**3*spTE**2*t_r + 9*b/2*c_r**3*spTE*t_r - 3*b/2*c_r**3*t_r + 8*b/2*c_r**2*c_t*spLE**2*t_r + 4*b/2*c_r**2*c_t*spLE**2*t_t - 4*b/2*c_r**2*c_t*spLE*t_r - 2*b/2*c_r**2*c_t*spLE*t_t - 8*b/2*c_r**2*c_t*spTE**2*t_r - 4*b/2*c_r**2*c_t*spTE**2*t_t + 12*b/2*c_r**2*c_t*spTE*t_r + 6*b/2*c_r**2*c_t*spTE*t_t - 4*b/2*c_r**2*c_t*t_r - 2*b/2*c_r**2*c_t*t_t - 20*b/2*c_r**2*sk*spLE**2 + 20*b/2*c_r**2*sk*spLE*t_r + 10*b/2*c_r**2*sk*spLE + 20*b/2*c_r**2*sk*spTE**2 - 20*b/2*c_r**2*sk*spTE*t_r - 30*b/2*c_r**2*sk*spTE + 10*b/2*c_r**2*sk*t_r + 10*b/2*c_r**2*sk - 20*b/2*c_r**2*spLE**2*spk + 20*b/2*c_r**2*spLE*spk*t_r + 10*b/2*c_r**2*spLE*spk + 20*b/2*c_r**2*spTE**2*spk - 20*b/2*c_r**2*spTE*spk*t_r - 30*b/2*c_r**2*spTE*spk + 10*b/2*c_r**2*spk*t_r + 10*b/2*c_r**2*spk + 6*b/2*c_r*c_t**2*spLE**2*t_r + 12*b/2*c_r*c_t**2*spLE**2*t_t - 3*b/2*c_r*c_t**2*spLE*t_r - 6*b/2*c_r*c_t**2*spLE*t_t - 6*b/2*c_r*c_t**2*spTE**2*t_r - 12*b/2*c_r*c_t**2*spTE**2*t_t + 9*b/2*c_r*c_t**2*spTE*t_r + 18*b/2*c_r*c_t**2*spTE*t_t - 3*b/2*c_r*c_t**2*t_r - 6*b/2*c_r*c_t**2*t_t - 40*b/2*c_r*c_t*sk*spLE**2 + 20*b/2*c_r*c_t*sk*spLE*t_r + 20*b/2*c_r*c_t*sk*spLE*t_t + 20*b/2*c_r*c_t*sk*spLE + 40*b/2*c_r*c_t*sk*spTE**2 - 20*b/2*c_r*c_t*sk*spTE*t_r - 20*b/2*c_r*c_t*sk*spTE*t_t - 60*b/2*c_r*c_t*sk*spTE + 10*b/2*c_r*c_t*sk*t_r + 10*b/2*c_r*c_t*sk*t_t + 20*b/2*c_r*c_t*sk - 40*b/2*c_r*c_t*spLE**2*spk + 20*b/2*c_r*c_t*spLE*spk*t_r + 20*b/2*c_r*c_t*spLE*spk*t_t + 20*b/2*c_r*c_t*spLE*spk + 40*b/2*c_r*c_t*spTE**2*spk - 20*b/2*c_r*c_t*spTE*spk*t_r - 20*b/2*c_r*c_t*spTE*spk*t_t - 60*b/2*c_r*c_t*spTE*spk + 10*b/2*c_r*c_t*spk*t_r + 10*b/2*c_r*c_t*spk*t_t + 20*b/2*c_r*c_t*spk - 80*b/2*c_r*sk**2*spLE + 80*b/2*c_r*sk**2*spTE - 40*b/2*c_r*sk**2 - 160*b/2*c_r*sk*spLE*spk + 160*b/2*c_r*sk*spTE*spk - 80*b/2*c_r*sk*spk - 80*b/2*c_r*spLE*spk**2 + 80*b/2*c_r*spTE*spk**2 - 40*b/2*c_r*spk**2 + 24*b/2*c_t**3*spLE**2*t_t - 12*b/2*c_t**3*spLE*t_t - 24*b/2*c_t**3*spTE**2*t_t + 36*b/2*c_t**3*spTE*t_t - 12*b/2*c_t**3*t_t - 60*b/2*c_t**2*sk*spLE**2 + 60*b/2*c_t**2*sk*spLE*t_t + 30*b/2*c_t**2*sk*spLE + 60*b/2*c_t**2*sk*spTE**2 - 60*b/2*c_t**2*sk*spTE*t_t - 90*b/2*c_t**2*sk*spTE + 30*b/2*c_t**2*sk*t_t + 30*b/2*c_t**2*sk - 60*b/2*c_t**2*spLE**2*spk + 60*b/2*c_t**2*spLE*spk*t_t + 30*b/2*c_t**2*spLE*spk + 60*b/2*c_t**2*spTE**2*spk - 60*b/2*c_t**2*spTE*spk*t_t - 90*b/2*c_t**2*spTE*spk + 30*b/2*c_t**2*spk*t_t + 30*b/2*c_t**2*spk - 160*b/2*c_t*sk**2*spLE + 160*b/2*c_t*sk**2*spTE - 80*b/2*c_t*sk**2 - 320*b/2*c_t*sk*spLE*spk + 320*b/2*c_t*sk*spTE*spk - 160*b/2*c_t*sk*spk - 160*b/2*c_t*spLE*spk**2 + 160*b/2*c_t*spTE*spk**2 - 80*b/2*c_t*spk**2 + 80*c_r**2*spLE*t_r*x_cg*y_cg + 80*c_r**2*spTE*t_r*x_cg*y_cg - 80*c_r**2*t_r*x_cg*y_cg + 40*c_r*c_t*spLE*t_r*x_cg*y_cg + 40*c_r*c_t*spLE*t_t*x_cg*y_cg + 40*c_r*c_t*spTE*t_r*x_cg*y_cg + 40*c_r*c_t*spTE*t_t*x_cg*y_cg - 40*c_r*c_t*t_r*x_cg*y_cg - 40*c_r*c_t*t_t*x_cg*y_cg - 240*c_r*sk*spLE*x_cg*y_cg - 240*c_r*sk*spTE*x_cg*y_cg + 240*c_r*sk*t_r*x_cg*y_cg + 240*c_r*sk*x_cg*y_cg - 240*c_r*spLE*spk*x_cg*y_cg - 240*c_r*spTE*spk*x_cg*y_cg + 240*c_r*spk*t_r*x_cg*y_cg + 240*c_r*spk*x_cg*y_cg + 80*c_t**2*spLE*t_t*x_cg*y_cg + 80*c_t**2*spTE*t_t*x_cg*y_cg - 80*c_t**2*t_t*x_cg*y_cg - 240*c_t*sk*spLE*x_cg*y_cg - 240*c_t*sk*spTE*x_cg*y_cg + 240*c_t*sk*t_t*x_cg*y_cg + 240*c_t*sk*x_cg*y_cg - 240*c_t*spLE*spk*x_cg*y_cg - 240*c_t*spTE*spk*x_cg*y_cg + 240*c_t*spk*t_t*x_cg*y_cg + 240*c_t*spk*x_cg*y_cg - 960*sk**2*x_cg*y_cg - 1920*sk*spk*x_cg*y_cg - 960*spk**2*x_cg*y_cg)/240 

# Ixz_new = b/2*p*x_cg*z_cg*(2*c_r**2*spLE*t_r + 2*c_r**2*spTE*t_r - 2*c_r**2*t_r + c_r*c_t*spLE*t_r + c_r*c_t*spLE*t_t + c_r*c_t*spTE*t_r + c_r*c_t*spTE*t_t - c_r*c_t*t_r - c_r*c_t*t_t - 6*c_r*sk*spLE - 6*c_r*sk*spTE + 6*c_r*sk*t_r + 6*c_r*sk - 6*c_r*spLE*spk - 6*c_r*spTE*spk + 6*c_r*spk*t_r + 6*c_r*spk + 2*c_t**2*spLE*t_t + 2*c_t**2*spTE*t_t - 2*c_t**2*t_t - 6*c_t*sk*spLE - 6*c_t*sk*spTE + 6*c_t*sk*t_t + 6*c_t*sk - 6*c_t*spLE*spk - 6*c_t*spTE*spk + 6*c_t*spk*t_t + 6*c_t*spk - 24*sk**2 - 48*sk*spk - 24*spk**2)/6 

# Iyz_new = b/2*p*y_cg*z_cg*(2*c_r**2*spLE*t_r + 2*c_r**2*spTE*t_r - 2*c_r**2*t_r + c_r*c_t*spLE*t_r + c_r*c_t*spLE*t_t + c_r*c_t*spTE*t_r + c_r*c_t*spTE*t_t - c_r*c_t*t_r - c_r*c_t*t_t - 6*c_r*sk*spLE - 6*c_r*sk*spTE + 6*c_r*sk*t_r + 6*c_r*sk - 6*c_r*spLE*spk - 6*c_r*spTE*spk + 6*c_r*spk*t_r + 6*c_r*spk + 2*c_t**2*spLE*t_t + 2*c_t**2*spTE*t_t - 2*c_t**2*t_t - 6*c_t*sk*spLE - 6*c_t*sk*spTE + 6*c_t*sk*t_t + 6*c_t*sk - 6*c_t*spLE*spk - 6*c_t*spTE*spk + 6*c_t*spk*t_t + 6*c_t*spk - 24*sk**2 - 48*sk*spk - 24*spk**2)/6 
















quit()

# replacement variables
Ixxs = np.zeros((len(names),len(types)))
Iyys = np.zeros((len(names),len(types)))
Izzs = np.zeros((len(names),len(types)))
Ixys = np.zeros((len(names),len(types)))
Ixzs = np.zeros((len(names),len(types)))
Iyzs = np.zeros((len(names),len(types)))
Ixxs[0,1] =  0.27104439671784700; Iyys[0,1] =  0.03075003481071670
Izzs[0,1] =  0.30039189780568200; Ixys[0,1] = -0.06745635886119230
Ixzs[0,1] = -0.00306973581152483; Iyzs[0,1] =  0.00380336762603344
Ixxs[1,1] =  0.27212442282588400; Iyys[1,1] =  0.03121298066761980
Izzs[1,1] =  0.29977679119786200; Ixys[1,1] = -0.06846910704295400
Ixzs[1,1] = -0.00396349754460123; Iyzs[1,1] = -0.00144473643314478
Ixxs[2,1] =  0.27054092901100300; Iyys[2,1] =  0.03107646267172250
Izzs[2,1] =  0.30122376701684600; Ixys[2,1] = -0.06821436097470010
Ixzs[2,1] =  0.00000000000000000; Iyzs[2,1] =  0.00000000000000000
Ixxs[3,1] =  0.27104434077205200; Iyys[3,1] =  0.01401005843227450
Izzs[3,1] =  0.28365186610306500; Ixys[3,1] =  0.00498199695406229
Ixzs[3,1] = -0.00205063436315037; Iyzs[3,1] =  0.00380316062659290

for i in range(len(names)):
    Ixxs[i,0] = Ixx_new.subs(vals[i]).subs(xcg,xcgs[i,0]).subs(ycg,\
        ycgs[i,0]).subs(zcg,zcgs[i,0])
    Iyys[i,0] = Iyy_new.subs(vals[i]).subs(xcg,xcgs[i,0]).subs(ycg,\
        ycgs[i,0]).subs(zcg,zcgs[i,0])
    Izzs[i,0] = Izz_new.subs(vals[i]).subs(xcg,xcgs[i,0]).subs(ycg,\
        ycgs[i,0]).subs(zcg,zcgs[i,0])
    Ixys[i,0] = Ixy_new.subs(vals[i]).subs(xcg,xcgs[i,0]).subs(ycg,\
        ycgs[i,0]).subs(zcg,zcgs[i,0])
    Ixzs[i,0] = Ixz_new.subs(vals[i]).subs(xcg,xcgs[i,0]).subs(ycg,\
        ycgs[i,0]).subs(zcg,zcgs[i,0])
    Iyzs[i,0] = Iyz_new.subs(vals[i]).subs(xcg,xcgs[i,0]).subs(ycg,\
        ycgs[i,0]).subs(zcg,zcgs[i,0])
    
    # percent difference
    if Ixxs[i,1] != 0.0:
        Ixxs[i,2] = (Ixxs[i,1]-Ixxs[i,0]) / Ixxs[i,1] * 100.0
    if Iyys[i,1] != 0.0:
        Iyys[i,2] = (Iyys[i,1]-Iyys[i,0]) / Iyys[i,1] * 100.0
    if Izzs[i,1] != 0.0:
        Izzs[i,2] = (Izzs[i,1]-Izzs[i,0]) / Izzs[i,1] * 100.0
    if Ixys[i,1] != 0.0:
        Ixys[i,2] = (Ixys[i,1]-Ixys[i,0]) / Ixys[i,1] * 100.0
    if Ixzs[i,1] != 0.0:
        Ixzs[i,2] = (Ixzs[i,1]-Ixzs[i,0]) / Ixzs[i,1] * 100.0
    if Iyzs[i,1] != 0.0:
        Iyzs[i,2] = (Iyzs[i,1]-Iyzs[i,0]) / Iyzs[i,1] * 100.0

print_table("Ixx",names,types,Ixxs)
print_table("Iyy",names,types,Iyys)
print_table("Izz",names,types,Izzs)
print_table("Ixy",names,types,Ixys)
print_table("Ixz",names,types,Ixzs)
print_table("Iyz",names,types,Iyzs)

# =xcg= -------all------- ---only mount---- ----no twist----- ----no sweep-----
#  me     -1.405723272365   -1.387081725555   -1.413822750967   -0.435312683560
#  SW     -1.405726540000   -1.387081730000   -1.413822750000   -0.435313130000
#  PD      0.000232451690    0.000000320449   -0.000000068380    0.000102555993
#             +++++                                                 +++++      
# =ycg= -------all------- ---only mount---- ----no twist----- ----no sweep-----
#  me      3.621621621622    3.621621621622    3.621621621622    3.621621621622
#  SW      3.621621150000    3.621621620000    3.621621620000    3.621620810000
#  PD     -0.000013022390   -0.000000044776   -0.000000044776   -0.000022410453
                                                                             
# =zcg= -------all------- ---only mount---- ----no twist----- ----no sweep-----
#  me      0.074666624280    0.151655891255    0.000000000000    0.074666624280
#  SW      0.074630280000    0.151655890000    0.000000000000    0.074657440000
#  PD     -0.048699107050   -0.000000827563    0.000000000000   -0.012301895094
#             -----                                                 -----  

# =Ixx= -------all------- ---only mount---- ----no twist----- ----no sweep-----
#  me      0.271045101355    0.272124087696    0.270540595784    0.271045101355
#  SW      0.271044396718    0.272124422826    0.270540929011    0.271044340772
#  PD     -0.000259971025    0.000123153268    0.000123170723   -0.000280611906
#             -----             +++++             +++++             -----      
# =Iyy= -------all------- ---only mount---- ----no twist----- ----no sweep-----
#  me      0.031970582878    0.033764820945    0.031076424337    0.014010195110
#  SW      0.030750034811    0.031212980668    0.031076462672    0.014010058432
#  PD     -3.969257513754   -8.175573825663    0.000123357352   -0.000975571836
#             -----             -----             +++++             -----      
# =Izz= -------all------- ---only mount---- ----no twist----- ----no sweep-----
#  me      0.301613049091    0.302328300817    0.301223396120    0.283652661323
#  SW      0.300391897806    0.299776791198    0.301223767017    0.283651866103
#  PD     -0.406519381510   -0.851136476923    0.000123129873   -0.000280350780
#             -----             -----             +++++             -----      
# =Ixy= -------all------- ---only mount---- ----no twist----- ----no sweep-----
#  me     -0.067456319176   -0.068469022626   -0.068214276922    0.004982079368
#  SW     -0.067456358861   -0.068469107043   -0.068214360975    0.004981996954
#  PD      0.000058831637    0.000123291304    0.000123218457   -0.001654237544
#                               +++++             +++++             -----      
# =Ixz= -------all------- ---only mount---- ----no twist----- ----no sweep-----
#  me      0.001702524337    0.003272719133    0.000000000000   -0.002050762914
#  SW     -0.003069735812   -0.003963497545    0.000000000000   -0.002050634363
#  PD    155.461591522120  182.571493900406    0.000000000000   -0.006268837119
#             +++++             +++++                               -----      
# =Iyz= -------all------- ---only mount---- ----no twist----- ----no sweep-----
#  me      0.003803223946   -0.001444734682    0.000000000000    0.003803223946
#  SW      0.003803367626   -0.001444736433    0.000000000000    0.003803160627
#  PD      0.003777694680    0.000121177306    0.000000000000   -0.001664926368
#             +++++             +++++                               -----   

quit()

# write out to text file as well as latex version
# print(sy.latex(Ixx_new.subs(w,sym("\\omega")).subs(L,sym("\\Lambda"))))




# # sympy
# replacement variables
# b_2 = 8.0
# cr = 2.0
# ct = 1.50
# tr = tt = 0.12
# Lr = 15.0
# p = 0.0175
# tanL = sy.tan(np.deg2rad(Lr))
# wr = np.deg2rad(2.0)
# w = np.deg2rad(18.0)
# xcg = -1.38708173
# ycg = 3.62162162
# zcg = 0.15165589

# Ixx_cpy = sy.Piecewise((b_2*p*(-1280*cr**2*tr*ycg**2 - 1280*cr**2*tr*zcg**2 - 640*cr*ct*tr*ycg**2 - 640*cr*ct*tr*zcg**2 - 640*cr*ct*tt*ycg**2 - 640*cr*ct*tt*zcg**2 - 1280*ct**2*tt*ycg**2 - 1280*ct**2*tt*zcg**2 + (-120*cr**4*tr**3*sy.sin(2*wr) + 120*cr**4*tr**3*sy.sin(2*w + 2*wr) + 210*cr**4*tr*sy.sin(2*wr) - 210*cr**4*tr*sy.sin(2*w + 2*wr) + 120*cr**3*ct*tr**3*sy.sin(2*wr) - 120*cr**3*ct*tr**3*sy.sin(2*w + 2*wr) + 360*cr**3*ct*tr**2*tt*sy.sin(2*wr) - 360*cr**3*ct*tr**2*tt*sy.sin(2*w + 2*wr) - 630*cr**3*ct*tr*sy.sin(2*wr) + 630*cr**3*ct*tr*sy.sin(2*w + 2*wr) - 210*cr**3*ct*tt*sy.sin(2*wr) + 210*cr**3*ct*tt*sy.sin(2*w + 2*wr) - 360*cr**2*ct**2*tr**2*tt*sy.sin(2*wr) + 360*cr**2*ct**2*tr**2*tt*sy.sin(2*w + 2*wr) - 360*cr**2*ct**2*tr*tt**2*sy.sin(2*wr) + 360*cr**2*ct**2*tr*tt**2*sy.sin(2*w + 2*wr) + 630*cr**2*ct**2*tr*sy.sin(2*wr) - 630*cr**2*ct**2*tr*sy.sin(2*w + 2*wr) + 630*cr**2*ct**2*tt*sy.sin(2*wr) - 630*cr**2*ct**2*tt*sy.sin(2*w + 2*wr) + 360*cr*ct**3*tr*tt**2*sy.sin(2*wr) - 360*cr*ct**3*tr*tt**2*sy.sin(2*w + 2*wr) - 210*cr*ct**3*tr*sy.sin(2*wr) + 210*cr*ct**3*tr*sy.sin(2*w + 2*wr) + 120*cr*ct**3*tt**3*sy.sin(2*wr) - 120*cr*ct**3*tt**3*sy.sin(2*w + 2*wr) - 630*cr*ct**3*tt*sy.sin(2*wr) + 630*cr*ct**3*tt*sy.sin(2*w + 2*wr) - 120*ct**4*tt**3*sy.sin(2*wr) + 120*ct**4*tt**3*sy.sin(2*w + 2*wr) + 210*ct**4*tt*sy.sin(2*wr) - 210*ct**4*tt*sy.sin(2*w + 2*wr) + 2*w**5*(64*b_2**2*cr**2*tr + 96*b_2**2*cr*ct*tr + 96*b_2**2*cr*ct*tt + 384*b_2**2*ct**2*tt + 16*cr**4*tr**3 + 28*cr**4*tr + 4*cr**3*ct*tr**3 + 12*cr**3*ct*tr**2*tt + 21*cr**3*ct*tr + 7*cr**3*ct*tt + 8*cr**2*ct**2*tr**2*tt + 8*cr**2*ct**2*tr*tt**2 + 14*cr**2*ct**2*tr + 14*cr**2*ct**2*tt + 12*cr*ct**3*tr*tt**2 + 7*cr*ct**3*tr + 4*cr*ct**3*tt**3 + 21*cr*ct**3*tt + 16*ct**4*tt**3 + 28*ct**4*tt) + 20*w**4*(-4*cr**4*tr**3*sy.sin(2*wr) + 7*cr**4*tr*sy.sin(2*wr) + 4*ct**4*tt**3*sy.sin(2*w + 2*wr) - 7*ct**4*tt*sy.sin(2*w + 2*wr)) + 10*w**3*(16*cr**4*tr**3*sy.cos(2*wr) - 28*cr**4*tr*sy.cos(2*wr) - 4*cr**3*ct*tr**3*sy.cos(2*wr) - 12*cr**3*ct*tr**2*tt*sy.cos(2*wr) + 21*cr**3*ct*tr*sy.cos(2*wr) + 7*cr**3*ct*tt*sy.cos(2*wr) - 12*cr*ct**3*tr*tt**2*sy.cos(2*w + 2*wr) + 7*cr*ct**3*tr*sy.cos(2*w + 2*wr) - 4*cr*ct**3*tt**3*sy.cos(2*w + 2*wr) + 21*cr*ct**3*tt*sy.cos(2*w + 2*wr) + 16*ct**4*tt**3*sy.cos(2*w + 2*wr) - 28*ct**4*tt*sy.cos(2*w + 2*wr)) + 30*w**2*(8*cr**4*tr**3*sy.sin(2*wr) - 14*cr**4*tr*sy.sin(2*wr) - 4*cr**3*ct*tr**3*sy.sin(2*wr) - 12*cr**3*ct*tr**2*tt*sy.sin(2*wr) + 21*cr**3*ct*tr*sy.sin(2*wr) + 7*cr**3*ct*tt*sy.sin(2*wr) + 4*cr**2*ct**2*tr**2*tt*sy.sin(2*wr) - 4*cr**2*ct**2*tr**2*tt*sy.sin(2*w + 2*wr) + 4*cr**2*ct**2*tr*tt**2*sy.sin(2*wr) - 4*cr**2*ct**2*tr*tt**2*sy.sin(2*w + 2*wr) - 7*cr**2*ct**2*tr*sy.sin(2*wr) + 7*cr**2*ct**2*tr*sy.sin(2*w + 2*wr) - 7*cr**2*ct**2*tt*sy.sin(2*wr) + 7*cr**2*ct**2*tt*sy.sin(2*w + 2*wr) + 12*cr*ct**3*tr*tt**2*sy.sin(2*w + 2*wr) - 7*cr*ct**3*tr*sy.sin(2*w + 2*wr) + 4*cr*ct**3*tt**3*sy.sin(2*w + 2*wr) - 21*cr*ct**3*tt*sy.sin(2*w + 2*wr) - 8*ct**4*tt**3*sy.sin(2*w + 2*wr) + 14*ct**4*tt*sy.sin(2*w + 2*wr)) + 15*w*(-16*cr**4*tr**3*sy.cos(2*wr) + 28*cr**4*tr*sy.cos(2*wr) + 12*cr**3*ct*tr**3*sy.cos(2*wr) + 4*cr**3*ct*tr**3*sy.cos(2*w + 2*wr) + 36*cr**3*ct*tr**2*tt*sy.cos(2*wr) + 12*cr**3*ct*tr**2*tt*sy.cos(2*w + 2*wr) - 63*cr**3*ct*tr*sy.cos(2*wr) - 21*cr**3*ct*tr*sy.cos(2*w + 2*wr) - 21*cr**3*ct*tt*sy.cos(2*wr) - 7*cr**3*ct*tt*sy.cos(2*w + 2*wr) - 24*cr**2*ct**2*tr**2*tt*sy.cos(2*wr) - 24*cr**2*ct**2*tr**2*tt*sy.cos(2*w + 2*wr) - 24*cr**2*ct**2*tr*tt**2*sy.cos(2*wr) - 24*cr**2*ct**2*tr*tt**2*sy.cos(2*w + 2*wr) + 42*cr**2*ct**2*tr*sy.cos(2*wr) + 42*cr**2*ct**2*tr*sy.cos(2*w + 2*wr) + 42*cr**2*ct**2*tt*sy.cos(2*wr) + 42*cr**2*ct**2*tt*sy.cos(2*w + 2*wr) + 12*cr*ct**3*tr*tt**2*sy.cos(2*wr) + 36*cr*ct**3*tr*tt**2*sy.cos(2*w + 2*wr) - 7*cr*ct**3*tr*sy.cos(2*wr) - 21*cr*ct**3*tr*sy.cos(2*w + 2*wr) + 4*cr*ct**3*tt**3*sy.cos(2*wr) + 12*cr*ct**3*tt**3*sy.cos(2*w + 2*wr) - 21*cr*ct**3*tt*sy.cos(2*wr) - 63*cr*ct**3*tt*sy.cos(2*w + 2*wr) - 16*ct**4*tt**3*sy.cos(2*w + 2*wr) + 28*ct**4*tt*sy.cos(2*w + 2*wr)))/w**5)/3840, (w > -sy.oo) & (w < sy.oo) & (w != 0)), (b_2*p*(32*b_2**2*cr**2*tr + 48*b_2**2*cr*ct*tr + 48*b_2**2*cr*ct*tt + 192*b_2**2*ct**2*tt + 16*cr**4*tr**3*sy.cos(wr)**2 - 28*cr**4*tr*sy.cos(wr)**2 + 28*cr**4*tr + 4*cr**3*ct*tr**3*sy.cos(wr)**2 + 12*cr**3*ct*tr**2*tt*sy.cos(wr)**2 - 21*cr**3*ct*tr*sy.cos(wr)**2 + 21*cr**3*ct*tr - 7*cr**3*ct*tt*sy.cos(wr)**2 + 7*cr**3*ct*tt + 8*cr**2*ct**2*tr**2*tt*sy.cos(wr)**2 + 8*cr**2*ct**2*tr*tt**2*sy.cos(wr)**2 - 14*cr**2*ct**2*tr*sy.cos(wr)**2 + 14*cr**2*ct**2*tr - 14*cr**2*ct**2*tt*sy.cos(wr)**2 + 14*cr**2*ct**2*tt - 320*cr**2*tr*ycg**2 - 320*cr**2*tr*zcg**2 + 12*cr*ct**3*tr*tt**2*sy.cos(wr)**2 - 7*cr*ct**3*tr*sy.cos(wr)**2 + 7*cr*ct**3*tr + 4*cr*ct**3*tt**3*sy.cos(wr)**2 - 21*cr*ct**3*tt*sy.cos(wr)**2 + 21*cr*ct**3*tt - 160*cr*ct*tr*ycg**2 - 160*cr*ct*tr*zcg**2 - 160*cr*ct*tt*ycg**2 - 160*cr*ct*tt*zcg**2 + 16*ct**4*tt**3*sy.cos(wr)**2 - 28*ct**4*tt*sy.cos(wr)**2 + 28*ct**4*tt - 320*ct**2*tt*ycg**2 - 320*ct**2*tt*zcg**2)/960, True))
# Iyy_cpy = b_2*p*(32*b_2**2*cr**2*tr*tanL**2 + 48*b_2**2*cr*ct*tr*tanL**2 + 48*b_2**2*cr*ct*tt*tanL**2 + 192*b_2**2*ct**2*tt*tanL**2 + 24*b_2*cr**3*tr*tanL + 32*b_2*cr**2*ct*tr*tanL + 16*b_2*cr**2*ct*tt*tanL + 24*b_2*cr*ct**2*tr*tanL + 48*b_2*cr*ct**2*tt*tanL + 96*b_2*ct**3*tt*tanL + 16*cr**4*tr**3 + 28*cr**4*tr + 4*cr**3*ct*tr**3 + 12*cr**3*ct*tr**2*tt + 21*cr**3*ct*tr + 7*cr**3*ct*tt + 8*cr**2*ct**2*tr**2*tt + 8*cr**2*ct**2*tr*tt**2 + 14*cr**2*ct**2*tr + 14*cr**2*ct**2*tt - 320*cr**2*tr*xcg**2 - 320*cr**2*tr*zcg**2 + 12*cr*ct**3*tr*tt**2 + 7*cr*ct**3*tr + 4*cr*ct**3*tt**3 + 21*cr*ct**3*tt - 160*cr*ct*tr*xcg**2 - 160*cr*ct*tr*zcg**2 - 160*cr*ct*tt*xcg**2 - 160*cr*ct*tt*zcg**2 + 16*ct**4*tt**3 + 28*ct**4*tt - 320*ct**2*tt*xcg**2 - 320*ct**2*tt*zcg**2)/960 
# Izz_cpy = sy.Piecewise((b_2*p*(128*b_2**2*cr**2*tr*tanL**2 + 192*b_2**2*cr*ct*tr*tanL**2 + 192*b_2**2*cr*ct*tt*tanL**2 + 768*b_2**2*ct**2*tt*tanL**2 + 96*b_2*cr**3*tr*tanL + 128*b_2*cr**2*ct*tr*tanL + 64*b_2*cr**2*ct*tt*tanL + 96*b_2*cr*ct**2*tr*tanL + 192*b_2*cr*ct**2*tt*tanL + 384*b_2*ct**3*tt*tanL - 1280*cr**2*tr*xcg**2 - 1280*cr**2*tr*ycg**2 - 640*cr*ct*tr*xcg**2 - 640*cr*ct*tr*ycg**2 - 640*cr*ct*tt*xcg**2 - 640*cr*ct*tt*ycg**2 - 1280*ct**2*tt*xcg**2 - 1280*ct**2*tt*ycg**2 + (120*cr**4*tr**3*sy.sin(2*wr) - 120*cr**4*tr**3*sy.sin(2*w + 2*wr) - 210*cr**4*tr*sy.sin(2*wr) + 210*cr**4*tr*sy.sin(2*w + 2*wr) - 120*cr**3*ct*tr**3*sy.sin(2*wr) + 120*cr**3*ct*tr**3*sy.sin(2*w + 2*wr) - 360*cr**3*ct*tr**2*tt*sy.sin(2*wr) + 360*cr**3*ct*tr**2*tt*sy.sin(2*w + 2*wr) + 630*cr**3*ct*tr*sy.sin(2*wr) - 630*cr**3*ct*tr*sy.sin(2*w + 2*wr) + 210*cr**3*ct*tt*sy.sin(2*wr) - 210*cr**3*ct*tt*sy.sin(2*w + 2*wr) + 360*cr**2*ct**2*tr**2*tt*sy.sin(2*wr) - 360*cr**2*ct**2*tr**2*tt*sy.sin(2*w + 2*wr) + 360*cr**2*ct**2*tr*tt**2*sy.sin(2*wr) - 360*cr**2*ct**2*tr*tt**2*sy.sin(2*w + 2*wr) - 630*cr**2*ct**2*tr*sy.sin(2*wr) + 630*cr**2*ct**2*tr*sy.sin(2*w + 2*wr) - 630*cr**2*ct**2*tt*sy.sin(2*wr) + 630*cr**2*ct**2*tt*sy.sin(2*w + 2*wr) - 360*cr*ct**3*tr*tt**2*sy.sin(2*wr) + 360*cr*ct**3*tr*tt**2*sy.sin(2*w + 2*wr) + 210*cr*ct**3*tr*sy.sin(2*wr) - 210*cr*ct**3*tr*sy.sin(2*w + 2*wr) - 120*cr*ct**3*tt**3*sy.sin(2*wr) + 120*cr*ct**3*tt**3*sy.sin(2*w + 2*wr) + 630*cr*ct**3*tt*sy.sin(2*wr) - 630*cr*ct**3*tt*sy.sin(2*w + 2*wr) + 120*ct**4*tt**3*sy.sin(2*wr) - 120*ct**4*tt**3*sy.sin(2*w + 2*wr) - 210*ct**4*tt*sy.sin(2*wr) + 210*ct**4*tt*sy.sin(2*w + 2*wr) + 2*w**5*(64*b_2**2*cr**2*tr + 96*b_2**2*cr*ct*tr + 96*b_2**2*cr*ct*tt + 384*b_2**2*ct**2*tt + 16*cr**4*tr**3 + 28*cr**4*tr + 4*cr**3*ct*tr**3 + 12*cr**3*ct*tr**2*tt + 21*cr**3*ct*tr + 7*cr**3*ct*tt + 8*cr**2*ct**2*tr**2*tt + 8*cr**2*ct**2*tr*tt**2 + 14*cr**2*ct**2*tr + 14*cr**2*ct**2*tt + 12*cr*ct**3*tr*tt**2 + 7*cr*ct**3*tr + 4*cr*ct**3*tt**3 + 21*cr*ct**3*tt + 16*ct**4*tt**3 + 28*ct**4*tt) + 20*w**4*(4*cr**4*tr**3*sy.sin(2*wr) - 7*cr**4*tr*sy.sin(2*wr) - 4*ct**4*tt**3*sy.sin(2*w + 2*wr) + 7*ct**4*tt*sy.sin(2*w + 2*wr)) + 10*w**3*(-16*cr**4*tr**3*sy.cos(2*wr) + 28*cr**4*tr*sy.cos(2*wr) + 4*cr**3*ct*tr**3*sy.cos(2*wr) + 12*cr**3*ct*tr**2*tt*sy.cos(2*wr) - 21*cr**3*ct*tr*sy.cos(2*wr) - 7*cr**3*ct*tt*sy.cos(2*wr) + 12*cr*ct**3*tr*tt**2*sy.cos(2*w + 2*wr) - 7*cr*ct**3*tr*sy.cos(2*w + 2*wr) + 4*cr*ct**3*tt**3*sy.cos(2*w + 2*wr) - 21*cr*ct**3*tt*sy.cos(2*w + 2*wr) - 16*ct**4*tt**3*sy.cos(2*w + 2*wr) + 28*ct**4*tt*sy.cos(2*w + 2*wr)) + 30*w**2*(-8*cr**4*tr**3*sy.sin(2*wr) + 14*cr**4*tr*sy.sin(2*wr) + 4*cr**3*ct*tr**3*sy.sin(2*wr) + 12*cr**3*ct*tr**2*tt*sy.sin(2*wr) - 21*cr**3*ct*tr*sy.sin(2*wr) - 7*cr**3*ct*tt*sy.sin(2*wr) - 4*cr**2*ct**2*tr**2*tt*sy.sin(2*wr) + 4*cr**2*ct**2*tr**2*tt*sy.sin(2*w + 2*wr) - 4*cr**2*ct**2*tr*tt**2*sy.sin(2*wr) + 4*cr**2*ct**2*tr*tt**2*sy.sin(2*w + 2*wr) + 7*cr**2*ct**2*tr*sy.sin(2*wr) - 7*cr**2*ct**2*tr*sy.sin(2*w + 2*wr) + 7*cr**2*ct**2*tt*sy.sin(2*wr) - 7*cr**2*ct**2*tt*sy.sin(2*w + 2*wr) - 12*cr*ct**3*tr*tt**2*sy.sin(2*w + 2*wr) + 7*cr*ct**3*tr*sy.sin(2*w + 2*wr) - 4*cr*ct**3*tt**3*sy.sin(2*w + 2*wr) + 21*cr*ct**3*tt*sy.sin(2*w + 2*wr) + 8*ct**4*tt**3*sy.sin(2*w + 2*wr) - 14*ct**4*tt*sy.sin(2*w + 2*wr)) + 15*w*(16*cr**4*tr**3*sy.cos(2*wr) - 28*cr**4*tr*sy.cos(2*wr) - 12*cr**3*ct*tr**3*sy.cos(2*wr) - 4*cr**3*ct*tr**3*sy.cos(2*w + 2*wr) - 36*cr**3*ct*tr**2*tt*sy.cos(2*wr) - 12*cr**3*ct*tr**2*tt*sy.cos(2*w + 2*wr) + 63*cr**3*ct*tr*sy.cos(2*wr) + 21*cr**3*ct*tr*sy.cos(2*w + 2*wr) + 21*cr**3*ct*tt*sy.cos(2*wr) + 7*cr**3*ct*tt*sy.cos(2*w + 2*wr) + 24*cr**2*ct**2*tr**2*tt*sy.cos(2*wr) + 24*cr**2*ct**2*tr**2*tt*sy.cos(2*w + 2*wr) + 24*cr**2*ct**2*tr*tt**2*sy.cos(2*wr) + 24*cr**2*ct**2*tr*tt**2*sy.cos(2*w + 2*wr) - 42*cr**2*ct**2*tr*sy.cos(2*wr) - 42*cr**2*ct**2*tr*sy.cos(2*w + 2*wr) - 42*cr**2*ct**2*tt*sy.cos(2*wr) - 42*cr**2*ct**2*tt*sy.cos(2*w + 2*wr) - 12*cr*ct**3*tr*tt**2*sy.cos(2*wr) - 36*cr*ct**3*tr*tt**2*sy.cos(2*w + 2*wr) + 7*cr*ct**3*tr*sy.cos(2*wr) + 21*cr*ct**3*tr*sy.cos(2*w + 2*wr) - 4*cr*ct**3*tt**3*sy.cos(2*wr) - 12*cr*ct**3*tt**3*sy.cos(2*w + 2*wr) + 21*cr*ct**3*tt*sy.cos(2*wr) + 63*cr*ct**3*tt*sy.cos(2*w + 2*wr) + 16*ct**4*tt**3*sy.cos(2*w + 2*wr) - 28*ct**4*tt*sy.cos(2*w + 2*wr)))/w**5)/3840, (w > -sy.oo) & (w < sy.oo) & (w != 0)), (b_2*p*(32*b_2**2*cr**2*tr/sy.cos(Lr)**2 + 48*b_2**2*cr*ct*tr/sy.cos(Lr)**2 + 48*b_2**2*cr*ct*tt/sy.cos(Lr)**2 + 192*b_2**2*ct**2*tt/sy.cos(Lr)**2 + 24*b_2*cr**3*tr*tanL + 32*b_2*cr**2*ct*tr*tanL + 16*b_2*cr**2*ct*tt*tanL + 24*b_2*cr*ct**2*tr*tanL + 48*b_2*cr*ct**2*tt*tanL + 96*b_2*ct**3*tt*tanL + 16*cr**4*tr**3*sy.sin(wr)**2 - 28*cr**4*tr*sy.sin(wr)**2 + 28*cr**4*tr + 4*cr**3*ct*tr**3*sy.sin(wr)**2 + 12*cr**3*ct*tr**2*tt*sy.sin(wr)**2 - 21*cr**3*ct*tr*sy.sin(wr)**2 + 21*cr**3*ct*tr - 7*cr**3*ct*tt*sy.sin(wr)**2 + 7*cr**3*ct*tt + 8*cr**2*ct**2*tr**2*tt*sy.sin(wr)**2 + 8*cr**2*ct**2*tr*tt**2*sy.sin(wr)**2 - 14*cr**2*ct**2*tr*sy.sin(wr)**2 + 14*cr**2*ct**2*tr - 14*cr**2*ct**2*tt*sy.sin(wr)**2 + 14*cr**2*ct**2*tt - 320*cr**2*tr*xcg**2 - 320*cr**2*tr*ycg**2 + 12*cr*ct**3*tr*tt**2*sy.sin(wr)**2 - 7*cr*ct**3*tr*sy.sin(wr)**2 + 7*cr*ct**3*tr + 4*cr*ct**3*tt**3*sy.sin(wr)**2 - 21*cr*ct**3*tt*sy.sin(wr)**2 + 21*cr*ct**3*tt - 160*cr*ct*tr*xcg**2 - 160*cr*ct*tr*ycg**2 - 160*cr*ct*tt*xcg**2 - 160*cr*ct*tt*ycg**2 + 16*ct**4*tt**3*sy.sin(wr)**2 - 28*ct**4*tt*sy.sin(wr)**2 + 28*ct**4*tt - 320*ct**2*tt*xcg**2 - 320*ct**2*tt*ycg**2)/960, True)) 
# Ixy_cpy = sy.Piecewise((b_2*p*(-2*b_2**2*cr**2*tr*tanL - 3*b_2**2*cr*ct*tr*tanL - 3*b_2**2*cr*ct*tt*tanL - 12*b_2**2*ct**2*tt*tanL + 15*b_2*(-24*cr**3*tr*sy.sin(wr) + 24*cr**3*tr*sy.sin(w + wr) + 48*cr**2*ct*tr*sy.sin(wr) - 48*cr**2*ct*tr*sy.sin(w + wr) + 24*cr**2*ct*tt*sy.sin(wr) - 24*cr**2*ct*tt*sy.sin(w + wr) - 24*cr*ct**2*tr*sy.sin(wr) + 24*cr*ct**2*tr*sy.sin(w + wr) - 48*cr*ct**2*tt*sy.sin(wr) + 48*cr*ct**2*tt*sy.sin(w + wr) - ct**3*tt*w**4*sy.sin(w + wr) + 24*ct**3*tt*sy.sin(wr) - 24*ct**3*tt*sy.sin(w + wr) + w**3*(cr**3*tr*sy.cos(wr) + cr*ct**2*tr*sy.cos(w + wr) + 2*cr*ct**2*tt*sy.cos(w + wr) - 4*ct**3*tt*sy.cos(w + wr)) + 2*w**2*(3*cr**3*tr*sy.sin(wr) - 2*cr**2*ct*tr*sy.sin(wr) + 2*cr**2*ct*tr*sy.sin(w + wr) - cr**2*ct*tt*sy.sin(wr) + cr**2*ct*tt*sy.sin(w + wr) - 3*cr*ct**2*tr*sy.sin(w + wr) - 6*cr*ct**2*tt*sy.sin(w + wr) + 6*ct**3*tt*sy.sin(w + wr)) + 6*w*(-3*cr**3*tr*sy.cos(wr) - cr**3*tr*sy.cos(w + wr) + 4*cr**2*ct*tr*sy.cos(wr) + 4*cr**2*ct*tr*sy.cos(w + wr) + 2*cr**2*ct*tt*sy.cos(wr) + 2*cr**2*ct*tt*sy.cos(w + wr) - cr*ct**2*tr*sy.cos(wr) - 3*cr*ct**2*tr*sy.cos(w + wr) - 2*cr*ct**2*tt*sy.cos(wr) - 6*cr*ct**2*tt*sy.cos(w + wr) + 4*ct**3*tt*sy.cos(w + wr)))/w**5 - 20*cr**2*tr*xcg*ycg - 10*cr*ct*tr*xcg*ycg - 10*cr*ct*tt*xcg*ycg - 20*ct**2*tt*xcg*ycg)/60, (w > -sy.oo) & (w < sy.oo) & (w != 0)), (-b_2*p*(8*b_2**2*cr**2*tr*tanL + 12*b_2**2*cr*ct*tr*tanL + 12*b_2**2*cr*ct*tt*tanL + 48*b_2**2*ct**2*tt*tanL + b_2*(3*cr**3*tr + 4*cr**2*ct*tr + 2*cr**2*ct*tt + 3*cr*ct**2*tr + 6*cr*ct**2*tt + 12*ct**3*tt)*sy.cos(wr) + 80*cr**2*tr*xcg*ycg + 40*cr*ct*tr*xcg*ycg + 40*cr*ct*tt*xcg*ycg + 80*ct**2*tt*xcg*ycg)/240, True)) 
# Ixz_cpy = sy.Piecewise((b_2*p*(-256*cr**2*tr*xcg*zcg - 128*cr*ct*tr*xcg*zcg - 128*cr*ct*tt*xcg*zcg - 256*ct**2*tt*xcg*zcg + (48*cr**4*tr**3*sy.sin(wr)**2 - 48*cr**4*tr**3*sy.sin(w + wr)**2 - 84*cr**4*tr*sy.sin(wr)**2 + 84*cr**4*tr*sy.sin(w + wr)**2 - 48*cr**3*ct*tr**3*sy.sin(wr)**2 + 48*cr**3*ct*tr**3*sy.sin(w + wr)**2 - 144*cr**3*ct*tr**2*tt*sy.sin(wr)**2 + 144*cr**3*ct*tr**2*tt*sy.sin(w + wr)**2 + 252*cr**3*ct*tr*sy.sin(wr)**2 - 252*cr**3*ct*tr*sy.sin(w + wr)**2 + 84*cr**3*ct*tt*sy.sin(wr)**2 - 84*cr**3*ct*tt*sy.sin(w + wr)**2 + 144*cr**2*ct**2*tr**2*tt*sy.sin(wr)**2 - 144*cr**2*ct**2*tr**2*tt*sy.sin(w + wr)**2 + 144*cr**2*ct**2*tr*tt**2*sy.sin(wr)**2 - 144*cr**2*ct**2*tr*tt**2*sy.sin(w + wr)**2 - 252*cr**2*ct**2*tr*sy.sin(wr)**2 + 252*cr**2*ct**2*tr*sy.sin(w + wr)**2 - 252*cr**2*ct**2*tt*sy.sin(wr)**2 + 252*cr**2*ct**2*tt*sy.sin(w + wr)**2 - 144*cr*ct**3*tr*tt**2*sy.sin(wr)**2 + 144*cr*ct**3*tr*tt**2*sy.sin(w + wr)**2 + 84*cr*ct**3*tr*sy.sin(wr)**2 - 84*cr*ct**3*tr*sy.sin(w + wr)**2 - 48*cr*ct**3*tt**3*sy.sin(wr)**2 + 48*cr*ct**3*tt**3*sy.sin(w + wr)**2 + 252*cr*ct**3*tt*sy.sin(wr)**2 - 252*cr*ct**3*tt*sy.sin(w + wr)**2 + 48*ct**4*tt**3*sy.sin(wr)**2 - 48*ct**4*tt**3*sy.sin(w + wr)**2 - 84*ct**4*tt*sy.sin(wr)**2 + 84*ct**4*tt*sy.sin(w + wr)**2 + 4*w**4*(8*cr**4*tr**3*sy.sin(wr)**2 - 4*cr**4*tr**3*sy.sin(w + wr)**2 - 4*cr**4*tr**3*sy.cos(w + wr)**2 - 14*cr**4*tr*sy.sin(wr)**2 + 7*cr**4*tr*sy.sin(w + wr)**2 + 7*cr**4*tr*sy.cos(w + wr)**2 - 4*ct**4*tt**3*sy.sin(w + wr)**2 + 4*ct**4*tt**3*sy.cos(w + wr)**2 + 7*ct**4*tt*sy.sin(w + wr)**2 - 7*ct**4*tt*sy.cos(w + wr)**2) + 2*w**3*(-16*cr**4*tr**3*sy.sin(2*wr) + 28*cr**4*tr*sy.sin(2*wr) + 4*cr**3*ct*tr**3*sy.sin(2*wr) + 12*cr**3*ct*tr**2*tt*sy.sin(2*wr) - 21*cr**3*ct*tr*sy.sin(2*wr) - 7*cr**3*ct*tt*sy.sin(2*wr) + 12*cr*ct**3*tr*tt**2*sy.sin(2*w + 2*wr) - 7*cr*ct**3*tr*sy.sin(2*w + 2*wr) + 4*cr*ct**3*tt**3*sy.sin(2*w + 2*wr) - 21*cr*ct**3*tt*sy.sin(2*w + 2*wr) - 16*ct**4*tt**3*sy.sin(2*w + 2*wr) + 28*ct**4*tt*sy.sin(2*w + 2*wr)) + 6*w**2*(-16*cr**4*tr**3*sy.sin(wr)**2 + 8*cr**4*tr**3*sy.sin(w + wr)**2 + 8*cr**4*tr**3*sy.cos(w + wr)**2 + 28*cr**4*tr*sy.sin(wr)**2 - 14*cr**4*tr*sy.sin(w + wr)**2 - 14*cr**4*tr*sy.cos(w + wr)**2 + 8*cr**3*ct*tr**3*sy.sin(wr)**2 - 4*cr**3*ct*tr**3*sy.sin(w + wr)**2 - 4*cr**3*ct*tr**3*sy.cos(w + wr)**2 + 24*cr**3*ct*tr**2*tt*sy.sin(wr)**2 - 12*cr**3*ct*tr**2*tt*sy.sin(w + wr)**2 - 12*cr**3*ct*tr**2*tt*sy.cos(w + wr)**2 - 42*cr**3*ct*tr*sy.sin(wr)**2 + 21*cr**3*ct*tr*sy.sin(w + wr)**2 + 21*cr**3*ct*tr*sy.cos(w + wr)**2 - 14*cr**3*ct*tt*sy.sin(wr)**2 + 7*cr**3*ct*tt*sy.sin(w + wr)**2 + 7*cr**3*ct*tt*sy.cos(w + wr)**2 - 8*cr**2*ct**2*tr**2*tt*sy.sin(wr)**2 + 8*cr**2*ct**2*tr**2*tt*sy.sin(w + wr)**2 - 8*cr**2*ct**2*tr*tt**2*sy.sin(wr)**2 + 8*cr**2*ct**2*tr*tt**2*sy.sin(w + wr)**2 + 14*cr**2*ct**2*tr*sy.sin(wr)**2 - 14*cr**2*ct**2*tr*sy.sin(w + wr)**2 + 14*cr**2*ct**2*tt*sy.sin(wr)**2 - 14*cr**2*ct**2*tt*sy.sin(w + wr)**2 - 12*cr*ct**3*tr*tt**2*sy.sin(w + wr)**2 + 12*cr*ct**3*tr*tt**2*sy.cos(w + wr)**2 + 7*cr*ct**3*tr*sy.sin(w + wr)**2 - 7*cr*ct**3*tr*sy.cos(w + wr)**2 - 4*cr*ct**3*tt**3*sy.sin(w + wr)**2 + 4*cr*ct**3*tt**3*sy.cos(w + wr)**2 + 21*cr*ct**3*tt*sy.sin(w + wr)**2 - 21*cr*ct**3*tt*sy.cos(w + wr)**2 + 8*ct**4*tt**3*sy.sin(w + wr)**2 - 8*ct**4*tt**3*sy.cos(w + wr)**2 - 14*ct**4*tt*sy.sin(w + wr)**2 + 14*ct**4*tt*sy.cos(w + wr)**2) + 3*w*(16*cr**4*tr**3*sy.sin(2*wr) - 28*cr**4*tr*sy.sin(2*wr) - 12*cr**3*ct*tr**3*sy.sin(2*wr) - 4*cr**3*ct*tr**3*sy.sin(2*w + 2*wr) - 36*cr**3*ct*tr**2*tt*sy.sin(2*wr) - 12*cr**3*ct*tr**2*tt*sy.sin(2*w + 2*wr) + 63*cr**3*ct*tr*sy.sin(2*wr) + 21*cr**3*ct*tr*sy.sin(2*w + 2*wr) + 21*cr**3*ct*tt*sy.sin(2*wr) + 7*cr**3*ct*tt*sy.sin(2*w + 2*wr) + 24*cr**2*ct**2*tr**2*tt*sy.sin(2*wr) + 24*cr**2*ct**2*tr**2*tt*sy.sin(2*w + 2*wr) + 24*cr**2*ct**2*tr*tt**2*sy.sin(2*wr) + 24*cr**2*ct**2*tr*tt**2*sy.sin(2*w + 2*wr) - 42*cr**2*ct**2*tr*sy.sin(2*wr) - 42*cr**2*ct**2*tr*sy.sin(2*w + 2*wr) - 42*cr**2*ct**2*tt*sy.sin(2*wr) - 42*cr**2*ct**2*tt*sy.sin(2*w + 2*wr) - 12*cr*ct**3*tr*tt**2*sy.sin(2*wr) - 36*cr*ct**3*tr*tt**2*sy.sin(2*w + 2*wr) + 7*cr*ct**3*tr*sy.sin(2*wr) + 21*cr*ct**3*tr*sy.sin(2*w + 2*wr) - 4*cr*ct**3*tt**3*sy.sin(2*wr) - 12*cr*ct**3*tt**3*sy.sin(2*w + 2*wr) + 21*cr*ct**3*tt*sy.sin(2*wr) + 63*cr*ct**3*tt*sy.sin(2*w + 2*wr) + 16*ct**4*tt**3*sy.sin(2*w + 2*wr) - 28*ct**4*tt*sy.sin(2*w + 2*wr)))/w**5)/768, (w > -sy.oo) & (w < sy.oo) & (w != 0)), (b_2*p*(-640*cr**2*tr*xcg*zcg - 320*cr*ct*tr*xcg*zcg - 320*cr*ct*tt*xcg*zcg - 640*ct**2*tt*xcg*zcg + (-16*cr**4*tr**3 + 28*cr**4*tr - 4*cr**3*ct*tr**3 - 12*cr**3*ct*tr**2*tt + 21*cr**3*ct*tr + 7*cr**3*ct*tt - 8*cr**2*ct**2*tr**2*tt - 8*cr**2*ct**2*tr*tt**2 + 14*cr**2*ct**2*tr + 14*cr**2*ct**2*tt - 12*cr*ct**3*tr*tt**2 + 7*cr*ct**3*tr - 4*cr*ct**3*tt**3 + 21*cr*ct**3*tt - 16*ct**4*tt**3 + 28*ct**4*tt)*sy.sin(2*wr))/1920, True)) 
# Iyz_cpy = sy.Piecewise((b_2*p*(3*b_2*(-24*cr**3*tr*sy.cos(wr) + 24*cr**3*tr*sy.cos(w + wr) + 48*cr**2*ct*tr*sy.cos(wr) - 48*cr**2*ct*tr*sy.cos(w + wr) + 24*cr**2*ct*tt*sy.cos(wr) - 24*cr**2*ct*tt*sy.cos(w + wr) - 24*cr*ct**2*tr*sy.cos(wr) + 24*cr*ct**2*tr*sy.cos(w + wr) - 48*cr*ct**2*tt*sy.cos(wr) + 48*cr*ct**2*tt*sy.cos(w + wr) - ct**3*tt*w**4*sy.cos(w + wr) + 24*ct**3*tt*sy.cos(wr) - 24*ct**3*tt*sy.cos(w + wr) + w**3*(-cr**3*tr*sy.sin(wr) - cr*ct**2*tr*sy.sin(w + wr) - 2*cr*ct**2*tt*sy.sin(w + wr) + 4*ct**3*tt*sy.sin(w + wr)) + 2*w**2*(3*cr**3*tr*sy.cos(wr) - 2*cr**2*ct*tr*sy.cos(wr) + 2*cr**2*ct*tr*sy.cos(w + wr) - cr**2*ct*tt*sy.cos(wr) + cr**2*ct*tt*sy.cos(w + wr) - 3*cr*ct**2*tr*sy.cos(w + wr) - 6*cr*ct**2*tt*sy.cos(w + wr) + 6*ct**3*tt*sy.cos(w + wr)) + 6*w*(3*cr**3*tr*sy.sin(wr) + cr**3*tr*sy.sin(w + wr) - 4*cr**2*ct*tr*sy.sin(wr) - 4*cr**2*ct*tr*sy.sin(w + wr) - 2*cr**2*ct*tt*sy.sin(wr) - 2*cr**2*ct*tt*sy.sin(w + wr) + cr*ct**2*tr*sy.sin(wr) + 3*cr*ct**2*tr*sy.sin(w + wr) + 2*cr*ct**2*tt*sy.sin(wr) + 6*cr*ct**2*tt*sy.sin(w + wr) - 4*ct**3*tt*sy.sin(w + wr)))/w**5 - 4*cr**2*tr*ycg*zcg - 2*cr*ct*tr*ycg*zcg - 2*cr*ct*tt*ycg*zcg - 4*ct**2*tt*ycg*zcg)/12, (w > -sy.oo) & (w < sy.oo) & (w != 0)), (b_2*p*(b_2*(3*cr**3*tr + 4*cr**2*ct*tr + 2*cr**2*ct*tt + 3*cr*ct**2*tr + 6*cr*ct**2*tt + 12*ct**3*tt)*sy.sin(wr) - 80*cr**2*tr*ycg*zcg - 40*cr*ct*tr*ycg*zcg - 40*cr*ct*tt*ycg*zcg - 80*ct**2*tt*ycg*zcg)/240, True)) 

# w = 0.0

# Ixx_w0 = sy.Piecewise((b_2*p*(-1280*cr**2*tr*ycg**2 - 1280*cr**2*tr*zcg**2 - 640*cr*ct*tr*ycg**2 - 640*cr*ct*tr*zcg**2 - 640*cr*ct*tt*ycg**2 - 640*cr*ct*tt*zcg**2 - 1280*ct**2*tt*ycg**2 - 1280*ct**2*tt*zcg**2 + (-120*cr**4*tr**3*sy.sin(2*wr) + 120*cr**4*tr**3*sy.sin(2*w + 2*wr) + 210*cr**4*tr*sy.sin(2*wr) - 210*cr**4*tr*sy.sin(2*w + 2*wr) + 120*cr**3*ct*tr**3*sy.sin(2*wr) - 120*cr**3*ct*tr**3*sy.sin(2*w + 2*wr) + 360*cr**3*ct*tr**2*tt*sy.sin(2*wr) - 360*cr**3*ct*tr**2*tt*sy.sin(2*w + 2*wr) - 630*cr**3*ct*tr*sy.sin(2*wr) + 630*cr**3*ct*tr*sy.sin(2*w + 2*wr) - 210*cr**3*ct*tt*sy.sin(2*wr) + 210*cr**3*ct*tt*sy.sin(2*w + 2*wr) - 360*cr**2*ct**2*tr**2*tt*sy.sin(2*wr) + 360*cr**2*ct**2*tr**2*tt*sy.sin(2*w + 2*wr) - 360*cr**2*ct**2*tr*tt**2*sy.sin(2*wr) + 360*cr**2*ct**2*tr*tt**2*sy.sin(2*w + 2*wr) + 630*cr**2*ct**2*tr*sy.sin(2*wr) - 630*cr**2*ct**2*tr*sy.sin(2*w + 2*wr) + 630*cr**2*ct**2*tt*sy.sin(2*wr) - 630*cr**2*ct**2*tt*sy.sin(2*w + 2*wr) + 360*cr*ct**3*tr*tt**2*sy.sin(2*wr) - 360*cr*ct**3*tr*tt**2*sy.sin(2*w + 2*wr) - 210*cr*ct**3*tr*sy.sin(2*wr) + 210*cr*ct**3*tr*sy.sin(2*w + 2*wr) + 120*cr*ct**3*tt**3*sy.sin(2*wr) - 120*cr*ct**3*tt**3*sy.sin(2*w + 2*wr) - 630*cr*ct**3*tt*sy.sin(2*wr) + 630*cr*ct**3*tt*sy.sin(2*w + 2*wr) - 120*ct**4*tt**3*sy.sin(2*wr) + 120*ct**4*tt**3*sy.sin(2*w + 2*wr) + 210*ct**4*tt*sy.sin(2*wr) - 210*ct**4*tt*sy.sin(2*w + 2*wr) + 2*w**5*(64*b_2**2*cr**2*tr + 96*b_2**2*cr*ct*tr + 96*b_2**2*cr*ct*tt + 384*b_2**2*ct**2*tt + 16*cr**4*tr**3 + 28*cr**4*tr + 4*cr**3*ct*tr**3 + 12*cr**3*ct*tr**2*tt + 21*cr**3*ct*tr + 7*cr**3*ct*tt + 8*cr**2*ct**2*tr**2*tt + 8*cr**2*ct**2*tr*tt**2 + 14*cr**2*ct**2*tr + 14*cr**2*ct**2*tt + 12*cr*ct**3*tr*tt**2 + 7*cr*ct**3*tr + 4*cr*ct**3*tt**3 + 21*cr*ct**3*tt + 16*ct**4*tt**3 + 28*ct**4*tt) + 20*w**4*(-4*cr**4*tr**3*sy.sin(2*wr) + 7*cr**4*tr*sy.sin(2*wr) + 4*ct**4*tt**3*sy.sin(2*w + 2*wr) - 7*ct**4*tt*sy.sin(2*w + 2*wr)) + 10*w**3*(16*cr**4*tr**3*sy.cos(2*wr) - 28*cr**4*tr*sy.cos(2*wr) - 4*cr**3*ct*tr**3*sy.cos(2*wr) - 12*cr**3*ct*tr**2*tt*sy.cos(2*wr) + 21*cr**3*ct*tr*sy.cos(2*wr) + 7*cr**3*ct*tt*sy.cos(2*wr) - 12*cr*ct**3*tr*tt**2*sy.cos(2*w + 2*wr) + 7*cr*ct**3*tr*sy.cos(2*w + 2*wr) - 4*cr*ct**3*tt**3*sy.cos(2*w + 2*wr) + 21*cr*ct**3*tt*sy.cos(2*w + 2*wr) + 16*ct**4*tt**3*sy.cos(2*w + 2*wr) - 28*ct**4*tt*sy.cos(2*w + 2*wr)) + 30*w**2*(8*cr**4*tr**3*sy.sin(2*wr) - 14*cr**4*tr*sy.sin(2*wr) - 4*cr**3*ct*tr**3*sy.sin(2*wr) - 12*cr**3*ct*tr**2*tt*sy.sin(2*wr) + 21*cr**3*ct*tr*sy.sin(2*wr) + 7*cr**3*ct*tt*sy.sin(2*wr) + 4*cr**2*ct**2*tr**2*tt*sy.sin(2*wr) - 4*cr**2*ct**2*tr**2*tt*sy.sin(2*w + 2*wr) + 4*cr**2*ct**2*tr*tt**2*sy.sin(2*wr) - 4*cr**2*ct**2*tr*tt**2*sy.sin(2*w + 2*wr) - 7*cr**2*ct**2*tr*sy.sin(2*wr) + 7*cr**2*ct**2*tr*sy.sin(2*w + 2*wr) - 7*cr**2*ct**2*tt*sy.sin(2*wr) + 7*cr**2*ct**2*tt*sy.sin(2*w + 2*wr) + 12*cr*ct**3*tr*tt**2*sy.sin(2*w + 2*wr) - 7*cr*ct**3*tr*sy.sin(2*w + 2*wr) + 4*cr*ct**3*tt**3*sy.sin(2*w + 2*wr) - 21*cr*ct**3*tt*sy.sin(2*w + 2*wr) - 8*ct**4*tt**3*sy.sin(2*w + 2*wr) + 14*ct**4*tt*sy.sin(2*w + 2*wr)) + 15*w*(-16*cr**4*tr**3*sy.cos(2*wr) + 28*cr**4*tr*sy.cos(2*wr) + 12*cr**3*ct*tr**3*sy.cos(2*wr) + 4*cr**3*ct*tr**3*sy.cos(2*w + 2*wr) + 36*cr**3*ct*tr**2*tt*sy.cos(2*wr) + 12*cr**3*ct*tr**2*tt*sy.cos(2*w + 2*wr) - 63*cr**3*ct*tr*sy.cos(2*wr) - 21*cr**3*ct*tr*sy.cos(2*w + 2*wr) - 21*cr**3*ct*tt*sy.cos(2*wr) - 7*cr**3*ct*tt*sy.cos(2*w + 2*wr) - 24*cr**2*ct**2*tr**2*tt*sy.cos(2*wr) - 24*cr**2*ct**2*tr**2*tt*sy.cos(2*w + 2*wr) - 24*cr**2*ct**2*tr*tt**2*sy.cos(2*wr) - 24*cr**2*ct**2*tr*tt**2*sy.cos(2*w + 2*wr) + 42*cr**2*ct**2*tr*sy.cos(2*wr) + 42*cr**2*ct**2*tr*sy.cos(2*w + 2*wr) + 42*cr**2*ct**2*tt*sy.cos(2*wr) + 42*cr**2*ct**2*tt*sy.cos(2*w + 2*wr) + 12*cr*ct**3*tr*tt**2*sy.cos(2*wr) + 36*cr*ct**3*tr*tt**2*sy.cos(2*w + 2*wr) - 7*cr*ct**3*tr*sy.cos(2*wr) - 21*cr*ct**3*tr*sy.cos(2*w + 2*wr) + 4*cr*ct**3*tt**3*sy.cos(2*wr) + 12*cr*ct**3*tt**3*sy.cos(2*w + 2*wr) - 21*cr*ct**3*tt*sy.cos(2*wr) - 63*cr*ct**3*tt*sy.cos(2*w + 2*wr) - 16*ct**4*tt**3*sy.cos(2*w + 2*wr) + 28*ct**4*tt*sy.cos(2*w + 2*wr)))/w**5)/3840, (w > -sy.oo) & (w < sy.oo) & (w != 0)), (b_2*p*(32*b_2**2*cr**2*tr + 48*b_2**2*cr*ct*tr + 48*b_2**2*cr*ct*tt + 192*b_2**2*ct**2*tt + 16*cr**4*tr**3*sy.cos(wr)**2 - 28*cr**4*tr*sy.cos(wr)**2 + 28*cr**4*tr + 4*cr**3*ct*tr**3*sy.cos(wr)**2 + 12*cr**3*ct*tr**2*tt*sy.cos(wr)**2 - 21*cr**3*ct*tr*sy.cos(wr)**2 + 21*cr**3*ct*tr - 7*cr**3*ct*tt*sy.cos(wr)**2 + 7*cr**3*ct*tt + 8*cr**2*ct**2*tr**2*tt*sy.cos(wr)**2 + 8*cr**2*ct**2*tr*tt**2*sy.cos(wr)**2 - 14*cr**2*ct**2*tr*sy.cos(wr)**2 + 14*cr**2*ct**2*tr - 14*cr**2*ct**2*tt*sy.cos(wr)**2 + 14*cr**2*ct**2*tt - 320*cr**2*tr*ycg**2 - 320*cr**2*tr*zcg**2 + 12*cr*ct**3*tr*tt**2*sy.cos(wr)**2 - 7*cr*ct**3*tr*sy.cos(wr)**2 + 7*cr*ct**3*tr + 4*cr*ct**3*tt**3*sy.cos(wr)**2 - 21*cr*ct**3*tt*sy.cos(wr)**2 + 21*cr*ct**3*tt - 160*cr*ct*tr*ycg**2 - 160*cr*ct*tr*zcg**2 - 160*cr*ct*tt*ycg**2 - 160*cr*ct*tt*zcg**2 + 16*ct**4*tt**3*sy.cos(wr)**2 - 28*ct**4*tt*sy.cos(wr)**2 + 28*ct**4*tt - 320*ct**2*tt*ycg**2 - 320*ct**2*tt*zcg**2)/960, True))
# Iyy_w0 = b_2*p*(32*b_2**2*cr**2*tr*tanL**2 + 48*b_2**2*cr*ct*tr*tanL**2 + 48*b_2**2*cr*ct*tt*tanL**2 + 192*b_2**2*ct**2*tt*tanL**2 + 24*b_2*cr**3*tr*tanL + 32*b_2*cr**2*ct*tr*tanL + 16*b_2*cr**2*ct*tt*tanL + 24*b_2*cr*ct**2*tr*tanL + 48*b_2*cr*ct**2*tt*tanL + 96*b_2*ct**3*tt*tanL + 16*cr**4*tr**3 + 28*cr**4*tr + 4*cr**3*ct*tr**3 + 12*cr**3*ct*tr**2*tt + 21*cr**3*ct*tr + 7*cr**3*ct*tt + 8*cr**2*ct**2*tr**2*tt + 8*cr**2*ct**2*tr*tt**2 + 14*cr**2*ct**2*tr + 14*cr**2*ct**2*tt - 320*cr**2*tr*xcg**2 - 320*cr**2*tr*zcg**2 + 12*cr*ct**3*tr*tt**2 + 7*cr*ct**3*tr + 4*cr*ct**3*tt**3 + 21*cr*ct**3*tt - 160*cr*ct*tr*xcg**2 - 160*cr*ct*tr*zcg**2 - 160*cr*ct*tt*xcg**2 - 160*cr*ct*tt*zcg**2 + 16*ct**4*tt**3 + 28*ct**4*tt - 320*ct**2*tt*xcg**2 - 320*ct**2*tt*zcg**2)/960 
# Izz_w0 = sy.Piecewise((b_2*p*(128*b_2**2*cr**2*tr*tanL**2 + 192*b_2**2*cr*ct*tr*tanL**2 + 192*b_2**2*cr*ct*tt*tanL**2 + 768*b_2**2*ct**2*tt*tanL**2 + 96*b_2*cr**3*tr*tanL + 128*b_2*cr**2*ct*tr*tanL + 64*b_2*cr**2*ct*tt*tanL + 96*b_2*cr*ct**2*tr*tanL + 192*b_2*cr*ct**2*tt*tanL + 384*b_2*ct**3*tt*tanL - 1280*cr**2*tr*xcg**2 - 1280*cr**2*tr*ycg**2 - 640*cr*ct*tr*xcg**2 - 640*cr*ct*tr*ycg**2 - 640*cr*ct*tt*xcg**2 - 640*cr*ct*tt*ycg**2 - 1280*ct**2*tt*xcg**2 - 1280*ct**2*tt*ycg**2 + (120*cr**4*tr**3*sy.sin(2*wr) - 120*cr**4*tr**3*sy.sin(2*w + 2*wr) - 210*cr**4*tr*sy.sin(2*wr) + 210*cr**4*tr*sy.sin(2*w + 2*wr) - 120*cr**3*ct*tr**3*sy.sin(2*wr) + 120*cr**3*ct*tr**3*sy.sin(2*w + 2*wr) - 360*cr**3*ct*tr**2*tt*sy.sin(2*wr) + 360*cr**3*ct*tr**2*tt*sy.sin(2*w + 2*wr) + 630*cr**3*ct*tr*sy.sin(2*wr) - 630*cr**3*ct*tr*sy.sin(2*w + 2*wr) + 210*cr**3*ct*tt*sy.sin(2*wr) - 210*cr**3*ct*tt*sy.sin(2*w + 2*wr) + 360*cr**2*ct**2*tr**2*tt*sy.sin(2*wr) - 360*cr**2*ct**2*tr**2*tt*sy.sin(2*w + 2*wr) + 360*cr**2*ct**2*tr*tt**2*sy.sin(2*wr) - 360*cr**2*ct**2*tr*tt**2*sy.sin(2*w + 2*wr) - 630*cr**2*ct**2*tr*sy.sin(2*wr) + 630*cr**2*ct**2*tr*sy.sin(2*w + 2*wr) - 630*cr**2*ct**2*tt*sy.sin(2*wr) + 630*cr**2*ct**2*tt*sy.sin(2*w + 2*wr) - 360*cr*ct**3*tr*tt**2*sy.sin(2*wr) + 360*cr*ct**3*tr*tt**2*sy.sin(2*w + 2*wr) + 210*cr*ct**3*tr*sy.sin(2*wr) - 210*cr*ct**3*tr*sy.sin(2*w + 2*wr) - 120*cr*ct**3*tt**3*sy.sin(2*wr) + 120*cr*ct**3*tt**3*sy.sin(2*w + 2*wr) + 630*cr*ct**3*tt*sy.sin(2*wr) - 630*cr*ct**3*tt*sy.sin(2*w + 2*wr) + 120*ct**4*tt**3*sy.sin(2*wr) - 120*ct**4*tt**3*sy.sin(2*w + 2*wr) - 210*ct**4*tt*sy.sin(2*wr) + 210*ct**4*tt*sy.sin(2*w + 2*wr) + 2*w**5*(64*b_2**2*cr**2*tr + 96*b_2**2*cr*ct*tr + 96*b_2**2*cr*ct*tt + 384*b_2**2*ct**2*tt + 16*cr**4*tr**3 + 28*cr**4*tr + 4*cr**3*ct*tr**3 + 12*cr**3*ct*tr**2*tt + 21*cr**3*ct*tr + 7*cr**3*ct*tt + 8*cr**2*ct**2*tr**2*tt + 8*cr**2*ct**2*tr*tt**2 + 14*cr**2*ct**2*tr + 14*cr**2*ct**2*tt + 12*cr*ct**3*tr*tt**2 + 7*cr*ct**3*tr + 4*cr*ct**3*tt**3 + 21*cr*ct**3*tt + 16*ct**4*tt**3 + 28*ct**4*tt) + 20*w**4*(4*cr**4*tr**3*sy.sin(2*wr) - 7*cr**4*tr*sy.sin(2*wr) - 4*ct**4*tt**3*sy.sin(2*w + 2*wr) + 7*ct**4*tt*sy.sin(2*w + 2*wr)) + 10*w**3*(-16*cr**4*tr**3*sy.cos(2*wr) + 28*cr**4*tr*sy.cos(2*wr) + 4*cr**3*ct*tr**3*sy.cos(2*wr) + 12*cr**3*ct*tr**2*tt*sy.cos(2*wr) - 21*cr**3*ct*tr*sy.cos(2*wr) - 7*cr**3*ct*tt*sy.cos(2*wr) + 12*cr*ct**3*tr*tt**2*sy.cos(2*w + 2*wr) - 7*cr*ct**3*tr*sy.cos(2*w + 2*wr) + 4*cr*ct**3*tt**3*sy.cos(2*w + 2*wr) - 21*cr*ct**3*tt*sy.cos(2*w + 2*wr) - 16*ct**4*tt**3*sy.cos(2*w + 2*wr) + 28*ct**4*tt*sy.cos(2*w + 2*wr)) + 30*w**2*(-8*cr**4*tr**3*sy.sin(2*wr) + 14*cr**4*tr*sy.sin(2*wr) + 4*cr**3*ct*tr**3*sy.sin(2*wr) + 12*cr**3*ct*tr**2*tt*sy.sin(2*wr) - 21*cr**3*ct*tr*sy.sin(2*wr) - 7*cr**3*ct*tt*sy.sin(2*wr) - 4*cr**2*ct**2*tr**2*tt*sy.sin(2*wr) + 4*cr**2*ct**2*tr**2*tt*sy.sin(2*w + 2*wr) - 4*cr**2*ct**2*tr*tt**2*sy.sin(2*wr) + 4*cr**2*ct**2*tr*tt**2*sy.sin(2*w + 2*wr) + 7*cr**2*ct**2*tr*sy.sin(2*wr) - 7*cr**2*ct**2*tr*sy.sin(2*w + 2*wr) + 7*cr**2*ct**2*tt*sy.sin(2*wr) - 7*cr**2*ct**2*tt*sy.sin(2*w + 2*wr) - 12*cr*ct**3*tr*tt**2*sy.sin(2*w + 2*wr) + 7*cr*ct**3*tr*sy.sin(2*w + 2*wr) - 4*cr*ct**3*tt**3*sy.sin(2*w + 2*wr) + 21*cr*ct**3*tt*sy.sin(2*w + 2*wr) + 8*ct**4*tt**3*sy.sin(2*w + 2*wr) - 14*ct**4*tt*sy.sin(2*w + 2*wr)) + 15*w*(16*cr**4*tr**3*sy.cos(2*wr) - 28*cr**4*tr*sy.cos(2*wr) - 12*cr**3*ct*tr**3*sy.cos(2*wr) - 4*cr**3*ct*tr**3*sy.cos(2*w + 2*wr) - 36*cr**3*ct*tr**2*tt*sy.cos(2*wr) - 12*cr**3*ct*tr**2*tt*sy.cos(2*w + 2*wr) + 63*cr**3*ct*tr*sy.cos(2*wr) + 21*cr**3*ct*tr*sy.cos(2*w + 2*wr) + 21*cr**3*ct*tt*sy.cos(2*wr) + 7*cr**3*ct*tt*sy.cos(2*w + 2*wr) + 24*cr**2*ct**2*tr**2*tt*sy.cos(2*wr) + 24*cr**2*ct**2*tr**2*tt*sy.cos(2*w + 2*wr) + 24*cr**2*ct**2*tr*tt**2*sy.cos(2*wr) + 24*cr**2*ct**2*tr*tt**2*sy.cos(2*w + 2*wr) - 42*cr**2*ct**2*tr*sy.cos(2*wr) - 42*cr**2*ct**2*tr*sy.cos(2*w + 2*wr) - 42*cr**2*ct**2*tt*sy.cos(2*wr) - 42*cr**2*ct**2*tt*sy.cos(2*w + 2*wr) - 12*cr*ct**3*tr*tt**2*sy.cos(2*wr) - 36*cr*ct**3*tr*tt**2*sy.cos(2*w + 2*wr) + 7*cr*ct**3*tr*sy.cos(2*wr) + 21*cr*ct**3*tr*sy.cos(2*w + 2*wr) - 4*cr*ct**3*tt**3*sy.cos(2*wr) - 12*cr*ct**3*tt**3*sy.cos(2*w + 2*wr) + 21*cr*ct**3*tt*sy.cos(2*wr) + 63*cr*ct**3*tt*sy.cos(2*w + 2*wr) + 16*ct**4*tt**3*sy.cos(2*w + 2*wr) - 28*ct**4*tt*sy.cos(2*w + 2*wr)))/w**5)/3840, (w > -sy.oo) & (w < sy.oo) & (w != 0)), (b_2*p*(32*b_2**2*cr**2*tr/sy.cos(Lr)**2 + 48*b_2**2*cr*ct*tr/sy.cos(Lr)**2 + 48*b_2**2*cr*ct*tt/sy.cos(Lr)**2 + 192*b_2**2*ct**2*tt/sy.cos(Lr)**2 + 24*b_2*cr**3*tr*tanL + 32*b_2*cr**2*ct*tr*tanL + 16*b_2*cr**2*ct*tt*tanL + 24*b_2*cr*ct**2*tr*tanL + 48*b_2*cr*ct**2*tt*tanL + 96*b_2*ct**3*tt*tanL + 16*cr**4*tr**3*sy.sin(wr)**2 - 28*cr**4*tr*sy.sin(wr)**2 + 28*cr**4*tr + 4*cr**3*ct*tr**3*sy.sin(wr)**2 + 12*cr**3*ct*tr**2*tt*sy.sin(wr)**2 - 21*cr**3*ct*tr*sy.sin(wr)**2 + 21*cr**3*ct*tr - 7*cr**3*ct*tt*sy.sin(wr)**2 + 7*cr**3*ct*tt + 8*cr**2*ct**2*tr**2*tt*sy.sin(wr)**2 + 8*cr**2*ct**2*tr*tt**2*sy.sin(wr)**2 - 14*cr**2*ct**2*tr*sy.sin(wr)**2 + 14*cr**2*ct**2*tr - 14*cr**2*ct**2*tt*sy.sin(wr)**2 + 14*cr**2*ct**2*tt - 320*cr**2*tr*xcg**2 - 320*cr**2*tr*ycg**2 + 12*cr*ct**3*tr*tt**2*sy.sin(wr)**2 - 7*cr*ct**3*tr*sy.sin(wr)**2 + 7*cr*ct**3*tr + 4*cr*ct**3*tt**3*sy.sin(wr)**2 - 21*cr*ct**3*tt*sy.sin(wr)**2 + 21*cr*ct**3*tt - 160*cr*ct*tr*xcg**2 - 160*cr*ct*tr*ycg**2 - 160*cr*ct*tt*xcg**2 - 160*cr*ct*tt*ycg**2 + 16*ct**4*tt**3*sy.sin(wr)**2 - 28*ct**4*tt*sy.sin(wr)**2 + 28*ct**4*tt - 320*ct**2*tt*xcg**2 - 320*ct**2*tt*ycg**2)/960, True)) 
# Ixy_w0 = sy.Piecewise((b_2*p*(-2*b_2**2*cr**2*tr*tanL - 3*b_2**2*cr*ct*tr*tanL - 3*b_2**2*cr*ct*tt*tanL - 12*b_2**2*ct**2*tt*tanL + 15*b_2*(-24*cr**3*tr*sy.sin(wr) + 24*cr**3*tr*sy.sin(w + wr) + 48*cr**2*ct*tr*sy.sin(wr) - 48*cr**2*ct*tr*sy.sin(w + wr) + 24*cr**2*ct*tt*sy.sin(wr) - 24*cr**2*ct*tt*sy.sin(w + wr) - 24*cr*ct**2*tr*sy.sin(wr) + 24*cr*ct**2*tr*sy.sin(w + wr) - 48*cr*ct**2*tt*sy.sin(wr) + 48*cr*ct**2*tt*sy.sin(w + wr) - ct**3*tt*w**4*sy.sin(w + wr) + 24*ct**3*tt*sy.sin(wr) - 24*ct**3*tt*sy.sin(w + wr) + w**3*(cr**3*tr*sy.cos(wr) + cr*ct**2*tr*sy.cos(w + wr) + 2*cr*ct**2*tt*sy.cos(w + wr) - 4*ct**3*tt*sy.cos(w + wr)) + 2*w**2*(3*cr**3*tr*sy.sin(wr) - 2*cr**2*ct*tr*sy.sin(wr) + 2*cr**2*ct*tr*sy.sin(w + wr) - cr**2*ct*tt*sy.sin(wr) + cr**2*ct*tt*sy.sin(w + wr) - 3*cr*ct**2*tr*sy.sin(w + wr) - 6*cr*ct**2*tt*sy.sin(w + wr) + 6*ct**3*tt*sy.sin(w + wr)) + 6*w*(-3*cr**3*tr*sy.cos(wr) - cr**3*tr*sy.cos(w + wr) + 4*cr**2*ct*tr*sy.cos(wr) + 4*cr**2*ct*tr*sy.cos(w + wr) + 2*cr**2*ct*tt*sy.cos(wr) + 2*cr**2*ct*tt*sy.cos(w + wr) - cr*ct**2*tr*sy.cos(wr) - 3*cr*ct**2*tr*sy.cos(w + wr) - 2*cr*ct**2*tt*sy.cos(wr) - 6*cr*ct**2*tt*sy.cos(w + wr) + 4*ct**3*tt*sy.cos(w + wr)))/w**5 - 20*cr**2*tr*xcg*ycg - 10*cr*ct*tr*xcg*ycg - 10*cr*ct*tt*xcg*ycg - 20*ct**2*tt*xcg*ycg)/60, (w > -sy.oo) & (w < sy.oo) & (w != 0)), (-b_2*p*(8*b_2**2*cr**2*tr*tanL + 12*b_2**2*cr*ct*tr*tanL + 12*b_2**2*cr*ct*tt*tanL + 48*b_2**2*ct**2*tt*tanL + b_2*(3*cr**3*tr + 4*cr**2*ct*tr + 2*cr**2*ct*tt + 3*cr*ct**2*tr + 6*cr*ct**2*tt + 12*ct**3*tt)*sy.cos(wr) + 80*cr**2*tr*xcg*ycg + 40*cr*ct*tr*xcg*ycg + 40*cr*ct*tt*xcg*ycg + 80*ct**2*tt*xcg*ycg)/240, True)) 
# Ixz_w0 = sy.Piecewise((b_2*p*(-256*cr**2*tr*xcg*zcg - 128*cr*ct*tr*xcg*zcg - 128*cr*ct*tt*xcg*zcg - 256*ct**2*tt*xcg*zcg + (48*cr**4*tr**3*sy.sin(wr)**2 - 48*cr**4*tr**3*sy.sin(w + wr)**2 - 84*cr**4*tr*sy.sin(wr)**2 + 84*cr**4*tr*sy.sin(w + wr)**2 - 48*cr**3*ct*tr**3*sy.sin(wr)**2 + 48*cr**3*ct*tr**3*sy.sin(w + wr)**2 - 144*cr**3*ct*tr**2*tt*sy.sin(wr)**2 + 144*cr**3*ct*tr**2*tt*sy.sin(w + wr)**2 + 252*cr**3*ct*tr*sy.sin(wr)**2 - 252*cr**3*ct*tr*sy.sin(w + wr)**2 + 84*cr**3*ct*tt*sy.sin(wr)**2 - 84*cr**3*ct*tt*sy.sin(w + wr)**2 + 144*cr**2*ct**2*tr**2*tt*sy.sin(wr)**2 - 144*cr**2*ct**2*tr**2*tt*sy.sin(w + wr)**2 + 144*cr**2*ct**2*tr*tt**2*sy.sin(wr)**2 - 144*cr**2*ct**2*tr*tt**2*sy.sin(w + wr)**2 - 252*cr**2*ct**2*tr*sy.sin(wr)**2 + 252*cr**2*ct**2*tr*sy.sin(w + wr)**2 - 252*cr**2*ct**2*tt*sy.sin(wr)**2 + 252*cr**2*ct**2*tt*sy.sin(w + wr)**2 - 144*cr*ct**3*tr*tt**2*sy.sin(wr)**2 + 144*cr*ct**3*tr*tt**2*sy.sin(w + wr)**2 + 84*cr*ct**3*tr*sy.sin(wr)**2 - 84*cr*ct**3*tr*sy.sin(w + wr)**2 - 48*cr*ct**3*tt**3*sy.sin(wr)**2 + 48*cr*ct**3*tt**3*sy.sin(w + wr)**2 + 252*cr*ct**3*tt*sy.sin(wr)**2 - 252*cr*ct**3*tt*sy.sin(w + wr)**2 + 48*ct**4*tt**3*sy.sin(wr)**2 - 48*ct**4*tt**3*sy.sin(w + wr)**2 - 84*ct**4*tt*sy.sin(wr)**2 + 84*ct**4*tt*sy.sin(w + wr)**2 + 4*w**4*(8*cr**4*tr**3*sy.sin(wr)**2 - 4*cr**4*tr**3*sy.sin(w + wr)**2 - 4*cr**4*tr**3*sy.cos(w + wr)**2 - 14*cr**4*tr*sy.sin(wr)**2 + 7*cr**4*tr*sy.sin(w + wr)**2 + 7*cr**4*tr*sy.cos(w + wr)**2 - 4*ct**4*tt**3*sy.sin(w + wr)**2 + 4*ct**4*tt**3*sy.cos(w + wr)**2 + 7*ct**4*tt*sy.sin(w + wr)**2 - 7*ct**4*tt*sy.cos(w + wr)**2) + 2*w**3*(-16*cr**4*tr**3*sy.sin(2*wr) + 28*cr**4*tr*sy.sin(2*wr) + 4*cr**3*ct*tr**3*sy.sin(2*wr) + 12*cr**3*ct*tr**2*tt*sy.sin(2*wr) - 21*cr**3*ct*tr*sy.sin(2*wr) - 7*cr**3*ct*tt*sy.sin(2*wr) + 12*cr*ct**3*tr*tt**2*sy.sin(2*w + 2*wr) - 7*cr*ct**3*tr*sy.sin(2*w + 2*wr) + 4*cr*ct**3*tt**3*sy.sin(2*w + 2*wr) - 21*cr*ct**3*tt*sy.sin(2*w + 2*wr) - 16*ct**4*tt**3*sy.sin(2*w + 2*wr) + 28*ct**4*tt*sy.sin(2*w + 2*wr)) + 6*w**2*(-16*cr**4*tr**3*sy.sin(wr)**2 + 8*cr**4*tr**3*sy.sin(w + wr)**2 + 8*cr**4*tr**3*sy.cos(w + wr)**2 + 28*cr**4*tr*sy.sin(wr)**2 - 14*cr**4*tr*sy.sin(w + wr)**2 - 14*cr**4*tr*sy.cos(w + wr)**2 + 8*cr**3*ct*tr**3*sy.sin(wr)**2 - 4*cr**3*ct*tr**3*sy.sin(w + wr)**2 - 4*cr**3*ct*tr**3*sy.cos(w + wr)**2 + 24*cr**3*ct*tr**2*tt*sy.sin(wr)**2 - 12*cr**3*ct*tr**2*tt*sy.sin(w + wr)**2 - 12*cr**3*ct*tr**2*tt*sy.cos(w + wr)**2 - 42*cr**3*ct*tr*sy.sin(wr)**2 + 21*cr**3*ct*tr*sy.sin(w + wr)**2 + 21*cr**3*ct*tr*sy.cos(w + wr)**2 - 14*cr**3*ct*tt*sy.sin(wr)**2 + 7*cr**3*ct*tt*sy.sin(w + wr)**2 + 7*cr**3*ct*tt*sy.cos(w + wr)**2 - 8*cr**2*ct**2*tr**2*tt*sy.sin(wr)**2 + 8*cr**2*ct**2*tr**2*tt*sy.sin(w + wr)**2 - 8*cr**2*ct**2*tr*tt**2*sy.sin(wr)**2 + 8*cr**2*ct**2*tr*tt**2*sy.sin(w + wr)**2 + 14*cr**2*ct**2*tr*sy.sin(wr)**2 - 14*cr**2*ct**2*tr*sy.sin(w + wr)**2 + 14*cr**2*ct**2*tt*sy.sin(wr)**2 - 14*cr**2*ct**2*tt*sy.sin(w + wr)**2 - 12*cr*ct**3*tr*tt**2*sy.sin(w + wr)**2 + 12*cr*ct**3*tr*tt**2*sy.cos(w + wr)**2 + 7*cr*ct**3*tr*sy.sin(w + wr)**2 - 7*cr*ct**3*tr*sy.cos(w + wr)**2 - 4*cr*ct**3*tt**3*sy.sin(w + wr)**2 + 4*cr*ct**3*tt**3*sy.cos(w + wr)**2 + 21*cr*ct**3*tt*sy.sin(w + wr)**2 - 21*cr*ct**3*tt*sy.cos(w + wr)**2 + 8*ct**4*tt**3*sy.sin(w + wr)**2 - 8*ct**4*tt**3*sy.cos(w + wr)**2 - 14*ct**4*tt*sy.sin(w + wr)**2 + 14*ct**4*tt*sy.cos(w + wr)**2) + 3*w*(16*cr**4*tr**3*sy.sin(2*wr) - 28*cr**4*tr*sy.sin(2*wr) - 12*cr**3*ct*tr**3*sy.sin(2*wr) - 4*cr**3*ct*tr**3*sy.sin(2*w + 2*wr) - 36*cr**3*ct*tr**2*tt*sy.sin(2*wr) - 12*cr**3*ct*tr**2*tt*sy.sin(2*w + 2*wr) + 63*cr**3*ct*tr*sy.sin(2*wr) + 21*cr**3*ct*tr*sy.sin(2*w + 2*wr) + 21*cr**3*ct*tt*sy.sin(2*wr) + 7*cr**3*ct*tt*sy.sin(2*w + 2*wr) + 24*cr**2*ct**2*tr**2*tt*sy.sin(2*wr) + 24*cr**2*ct**2*tr**2*tt*sy.sin(2*w + 2*wr) + 24*cr**2*ct**2*tr*tt**2*sy.sin(2*wr) + 24*cr**2*ct**2*tr*tt**2*sy.sin(2*w + 2*wr) - 42*cr**2*ct**2*tr*sy.sin(2*wr) - 42*cr**2*ct**2*tr*sy.sin(2*w + 2*wr) - 42*cr**2*ct**2*tt*sy.sin(2*wr) - 42*cr**2*ct**2*tt*sy.sin(2*w + 2*wr) - 12*cr*ct**3*tr*tt**2*sy.sin(2*wr) - 36*cr*ct**3*tr*tt**2*sy.sin(2*w + 2*wr) + 7*cr*ct**3*tr*sy.sin(2*wr) + 21*cr*ct**3*tr*sy.sin(2*w + 2*wr) - 4*cr*ct**3*tt**3*sy.sin(2*wr) - 12*cr*ct**3*tt**3*sy.sin(2*w + 2*wr) + 21*cr*ct**3*tt*sy.sin(2*wr) + 63*cr*ct**3*tt*sy.sin(2*w + 2*wr) + 16*ct**4*tt**3*sy.sin(2*w + 2*wr) - 28*ct**4*tt*sy.sin(2*w + 2*wr)))/w**5)/768, (w > -sy.oo) & (w < sy.oo) & (w != 0)), (b_2*p*(-640*cr**2*tr*xcg*zcg - 320*cr*ct*tr*xcg*zcg - 320*cr*ct*tt*xcg*zcg - 640*ct**2*tt*xcg*zcg + (-16*cr**4*tr**3 + 28*cr**4*tr - 4*cr**3*ct*tr**3 - 12*cr**3*ct*tr**2*tt + 21*cr**3*ct*tr + 7*cr**3*ct*tt - 8*cr**2*ct**2*tr**2*tt - 8*cr**2*ct**2*tr*tt**2 + 14*cr**2*ct**2*tr + 14*cr**2*ct**2*tt - 12*cr*ct**3*tr*tt**2 + 7*cr*ct**3*tr - 4*cr*ct**3*tt**3 + 21*cr*ct**3*tt - 16*ct**4*tt**3 + 28*ct**4*tt)*sy.sin(2*wr))/1920, True)) 
# Iyz_w0 = sy.Piecewise((b_2*p*(3*b_2*(-24*cr**3*tr*sy.cos(wr) + 24*cr**3*tr*sy.cos(w + wr) + 48*cr**2*ct*tr*sy.cos(wr) - 48*cr**2*ct*tr*sy.cos(w + wr) + 24*cr**2*ct*tt*sy.cos(wr) - 24*cr**2*ct*tt*sy.cos(w + wr) - 24*cr*ct**2*tr*sy.cos(wr) + 24*cr*ct**2*tr*sy.cos(w + wr) - 48*cr*ct**2*tt*sy.cos(wr) + 48*cr*ct**2*tt*sy.cos(w + wr) - ct**3*tt*w**4*sy.cos(w + wr) + 24*ct**3*tt*sy.cos(wr) - 24*ct**3*tt*sy.cos(w + wr) + w**3*(-cr**3*tr*sy.sin(wr) - cr*ct**2*tr*sy.sin(w + wr) - 2*cr*ct**2*tt*sy.sin(w + wr) + 4*ct**3*tt*sy.sin(w + wr)) + 2*w**2*(3*cr**3*tr*sy.cos(wr) - 2*cr**2*ct*tr*sy.cos(wr) + 2*cr**2*ct*tr*sy.cos(w + wr) - cr**2*ct*tt*sy.cos(wr) + cr**2*ct*tt*sy.cos(w + wr) - 3*cr*ct**2*tr*sy.cos(w + wr) - 6*cr*ct**2*tt*sy.cos(w + wr) + 6*ct**3*tt*sy.cos(w + wr)) + 6*w*(3*cr**3*tr*sy.sin(wr) + cr**3*tr*sy.sin(w + wr) - 4*cr**2*ct*tr*sy.sin(wr) - 4*cr**2*ct*tr*sy.sin(w + wr) - 2*cr**2*ct*tt*sy.sin(wr) - 2*cr**2*ct*tt*sy.sin(w + wr) + cr*ct**2*tr*sy.sin(wr) + 3*cr*ct**2*tr*sy.sin(w + wr) + 2*cr*ct**2*tt*sy.sin(wr) + 6*cr*ct**2*tt*sy.sin(w + wr) - 4*ct**3*tt*sy.sin(w + wr)))/w**5 - 4*cr**2*tr*ycg*zcg - 2*cr*ct*tr*ycg*zcg - 2*cr*ct*tt*ycg*zcg - 4*ct**2*tt*ycg*zcg)/12, (w > -sy.oo) & (w < sy.oo) & (w != 0)), (b_2*p*(b_2*(3*cr**3*tr + 4*cr**2*ct*tr + 2*cr**2*ct*tt + 3*cr*ct**2*tr + 6*cr*ct**2*tt + 12*ct**3*tt)*sy.sin(wr) - 80*cr**2*tr*ycg*zcg - 40*cr*ct*tr*ycg*zcg - 40*cr*ct*tt*ycg*zcg - 80*ct**2*tt*ycg*zcg)/240, True)) 

# Ixx_w0 = ksubs(Ixx_w0)
# Iyy_w0 = ksubs(Iyy_w0)
# Izz_w0 = ksubs(Izz_w0)
# Ixy_w0 = ksubs(Ixy_w0)
# Ixz_w0 = ksubs(Ixz_w0)
# Iyz_w0 = ksubs(Iyz_w0)

# print("Ixx diff =", simp(Ixx_new - Ixx_cpy))
# print("Iyy diff =", simp(Iyy_new - Iyy_cpy))
# print("Izz diff =", simp(Izz_new - Izz_cpy))
# print("Ixy diff =", simp(Ixy_new - Ixy_cpy))
# print("Ixz diff =", simp(Ixz_new - Ixz_cpy))
# print("Iyz diff =", simp(Iyz_new - Iyz_cpy))

# replacements
# Piecewise ==> sy.Piecewise
# b/2 ==> b_2
# c_r ==> cr
# c_t ==> ct
# tan(Λ_r) ==> tanL
# Λ_r ==> Lr
# t_r ==> tr
# t_t ==> tt
# ω_r ==> wr
# ω ==> w
# x_cg ==> xcg
# y_cg ==> ycg
# z_cg ==> zcg
# sin ==> sy.sin
# cos ==> sy.cos
# oo ==> sy.oo
# Ne ==> fix these to !=

# # output
# Ixx_cpy  = Piecewise((b/2*p*(-1280*c_r**2*t_r*y_cg**2 - 1280*c_r**2*t_r*z_cg**2 - 640*c_r*c_t*t_r*y_cg**2 - 640*c_r*c_t*t_r*z_cg**2 - 640*c_r*c_t*t_t*y_cg**2 - 640*c_r*c_t*t_t*z_cg**2 - 1280*c_t**2*t_t*y_cg**2 - 1280*c_t**2*t_t*z_cg**2 + (-120*c_r**4*t_r**3*sin(2*ω_r) + 120*c_r**4*t_r**3*sin(2*ω + 2*ω_r) + 210*c_r**4*t_r*sin(2*ω_r) - 210*c_r**4*t_r*sin(2*ω + 2*ω_r) + 120*c_r**3*c_t*t_r**3*sin(2*ω_r) - 120*c_r**3*c_t*t_r**3*sin(2*ω + 2*ω_r) + 360*c_r**3*c_t*t_r**2*t_t*sin(2*ω_r) - 360*c_r**3*c_t*t_r**2*t_t*sin(2*ω + 2*ω_r) - 630*c_r**3*c_t*t_r*sin(2*ω_r) + 630*c_r**3*c_t*t_r*sin(2*ω + 2*ω_r) - 210*c_r**3*c_t*t_t*sin(2*ω_r) + 210*c_r**3*c_t*t_t*sin(2*ω + 2*ω_r) - 360*c_r**2*c_t**2*t_r**2*t_t*sin(2*ω_r) + 360*c_r**2*c_t**2*t_r**2*t_t*sin(2*ω + 2*ω_r) - 360*c_r**2*c_t**2*t_r*t_t**2*sin(2*ω_r) + 360*c_r**2*c_t**2*t_r*t_t**2*sin(2*ω + 2*ω_r) + 630*c_r**2*c_t**2*t_r*sin(2*ω_r) - 630*c_r**2*c_t**2*t_r*sin(2*ω + 2*ω_r) + 630*c_r**2*c_t**2*t_t*sin(2*ω_r) - 630*c_r**2*c_t**2*t_t*sin(2*ω + 2*ω_r) + 360*c_r*c_t**3*t_r*t_t**2*sin(2*ω_r) - 360*c_r*c_t**3*t_r*t_t**2*sin(2*ω + 2*ω_r) - 210*c_r*c_t**3*t_r*sin(2*ω_r) + 210*c_r*c_t**3*t_r*sin(2*ω + 2*ω_r) + 120*c_r*c_t**3*t_t**3*sin(2*ω_r) - 120*c_r*c_t**3*t_t**3*sin(2*ω + 2*ω_r) - 630*c_r*c_t**3*t_t*sin(2*ω_r) + 630*c_r*c_t**3*t_t*sin(2*ω + 2*ω_r) - 120*c_t**4*t_t**3*sin(2*ω_r) + 120*c_t**4*t_t**3*sin(2*ω + 2*ω_r) + 210*c_t**4*t_t*sin(2*ω_r) - 210*c_t**4*t_t*sin(2*ω + 2*ω_r) + 2*ω**5*(64*b/2**2*c_r**2*t_r + 96*b/2**2*c_r*c_t*t_r + 96*b/2**2*c_r*c_t*t_t + 384*b/2**2*c_t**2*t_t + 16*c_r**4*t_r**3 + 28*c_r**4*t_r + 4*c_r**3*c_t*t_r**3 + 12*c_r**3*c_t*t_r**2*t_t + 21*c_r**3*c_t*t_r + 7*c_r**3*c_t*t_t + 8*c_r**2*c_t**2*t_r**2*t_t + 8*c_r**2*c_t**2*t_r*t_t**2 + 14*c_r**2*c_t**2*t_r + 14*c_r**2*c_t**2*t_t + 12*c_r*c_t**3*t_r*t_t**2 + 7*c_r*c_t**3*t_r + 4*c_r*c_t**3*t_t**3 + 21*c_r*c_t**3*t_t + 16*c_t**4*t_t**3 + 28*c_t**4*t_t) + 20*ω**4*(-4*c_r**4*t_r**3*sin(2*ω_r) + 7*c_r**4*t_r*sin(2*ω_r) + 4*c_t**4*t_t**3*sin(2*ω + 2*ω_r) - 7*c_t**4*t_t*sin(2*ω + 2*ω_r)) + 10*ω**3*(16*c_r**4*t_r**3*cos(2*ω_r) - 28*c_r**4*t_r*cos(2*ω_r) - 4*c_r**3*c_t*t_r**3*cos(2*ω_r) - 12*c_r**3*c_t*t_r**2*t_t*cos(2*ω_r) + 21*c_r**3*c_t*t_r*cos(2*ω_r) + 7*c_r**3*c_t*t_t*cos(2*ω_r) - 12*c_r*c_t**3*t_r*t_t**2*cos(2*ω + 2*ω_r) + 7*c_r*c_t**3*t_r*cos(2*ω + 2*ω_r) - 4*c_r*c_t**3*t_t**3*cos(2*ω + 2*ω_r) + 21*c_r*c_t**3*t_t*cos(2*ω + 2*ω_r) + 16*c_t**4*t_t**3*cos(2*ω + 2*ω_r) - 28*c_t**4*t_t*cos(2*ω + 2*ω_r)) + 30*ω**2*(8*c_r**4*t_r**3*sin(2*ω_r) - 14*c_r**4*t_r*sin(2*ω_r) - 4*c_r**3*c_t*t_r**3*sin(2*ω_r) - 12*c_r**3*c_t*t_r**2*t_t*sin(2*ω_r) + 21*c_r**3*c_t*t_r*sin(2*ω_r) + 7*c_r**3*c_t*t_t*sin(2*ω_r) + 4*c_r**2*c_t**2*t_r**2*t_t*sin(2*ω_r) - 4*c_r**2*c_t**2*t_r**2*t_t*sin(2*ω + 2*ω_r) + 4*c_r**2*c_t**2*t_r*t_t**2*sin(2*ω_r) - 4*c_r**2*c_t**2*t_r*t_t**2*sin(2*ω + 2*ω_r) - 7*c_r**2*c_t**2*t_r*sin(2*ω_r) + 7*c_r**2*c_t**2*t_r*sin(2*ω + 2*ω_r) - 7*c_r**2*c_t**2*t_t*sin(2*ω_r) + 7*c_r**2*c_t**2*t_t*sin(2*ω + 2*ω_r) + 12*c_r*c_t**3*t_r*t_t**2*sin(2*ω + 2*ω_r) - 7*c_r*c_t**3*t_r*sin(2*ω + 2*ω_r) + 4*c_r*c_t**3*t_t**3*sin(2*ω + 2*ω_r) - 21*c_r*c_t**3*t_t*sin(2*ω + 2*ω_r) - 8*c_t**4*t_t**3*sin(2*ω + 2*ω_r) + 14*c_t**4*t_t*sin(2*ω + 2*ω_r)) + 15*ω*(-16*c_r**4*t_r**3*cos(2*ω_r) + 28*c_r**4*t_r*cos(2*ω_r) + 12*c_r**3*c_t*t_r**3*cos(2*ω_r) + 4*c_r**3*c_t*t_r**3*cos(2*ω + 2*ω_r) + 36*c_r**3*c_t*t_r**2*t_t*cos(2*ω_r) + 12*c_r**3*c_t*t_r**2*t_t*cos(2*ω + 2*ω_r) - 63*c_r**3*c_t*t_r*cos(2*ω_r) - 21*c_r**3*c_t*t_r*cos(2*ω + 2*ω_r) - 21*c_r**3*c_t*t_t*cos(2*ω_r) - 7*c_r**3*c_t*t_t*cos(2*ω + 2*ω_r) - 24*c_r**2*c_t**2*t_r**2*t_t*cos(2*ω_r) - 24*c_r**2*c_t**2*t_r**2*t_t*cos(2*ω + 2*ω_r) - 24*c_r**2*c_t**2*t_r*t_t**2*cos(2*ω_r) - 24*c_r**2*c_t**2*t_r*t_t**2*cos(2*ω + 2*ω_r) + 42*c_r**2*c_t**2*t_r*cos(2*ω_r) + 42*c_r**2*c_t**2*t_r*cos(2*ω + 2*ω_r) + 42*c_r**2*c_t**2*t_t*cos(2*ω_r) + 42*c_r**2*c_t**2*t_t*cos(2*ω + 2*ω_r) + 12*c_r*c_t**3*t_r*t_t**2*cos(2*ω_r) + 36*c_r*c_t**3*t_r*t_t**2*cos(2*ω + 2*ω_r) - 7*c_r*c_t**3*t_r*cos(2*ω_r) - 21*c_r*c_t**3*t_r*cos(2*ω + 2*ω_r) + 4*c_r*c_t**3*t_t**3*cos(2*ω_r) + 12*c_r*c_t**3*t_t**3*cos(2*ω + 2*ω_r) - 21*c_r*c_t**3*t_t*cos(2*ω_r) - 63*c_r*c_t**3*t_t*cos(2*ω + 2*ω_r) - 16*c_t**4*t_t**3*cos(2*ω + 2*ω_r) + 28*c_t**4*t_t*cos(2*ω + 2*ω_r)))/ω**5)/3840, (ω >= -oo) & (ω < oo)), (b/2*p*(32*b/2**2*c_r**2*t_r + 48*b/2**2*c_r*c_t*t_r + 48*b/2**2*c_r*c_t*t_t + 192*b/2**2*c_t**2*t_t + 16*c_r**4*t_r**3*cos(ω_r)**2 - 28*c_r**4*t_r*cos(ω_r)**2 + 28*c_r**4*t_r + 4*c_r**3*c_t*t_r**3*cos(ω_r)**2 + 12*c_r**3*c_t*t_r**2*t_t*cos(ω_r)**2 - 21*c_r**3*c_t*t_r*cos(ω_r)**2 + 21*c_r**3*c_t*t_r - 7*c_r**3*c_t*t_t*cos(ω_r)**2 + 7*c_r**3*c_t*t_t + 8*c_r**2*c_t**2*t_r**2*t_t*cos(ω_r)**2 + 8*c_r**2*c_t**2*t_r*t_t**2*cos(ω_r)**2 - 14*c_r**2*c_t**2*t_r*cos(ω_r)**2 + 14*c_r**2*c_t**2*t_r - 14*c_r**2*c_t**2*t_t*cos(ω_r)**2 + 14*c_r**2*c_t**2*t_t - 320*c_r**2*t_r*y_cg**2 - 320*c_r**2*t_r*z_cg**2 + 12*c_r*c_t**3*t_r*t_t**2*cos(ω_r)**2 - 7*c_r*c_t**3*t_r*cos(ω_r)**2 + 7*c_r*c_t**3*t_r + 4*c_r*c_t**3*t_t**3*cos(ω_r)**2 - 21*c_r*c_t**3*t_t*cos(ω_r)**2 + 21*c_r*c_t**3*t_t - 160*c_r*c_t*t_r*y_cg**2 - 160*c_r*c_t*t_r*z_cg**2 - 160*c_r*c_t*t_t*y_cg**2 - 160*c_r*c_t*t_t*z_cg**2 + 16*c_t**4*t_t**3*cos(ω_r)**2 - 28*c_t**4*t_t*cos(ω_r)**2 + 28*c_t**4*t_t - 320*c_t**2*t_t*y_cg**2 - 320*c_t**2*t_t*z_cg**2)/960, True)) 

# Iyy_cpy  = b/2*p*(32*b/2**2*c_r**2*t_r*tan(Λ_r)**2 + 48*b/2**2*c_r*c_t*t_r*tan(Λ_r)**2 + 48*b/2**2*c_r*c_t*t_t*tan(Λ_r)**2 + 192*b/2**2*c_t**2*t_t*tan(Λ_r)**2 + 24*b/2*c_r**3*t_r*tan(Λ_r) + 32*b/2*c_r**2*c_t*t_r*tan(Λ_r) + 16*b/2*c_r**2*c_t*t_t*tan(Λ_r) + 24*b/2*c_r*c_t**2*t_r*tan(Λ_r) + 48*b/2*c_r*c_t**2*t_t*tan(Λ_r) + 96*b/2*c_t**3*t_t*tan(Λ_r) + 16*c_r**4*t_r**3 + 28*c_r**4*t_r + 4*c_r**3*c_t*t_r**3 + 12*c_r**3*c_t*t_r**2*t_t + 21*c_r**3*c_t*t_r + 7*c_r**3*c_t*t_t + 8*c_r**2*c_t**2*t_r**2*t_t + 8*c_r**2*c_t**2*t_r*t_t**2 + 14*c_r**2*c_t**2*t_r + 14*c_r**2*c_t**2*t_t - 320*c_r**2*t_r*x_cg**2 - 320*c_r**2*t_r*z_cg**2 + 12*c_r*c_t**3*t_r*t_t**2 + 7*c_r*c_t**3*t_r + 4*c_r*c_t**3*t_t**3 + 21*c_r*c_t**3*t_t - 160*c_r*c_t*t_r*x_cg**2 - 160*c_r*c_t*t_r*z_cg**2 - 160*c_r*c_t*t_t*x_cg**2 - 160*c_r*c_t*t_t*z_cg**2 + 16*c_t**4*t_t**3 + 28*c_t**4*t_t - 320*c_t**2*t_t*x_cg**2 - 320*c_t**2*t_t*z_cg**2)/960 

# Izz_cpy  = Piecewise((b/2*p*(128*b/2**2*c_r**2*t_r*tan(Λ_r)**2 + 192*b/2**2*c_r*c_t*t_r*tan(Λ_r)**2 + 192*b/2**2*c_r*c_t*t_t*tan(Λ_r)**2 + 768*b/2**2*c_t**2*t_t*tan(Λ_r)**2 + 96*b/2*c_r**3*t_r*tan(Λ_r) + 128*b/2*c_r**2*c_t*t_r*tan(Λ_r) + 64*b/2*c_r**2*c_t*t_t*tan(Λ_r) + 96*b/2*c_r*c_t**2*t_r*tan(Λ_r) + 192*b/2*c_r*c_t**2*t_t*tan(Λ_r) + 384*b/2*c_t**3*t_t*tan(Λ_r) - 1280*c_r**2*t_r*x_cg**2 - 1280*c_r**2*t_r*y_cg**2 - 640*c_r*c_t*t_r*x_cg**2 - 640*c_r*c_t*t_r*y_cg**2 - 640*c_r*c_t*t_t*x_cg**2 - 640*c_r*c_t*t_t*y_cg**2 - 1280*c_t**2*t_t*x_cg**2 - 1280*c_t**2*t_t*y_cg**2 + (120*c_r**4*t_r**3*sin(2*ω_r) - 120*c_r**4*t_r**3*sin(2*ω + 2*ω_r) - 210*c_r**4*t_r*sin(2*ω_r) + 210*c_r**4*t_r*sin(2*ω + 2*ω_r) - 120*c_r**3*c_t*t_r**3*sin(2*ω_r) + 120*c_r**3*c_t*t_r**3*sin(2*ω + 2*ω_r) - 360*c_r**3*c_t*t_r**2*t_t*sin(2*ω_r) + 360*c_r**3*c_t*t_r**2*t_t*sin(2*ω + 2*ω_r) + 630*c_r**3*c_t*t_r*sin(2*ω_r) - 630*c_r**3*c_t*t_r*sin(2*ω + 2*ω_r) + 210*c_r**3*c_t*t_t*sin(2*ω_r) - 210*c_r**3*c_t*t_t*sin(2*ω + 2*ω_r) + 360*c_r**2*c_t**2*t_r**2*t_t*sin(2*ω_r) - 360*c_r**2*c_t**2*t_r**2*t_t*sin(2*ω + 2*ω_r) + 360*c_r**2*c_t**2*t_r*t_t**2*sin(2*ω_r) - 360*c_r**2*c_t**2*t_r*t_t**2*sin(2*ω + 2*ω_r) - 630*c_r**2*c_t**2*t_r*sin(2*ω_r) + 630*c_r**2*c_t**2*t_r*sin(2*ω + 2*ω_r) - 630*c_r**2*c_t**2*t_t*sin(2*ω_r) + 630*c_r**2*c_t**2*t_t*sin(2*ω + 2*ω_r) - 360*c_r*c_t**3*t_r*t_t**2*sin(2*ω_r) + 360*c_r*c_t**3*t_r*t_t**2*sin(2*ω + 2*ω_r) + 210*c_r*c_t**3*t_r*sin(2*ω_r) - 210*c_r*c_t**3*t_r*sin(2*ω + 2*ω_r) - 120*c_r*c_t**3*t_t**3*sin(2*ω_r) + 120*c_r*c_t**3*t_t**3*sin(2*ω + 2*ω_r) + 630*c_r*c_t**3*t_t*sin(2*ω_r) - 630*c_r*c_t**3*t_t*sin(2*ω + 2*ω_r) + 120*c_t**4*t_t**3*sin(2*ω_r) - 120*c_t**4*t_t**3*sin(2*ω + 2*ω_r) - 210*c_t**4*t_t*sin(2*ω_r) + 210*c_t**4*t_t*sin(2*ω + 2*ω_r) + 2*ω**5*(64*b/2**2*c_r**2*t_r + 96*b/2**2*c_r*c_t*t_r + 96*b/2**2*c_r*c_t*t_t + 384*b/2**2*c_t**2*t_t + 16*c_r**4*t_r**3 + 28*c_r**4*t_r + 4*c_r**3*c_t*t_r**3 + 12*c_r**3*c_t*t_r**2*t_t + 21*c_r**3*c_t*t_r + 7*c_r**3*c_t*t_t + 8*c_r**2*c_t**2*t_r**2*t_t + 8*c_r**2*c_t**2*t_r*t_t**2 + 14*c_r**2*c_t**2*t_r + 14*c_r**2*c_t**2*t_t + 12*c_r*c_t**3*t_r*t_t**2 + 7*c_r*c_t**3*t_r + 4*c_r*c_t**3*t_t**3 + 21*c_r*c_t**3*t_t + 16*c_t**4*t_t**3 + 28*c_t**4*t_t) + 20*ω**4*(4*c_r**4*t_r**3*sin(2*ω_r) - 7*c_r**4*t_r*sin(2*ω_r) - 4*c_t**4*t_t**3*sin(2*ω + 2*ω_r) + 7*c_t**4*t_t*sin(2*ω + 2*ω_r)) + 10*ω**3*(-16*c_r**4*t_r**3*cos(2*ω_r) + 28*c_r**4*t_r*cos(2*ω_r) + 4*c_r**3*c_t*t_r**3*cos(2*ω_r) + 12*c_r**3*c_t*t_r**2*t_t*cos(2*ω_r) - 21*c_r**3*c_t*t_r*cos(2*ω_r) - 7*c_r**3*c_t*t_t*cos(2*ω_r) + 12*c_r*c_t**3*t_r*t_t**2*cos(2*ω + 2*ω_r) - 7*c_r*c_t**3*t_r*cos(2*ω + 2*ω_r) + 4*c_r*c_t**3*t_t**3*cos(2*ω + 2*ω_r) - 21*c_r*c_t**3*t_t*cos(2*ω + 2*ω_r) - 16*c_t**4*t_t**3*cos(2*ω + 2*ω_r) + 28*c_t**4*t_t*cos(2*ω + 2*ω_r)) + 30*ω**2*(-8*c_r**4*t_r**3*sin(2*ω_r) + 14*c_r**4*t_r*sin(2*ω_r) + 4*c_r**3*c_t*t_r**3*sin(2*ω_r) + 12*c_r**3*c_t*t_r**2*t_t*sin(2*ω_r) - 21*c_r**3*c_t*t_r*sin(2*ω_r) - 7*c_r**3*c_t*t_t*sin(2*ω_r) - 4*c_r**2*c_t**2*t_r**2*t_t*sin(2*ω_r) + 4*c_r**2*c_t**2*t_r**2*t_t*sin(2*ω + 2*ω_r) - 4*c_r**2*c_t**2*t_r*t_t**2*sin(2*ω_r) + 4*c_r**2*c_t**2*t_r*t_t**2*sin(2*ω + 2*ω_r) + 7*c_r**2*c_t**2*t_r*sin(2*ω_r) - 7*c_r**2*c_t**2*t_r*sin(2*ω + 2*ω_r) + 7*c_r**2*c_t**2*t_t*sin(2*ω_r) - 7*c_r**2*c_t**2*t_t*sin(2*ω + 2*ω_r) - 12*c_r*c_t**3*t_r*t_t**2*sin(2*ω + 2*ω_r) + 7*c_r*c_t**3*t_r*sin(2*ω + 2*ω_r) - 4*c_r*c_t**3*t_t**3*sin(2*ω + 2*ω_r) + 21*c_r*c_t**3*t_t*sin(2*ω + 2*ω_r) + 8*c_t**4*t_t**3*sin(2*ω + 2*ω_r) - 14*c_t**4*t_t*sin(2*ω + 2*ω_r)) + 15*ω*(16*c_r**4*t_r**3*cos(2*ω_r) - 28*c_r**4*t_r*cos(2*ω_r) - 12*c_r**3*c_t*t_r**3*cos(2*ω_r) - 4*c_r**3*c_t*t_r**3*cos(2*ω + 2*ω_r) - 36*c_r**3*c_t*t_r**2*t_t*cos(2*ω_r) - 12*c_r**3*c_t*t_r**2*t_t*cos(2*ω + 2*ω_r) + 63*c_r**3*c_t*t_r*cos(2*ω_r) + 21*c_r**3*c_t*t_r*cos(2*ω + 2*ω_r) + 21*c_r**3*c_t*t_t*cos(2*ω_r) + 7*c_r**3*c_t*t_t*cos(2*ω + 2*ω_r) + 24*c_r**2*c_t**2*t_r**2*t_t*cos(2*ω_r) + 24*c_r**2*c_t**2*t_r**2*t_t*cos(2*ω + 2*ω_r) + 24*c_r**2*c_t**2*t_r*t_t**2*cos(2*ω_r) + 24*c_r**2*c_t**2*t_r*t_t**2*cos(2*ω + 2*ω_r) - 42*c_r**2*c_t**2*t_r*cos(2*ω_r) - 42*c_r**2*c_t**2*t_r*cos(2*ω + 2*ω_r) - 42*c_r**2*c_t**2*t_t*cos(2*ω_r) - 42*c_r**2*c_t**2*t_t*cos(2*ω + 2*ω_r) - 12*c_r*c_t**3*t_r*t_t**2*cos(2*ω_r) - 36*c_r*c_t**3*t_r*t_t**2*cos(2*ω + 2*ω_r) + 7*c_r*c_t**3*t_r*cos(2*ω_r) + 21*c_r*c_t**3*t_r*cos(2*ω + 2*ω_r) - 4*c_r*c_t**3*t_t**3*cos(2*ω_r) - 12*c_r*c_t**3*t_t**3*cos(2*ω + 2*ω_r) + 21*c_r*c_t**3*t_t*cos(2*ω_r) + 63*c_r*c_t**3*t_t*cos(2*ω + 2*ω_r) + 16*c_t**4*t_t**3*cos(2*ω + 2*ω_r) - 28*c_t**4*t_t*cos(2*ω + 2*ω_r)))/ω**5)/3840, (ω >= -oo) & (ω < oo)), (b/2*p*(32*b/2**2*c_r**2*t_r/cos(Λ_r)**2 + 48*b/2**2*c_r*c_t*t_r/cos(Λ_r)**2 + 48*b/2**2*c_r*c_t*t_t/cos(Λ_r)**2 + 192*b/2**2*c_t**2*t_t/cos(Λ_r)**2 + 24*b/2*c_r**3*t_r*tan(Λ_r) + 32*b/2*c_r**2*c_t*t_r*tan(Λ_r) + 16*b/2*c_r**2*c_t*t_t*tan(Λ_r) + 24*b/2*c_r*c_t**2*t_r*tan(Λ_r) + 48*b/2*c_r*c_t**2*t_t*tan(Λ_r) + 96*b/2*c_t**3*t_t*tan(Λ_r) + 16*c_r**4*t_r**3*sin(ω_r)**2 - 28*c_r**4*t_r*sin(ω_r)**2 + 28*c_r**4*t_r + 4*c_r**3*c_t*t_r**3*sin(ω_r)**2 + 12*c_r**3*c_t*t_r**2*t_t*sin(ω_r)**2 - 21*c_r**3*c_t*t_r*sin(ω_r)**2 + 21*c_r**3*c_t*t_r - 7*c_r**3*c_t*t_t*sin(ω_r)**2 + 7*c_r**3*c_t*t_t + 8*c_r**2*c_t**2*t_r**2*t_t*sin(ω_r)**2 + 8*c_r**2*c_t**2*t_r*t_t**2*sin(ω_r)**2 - 14*c_r**2*c_t**2*t_r*sin(ω_r)**2 + 14*c_r**2*c_t**2*t_r - 14*c_r**2*c_t**2*t_t*sin(ω_r)**2 + 14*c_r**2*c_t**2*t_t - 320*c_r**2*t_r*x_cg**2 - 320*c_r**2*t_r*y_cg**2 + 12*c_r*c_t**3*t_r*t_t**2*sin(ω_r)**2 - 7*c_r*c_t**3*t_r*sin(ω_r)**2 + 7*c_r*c_t**3*t_r + 4*c_r*c_t**3*t_t**3*sin(ω_r)**2 - 21*c_r*c_t**3*t_t*sin(ω_r)**2 + 21*c_r*c_t**3*t_t - 160*c_r*c_t*t_r*x_cg**2 - 160*c_r*c_t*t_r*y_cg**2 - 160*c_r*c_t*t_t*x_cg**2 - 160*c_r*c_t*t_t*y_cg**2 + 16*c_t**4*t_t**3*sin(ω_r)**2 - 28*c_t**4*t_t*sin(ω_r)**2 + 28*c_t**4*t_t - 320*c_t**2*t_t*x_cg**2 - 320*c_t**2*t_t*y_cg**2)/960, True)) 

# Ixy_cpy  = Piecewise((b/2*p*(-2*b/2**2*c_r**2*t_r*tan(Λ_r) - 3*b/2**2*c_r*c_t*t_r*tan(Λ_r) - 3*b/2**2*c_r*c_t*t_t*tan(Λ_r) - 12*b/2**2*c_t**2*t_t*tan(Λ_r) + 15*b/2*(-24*c_r**3*t_r*sin(ω_r) + 24*c_r**3*t_r*sin(ω + ω_r) + 48*c_r**2*c_t*t_r*sin(ω_r) - 48*c_r**2*c_t*t_r*sin(ω + ω_r) + 24*c_r**2*c_t*t_t*sin(ω_r) - 24*c_r**2*c_t*t_t*sin(ω + ω_r) - 24*c_r*c_t**2*t_r*sin(ω_r) + 24*c_r*c_t**2*t_r*sin(ω + ω_r) - 48*c_r*c_t**2*t_t*sin(ω_r) + 48*c_r*c_t**2*t_t*sin(ω + ω_r) - c_t**3*t_t*ω**4*sin(ω + ω_r) + 24*c_t**3*t_t*sin(ω_r) - 24*c_t**3*t_t*sin(ω + ω_r) + ω**3*(c_r**3*t_r*cos(ω_r) + c_r*c_t**2*t_r*cos(ω + ω_r) + 2*c_r*c_t**2*t_t*cos(ω + ω_r) - 4*c_t**3*t_t*cos(ω + ω_r)) + 2*ω**2*(3*c_r**3*t_r*sin(ω_r) - 2*c_r**2*c_t*t_r*sin(ω_r) + 2*c_r**2*c_t*t_r*sin(ω + ω_r) - c_r**2*c_t*t_t*sin(ω_r) + c_r**2*c_t*t_t*sin(ω + ω_r) - 3*c_r*c_t**2*t_r*sin(ω + ω_r) - 6*c_r*c_t**2*t_t*sin(ω + ω_r) + 6*c_t**3*t_t*sin(ω + ω_r)) + 6*ω*(-3*c_r**3*t_r*cos(ω_r) - c_r**3*t_r*cos(ω + ω_r) + 4*c_r**2*c_t*t_r*cos(ω_r) + 4*c_r**2*c_t*t_r*cos(ω + ω_r) + 2*c_r**2*c_t*t_t*cos(ω_r) + 2*c_r**2*c_t*t_t*cos(ω + ω_r) - c_r*c_t**2*t_r*cos(ω_r) - 3*c_r*c_t**2*t_r*cos(ω + ω_r) - 2*c_r*c_t**2*t_t*cos(ω_r) - 6*c_r*c_t**2*t_t*cos(ω + ω_r) + 4*c_t**3*t_t*cos(ω + ω_r)))/ω**5 - 20*c_r**2*t_r*x_cg*y_cg - 10*c_r*c_t*t_r*x_cg*y_cg - 10*c_r*c_t*t_t*x_cg*y_cg - 20*c_t**2*t_t*x_cg*y_cg)/60, (ω >= -oo) & (ω < oo)), (-b/2*p*(8*b/2**2*c_r**2*t_r*tan(Λ_r) + 12*b/2**2*c_r*c_t*t_r*tan(Λ_r) + 12*b/2**2*c_r*c_t*t_t*tan(Λ_r) + 48*b/2**2*c_t**2*t_t*tan(Λ_r) + b/2*(3*c_r**3*t_r + 4*c_r**2*c_t*t_r + 2*c_r**2*c_t*t_t + 3*c_r*c_t**2*t_r + 6*c_r*c_t**2*t_t + 12*c_t**3*t_t)*cos(ω_r) + 80*c_r**2*t_r*x_cg*y_cg + 40*c_r*c_t*t_r*x_cg*y_cg + 40*c_r*c_t*t_t*x_cg*y_cg + 80*c_t**2*t_t*x_cg*y_cg)/240, True)) 

# Ixz_cpy  = Piecewise((b/2*p*(-256*c_r**2*t_r*x_cg*z_cg - 128*c_r*c_t*t_r*x_cg*z_cg - 128*c_r*c_t*t_t*x_cg*z_cg - 256*c_t**2*t_t*x_cg*z_cg + (48*c_r**4*t_r**3*sin(ω_r)**2 - 48*c_r**4*t_r**3*sin(ω + ω_r)**2 - 84*c_r**4*t_r*sin(ω_r)**2 + 84*c_r**4*t_r*sin(ω + ω_r)**2 - 48*c_r**3*c_t*t_r**3*sin(ω_r)**2 + 48*c_r**3*c_t*t_r**3*sin(ω + ω_r)**2 - 144*c_r**3*c_t*t_r**2*t_t*sin(ω_r)**2 + 144*c_r**3*c_t*t_r**2*t_t*sin(ω + ω_r)**2 + 252*c_r**3*c_t*t_r*sin(ω_r)**2 - 252*c_r**3*c_t*t_r*sin(ω + ω_r)**2 + 84*c_r**3*c_t*t_t*sin(ω_r)**2 - 84*c_r**3*c_t*t_t*sin(ω + ω_r)**2 + 144*c_r**2*c_t**2*t_r**2*t_t*sin(ω_r)**2 - 144*c_r**2*c_t**2*t_r**2*t_t*sin(ω + ω_r)**2 + 144*c_r**2*c_t**2*t_r*t_t**2*sin(ω_r)**2 - 144*c_r**2*c_t**2*t_r*t_t**2*sin(ω + ω_r)**2 - 252*c_r**2*c_t**2*t_r*sin(ω_r)**2 + 252*c_r**2*c_t**2*t_r*sin(ω + ω_r)**2 - 252*c_r**2*c_t**2*t_t*sin(ω_r)**2 + 252*c_r**2*c_t**2*t_t*sin(ω + ω_r)**2 - 144*c_r*c_t**3*t_r*t_t**2*sin(ω_r)**2 + 144*c_r*c_t**3*t_r*t_t**2*sin(ω + ω_r)**2 + 84*c_r*c_t**3*t_r*sin(ω_r)**2 - 84*c_r*c_t**3*t_r*sin(ω + ω_r)**2 - 48*c_r*c_t**3*t_t**3*sin(ω_r)**2 + 48*c_r*c_t**3*t_t**3*sin(ω + ω_r)**2 + 252*c_r*c_t**3*t_t*sin(ω_r)**2 - 252*c_r*c_t**3*t_t*sin(ω + ω_r)**2 + 48*c_t**4*t_t**3*sin(ω_r)**2 - 48*c_t**4*t_t**3*sin(ω + ω_r)**2 - 84*c_t**4*t_t*sin(ω_r)**2 + 84*c_t**4*t_t*sin(ω + ω_r)**2 + 4*ω**4*(8*c_r**4*t_r**3*sin(ω_r)**2 - 4*c_r**4*t_r**3*sin(ω + ω_r)**2 - 4*c_r**4*t_r**3*cos(ω + ω_r)**2 - 14*c_r**4*t_r*sin(ω_r)**2 + 7*c_r**4*t_r*sin(ω + ω_r)**2 + 7*c_r**4*t_r*cos(ω + ω_r)**2 - 4*c_t**4*t_t**3*sin(ω + ω_r)**2 + 4*c_t**4*t_t**3*cos(ω + ω_r)**2 + 7*c_t**4*t_t*sin(ω + ω_r)**2 - 7*c_t**4*t_t*cos(ω + ω_r)**2) + 2*ω**3*(-16*c_r**4*t_r**3*sin(2*ω_r) + 28*c_r**4*t_r*sin(2*ω_r) + 4*c_r**3*c_t*t_r**3*sin(2*ω_r) + 12*c_r**3*c_t*t_r**2*t_t*sin(2*ω_r) - 21*c_r**3*c_t*t_r*sin(2*ω_r) - 7*c_r**3*c_t*t_t*sin(2*ω_r) + 12*c_r*c_t**3*t_r*t_t**2*sin(2*ω + 2*ω_r) - 7*c_r*c_t**3*t_r*sin(2*ω + 2*ω_r) + 4*c_r*c_t**3*t_t**3*sin(2*ω + 2*ω_r) - 21*c_r*c_t**3*t_t*sin(2*ω + 2*ω_r) - 16*c_t**4*t_t**3*sin(2*ω + 2*ω_r) + 28*c_t**4*t_t*sin(2*ω + 2*ω_r)) + 6*ω**2*(-16*c_r**4*t_r**3*sin(ω_r)**2 + 8*c_r**4*t_r**3*sin(ω + ω_r)**2 + 8*c_r**4*t_r**3*cos(ω + ω_r)**2 + 28*c_r**4*t_r*sin(ω_r)**2 - 14*c_r**4*t_r*sin(ω + ω_r)**2 - 14*c_r**4*t_r*cos(ω + ω_r)**2 + 8*c_r**3*c_t*t_r**3*sin(ω_r)**2 - 4*c_r**3*c_t*t_r**3*sin(ω + ω_r)**2 - 4*c_r**3*c_t*t_r**3*cos(ω + ω_r)**2 + 24*c_r**3*c_t*t_r**2*t_t*sin(ω_r)**2 - 12*c_r**3*c_t*t_r**2*t_t*sin(ω + ω_r)**2 - 12*c_r**3*c_t*t_r**2*t_t*cos(ω + ω_r)**2 - 42*c_r**3*c_t*t_r*sin(ω_r)**2 + 21*c_r**3*c_t*t_r*sin(ω + ω_r)**2 + 21*c_r**3*c_t*t_r*cos(ω + ω_r)**2 - 14*c_r**3*c_t*t_t*sin(ω_r)**2 + 7*c_r**3*c_t*t_t*sin(ω + ω_r)**2 + 7*c_r**3*c_t*t_t*cos(ω + ω_r)**2 - 8*c_r**2*c_t**2*t_r**2*t_t*sin(ω_r)**2 + 8*c_r**2*c_t**2*t_r**2*t_t*sin(ω + ω_r)**2 - 8*c_r**2*c_t**2*t_r*t_t**2*sin(ω_r)**2 + 8*c_r**2*c_t**2*t_r*t_t**2*sin(ω + ω_r)**2 + 14*c_r**2*c_t**2*t_r*sin(ω_r)**2 - 14*c_r**2*c_t**2*t_r*sin(ω + ω_r)**2 + 14*c_r**2*c_t**2*t_t*sin(ω_r)**2 - 14*c_r**2*c_t**2*t_t*sin(ω + ω_r)**2 - 12*c_r*c_t**3*t_r*t_t**2*sin(ω + ω_r)**2 + 12*c_r*c_t**3*t_r*t_t**2*cos(ω + ω_r)**2 + 7*c_r*c_t**3*t_r*sin(ω + ω_r)**2 - 7*c_r*c_t**3*t_r*cos(ω + ω_r)**2 - 4*c_r*c_t**3*t_t**3*sin(ω + ω_r)**2 + 4*c_r*c_t**3*t_t**3*cos(ω + ω_r)**2 + 21*c_r*c_t**3*t_t*sin(ω + ω_r)**2 - 21*c_r*c_t**3*t_t*cos(ω + ω_r)**2 + 8*c_t**4*t_t**3*sin(ω + ω_r)**2 - 8*c_t**4*t_t**3*cos(ω + ω_r)**2 - 14*c_t**4*t_t*sin(ω + ω_r)**2 + 14*c_t**4*t_t*cos(ω + ω_r)**2) + 3*ω*(16*c_r**4*t_r**3*sin(2*ω_r) - 28*c_r**4*t_r*sin(2*ω_r) - 12*c_r**3*c_t*t_r**3*sin(2*ω_r) - 4*c_r**3*c_t*t_r**3*sin(2*ω + 2*ω_r) - 36*c_r**3*c_t*t_r**2*t_t*sin(2*ω_r) - 12*c_r**3*c_t*t_r**2*t_t*sin(2*ω + 2*ω_r) + 63*c_r**3*c_t*t_r*sin(2*ω_r) + 21*c_r**3*c_t*t_r*sin(2*ω + 2*ω_r) + 21*c_r**3*c_t*t_t*sin(2*ω_r) + 7*c_r**3*c_t*t_t*sin(2*ω + 2*ω_r) + 24*c_r**2*c_t**2*t_r**2*t_t*sin(2*ω_r) + 24*c_r**2*c_t**2*t_r**2*t_t*sin(2*ω + 2*ω_r) + 24*c_r**2*c_t**2*t_r*t_t**2*sin(2*ω_r) + 24*c_r**2*c_t**2*t_r*t_t**2*sin(2*ω + 2*ω_r) - 42*c_r**2*c_t**2*t_r*sin(2*ω_r) - 42*c_r**2*c_t**2*t_r*sin(2*ω + 2*ω_r) - 42*c_r**2*c_t**2*t_t*sin(2*ω_r) - 42*c_r**2*c_t**2*t_t*sin(2*ω + 2*ω_r) - 12*c_r*c_t**3*t_r*t_t**2*sin(2*ω_r) - 36*c_r*c_t**3*t_r*t_t**2*sin(2*ω + 2*ω_r) + 7*c_r*c_t**3*t_r*sin(2*ω_r) + 21*c_r*c_t**3*t_r*sin(2*ω + 2*ω_r) - 4*c_r*c_t**3*t_t**3*sin(2*ω_r) - 12*c_r*c_t**3*t_t**3*sin(2*ω + 2*ω_r) + 21*c_r*c_t**3*t_t*sin(2*ω_r) + 63*c_r*c_t**3*t_t*sin(2*ω + 2*ω_r) + 16*c_t**4*t_t**3*sin(2*ω + 2*ω_r) - 28*c_t**4*t_t*sin(2*ω + 2*ω_r)))/ω**5)/768, (ω >= -oo) & (ω < oo)), (b/2*p*(-640*c_r**2*t_r*x_cg*z_cg - 320*c_r*c_t*t_r*x_cg*z_cg - 320*c_r*c_t*t_t*x_cg*z_cg - 640*c_t**2*t_t*x_cg*z_cg + (-16*c_r**4*t_r**3 + 28*c_r**4*t_r - 4*c_r**3*c_t*t_r**3 - 12*c_r**3*c_t*t_r**2*t_t + 21*c_r**3*c_t*t_r + 7*c_r**3*c_t*t_t - 8*c_r**2*c_t**2*t_r**2*t_t - 8*c_r**2*c_t**2*t_r*t_t**2 + 14*c_r**2*c_t**2*t_r + 14*c_r**2*c_t**2*t_t - 12*c_r*c_t**3*t_r*t_t**2 + 7*c_r*c_t**3*t_r - 4*c_r*c_t**3*t_t**3 + 21*c_r*c_t**3*t_t - 16*c_t**4*t_t**3 + 28*c_t**4*t_t)*sin(2*ω_r))/1920, True)) 

# Iyz_cpy  = Piecewise((b/2*p*(3*b/2*(-24*c_r**3*t_r*cos(ω_r) + 24*c_r**3*t_r*cos(ω + ω_r) + 48*c_r**2*c_t*t_r*cos(ω_r) - 48*c_r**2*c_t*t_r*cos(ω + ω_r) + 24*c_r**2*c_t*t_t*cos(ω_r) - 24*c_r**2*c_t*t_t*cos(ω + ω_r) - 24*c_r*c_t**2*t_r*cos(ω_r) + 24*c_r*c_t**2*t_r*cos(ω + ω_r) - 48*c_r*c_t**2*t_t*cos(ω_r) + 48*c_r*c_t**2*t_t*cos(ω + ω_r) - c_t**3*t_t*ω**4*cos(ω + ω_r) + 24*c_t**3*t_t*cos(ω_r) - 24*c_t**3*t_t*cos(ω + ω_r) + ω**3*(-c_r**3*t_r*sin(ω_r) - c_r*c_t**2*t_r*sin(ω + ω_r) - 2*c_r*c_t**2*t_t*sin(ω + ω_r) + 4*c_t**3*t_t*sin(ω + ω_r)) + 2*ω**2*(3*c_r**3*t_r*cos(ω_r) - 2*c_r**2*c_t*t_r*cos(ω_r) + 2*c_r**2*c_t*t_r*cos(ω + ω_r) - c_r**2*c_t*t_t*cos(ω_r) + c_r**2*c_t*t_t*cos(ω + ω_r) - 3*c_r*c_t**2*t_r*cos(ω + ω_r) - 6*c_r*c_t**2*t_t*cos(ω + ω_r) + 6*c_t**3*t_t*cos(ω + ω_r)) + 6*ω*(3*c_r**3*t_r*sin(ω_r) + c_r**3*t_r*sin(ω + ω_r) - 4*c_r**2*c_t*t_r*sin(ω_r) - 4*c_r**2*c_t*t_r*sin(ω + ω_r) - 2*c_r**2*c_t*t_t*sin(ω_r) - 2*c_r**2*c_t*t_t*sin(ω + ω_r) + c_r*c_t**2*t_r*sin(ω_r) + 3*c_r*c_t**2*t_r*sin(ω + ω_r) + 2*c_r*c_t**2*t_t*sin(ω_r) + 6*c_r*c_t**2*t_t*sin(ω + ω_r) - 4*c_t**3*t_t*sin(ω + ω_r)))/ω**5 - 4*c_r**2*t_r*y_cg*z_cg - 2*c_r*c_t*t_r*y_cg*z_cg - 2*c_r*c_t*t_t*y_cg*z_cg - 4*c_t**2*t_t*y_cg*z_cg)/12, (ω >= -oo) & (ω < oo)), (b/2*p*(b/2*(3*c_r**3*t_r + 4*c_r**2*c_t*t_r + 2*c_r**2*c_t*t_t + 3*c_r*c_t**2*t_r + 6*c_r*c_t**2*t_t + 12*c_t**3*t_t)*sin(ω_r) - 80*c_r**2*t_r*y_cg*z_cg - 40*c_r*c_t*t_r*y_cg*z_cg - 40*c_r*c_t*t_t*y_cg*z_cg - 80*c_t**2*t_t*y_cg*z_cg)/240, True)) 

# Ixx_w0 = b/2*p*(32*b/2**2*c_r**2*t_r + 48*b/2**2*c_r*c_t*t_r + 48*b/2**2*c_r*c_t*t_t + 192*b/2**2*c_t**2*t_t + 16*c_r**4*t_r**3*cos(ω_r)**2 - 28*c_r**4*t_r*cos(ω_r)**2 + 28*c_r**4*t_r + 4*c_r**3*c_t*t_r**3*cos(ω_r)**2 + 12*c_r**3*c_t*t_r**2*t_t*cos(ω_r)**2 - 21*c_r**3*c_t*t_r*cos(ω_r)**2 + 21*c_r**3*c_t*t_r - 7*c_r**3*c_t*t_t*cos(ω_r)**2 + 7*c_r**3*c_t*t_t + 8*c_r**2*c_t**2*t_r**2*t_t*cos(ω_r)**2 + 8*c_r**2*c_t**2*t_r*t_t**2*cos(ω_r)**2 - 14*c_r**2*c_t**2*t_r*cos(ω_r)**2 + 14*c_r**2*c_t**2*t_r - 14*c_r**2*c_t**2*t_t*cos(ω_r)**2 + 14*c_r**2*c_t**2*t_t - 320*c_r**2*t_r*y_cg**2 - 320*c_r**2*t_r*z_cg**2 + 12*c_r*c_t**3*t_r*t_t**2*cos(ω_r)**2 - 7*c_r*c_t**3*t_r*cos(ω_r)**2 + 7*c_r*c_t**3*t_r + 4*c_r*c_t**3*t_t**3*cos(ω_r)**2 - 21*c_r*c_t**3*t_t*cos(ω_r)**2 + 21*c_r*c_t**3*t_t - 160*c_r*c_t*t_r*y_cg**2 - 160*c_r*c_t*t_r*z_cg**2 - 160*c_r*c_t*t_t*y_cg**2 - 160*c_r*c_t*t_t*z_cg**2 + 16*c_t**4*t_t**3*cos(ω_r)**2 - 28*c_t**4*t_t*cos(ω_r)**2 + 28*c_t**4*t_t - 320*c_t**2*t_t*y_cg**2 - 320*c_t**2*t_t*z_cg**2)/960
# Iyy_w0 = b/2*p*(32*b/2**2*c_r**2*t_r*tan(Λ_r)**2 + 48*b/2**2*c_r*c_t*t_r*tan(Λ_r)**2 + 48*b/2**2*c_r*c_t*t_t*tan(Λ_r)**2 + 192*b/2**2*c_t**2*t_t*tan(Λ_r)**2 + 24*b/2*c_r**3*t_r*tan(Λ_r) + 32*b/2*c_r**2*c_t*t_r*tan(Λ_r) + 16*b/2*c_r**2*c_t*t_t*tan(Λ_r) + 24*b/2*c_r*c_t**2*t_r*tan(Λ_r) + 48*b/2*c_r*c_t**2*t_t*tan(Λ_r) + 96*b/2*c_t**3*t_t*tan(Λ_r) + 16*c_r**4*t_r**3 + 28*c_r**4*t_r + 4*c_r**3*c_t*t_r**3 + 12*c_r**3*c_t*t_r**2*t_t + 21*c_r**3*c_t*t_r + 7*c_r**3*c_t*t_t + 8*c_r**2*c_t**2*t_r**2*t_t + 8*c_r**2*c_t**2*t_r*t_t**2 + 14*c_r**2*c_t**2*t_r + 14*c_r**2*c_t**2*t_t - 320*c_r**2*t_r*x_cg**2 - 320*c_r**2*t_r*z_cg**2 + 12*c_r*c_t**3*t_r*t_t**2 + 7*c_r*c_t**3*t_r + 4*c_r*c_t**3*t_t**3 + 21*c_r*c_t**3*t_t - 160*c_r*c_t*t_r*x_cg**2 - 160*c_r*c_t*t_r*z_cg**2 - 160*c_r*c_t*t_t*x_cg**2 - 160*c_r*c_t*t_t*z_cg**2 + 16*c_t**4*t_t**3 + 28*c_t**4*t_t - 320*c_t**2*t_t*x_cg**2 - 320*c_t**2*t_t*z_cg**2)/960
# Izz_w0 = b/2*p*(32*b/2**2*c_r**2*t_r/cos(Λ_r)**2 + 48*b/2**2*c_r*c_t*t_r/cos(Λ_r)**2 + 48*b/2**2*c_r*c_t*t_t/cos(Λ_r)**2 + 192*b/2**2*c_t**2*t_t/cos(Λ_r)**2 + 24*b/2*c_r**3*t_r*tan(Λ_r) + 32*b/2*c_r**2*c_t*t_r*tan(Λ_r) + 16*b/2*c_r**2*c_t*t_t*tan(Λ_r) + 24*b/2*c_r*c_t**2*t_r*tan(Λ_r) + 48*b/2*c_r*c_t**2*t_t*tan(Λ_r) + 96*b/2*c_t**3*t_t*tan(Λ_r) + 16*c_r**4*t_r**3*sin(ω_r)**2 - 28*c_r**4*t_r*sin(ω_r)**2 + 28*c_r**4*t_r + 4*c_r**3*c_t*t_r**3*sin(ω_r)**2 + 12*c_r**3*c_t*t_r**2*t_t*sin(ω_r)**2 - 21*c_r**3*c_t*t_r*sin(ω_r)**2 + 21*c_r**3*c_t*t_r - 7*c_r**3*c_t*t_t*sin(ω_r)**2 + 7*c_r**3*c_t*t_t + 8*c_r**2*c_t**2*t_r**2*t_t*sin(ω_r)**2 + 8*c_r**2*c_t**2*t_r*t_t**2*sin(ω_r)**2 - 14*c_r**2*c_t**2*t_r*sin(ω_r)**2 + 14*c_r**2*c_t**2*t_r - 14*c_r**2*c_t**2*t_t*sin(ω_r)**2 + 14*c_r**2*c_t**2*t_t - 320*c_r**2*t_r*x_cg**2 - 320*c_r**2*t_r*y_cg**2 + 12*c_r*c_t**3*t_r*t_t**2*sin(ω_r)**2 - 7*c_r*c_t**3*t_r*sin(ω_r)**2 + 7*c_r*c_t**3*t_r + 4*c_r*c_t**3*t_t**3*sin(ω_r)**2 - 21*c_r*c_t**3*t_t*sin(ω_r)**2 + 21*c_r*c_t**3*t_t - 160*c_r*c_t*t_r*x_cg**2 - 160*c_r*c_t*t_r*y_cg**2 - 160*c_r*c_t*t_t*x_cg**2 - 160*c_r*c_t*t_t*y_cg**2 + 16*c_t**4*t_t**3*sin(ω_r)**2 - 28*c_t**4*t_t*sin(ω_r)**2 + 28*c_t**4*t_t - 320*c_t**2*t_t*x_cg**2 - 320*c_t**2*t_t*y_cg**2)/960
# Ixy_w0 = b/2*p*(-8*b/2**2*c_r**2*t_r*tan(Λ_r) - 12*b/2**2*c_r*c_t*t_r*tan(Λ_r) - 12*b/2**2*c_r*c_t*t_t*tan(Λ_r) - 48*b/2**2*c_t**2*t_t*tan(Λ_r) + 3*b/2*c_r**3*t_r*cos(ω_r) + 4*b/2*c_r**2*c_t*t_r*cos(ω_r) + 2*b/2*c_r**2*c_t*t_t*cos(ω_r) + 3*b/2*c_r*c_t**2*t_r*cos(ω_r) + 6*b/2*c_r*c_t**2*t_t*cos(ω_r) + 12*b/2*c_t**3*t_t*cos(ω_r) - 80*c_r**2*t_r*x_cg*y_cg - 40*c_r*c_t*t_r*x_cg*y_cg - 40*c_r*c_t*t_t*x_cg*y_cg - 80*c_t**2*t_t*x_cg*y_cg)/240
# Ixz_w0 = b/2*p*(-16*c_r**4*t_r**3*sin(2*ω_r) + 28*c_r**4*t_r*sin(2*ω_r) - 4*c_r**3*c_t*t_r**3*sin(2*ω_r) - 12*c_r**3*c_t*t_r**2*t_t*sin(2*ω_r) + 21*c_r**3*c_t*t_r*sin(2*ω_r) + 7*c_r**3*c_t*t_t*sin(2*ω_r) - 8*c_r**2*c_t**2*t_r**2*t_t*sin(2*ω_r) - 8*c_r**2*c_t**2*t_r*t_t**2*sin(2*ω_r) + 14*c_r**2*c_t**2*t_r*sin(2*ω_r) + 14*c_r**2*c_t**2*t_t*sin(2*ω_r) - 640*c_r**2*t_r*x_cg*z_cg - 12*c_r*c_t**3*t_r*t_t**2*sin(2*ω_r) + 7*c_r*c_t**3*t_r*sin(2*ω_r) - 4*c_r*c_t**3*t_t**3*sin(2*ω_r) + 21*c_r*c_t**3*t_t*sin(2*ω_r) - 320*c_r*c_t*t_r*x_cg*z_cg - 320*c_r*c_t*t_t*x_cg*z_cg - 16*c_t**4*t_t**3*sin(2*ω_r) + 28*c_t**4*t_t*sin(2*ω_r) - 640*c_t**2*t_t*x_cg*z_cg)/1920
# Iyz_w0 = -b/2*p*(3*b/2*c_r**3*t_r*sin(ω_r) + 4*b/2*c_r**2*c_t*t_r*sin(ω_r) + 2*b/2*c_r**2*c_t*t_t*sin(ω_r) + 3*b/2*c_r*c_t**2*t_r*sin(ω_r) + 6*b/2*c_r*c_t**2*t_t*sin(ω_r) + 12*b/2*c_t**3*t_t*sin(ω_r) + 80*c_r**2*t_r*y_cg*z_cg + 40*c_r*c_t*t_r*y_cg*z_cg + 40*c_r*c_t*t_t*y_cg*z_cg + 80*c_t**2*t_t*y_cg*z_cg)/240

Ixx_cp1 = 0.272124422825884
Iyy_cp1 = 0.0312129806676198
Izz_cp1 = 0.299776791197862
Ixy_cp1 = -0.068469107042954 # SEE THE SIGN CHANGES HERE!!!
Ixz_cp1 = -0.00396349754460123
Iyz_cp1 = -0.00144473643314478 # SEE THE SIGN CHANGES HERE!!!

Ixx_lef = 0.272124422825884
Iyy_lef = 0.0312129806676198
Izz_lef = 0.299776791197862
Ixy_lef = 0.068469107042954 # SEE THE SIGN CHANGES HERE!!!
Ixz_lef = -0.00396349754460123
Iyz_lef = 0.00144473643314478 # SEE THE SIGN CHANGES HERE!!!

# print("Ixx_cpy  =", Ixx_cpy)
# print("Iyy_cpy  =", Iyy_cpy)
# print("Izz_cpy  =", Izz_cpy)
# print("Ixy_cpy  =", Ixy_cpy)
# print("Ixz_cpy  =", Ixz_cpy)
# print("Iyz_cpy  =", Iyz_cpy)

# print("Ixx_w0   =", Ixx_w0, "\n")
# print("Iyy_w0   =", Iyy_w0, "\n")
# print("Izz_w0   =", Izz_w0, "\n")
# print("Ixy_w0   =", Ixy_w0, "\n")
# print("Ixz_w0   =", Ixz_w0, "\n")
# print("Iyz_w0   =", Iyz_w0, "\n")

print("Ixx_cp1 =", Ixx_cp1)
print("Iyy_cp1 =", Iyy_cp1)
print("Izz_cp1 =", Izz_cp1)
print("Ixy_cp1 =", Ixy_cp1)
print("Ixz_cp1 =", Ixz_cp1)
print("Iyz_cp1 =", Iyz_cp1)
# print()
# print("Ixx % diff =", (Ixx_cpy - Ixx_cp1)/Ixx_cp1 * 100.0)
# print("Iyy % diff =", (Iyy_cpy - Iyy_cp1)/Iyy_cp1 * 100.0)
# print("Izz % diff =", (Izz_cpy - Izz_cp1)/Izz_cp1 * 100.0)
# print("Ixy % diff =", (Ixy_cpy - Ixy_cp1)/Ixy_cp1 * 100.0)
# print("Ixz % diff =", (Ixz_cpy - Ixz_cp1)/Ixz_cp1 * 100.0)
# print("Iyz % diff =", (Iyz_cpy - Iyz_cp1)/Iyz_cp1 * 100.0)



# ################ just mounting angle
# declaring variables...
# creating bounds...
# integrating...
#          mass...
#          moments...
#          sub inertias...
#          full inertias...
#                  Ixx...
#                  Iyy...
#                  Izz...
#                  Ixy...
#                  Ixz...
#                  Iyz...
# simplifying...
#          mass and moments...
#          moments of inertia...
#          products of inertia...
# reporting...


# m = b/2*k_a*p/6

# Myz = -b/2*p*(4*b/2*c_r**2*t_r*tan(Λ_r) + 4*b/2*c_r*c_t*t_r*tan(Λ_r) + 4*b/2*c_r*c_t*t_t*tan(Λ_r) + 12*b/2*c_t**2*t_t*tan(Λ_r) + 3*c_r**3*t_r + 2*c_r**2*c_t*t_r + c_r**2*c_t*t_t + c_r*c_t**2*t_r + 2*c_r*c_t**2*t_t + 3*c_t**3*t_t)/48

# Mxz = b/2**2*p*(c_r**2*t_r + c_r*c_t*t_r + c_r*c_t*t_t + 3*c_t**2*t_t)/12

# Mxy = 0

# Ixx_new = b/2*p*(32*b/2**2*c_r**2*t_r + 48*b/2**2*c_r*c_t*t_r + 48*b/2**2*c_r*c_t*t_t + 192*b/2**2*c_t**2*t_t + 16*c_r**4*t_r**3*cos(ω_r)**2 - 28*c_r**4*t_r*cos(ω_r)**2 + 28*c_r**4*t_r + 4*c_r**3*c_t*t_r**3*cos(ω_r)**2 + 12*c_r**3*c_t*t_r**2*t_t*cos(ω_r)**2 - 21*c_r**3*c_t*t_r*cos(ω_r)**2 + 21*c_r**3*c_t*t_r - 7*c_r**3*c_t*t_t*cos(ω_r)**2 + 7*c_r**3*c_t*t_t + 8*c_r**2*c_t**2*t_r**2*t_t*cos(ω_r)**2 + 8*c_r**2*c_t**2*t_r*t_t**2*cos(ω_r)**2 - 14*c_r**2*c_t**2*t_r*cos(ω_r)**2 + 14*c_r**2*c_t**2*t_r - 14*c_r**2*c_t**2*t_t*cos(ω_r)**2 + 14*c_r**2*c_t**2*t_t - 320*c_r**2*t_r*y_cg**2 - 320*c_r**2*t_r*z_cg**2 + 12*c_r*c_t**3*t_r*t_t**2*cos(ω_r)**2 - 7*c_r*c_t**3*t_r*cos(ω_r)**2 + 7*c_r*c_t**3*t_r + 4*c_r*c_t**3*t_t**3*cos(ω_r)**2 - 21*c_r*c_t**3*t_t*cos(ω_r)**2 + 21*c_r*c_t**3*t_t - 160*c_r*c_t*t_r*y_cg**2 - 160*c_r*c_t*t_r*z_cg**2 - 160*c_r*c_t*t_t*y_cg**2 - 160*c_r*c_t*t_t*z_cg**2 + 16*c_t**4*t_t**3*cos(ω_r)**2 - 28*c_t**4*t_t*cos(ω_r)**2 + 28*c_t**4*t_t - 320*c_t**2*t_t*y_cg**2 - 320*c_t**2*t_t*z_cg**2)/960
# Iyy_new = b/2*p*(32*b/2**2*c_r**2*t_r*tan(Λ_r)**2 + 48*b/2**2*c_r*c_t*t_r*tan(Λ_r)**2 + 48*b/2**2*c_r*c_t*t_t*tan(Λ_r)**2 + 192*b/2**2*c_t**2*t_t*tan(Λ_r)**2 + 24*b/2*c_r**3*t_r*tan(Λ_r) + 32*b/2*c_r**2*c_t*t_r*tan(Λ_r) + 16*b/2*c_r**2*c_t*t_t*tan(Λ_r) + 24*b/2*c_r*c_t**2*t_r*tan(Λ_r) + 48*b/2*c_r*c_t**2*t_t*tan(Λ_r) + 96*b/2*c_t**3*t_t*tan(Λ_r) + 16*c_r**4*t_r**3 + 28*c_r**4*t_r + 4*c_r**3*c_t*t_r**3 + 12*c_r**3*c_t*t_r**2*t_t + 21*c_r**3*c_t*t_r + 7*c_r**3*c_t*t_t + 8*c_r**2*c_t**2*t_r**2*t_t + 8*c_r**2*c_t**2*t_r*t_t**2 + 14*c_r**2*c_t**2*t_r + 14*c_r**2*c_t**2*t_t - 320*c_r**2*t_r*x_cg**2 - 320*c_r**2*t_r*z_cg**2 + 12*c_r*c_t**3*t_r*t_t**2 + 7*c_r*c_t**3*t_r + 4*c_r*c_t**3*t_t**3 + 21*c_r*c_t**3*t_t - 160*c_r*c_t*t_r*x_cg**2 - 160*c_r*c_t*t_r*z_cg**2 - 160*c_r*c_t*t_t*x_cg**2 - 160*c_r*c_t*t_t*z_cg**2 + 16*c_t**4*t_t**3 + 28*c_t**4*t_t - 320*c_t**2*t_t*x_cg**2 - 320*c_t**2*t_t*z_cg**2)/960
# Izz_new = b/2*p*(32*b/2**2*c_r**2*t_r/cos(Λ_r)**2 + 48*b/2**2*c_r*c_t*t_r/cos(Λ_r)**2 + 48*b/2**2*c_r*c_t*t_t/cos(Λ_r)**2 + 192*b/2**2*c_t**2*t_t/cos(Λ_r)**2 + 24*b/2*c_r**3*t_r*tan(Λ_r) + 32*b/2*c_r**2*c_t*t_r*tan(Λ_r) + 16*b/2*c_r**2*c_t*t_t*tan(Λ_r) + 24*b/2*c_r*c_t**2*t_r*tan(Λ_r) + 48*b/2*c_r*c_t**2*t_t*tan(Λ_r) + 96*b/2*c_t**3*t_t*tan(Λ_r) + 16*c_r**4*t_r**3*sin(ω_r)**2 - 28*c_r**4*t_r*sin(ω_r)**2 + 28*c_r**4*t_r + 4*c_r**3*c_t*t_r**3*sin(ω_r)**2 + 12*c_r**3*c_t*t_r**2*t_t*sin(ω_r)**2 - 21*c_r**3*c_t*t_r*sin(ω_r)**2 + 21*c_r**3*c_t*t_r - 7*c_r**3*c_t*t_t*sin(ω_r)**2 + 7*c_r**3*c_t*t_t + 8*c_r**2*c_t**2*t_r**2*t_t*sin(ω_r)**2 + 8*c_r**2*c_t**2*t_r*t_t**2*sin(ω_r)**2 - 14*c_r**2*c_t**2*t_r*sin(ω_r)**2 + 14*c_r**2*c_t**2*t_r - 14*c_r**2*c_t**2*t_t*sin(ω_r)**2 + 14*c_r**2*c_t**2*t_t - 320*c_r**2*t_r*x_cg**2 - 320*c_r**2*t_r*y_cg**2 + 12*c_r*c_t**3*t_r*t_t**2*sin(ω_r)**2 - 7*c_r*c_t**3*t_r*sin(ω_r)**2 + 7*c_r*c_t**3*t_r + 4*c_r*c_t**3*t_t**3*sin(ω_r)**2 - 21*c_r*c_t**3*t_t*sin(ω_r)**2 + 21*c_r*c_t**3*t_t - 160*c_r*c_t*t_r*x_cg**2 - 160*c_r*c_t*t_r*y_cg**2 - 160*c_r*c_t*t_t*x_cg**2 - 160*c_r*c_t*t_t*y_cg**2 + 16*c_t**4*t_t**3*sin(ω_r)**2 - 28*c_t**4*t_t*sin(ω_r)**2 + 28*c_t**4*t_t - 320*c_t**2*t_t*x_cg**2 - 320*c_t**2*t_t*y_cg**2)/960
# Ixy_new = b/2*p*(-8*b/2**2*c_r**2*t_r*tan(Λ_r) - 12*b/2**2*c_r*c_t*t_r*tan(Λ_r) - 12*b/2**2*c_r*c_t*t_t*tan(Λ_r) - 48*b/2**2*c_t**2*t_t*tan(Λ_r) + 3*b/2*c_r**3*t_r*cos(ω_r) + 4*b/2*c_r**2*c_t*t_r*cos(ω_r) + 2*b/2*c_r**2*c_t*t_t*cos(ω_r) + 3*b/2*c_r*c_t**2*t_r*cos(ω_r) + 6*b/2*c_r*c_t**2*t_t*cos(ω_r) + 12*b/2*c_t**3*t_t*cos(ω_r) - 80*c_r**2*t_r*x_cg*y_cg - 40*c_r*c_t*t_r*x_cg*y_cg - 40*c_r*c_t*t_t*x_cg*y_cg - 80*c_t**2*t_t*x_cg*y_cg)/240
# Ixz_new = b/2*p*(-16*c_r**4*t_r**3*sin(2*ω_r) + 28*c_r**4*t_r*sin(2*ω_r) - 4*c_r**3*c_t*t_r**3*sin(2*ω_r) - 12*c_r**3*c_t*t_r**2*t_t*sin(2*ω_r) + 21*c_r**3*c_t*t_r*sin(2*ω_r) + 7*c_r**3*c_t*t_t*sin(2*ω_r) - 8*c_r**2*c_t**2*t_r**2*t_t*sin(2*ω_r) - 8*c_r**2*c_t**2*t_r*t_t**2*sin(2*ω_r) + 14*c_r**2*c_t**2*t_r*sin(2*ω_r) + 14*c_r**2*c_t**2*t_t*sin(2*ω_r) - 640*c_r**2*t_r*x_cg*z_cg - 12*c_r*c_t**3*t_r*t_t**2*sin(2*ω_r) + 7*c_r*c_t**3*t_r*sin(2*ω_r) - 4*c_r*c_t**3*t_t**3*sin(2*ω_r) + 21*c_r*c_t**3*t_t*sin(2*ω_r) - 320*c_r*c_t*t_r*x_cg*z_cg - 320*c_r*c_t*t_t*x_cg*z_cg - 16*c_t**4*t_t**3*sin(2*ω_r) + 28*c_t**4*t_t*sin(2*ω_r) - 640*c_t**2*t_t*x_cg*z_cg)/1920
# Iyz_new = -b/2*p*(3*b/2*c_r**3*t_r*sin(ω_r) + 4*b/2*c_r**2*c_t*t_r*sin(ω_r) + 2*b/2*c_r**2*c_t*t_t*sin(ω_r) + 3*b/2*c_r*c_t**2*t_r*sin(ω_r) + 6*b/2*c_r*c_t**2*t_t*sin(ω_r) + 12*b/2*c_t**3*t_t*sin(ω_r) + 80*c_r**2*t_r*y_cg*z_cg + 40*c_r*c_t*t_r*y_cg*z_cg + 40*c_r*c_t*t_t*y_cg*z_cg + 80*c_t**2*t_t*y_cg*z_cg)/240


# replacement variables
b_2 = 8.0
cr = 2.0
ct = 1.50
tr = tt = 0.12
Lr = np.deg2rad(15.0)
p = 0.0175
tanL = sy.tan(Lr)
wr = np.deg2rad(20.0)
w = np.deg2rad(18.0)
xcg = -1.38708173
ycg = 3.62162162
zcg = 0.15165589

# # reformatted coefficients, and some changed
ka = tr*cr*(2*cr+ct) + tt*ct*(cr+2*ct)
kb = tr*cr*(cr+ct) + tt*ct*(cr+3*ct)
kc = tr*cr*(3*cr**2+2*cr*ct+ct**2) + tt*ct*(cr**2+2*cr*ct+3*ct**2)
kd1 = tr*cr*(2*cr+3*ct) + tt*ct*(3*cr+12*ct)
ke = tr*cr*(3*cr**2+4*cr*ct+3*ct**2) + 2*tt*ct*(cr**2+3*cr*ct+6*ct**2)
kf1 = tr*cr*(4*cr**3+3*cr**2*ct+2*cr*ct**2+ct**3) + \
    tt*ct*(cr**3+2*cr**2*ct+3*cr*ct**2+4*ct**3)
kg = tr**3*cr**3*(4*cr+ct) + tr**2*cr**2*tt*ct*(3*cr+2*ct) + \
    tr*cr*tt**2*ct**2*(2*cr+3*ct) + tt**3*ct**3*(cr+4*ct)
m = b_2*ka*p/6
minv = 6/b_2/ka/p

# simplifying equations
Ixx_new = m * ( 1/160 * (16*b_2**2*kd1 + sy.cos(wr)**2*(4*kg - 7*kf1) + 7*kf1 )/ka - ycg**2 - zcg**2 )
Iyy_new = m * ( 1/160 * (4*b_2*(b_2*sy.tan(Lr)**2*kd1 + 2*sy.tan(Lr)*ke) + 4*kg + 7*kf1 )/ka  - xcg**2 - zcg**2)
Izz_new = m * ( 1/160 * (4*b_2*(b_2*(sy.tan(Lr)**2 + 1)*kd1 + 2*sy.tan(Lr)*ke) + sy.sin(wr)**2*(4*kg - 7*kf1) + 7*kf1)/ka - xcg**2 - ycg**2 )
Ixy_new =-m * ( 1/ 40 * b_2*(4*b_2*kd1*sy.tan(Lr) - ke*sy.cos(wr) )/ka + xcg*ycg ) # test out left twist wing
# add delta function
Ixz_new =-m * ( 1/320 * (4*kg*sy.sin(2*wr) - 7*kf1*sy.sin(2*wr))/ka + xcg*zcg )
Iyz_new =-m * ( 1/ 40 * ( b_2*ke*sy.sin(wr) )/ka + ycg*zcg )

print()
print("Ixx_new =", Ixx_new)
print("Iyy_new =", Iyy_new)
print("Izz_new =", Izz_new)
print("Ixy_new =", Ixy_new)
print("Ixz_new =", Ixz_new)
print("Iyz_new =", Iyz_new)
print()
print("Ixx % diff =", (Ixx_new - Ixx_cp1)/Ixx_cp1 * 100.0)
print("Iyy % diff =", (Iyy_new - Iyy_cp1)/Iyy_cp1 * 100.0)
print("Izz % diff =", (Izz_new - Izz_cp1)/Izz_cp1 * 100.0)
print("Ixy % diff =", (Ixy_new - Ixy_cp1)/Ixy_cp1 * 100.0)
print("Ixz % diff =", (Ixz_new - Ixz_cp1)/Ixz_cp1 * 100.0)
print("Iyz % diff =", (Iyz_new - Iyz_cp1)/Iyz_cp1 * 100.0)

