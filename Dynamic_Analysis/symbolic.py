import sympy as sy


sym = sy.Symbol
simp = sy.simplify
exp = sy.expand
sqrt = sy.sqrt
fact = sy.factor

g = sym("g")
Vo = sym("Vo")
Kzmu = sym("Kzmu")
Kxmu = sym("Kxmu")
Kmal = sym("Kmal")
Kzal = sym("Kzal")
Kmqb = sym("Kmqb")
Kxal = sym("Kxal")
Kyb = sym("Kyb")
Kyp = sym("Kyp")
Kyr = sym("Kyr")
Klb = sym("Klb")
Klp = sym("Klp")
Klr = sym("Klr")
Knb = sym("Knb")
Knp = sym("Knp")
Knr = sym("Knr")

# ratios
RPHs = Kmal / (Kmal - Kzal * Kmqb)
RPHd = Kxal * Kmqb / (Kmal - Kzal * Kmqb)
RPHp = RPHs * (Kzal + Kmqb) / (Kmal - Kzal * Kmqb)
RDs = (Klb*(1-(1-Kyr)*Knp) - Kyb*Klr*Knp) / Klp
RDs = (Klb - Kyb*Klr*Knp) / Klp #### simplified
RDc = Klr * Knp / Klp
RDp = (Klr*Knb-Klb*Knr)/Klp/(Knb+Kyb*Knr) - RDs/Klp
RDp = - RDs/Klp #### simplified


# sigma, omega_d # phugiod
S = - g/Vo *( -Kzmu/2 * ( -Kxmu/Kzmu + RPHd - RPHp ) )
W =   g/Vo *( -Kzmu/2 * sqrt( 4/Kzmu*RPHs - ( Kxmu/( -Kzmu + RPHd ) )**2 ) )
# sigma, omega_d # dutch roll
S = -1/2*g/Vo * (Kyb + Knr - RDc + RDp)
W = g/Vo * sy.sqrt( (1-Kyr)*Knb + Kyb*Knr + RDs - ((Kyb+Knr)/2)**2 )

Wn = sqrt( S**2 + W**2 )

# print(Wn)
# print()
# print(exp(Wn))
# print()
# print(simp(Wn))
# print()

Z = S / Wn

Do = sym("Do")
W = sym("W")
Lo = sym("Lo")
Da = sym("Da")
La = sym("La")
ma = sym("ma")
mq = sym("mq")
Iyy = sym("Iyy")
Yb = sym("Yb")
Yp = sym("Yp")
Yr = sym("Yr")
lb = sym("lb")
lp = sym("lp")
lr = sym("lr")
nb = sym("nb")
np = sym("np")
nr = sym("nr")
Ixx = sym("Ixx")
Izz = sym("Izz")

replacers = [
    [Kxmu,-2*Do/W],
    [Kzmu,2*Lo/W],
    [Kxal,(Lo-Da)/W],
    [Kzal,(-Do-La)/W],
    [Kmal,Vo**2*ma/g**2/Iyy],
    [Kmqb,Vo*mq/g/Iyy],
    [Kyb,Yb/W],
    [Kyp,g*Yp/Vo/W],
    [Kyr,g*Yr/Vo/W],
    [Klb,Vo**2*lb/g**2/Ixx],
    [Klp,Vo*lp/g/Ixx],
    [Klr,Vo*lr/g/Ixx],
    [Knb,Vo**2*nb/g**2/Ixx],
    [Knp,Vo*np/g/Ixx],
    [Knr,Vo*nr/g/Ixx]
]
Wn = Wn.subs(replacers)
print(Wn)
print()
print(exp(Wn))
print()
print(simp(Wn))
# Z = Z.subs(replacers)

# print(Z)
# print()
# # print(exp(Z))
# print()
# print(simp(Z))


# a = sym("a")
# b = sym("b")
# c = sym("c")
# d = sym("d")

# one = (a+b-c+d)**2 - (a+b)**2
# two = 2*(a+b)*(d-c) + (d-c)**2
# print(simp(one-two))



