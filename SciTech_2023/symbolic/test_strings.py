import numpy as np
import sympy as sy
from sympy.parsing.sympy_parser import parse_expr

a = r"35 c_r^4 \tau_r^3 + 15 c_r^4 \tau_r^2 \tau_t + 5 c_r^4 \tau_r \tau_t^2 + c_r^4 \tau_t^3 + 20 c_r^3 c_t \tau_r^3 + 20 c_r^3 c_t \tau_r^2 \tau_t + 12 c_r^3 c_t \tau_r \tau_t^2 + 4 c_r^3 c_t \tau_t^3 + 10 c_r^2 c_t^2 \tau_r^3 + 18 c_r^2 c_t^2 \tau_r^2 \tau_t + 18 c_r^2 c_t^2 \tau_r \tau_t^2 + 10 c_r^2 c_t^2 \tau_t^3 + 4 c_r c_t^3 \tau_r^3 + 12 c_r c_t^3 \tau_r^2 \tau_t + 20 c_r c_t^3 \tau_r \tau_t^2 + 20 c_r c_t^3 \tau_t^3 + c_t^4 \tau_r^3 + 5 c_t^4 \tau_r^2 \tau_t + 15 c_t^4 \tau_r \tau_t^2 + 35 c_t^4 \tau_t^3"

b = r"104458200 a_4^3 c_r^4 \tau_r^3 + 44767800 a_4^3 c_r^4 \tau_r^2 \tau_t + 14922600 a_4^3 c_r^4 \tau_r \tau_t^2 + 2984520 a_4^3 c_r^4 \tau_t^3 + 59690400 a_4^3 c_r^3 c_t \tau_r^3 + 59690400 a_4^3 c_r^3 c_t \tau_r^2 \tau_t + 35814240 a_4^3 c_r^3 c_t \tau_r \tau_t^2 + 11938080 a_4^3 c_r^3 c_t \tau_t^3 + 29845200 a_4^3 c_r^2 c_t^2 \tau_r^3 + 53721360 a_4^3 c_r^2 c_t^2 \tau_r^2 \tau_t + 53721360 a_4^3 c_r^2 c_t^2 \tau_r \tau_t^2 + 29845200 a_4^3 c_r^2 c_t^2 \tau_t^3 + 11938080 a_4^3 c_r c_t^3 \tau_r^3 + 35814240 a_4^3 c_r c_t^3 \tau_r^2 \tau_t + 59690400 a_4^3 c_r c_t^3 \tau_r \tau_t^2 + 59690400 a_4^3 c_r c_t^3 \tau_t^3 + 2984520 a_4^3 c_t^4 \tau_r^3 + 14922600 a_4^3 c_t^4 \tau_r^2 \tau_t + 44767800 a_4^3 c_t^4 \tau_r \tau_t^2 + 104458200 a_4^3 c_t^4 \tau_t^3"


# determine part to remove
get_rem = b.split()
i = 0
rem = ""
end_loop = False
while not end_loop:
    if not get_rem[i].isnumeric():
        if get_rem[i][0] == "a":
            rem += get_rem[i] + " "
        else:
            end_loop = not end_loop
    i += 1
    
b = b.replace(rem,"")
a_split = a.split()
b_split = b.split()

a_num = int(a_split[0])
b_num = int(b_split[0])

c_num = int(b_num / a_num)

c = ""

# np.gcd.reduce(array) greatest commond divisor

for i in range(len(b_split)):

    if b_split[i].isnumeric():
        b_split[i] = str( int( int(b_split[i]) / c_num) )
    
    c += b_split[i]
    if i != len(b_split)-1:
        c += " "

c = c.replace(" 1 "," ")

print()
print(str(c_num) + " " + rem)
print()
print(c==a)
# print(a)
# print(c)

# d = "15519504 a_0^3 + 38798760 a_0^2 a_1 + 29099070 a_0^2 a_2 + 23279256 a_0^2 a_3 + 19399380 a_0^2 a_4 + 33256080 a_0 a_1^2 + 51731680 a_0 a_1 a_2 + 42325920 a_0 a_1 a_3 + 35814240 a_0 a_1 a_4 + 21162960 a_0 a_2^2 + 35814240 a_0 a_2 a_3 + 31039008 a_0 a_2 a_4 + 15519504 a_0 a_3^2 + 27387360 a_0 a_3 a_4 + 12252240 a_0 a_4^2"
# d = "9699690 a_1^3 + 23279256 a_1^2 a_2 + 19399380 a_1^2 a_3 + 16628040 a_1^2 a_4 + 19399380 a_1 a_2^2 + 33256080 a_1 a_2 a_3 + 29099070 a_1 a_2 a_4 + 14549535 a_1 a_3^2 + 25865840 a_1 a_3 a_4 + 11639628 a_1 a_4^2"
# d = "5542680 a_2^3 + 14549535 a_2^2 a_3 + 12932920 a_2^2 a_4 + 12932920 a_2 a_3^2 + 23279256 a_2 a_3 a_4 + 10581480 a_2 a_4^2"
# d = "3879876 a_3^3 + 10581480 a_3^2 a_4 + 9699690 a_3 a_4^2"
# d = "2984520 a_4^3"
d = "15519504 a_0^3 + 38798760 a_0^2 a_1 + 29099070 a_0^2 a_2 + 23279256 a_0^2 a_3 + 19399380 a_0^2 a_4 + 33256080 a_0 a_1^2 + 51731680 a_0 a_1 a_2 + 42325920 a_0 a_1 a_3 + 35814240 a_0 a_1 a_4 + 21162960 a_0 a_2^2 + 35814240 a_0 a_2 a_3 + 31039008 a_0 a_2 a_4 + 15519504 a_0 a_3^2 + 27387360 a_0 a_3 a_4 + 12252240 a_0 a_4^2 + 9699690 a_1^3 + 23279256 a_1^2 a_2 + 19399380 a_1^2 a_3 + 16628040 a_1^2 a_4 + 19399380 a_1 a_2^2 + 33256080 a_1 a_2 a_3 + 29099070 a_1 a_2 a_4 + 14549535 a_1 a_3^2 + 25865840 a_1 a_3 a_4 + 11639628 a_1 a_4^2 + 5542680 a_2^3 + 14549535 a_2^2 a_3 + 12932920 a_2^2 a_4 + 12932920 a_2 a_3^2 + 23279256 a_2 a_3 a_4 + 10581480 a_2 a_4^2 + 3879876 a_3^3 + 10581480 a_3^2 a_4 + 9699690 a_3 a_4^2 + 2984520 a_4^3"

d_split = d.split()

factors = []

for i in range(len(d_split)):

    if d_split[i].isnumeric():
        factors.append( int(d_split[i]) )

        # print(d_split[i+1], factors[-1] / 38798760.0)


divisor = np.gcd.reduce(factors)
print(divisor)
print(np.gcd.reduce([856,770,644,553,484]))

a0,a1,a2,a3,a4,a5 = sy.symbols("a_0 a_1 a_2 a_3 a_4 a_5")
A,B,C,D,E,F,G,H,I,J,K,L,M = sy.symbols("A B C D E F G H I J K L M")

eq = (A*a0 + B*a1 + C*a2 + D*a3 + E*a4)**3# + (G*a0 + H*a1 + I*a2 + J*a3 + K*a4)
print(sy.simplify(sy.expand(eq)))

## 2 a_0
## 46189 a_1
## 4199 a_2
## 176358 a_3

# print(factors[0]**(1./3.))

dsy = d.replace(" a"," * a")
dsy = dsy.replace("^","**")
d_sym = parse_expr(dsy)

total_divisor = "38798760"

e = ""

for i in range(len(d_split)):
    if d_split[i].isnumeric():
        inputs = d_split[i] + "/" + total_divisor
        output = sy.simplify(parse_expr(inputs))
        # print(inputs,output)

        d_split[i] = str(output)

    e += d_split[i]
    if i != len(d_split)-1:
        e += " "

print(e)


