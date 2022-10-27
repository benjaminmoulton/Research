
a = r"35 c_r^4 \tau_r^3 + 15 c_r^4 \tau_r^2 \tau_t + 5 c_r^4 \tau_r \tau_t^2 + c_r^4 \tau_t^3 + 20 c_r^3 c_t \tau_r^3 + 20 c_r^3 c_t \tau_r^2 \tau_t + 12 c_r^3 c_t \tau_r \tau_t^2 + 4 c_r^3 c_t \tau_t^3 + 10 c_r^2 c_t^2 \tau_r^3 + 18 c_r^2 c_t^2 \tau_r^2 \tau_t + 18 c_r^2 c_t^2 \tau_r \tau_t^2 + 10 c_r^2 c_t^2 \tau_t^3 + 4 c_r c_t^3 \tau_r^3 + 12 c_r c_t^3 \tau_r^2 \tau_t + 20 c_r c_t^3 \tau_r \tau_t^2 + 20 c_r c_t^3 \tau_t^3 + c_t^4 \tau_r^3 + 5 c_t^4 \tau_r^2 \tau_t + 15 c_t^4 \tau_r \tau_t^2 + 35 c_t^4 \tau_t^3"

b = r"543182640 a_0 a_3^2 c_r^4 \tau_r^3 + 232792560 a_0 a_3^2 c_r^4 \tau_r^2 \tau_t + 77597520 a_0 a_3^2 c_r^4 \tau_r \tau_t^2 + 15519504 a_0 a_3^2 c_r^4 \tau_t^3 + 310390080 a_0 a_3^2 c_r^3 c_t \tau_r^3 + 310390080 a_0 a_3^2 c_r^3 c_t \tau_r^2 \tau_t + 186234048 a_0 a_3^2 c_r^3 c_t \tau_r \tau_t^2 + 62078016 a_0 a_3^2 c_r^3 c_t \tau_t^3 + 155195040 a_0 a_3^2 c_r^2 c_t^2 \tau_r^3 + 279351072 a_0 a_3^2 c_r^2 c_t^2 \tau_r^2 \tau_t + 279351072 a_0 a_3^2 c_r^2 c_t^2 \tau_r \tau_t^2 + 155195040 a_0 a_3^2 c_r^2 c_t^2 \tau_t^3 + 62078016 a_0 a_3^2 c_r c_t^3 \tau_r^3 + 186234048 a_0 a_3^2 c_r c_t^3 \tau_r^2 \tau_t + 310390080 a_0 a_3^2 c_r c_t^3 \tau_r \tau_t^2 + 310390080 a_0 a_3^2 c_r c_t^3 \tau_t^3 + 15519504 a_0 a_3^2 c_t^4 \tau_r^3 + 77597520 a_0 a_3^2 c_t^4 \tau_r^2 \tau_t + 232792560 a_0 a_3^2 c_t^4 \tau_r \tau_t^2 + 543182640 a_0 a_3^2 c_t^4 \tau_t^3"


rem = "a_0 a_3^2 "
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

print(c_num,c==a)
# print(a)
# print(c)