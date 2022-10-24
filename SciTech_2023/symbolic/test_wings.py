import numpy as np
from matplotlib import pyplot as plt

n = 201
theta = np.linspace(0.0,np.pi,n)
z_b = -0.5 * np.cos(theta)

y_circle = (0.5**2. - z_b**2.)**0.5 * 2.0
c_up =   y_circle * 0.5
c_lo = - y_circle * 0.5

p_up = y_circle *  0.25
p_lo = y_circle * -0.75

RA = 8.0
e_c = 4. / np.pi / RA * np.sin(theta)
e_up = e_c * 0.25
e_lo = e_c * -0.75


# approx elliptic wing
n_half = (n-1)/2
n_half = int(n_half)

c_root_ell = e_c[n_half]
c_root_ELL = y_circle[n_half]

a_c = y_circle * c_root_ell / c_root_ELL

a_up = a_c * 0.25
a_lo = a_c * -0.75

plt.plot(z_b,c_up,"k")
plt.plot(z_b,c_lo,"k")
plt.plot(z_b,c_up*0.5,"k")
plt.plot(z_b,p_up,"b")
plt.plot(z_b,p_lo,"b")
plt.plot(z_b,a_up,"r")
plt.plot(z_b,a_lo,"r")
plt.plot(z_b,e_up,"g")
plt.plot(z_b,e_lo,"g")
plt.axis("equal")
plt.show()
