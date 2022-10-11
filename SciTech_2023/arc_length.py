import numpy as np
from matplotlib import pyplot as plt

a0 = 2.980
a1 = -1.320
a2 = -3.286
a3 = 2.441
a4 = -0.815

def y(x,tau):
    return tau * (a0*x**0.5+a1*x+a2*x**2.+a3*x**3.+a4*x**4.)

ts = np.linspace(0.01,0.4,100)
ths = np.linspace(0.0,np.pi,200)
xs = 0.5 * ( 1. - np.cos(ths) )

xaf = np.append(np.flip(xs)[:-1],xs)
xdfsq = (xaf[1:] - xaf[:-1])**2.

ls = np.zeros(ts.shape[0])

for i in range(ls.shape[0]):

    ys = y(xs,ts[i])
    yaf = np.append(np.flip(ys)[:-1],-ys)
    ydfsq = (yaf[1:] - yaf[:-1])**2.

    ls[i] = np.sum( (xdfsq + ydfsq)**0.5 )
    print("{:< 19.16f}\t{:< 19.16f}".format(ts[i],ls[i]))


plt.plot(ts,ls,"b.")
plt.show()