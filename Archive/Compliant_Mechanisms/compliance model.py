#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 10 10:08:09 2019

@author: benjamin
"""
from aerolab import gcode as gc
from scipy.interpolate import interp1d as interp

fx, fy = gc.NACAFileReader('NACA 0009.txt')

percent = .35 # % along chord length at which the mechanism begins.
n = 4 # number of mechanisms
port = (1 - percent) / (n+1) # portion each piece takes up
top = interp(fx[:len(fx)//2],fy[:len(fy)//2])

h = [] # changing height for each pt of the mechanism
for j in range(n+2):
    xval = percent + j*port
    yval = top(xval)
    h.append(float(yval))

print(h)

#for i in range(len(fx)):
#    fx[i] -= 0.35
#    fy[i] += 0.09/2





import math as m
change = -5
# eventually input list of lenght of each peice
# model NACA 0016
t = 1*.09
print(port)
th1 = 0;
th2 = 86.5 + change;
r1 = t
r2 = m.sqrt(port**2 + (h[0] - h[1])**2)
print(r2)
r4 = r2
r3 = r1 - 2* (h[0] - h[1])
print(r3)

xt = []
yt = []
xb = []
yb = []
xt.append(0)
yt.append(0)
xb.append(r1)
yb.append(0)

for i in range(n):
    th2p = th2 - th1
    delt = (r1**2 + r2**2 - 2*r1*r2*m.cos(m.radians(th2p)))**0.5
    beta = m.acos((r1**2 + delt**2 - r2**2) / (2*r1*delt))
    psi = m.acos((r3**2 + delt**2 - r4**2) / (2*r3*delt))
    lambd = m.acos((r4**2 + delt**2 - r3**2) / (2*r4*delt))
    th3 = psi - (beta - m.radians(th1))
    th4 = m.pi - lambd + (beta + m.radians(th1))
    
    th3 = m.degrees(th3)
    th4 = m.degrees(th4)    
    
    xtop = r2*m.cos(m.radians(th2)) + xt[i]
    ytop = r2*m.sin(m.radians(th2)) + yt[i]
    xbot = r2*m.cos(m.radians(th2)) + r3*m.cos(m.radians(th3)) + xt[i]
    ybot = r2*m.sin(m.radians(th2)) + r3*m.sin(m.radians(th3)) + yb[i]
    xt.append(xtop)
    yt.append(ytop)
    xb.append(xbot)
    yb.append(ybot)
    
    r1 = r3
    th1 = th3
    th2 -= th1
    r2 = m.sqrt(port**2 + (h[i+1] - h[i+2])**2)
    r4 = r2 #- .02
    r3 = r1 - 2* (h[i+1] - h[i+2])
    
    if(i == n-1):
        print('ha!')
        a3 = r3/2
        b3 = port
        xp = r2*m.cos(m.radians(th2)) + a3*m.cos(m.radians(th3)) - b3*m.sin(m.radians(th3)) + xt[i]
        yp = r2*m.sin(m.radians(th2)) + a3*m.sin(m.radians(th3)) + b3*m.cos(m.radians(th3)) + yt[i]
        xt.append(xp)
        yt.append(yp)
        xb.append(xp)
        yb.append(yp)
        

from matplotlib import pyplot as plt
plt.axis('equal')
plt.plot(xt,yt)
plt.plot(xb,yb)
plt.plot(fy,fx)
for g in range(n):
    plt.plot([xt[g+1],xb[g+1]],[yt[g+1],yb[g+1]])
plt.plot()


lent = m.sqrt((xt[1]-xt[0])**2 + (yt[1]-yt[0])**2)
lenb = m.sqrt((xb[1]-xb[0])**2 + (yb[1]-yb[0])**2)
print(lent - lenb)

lent = m.sqrt((xt[2]-xt[1])**2 + (yt[2]-yt[1])**2)
lenb = m.sqrt((xb[2]-xb[1])**2 + (yb[2]-yb[1])**2)
print(lent - lenb)












    