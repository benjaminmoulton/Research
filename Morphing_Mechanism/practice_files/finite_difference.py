#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 22 14:53:08 2019

@author: benjamin
"""

x_max = 0.5

points = 100

T = 0.12

A = 0.1

from math import cos,pi,sin

f = 2*pi / T

x = []
for i in range(points):
    x.append( i * x_max / points )

y = []

for j in range(points):
    y.append( A * cos( f * x[j] ) )

def FiniteDifference(x,y):
   
    f_prime = []
    
    for i in range(2, len(x) - 3):
        
        derivative = ( - y[i+2] + 8 * y[i+1] - 8 * y[i-1] + y[i-2] ) / 12
        f_prime.append( derivative )
    
    return f_prime
 

slope = []
for i in range(len(x)):
    slope_point = -A * f * sin( f * x[i])
    
    if(slope_point == 0):
        slope_point = 0.000000000000000000000000000000000001
    
    slope.append(slope_point)



def Shifter(x,y,deriv,thickness):
    
    tx = []
    ty = []
    bx = []
    by = []
    
    for i in range(2, len(x) - 3):
        invslope = -1 / deriv[i-2]
        xmid = ( x[i+1] - x[i] ) / 2  + x[i]
        ymid = ( y[i+1] - y[i] ) / 2  + y[i]
        xchange = ( thickness**2 / (1 + (invslope**2)) )**(1/2)
        ychange = ( thickness**2 / (1 + (invslope**-2)) )**(1/2)
        xtnew = xmid - (xchange if invslope < 0 else -xchange)
        ytnew = ymid + ychange
        ty.append(ytnew)
        tx.append(xtnew)
        xbnew = xmid + (xchange if invslope < 0 else -xchange)
        ybnew = ymid - ychange
        bx.append(xbnew)
        by.append(ybnew)
    
    
    return tx,ty,bx,by


thickness = 0.01
    
f_prime = FiniteDifference(x,y)

tx,ty,bx,by = Shifter(x,y,f_prime,thickness)
ntx,nty,nbx,nby = Shifter(x,y,slope,thickness)

from matplotlib import pyplot as plt

plt.axis('equal')
plt.plot(x,y)
#plt.plot(x[2:-3],f_prime)
plt.plot(tx,ty)
plt.plot(bx,by)
plt.show()



plt.axis('equal')
plt.plot(x,y)
#plt.plot(x[2:-3],f_prime)
plt.plot(ntx,nty)
plt.plot(nbx,nby)
plt.show()




