#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 18 11:32:13 2019

@author: benjamin
"""


def F(x,A,f):
    from math import cos
    F = 0
#    F = A * cos(f * x) + f
#    print(cos(f*x))
    return F

def A(x,a,f,rangeA,sin,pi):
    amp = rangeA / 4
    A = 0
    
    A = amp * sin( f*(x)  + pi/2)
    
    return A


n = 5
points_per_period = 30
T = 0.5

#print(vp)

z = int(n * points_per_period)
#print(z)


x_max = 1.2

x = []

for i in range(z):
    x.append(x_max * i / (z) )

y = []
yreg = []
yamp = []


from math import cos,pi,acos,sin,atan,degrees

f = 2*pi / T

A_f_change = 0.05*f

spine_width = 0.05


x_zero = acos(0) / f
#print(x_zero)
slope = -0.25*f*sin(f*x_zero)
theta = atan(slope)

Amax = 0.5 * spine_width / cos(theta)
print(Amax)
Amin = spine_width / 2
rangeA = Amax - Amin



for i in range(len(x)):
    y.append( 0.25 * cos( f  * x[i]) )
    F.new = F(x[i],A_f_change,f )
#    print(F.new)
    yreg.append( 0.25 * cos( f  * x[i]) ) #(  F.new * x[i]) )# + spine_width)
    
    
    
    yamp.append(( 0.25 ) * cos( f * x[i]) + A(x[i],0.25,f,rangeA,sin,pi) + spine_width/2 )

from matplotlib import pyplot as p

p.axis('equal')
p.plot(x,y) #,'bo')
#p.plot(x,yreg)
p.plot(x,yamp)
p.show()

