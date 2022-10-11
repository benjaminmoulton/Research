#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 25 12:44:40 2019

@author: benjamin
"""

x_max = 0.6

points = 200

T = 0.12

A = 0.1

spine_width = 0.05
G = 0.5

from math import cos,pi,sin,acos,atan

f = 2*pi / T

x = []
for i in range(points):
    x.append( i * x_max / points )

y = []


y_s = []
y_sb = []
y_st = []

for j in range(points):
    
    yp = A * cos( f * x[j] )
    y.append(yp)
    
    multiplier = 0
    
    sign = -1
    
    if(yp <= 0):
        multiplier = cos( f * x[j] + pi )
        sign = -1
        
    else:
        multiplier = cos( f * x[j] )
        sign = 1
    
    newpt = A * cos( f * x[j] ) * multiplier + spine_width/2
    # A * cos( f * x[j] ) + spine_width + ( cos( 2*f * x[j] ) ) * spine_width / 4   #+ spine_width / 2
    
    
    newtoppt = A * cos( f * x[j] ) + ( -G * cos( 2 * f * x[j] ) + (1 + (G)) )  * spine_width / 2
    
    newbotpt = A * cos( f * x[j] ) - ( -G * cos( 2 * f * x[j] ) + (1 + (G)) )  * spine_width / 2
    
    
    
    y_s.append(newpt)
    y_st.append(newtoppt)
    y_sb.append(newbotpt)
    




x_at_zero_y = acos(0) / f
slope = - abs( - f* sin(f * x_at_zero_y ) )
theta = atan(slope)
h = ( spine_width / 2 ) / cos( theta )
print(h)






from matplotlib import pyplot as p

p.axis('equal')
p.plot(x,y) #,'bo')
#p.plot(x,y_s)
p.plot(x,y_st)
p.plot(x,y_sb)
p.show()



