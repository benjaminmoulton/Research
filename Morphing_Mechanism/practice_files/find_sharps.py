#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 21 14:57:41 2019

@author: benjamin
"""
class point:
    x = 0
    y = 0
    coord = '( %.4f, %.4f )' % (x,y)

P1 = (2,50)

#class point(float,float):
#    print('banna')



import gcode as gc
from matplotlib import pyplot as p
from math import acos, degrees

x,y = gc.NACAFileReader('/home/benjamin/Desktop/PythonCode/aerolab/sharp_curve_points.txt')
u,v = gc.NACAFileReader('/home/benjamin/Desktop/PythonCode/aerolab/sharp_curve_points_a.txt')
p.axis('equal')
p.plot(x,y)
p.plot(u,v)

sharp = []

for i in range(len(x)-3):
    p2 = ( x[i],   y[i]   )
    p1 = ( x[i+1], y[i+1] )
    p3 = ( x[i+2], y[i+2] )
    
#    p2p1slope = ( p1[1] - p2[1] ) / ( p1[0] - p2[0] )
#    p1p3slope = ( p3[1] - p2[1] ) / ( p3[0] - p2[0] )
    
    a = (p1[0] - p2[0], p1[1] - p2[1])
    b = (p1[0] - p3[0], p1[1] - p3[1])
    
    a_magnitude = ( a[0]**2 + a[1]**2 )**0.5
    b_magnitude = ( b[0]**2 + b[1]**2 )**0.5
    
#    print(b_magnitude)
    
    a_dot_b = a[0] * b[0] + a[1] * b[1]
    
    
    theta_rad = acos( a_dot_b / ( a_magnitude * b_magnitude ) )
    theta = degrees( theta_rad )
    if(theta < 90):
        sharp.append(p1)



for j in range(len(sharp)):
    p.plot(sharp[j][0],sharp[j][1],'bo')

sharp2 = []

for i in range(len(x)-3):
    p2 = ( u[i],   v[i]   )
    p1 = ( u[i+1], v[i+1] )
    p3 = ( u[i+2], v[i+2] )
    
#    p2p1slope = ( p1[1] - p2[1] ) / ( p1[0] - p2[0] )
#    p1p3slope = ( p3[1] - p2[1] ) / ( p3[0] - p2[0] )
    
    a = (p1[0] - p2[0], p1[1] - p2[1])
    b = (p1[0] - p3[0], p1[1] - p3[1])
    
    a_magnitude = ( a[0]**2 + a[1]**2 )**0.5
    b_magnitude = ( b[0]**2 + b[1]**2 )**0.5
    
#    print(b_magnitude)
    
    a_dot_b = a[0] * b[0] + a[1] * b[1]
    
    
    theta_rad = acos( a_dot_b / ( a_magnitude * b_magnitude ) )
    theta = degrees( theta_rad )
    if(theta < 90):
        sharp2.append(p1)


for j in range(len(sharp)):
    p.plot(sharp2[j][0],sharp2[j][1],'bo')

p.show()



thck = 0.24
spineWidth = 0.1

d = spineWidth * thck / 2








































P1 = (2,50)


g = point()
g.x = 1
g.y = 4
g.coord
print(g.coord)















