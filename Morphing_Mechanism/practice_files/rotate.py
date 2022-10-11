#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 10 11:12:00 2019

@author: benjamin
"""

import NACApoints as nc
from matplotlib import pyplot as p
from math import sin,cos,radians,degrees,atan

L = 5

W = 5

t = 0.16

theta = degrees(atan(L/W))
print("theta is %.3f deg" % theta)

num = 150

xtop = []

for i in range(num):
    xval = i/num
    xtop.append(xval)

ytop = nc.getYht(t,xtop)[::-1]

ybot = ytop[::-1].copy()

for i in range(len(ybot)):
    ybot[i] *= -1

y = ytop + ybot
x = xtop + xtop[::-1]

#p.subplot(211)
#p.axis('equal')
#p.plot(x,y)
#p.title('regular')
#p.show()

angle = -radians(theta)

xr = []
yr = []

for i in range(len(x)):
    ptx = x[i]
    pty = y[i]
    
    ptxnew = ptx*cos(angle) + pty*sin(angle) # rotate by given angle
    ptynew = pty*cos(angle) - ptx*sin(angle)
    
    ptx = ptxnew
    pty = ptynew
    
    xr.append(ptx)
    yr.append(pty)


p.axis('equal')
p.plot(xr,yr)
p.title('rotated')
p.show()

max_xr = max(xr)
max_yr = max(yr)
print("max x at c=1 : %f" % max_xr)
print("max y at c=1 : %f" % max_yr)
#print(len(xr))
#print(xr[149:151])
#print(yr[149:151])

#for i in range(len(xr)):
#    xr[i] *= 2
#    yr[i] *= 2
#
#p.axis('equal')
#p.plot(xr,yr)
#p.title('rotated')
#p.show()

max_Twoxr = max(xr)
max_Twoyr = max(yr)

x_factor = W/max_xr
y_factor = L/max_yr

print("x factor %f" % x_factor)
print("y factor %f" % y_factor)


xx = []
xy = []
yx = []
yy = []

for i in range(len(xr)):
    xx.append(xr[i] * x_factor)
    xy.append(yr[i] * x_factor)
    
    yx.append(xr[i] * y_factor)
    yy.append(yr[i] * y_factor)

print(max(xx),max(xy))
print(max(yx),max(yy))


xnew = []
ynew = []

for i in range(len(xr)):
    xnew.append(xr[i] * x_factor)
    ynew.append(yr[i] * y_factor)

print("max xnew %f" % max(xnew))
print("max ynew %f" % max(ynew))

p.axis('equal')
p.plot(xnew,ynew)
p.plot([0,0,W,W,0],[0,L,L,0,0])
p.show()


xzero = []
yzero = []

for i in range(len(xnew)):
    ptx = xnew[i]
    pty = ynew[i]
    
    ptxnew = ptx*cos(-angle) + pty*sin(-angle) # rotate by given angle
    ptynew = pty*cos(-angle) - ptx*sin(-angle)
    
    ptx = ptxnew
    pty = ptynew
    
    xzero.append(ptx)
    yzero.append(pty)

print(xzero[149:151])
print(yzero[149:151])

midpt = (len(xzero)//2) -1
xpoint = ( xzero[midpt] + xzero[midpt + 1] ) / 2
ypoint = ( yzero[midpt] + yzero[midpt + 1] ) / 2

chord = ( xpoint**2 + ypoint**2 )**0.5

print("chord should be %f" % chord)

print("max xzero %f" % max(xzero))

p.axis('equal')
p.plot(xzero,yzero)
p.plot([0,max(xzero)],[0,0])
p.show()

print("chord is 6.5")





















