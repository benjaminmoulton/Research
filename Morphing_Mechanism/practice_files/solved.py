#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May  3 11:29:56 2019

@author: benjamin
"""

from math import sin,cos,pi
from matplotlib import pyplot as p


n = 150 # points
A = 0.2 # amplitude
t = 0.12 # period
f = 2 * pi / t # frequency
spine_width = 0.05
d = spine_width / 2

xmin_new = 0
xmax_new = 0.5



x_new = []
y_new = []
xb_new = []
yb_new = []
xt_new = []
yt_new = []

for i in range(n):
    
    xval_new = ( i / n ) * (xmax_new - xmin_new) + xmin_new
    yval_new = A * cos(f * xval_new)
    x_new.append(xval_new)
    y_new.append(yval_new)
    
    xbshift = ( d * (-f * A * sin(f * xval_new)) / ( ( 1 + ( -f*A*sin(f * xval_new) )**2 )**0.5 )  )
    ybshift = (  d / ( ( 1 + ( -f*A*sin(f * xval_new) )**2 )**0.5 )  )
    
    xdval_new = xval_new + xbshift
    ydval_new = yval_new - ybshift
    xb_new.append(xdval_new)
    yb_new.append(ydval_new)
        
    xtval_new = xval_new - xbshift
    ytval_new = yval_new + ybshift
    xt_new.append(xtval_new)
    yt_new.append(ytval_new)



num = int( n / ( 2 * (xmax_new - xmin_new) / t ) ) # num points per half period

b_ind_to_rem = [] # list initialized for the indices to be removed from the bottom list of points
t_ind_to_rem = []

for a in range(n):
    for b in range(-num,num):
        if(( a + b >= 0 ) and ( a + b < n )):
            dist = ( ( xb_new[a] - x_new[a+b] )**2 + ( yb_new[a] - y_new[a+b] )**2 )**0.5
            if(dist < d*0.99):
                b_ind_to_rem.append(a+b)
                
            dist = ( ( xt_new[a] - x_new[a+b] )**2 + ( yt_new[a] - y_new[a+b] )**2 )**0.5
            if(dist < d*0.99):
                t_ind_to_rem.append(a)

b_indices = list(set(b_ind_to_rem))
b_indices.sort()

t_indices = list(set(t_ind_to_rem))
t_indices.sort()

bx,by,tx,ty = [],[],[],[]

for i in range(len(b_indices)-1,-1,-1):
    bx.append(xb_new[b_indices[i]])
    by.append(yb_new[b_indices[i]])
    del xb_new[b_indices[i]]
    del yb_new[b_indices[i]]
for i in range(len(t_indices)-1,-1,-1):
    tx.append(xt_new[t_indices[i]])
    ty.append(yt_new[t_indices[i]])
    del xt_new[t_indices[i]]
    del yt_new[t_indices[i]]    




p.axis('equal')
p.plot(x_new,y_new)
p.plot(xb_new,yb_new)
p.plot(xt_new,yt_new)
p.plot(bx,by,'bo')
p.plot(tx,ty,'ro')
p.show()





