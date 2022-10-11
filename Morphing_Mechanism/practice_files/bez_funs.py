#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 29 15:00:17 2019

@author: benjamin
"""

def addSpine(x,y,dxlist,dylist,d):
    
    # the spine coordinate lists are initialized
    spinex = []
    
    spiney = []
    
    # number of points for iteration is determined
    
    n = len(x)
    
    
    # d values will be set to a list for use when determining self intersecting points
    
    d_vals = []
    
    # a for loop runs through each x y coordinated given, and determines
    #  the spined points on either side of the middle wave
    
    for i in range(n):
        
        # the individual points are set for calculation
        
        xval = x[i]
        
        yval = y[i]
    
        
        
        # if is T Local calculation, d is set for each calculation
        
#        if(isT):
#            
#            d = sw * Yht(t,xval)
        
        # for use in later calculation, the d values are stored to a list
        
        d_vals.append(abs(d))
        
        # the dx and dy values are determined for future use
        
        dx = dxlist[i]
        
        dy = dylist[i]
        
        # the shift values are found for use in determining the d coordinates
        
        xshift = d * dy / ( ( dx**2 + dy**2 )**0.5 )
        
        yshift = d * dx / ( ( dx**2 + dy**2 )**0.5 )
        
        # the d coordinates are found
        
        xd = xval + xshift
        
        yd = yval - yshift
        
        # the coordinates are added to the spinex and spiney lists
        
        spinex.append(xd)
        
        spiney.append(yd)
    
    # 2 lists are initialized to hold the indices where self intersection occurs
    #  on the spine wave
    
#    spinex, spiney = deintersect(x,y,spinex,spiney,n,pphp,d,d_vals,minus)
    
#    ##############################################################################3333
##    print(max(spinex))
#    
#    count = 0
#    for i in range(len(spinex)):
#        if(spinex[i - count] <= min(x)):
#            spinex.pop(i - count)
#            spiney.pop(i - count)
#            count += 1
#        elif(spinex[i - count] >= max(x)):
#            spinex.pop(i - count)
#            spiney.pop(i - count)
#            count += 1
#    ######################################################333333333333333333333333333
#    # the point is added to the start and end of the spine so that the curve meets correctly

    
#    x_begin = x[0]
#    
#    y_begin = y[0] - d
#    
#    spinex.insert(0,x_begin)
#    
#    spiney.insert(0,y_begin)
#    
#    x_end = x[len(x) - 1]
#    
#    y_end = y[len(y) - 1] - d
#    
#    spinex.append(x_end)
#    
#    spiney.append(y_end)
#    
#    # the spinex and spiney lists are then returned to the user
    
    return spinex, spiney



















a = [
     [1,2,3],
     [2,3,4]     
     ]

x,y = a
print(x)
print(y)

z = x,y

print(z)

a = [
     (1,2),
     (2,3)
     ]

x,y = a
print(x)
print(y)

z = x,y

print(z)
















