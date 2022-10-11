#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 17 14:28:31 2019

@author: benjamin

this document will make an airfoil skin shape.
given a desired LE%, TE %, and a skin thickness,
it will determine the middle points to be cut out for the hot wire cutter
such that the foil shape will resemble a 
"""

from matplotlib import pyplot as plt
import gcode as gc



x,y = gc.NACAFileReader('/home/benjamin/Desktop/foilshapes/NACA 0016.txt')

plt.subplot(211)
plt.axis('equal')
plt.plot(x,y)
plt.show()

def skinFoil(x,y,c,t,LE,TE,skin_thickness,isFlat,gap = 0.0):
    # this function assumes the airfoil has no camber...
    #  skin_thickness is a % of the thickness
    from scipy.interpolate import interp1d as itp
    
    sk = skin_thickness
    
    x_inner = []
    
    y_inner = []
    
    x_dechorded = []
    
    y_dechorded = []
    
    # remove chord
    
    for i in range(len(x)):
        
        xval = x[i] / c
        
        yval = y[i] / c
        
        x_dechorded.append(xval)
        
        y_dechorded.append(yval)
    
    a = x_dechorded[:len(x)//2]
    
    b = y_dechorded[:len(x)//2]
    
    Outer = itp(a,b)
    
    TE_works = False
    
    if(Outer(TE) < sk*t):
        
        TE_works = False
        
    else:
        
        TE_works = True
    
    LE_works = False
    
    if(Outer(LE) < sk*t):
        
        LE_works = False
        
    else:
        
        LE_works = True
    
    # find starting locations if LE and or TE don't work
      # maybe split around tmax
    
#    if(not TE_works):
#        
    
    # add beginning values
    
    x_inner.append(TE)
    y_inner.append(0)
    
    if( not isFlat ):
        
        xval = TE
        yval = -Outer(TE) + sk * t
        x_inner.append(xval)
        y_inner.append(yval)
        
        for i in range(len(x_dechorded)-1,-1,-1):
            
            xval = x_dechorded[i]
            
            yval = y_dechorded[i]
            
            if(yval >= 0):
                
                if(xval <= TE and xval >= LE):
                    ixval = xval
                    iyval = yval - sk * t
                    
                    x_inner.append(ixval)
                    
                    y_inner.append(iyval)
            
            elif(yval < 0):
                
                if(xval <= TE and xval >= LE):
                    ixval = xval
                    iyval = yval + sk * t
                    
                    x_inner.append(ixval)
                    
                    y_inner.append(iyval)
        
        xval = TE
        yval = Outer(TE) - sk * t
        x_inner.append(xval)
        y_inner.append(yval)
    
    elif( isFlat ):
        
        x0 = TE
        y0 = gap*t/2
        x1 = LE
        y1 = gap*t/2
        x2 = LE
        y2 = -gap*t/2
        x3 = TE
        y3 = -gap*t/2
        
        x_inner.append(x0)
        
        y_inner.append(y0)
        
        x_inner.append(x1)
        
        y_inner.append(y1)
        
        x_inner.append(x2)
        
        y_inner.append(y2)
        
        x_inner.append(x3)
        
        y_inner.append(y3)
    
    # add ending values
    
    x_inner.append(TE)
    y_inner.append(0)
    
    
    # merge lists
    
    skx = x + x_inner
    
    sky = y + y_inner
        
        
    
    
    
    return skx,sky

ch = 1
th = 0.16
le = 0.2
te = 0.6
sk = 0.3
isFlat = False
gap = 0.5

sx,sy =skinFoil(x,y,ch,th,le,te,sk,isFlat,gap)

plt.subplot(212)
plt.axis('equal')
plt.plot(sx,sy)
plt.show()



























