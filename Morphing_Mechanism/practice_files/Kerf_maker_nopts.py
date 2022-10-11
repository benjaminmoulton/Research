#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 18 16:48:03 2019

@author: benjamin
"""

import point as pt

#this function determines the outline of an airfoil shape (even fishbone) given 
# its data points and kerf value. shifts the foil shape to include this kerf value 
# for the whole geometry
def kerfMaker(airfoil,kerf):  
    
    """ """
    x = airfoil.x(); y = airfoil.y()
    
    """ """
    
    xone = x[:(len(x)//2) + 1]
    xtwo = x[(len(x)//2) + 1:]
    yone = y[:(len(x)//2) + 1]
    ytwo = y[(len(x)//2) + 1:]
    
    kerfedy = []
    kerfedx = []
    
    for j in range(len(yone)-1):
        
        rise = (yone[j+1] - yone[j])
        run = (xone[j+1] - xone[j])
        
        if(run == 0):
            
            if(rise < 0):
                
                yint = yone[j] + kerf
                xint = xone[j] - kerf
                kerfedy.append(yint)
                kerfedx.append(xint)
                nexx = xone[j+1] - kerf
                nexy = yone[j+1] + kerf
                kerfedy.append(nexy)
                kerfedx.append(nexx) 
                
            else: # rise > 0
                
                nexx = xone[j] + kerf
                nexy = yone[j] + kerf
                kerfedy.append(nexy)
                kerfedx.append(nexx)
                yint = yone[j+1] + kerf
                xint = xone[j+1] + kerf
                kerfedy.append(yint)
                kerfedx.append(xint)
        elif(rise != 0 and run != 0):
            
            slope = rise/run
            inslope = 1 / (-slope)
            xmid = ((xone[j+1] - xone[j]) / 2.0 ) + xone[j]
            ymid = ((yone[j+1] - yone[j]) / 2.0 ) + yone[j]
            xchange = ( kerf**2.0 / (1.0 + (inslope**2.0)) )**(0.5)
            ychange = ( kerf**2.0 / (1.0 + (inslope**(-2.0))) )**(0.5)
        
            xtnew = xmid - (xchange if inslope < 0.0 else -xchange)
            ytnew = ymid + ychange
            kerfedy.append(ytnew)
            kerfedx.append(xtnew)
    
    for j in range(len(ytwo)-1):
        
        rise = (ytwo[j+1] - ytwo[j])
        run = (xtwo[j+1] - xtwo[j])
        
        if(run == 0):
            
            if(rise > 0):
                
                nexx = xtwo[j] + kerf
                nexy = ytwo[j] - kerf
                kerfedy.append(nexy)
                kerfedx.append(nexx)
                yint = ytwo[j+1] - kerf
                xint = xtwo[j+1] + kerf
                kerfedy.append(yint)
                kerfedx.append(xint)
                
            else: # rise < 0
                yint = ytwo[j] - kerf
                xint = xtwo[j] - kerf
                kerfedy.append(yint)
                kerfedx.append(xint)
                nexx = xtwo[j+1] - kerf
                nexy = ytwo[j+1] - kerf
                kerfedy.append(nexy)
                kerfedx.append(nexx) 
                
        elif(rise != 0 and run != 0):
            slope = rise/run
            inslope = 1 / (-slope)
            xmid = ((xtwo[j+1] - xtwo[j]) / 2.0 ) + xtwo[j]
            ymid = ((ytwo[j+1] - ytwo[j]) / 2.0 ) + ytwo[j]
            xchange = ( kerf**2.0 / (1.0 + (inslope**2.0)) )**(0.5)
            ychange = ( kerf**2.0 / (1.0 + (inslope**(-2.0))) )**(0.5)
            
            xtnew = xmid + (xchange if inslope < 0 else -xchange)
            ytnew = ymid - ychange
            kerfedy.append(ytnew)
            kerfedx.append(xtnew)
    
    firstslope = (kerfedy[0] - kerfedy[1])  /  (kerfedx[0] - kerfedx[1])
    firstrun = - kerfedy[0] / firstslope
    firstx = kerfedx[0] + firstrun
    firsty = 0.0
    
    lastslope = (kerfedy[len(kerfedy)-1] - kerfedy[len(kerfedy)-2])  /  (kerfedx[len(kerfedy)-1] - kerfedx[len(kerfedy)-2])
    lastrun = - kerfedy[len(kerfedy)-1] / lastslope
    lastx = kerfedx[len(kerfedx)-1] + lastrun
    lasty = 0.0
    
    kerfedx.insert(0,firstx)
    kerfedy.insert(0,firsty)
    
    kerfedx.append(lastx)
    kerfedy.append(lasty)
    
    kerfed = pt.points(kerfedx,kerfedy)
    return kerfed
