#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 10 16:24:26 2018

@author: benjamin
"""
from math import radians, cos, sin,tan



##this function creates a square
#def square (height = 4):
#    
#    # the x values are created
#    x = [0,0,height,height,0]
#    
#    # the y values are created
#    y = [0,height,height,0,0]
#    
#    # the function returns the x and y values for the square to the user
#    return x,y



##this function creates an xy circle centered at (x,y) = (r,r)
#def circle (radius = 4,n = 100): # defines circle who starts at rightmost of circle, moves ccw
#    x = []
#    y = []
#    
#    # circle begins at the most right
#    x.append(2*radius)
#    
#    # the circle begins with the y value at the far right of the circle
#    y.append(radius)
#    
#    # a value is determined for the given set of points to be how far apart the points are interspersed, in degrees
#    degStep = 360 / n
#    
#    # begins at zero
#    Theta = 0
#    
#    # an iterative process to find the x,y coordinates for each point along the circle
#    for i in range(n):
#        x.append(radius + radius*cos(radians(Theta)))
#        y.append(radius + radius*sin(radians(Theta)))
#        Theta += degStep
#    
#    # the x,y points are returned to the user
#    return x,y



#a function that has not yet been developed
def CenterOfMass():
    return 0,0



#a function that is still to be developed, but is in use.
def NACAcg(xvals,chord): # finds the cg of a given NACA 4 digit airfoil
    x = chord / 4
    x += min(xvals)
    #some code that finds the exact position of y of cg on wing using NACAfour number
    y = 0
    return x,y

#def getPolar(a,b): # given two points, (x,y) and (a,b), it finds the polar position of (a,b) with respect to (x,y) (theta in rads)
#    radius = sqrt(pow(a,2) + pow(b,2))
#    #write code that finds theta from a specific point.
#    if a != 0 and b != 0:
#        theta = atan(b/a)
#        if b > 0 and a < 0:
#            theta += pi / 2
#    elif a == 0 and b > 0:
#        theta = pi / 2
#    elif a == 0 and b < 0:
#        theta = (3 * pi) / 2
#    elif b == 0 and a < 0:
#        theta = pi
#    elif b == 0 and a > 0:
#        theta = 0
#    return radius, theta
#
#def getCartesian(radius,theta):
#    a = radius * cos(theta)
#    b = radius * sin(theta)
#    return a,b



   
#function that creates washout on a foilshape, given the x and y coordinates and the angle to be twisted
def giveWashout(x,y,twistangle = 10): 
    xc,yc = NACAcg(x,getChord(x))
    
    # the angle is made negative and changed to radians for mathematical purposes
    twist = radians(-twistangle)
    a = []
    b = []
    
    # an iterative process to take each point in x,y and rotate it about the center of chord
    for i in range(len(x)):
        pta = x[i] - xc # set about origin
        ptb = y[i] - yc
        
        ptanew = pta*cos(twist) + ptb*sin(twist) # rotate by given angle
        ptbnew = ptb*cos(twist) - pta*sin(twist)
        
        pta = xc + ptanew # return to position
        ptb = yc + ptbnew
        
        a.append(pta) # append to x and y lists respectively
        b.append(ptb)
    return a,b
 
    


#this function takes in a list, and rounds each value in the list for Gcode manipulation
def roundValuesG(x):
    for q in range(len(x)):
        x[q] = round(x[q],3)
    return x
    



def giveSweep(u,sweepAngle=0,span=1): # if there is sweep, this function returns the foilshape of the UV side shifted for the sweep
    #assuming degree is from front of wing to parallel from old wing root ##### we dont need y vals, right?
   
    sweptu = []
    shiftValue = tan(radians(sweepAngle)) * span # the shift is calculated
    for i in range(len(u)):
        sweptu.append(shiftValue + u[i]) # and implemented for all x values where + is back toward the tail
    return sweptu  # returns the airfoil that was swept



def giveDihedral(v,dihedralAngle=0,span=1): # if there is dihedral, this function returns the foilshape coords of the uv side shifted
    
    dihedralv = []
    shiftValue = cos(radians(dihedralAngle)) * span
    for i in range(len(v)):
        dihedralv.append(shiftValue + v[i])
    return dihedralv



#aligns tip quarter chord with root quarter chord
def alignQuarterChord(x,u):
    rchord = getChord(x)
    tchord = getChord(u)
    shiftval = (rchord/4) - (tchord/4)
    for e in range(len(u)):
        u[e] += shiftval
    return u



def getChord(x): # given the x values for a set of data points (supposing no twist) returns the chord length of the airfoil
    maxi = max(x)
    mini = min(x)
    chord = round(maxi - mini,2)
    return chord


#determines the thickness of a foil given its y values
def getThickness(x,y):
    from scipy.interpolate import interp1d as inter
    f = inter(x[int((len(y)/2)+1):],y[int((len(y)/2)+1):])
    ymax = max(y)
    minAtYmax = f(x[y.index(ymax)])
    t = round(max(y)-minAtYmax,2)
    return t


#this function takes a given set of points for the x,y and u,v shapes of
# a foam wing, and determines what points the CNC hotwire machine should follow
# in cutting out the wing given different shapes,offsets, and the wire length
def fixAxisPoints(x,y,u,v,span,offset,wireLength):
    #offset on the left side of a right wing
    rightSpace = wireLength - span - offset
    
    Ax = []
    Ay = []
    Au = []
    Av = []
    
    # an iterative procedure is followed to determine where the points should be
    # based on the slope in the x and y directions between the two sets of points
    for r in range(len(x)):
        slopex = (u[r] - x[r]) / span
        slopey = (v[r] - y[r]) / span
        Ax.append(x[r] - slopex*offset)
        Ay.append(y[r] - slopey*offset)
        Au.append(u[r] + slopex*rightSpace)
        Av.append(v[r] + slopey*rightSpace)
    
    return Ax,Ay,Au,Av

#this function shifts the y and v data up the axis so it will fit in the machine constraints
def shiftUp(y,v,height):
    maxheight = max(max(y),max(v))
    shift = height - maxheight - 0.125
    
    for i in range(len(y)):
        y[i] += shift
        v[i] += shift
    return y,v

#this function shifts the x and u data right on axis so it fits in the machine
def shiftRight(x,u):
    shift = 1#inch
    
    for j in range(len(x)):
        x[j] += shift
        u[j] += shift
    return x,u
    


















