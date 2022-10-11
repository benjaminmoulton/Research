#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 18 15:27:29 2018

@author: benjamin
G Code Creator

"""
# receive params
# root airfoil, tip airfoil, twist, root chord, tip chord, sweep, root kerf, tip kerf, span, invert?
import math
import numpy

def main():
#    prms = getParams()
    
    # print(prms[3]) check to make sure can access prms
#    if(prms[9] == False):
#        print("It is false")
#    else:
#        print("It is true")
    
    rtairfoil,tpairfoil,twist,rtchord,tpchord,sweep,rtkerf,tpkerf,span,fishbone,isLeftWing = getParams()
    
    _,_,root = rtairfoil.partition(' ')
    _,_,tip = tpairfoil.partition(' ')
    rtfoil = int(root)
    tpfoil = int(tip)
    
    print(rtfoil)
    
    c = rtchord
    n = 20
    
    
    m,p,xx,decimal_along_camber,yvals = getCosineClustering(c,n,rtfoil)
    
    print(decimal_along_camber)
    
    
    # function to create file and get shapes of airfoils
    # output as seen below
    
    # if rtfoil == tpfoil, then tfoil = same as rfoil. else tfoil = own data set
    
    rfoil = [[0,1,2,3,3,2,1,0,0],
             [0,0,0,0,1,1,1,1,0]]
    
    tfoil = [[0,1,2,3,3,2,1,2,2],
             [0,0,0,0,1,1,1,1,0]]
    
    # if this wing has sweep, call defineSweep function, which will add sweep increments to the tip foil
    if(sweep != 0): 
        tfoil = defineSweep(sweep,tfoil,span)
        
    # if this wing is the left wing, call swapWings function, and cut out "oppositely"
    if(isLeftWing):
        rfoil, tfoil = swapWings(rfoil,tfoil)
        
    
    
    
    return 0



def getParams():
    rtairfoil = "NACA 2410" # 0
    tpairfoil = "NACA 2410" # 1
    twist = 0               # 2
    rtchord = 1             # 3
    tpchord = 1             # 4
    sweep = 0               # 5 # This is the degree value, where + is CW(right wing) toward the tail.
    rtkerf = 0              # 6
    tpkerf = 0              # 7
    span = 1                # 8 # semispan of wing (ie. left or right only)
    fishbone = False        # 9
    isLeftWing = False # Whether the wing to cut is the left wing or not(right)
    params = [rtairfoil,tpairfoil,twist,rtchord,tpchord,sweep,rtkerf,tpkerf,span,fishbone,isLeftWing]
    return params

#def defineSweep(sweepAngle,airfoil,span): # if there is sweep, this function returns the foilshape of the ZA side shifted for the sweep
#    
#    shiftValue = math.tan(math.radians(sweepAngle)) * span # the shift is calculated
#    for i in range(len(airfoil[0])):
#        airfoil[0][i] += shiftValue # and implemented for all x values where + is back toward the tail
#        
#    return airfoil # returns the airfoil that was swept

def swapWings(xy,za):
    temp = xy
    xy = za
    za = temp
    return xy,za

#def GcodeFileWriter(airfoilData):  # writes a given data file into Gcode for CNC hotwiring
#    f = open("gcode.txt","w+")
#
#    for i in range(len(airfoilData[0])):
#        f.write("G01 X" + str(airfoilData[0][i]) + " Y" + str(airfoilData[1][i]) 
#                + " Z" + str(airfoilData[0][i]) + " A" + str(airfoilData[1][i])+ "\n")
#    f.write("G2")
#    f.close()
#
#    return 0

def getNACAfour (nacaFour): # used to find values for NACA four digit airfoil
    #               and returns the values to the user for computation
#    nacaFour = input ("What NACA 4-digit airfoil would you like to use? ")
#    nacaFour = int(nacaFour)
    #returns values for camber, location, and chord to length ratio
    camber = int(nacaFour / 1000 )
    camber *= 100 #max camber returns decimal value of % camber
    
    camberLoc = 10 * int((nacaFour -(camber * 10)) / 100)
    # location of maximum camber (decimal value of % along wing)
    
    chordToLengthRatio = int(nacaFour - (camber * 10) - (camberLoc * 10))
    # ratio of chord to length ie, it is xx% thick as it is long
    m = camber * 0.0001
    p = camberLoc / 100
    xx = chordToLengthRatio / 100
    
    return m,p,xx

def getCosineClustering(c,n,foil):
    m,p,xx = getNACAfour(foil)
    decimal_along_camber = []
    yvals = []
    dtheta = 2* numpy.pi / (n - 1)
    for i in range(1,int(n/2)+1):
        decimal_along_camber.append( 0.5 * (1. - math.cos((i - 0.5) * dtheta)) )
        yvals.append(meanCamberLine(decimal_along_camber[i-1],c,m,p))
#    print(decimal_along_camber)
    return m,p,xx,decimal_along_camber,yvals

# this function returns the y value for the mean camber line for the airfoil created 
def meanCamberLine (x,c,m,p):# where x is on x axis and c is chord length
                            # ie x/c returns value from 0 to 1
    yc = x
    if x >= 0 and x < p*c:
        yc = (m/pow(p,2)) * (2*p*(x/c) - pow((x/c),2))
    elif x >= p*c and x <= c:
        yc = (m / pow((1-p),2)) * ((1 - 2*p) + 2*p*(x/c) - pow((x/c),2))
    else:
        exit
    return yc


main()



"""

Do we want this to be able to do dihedral? twist?

I think kerf is something we will have to calculate with cosine clustering

figure out kerf,



Read in multiple file types, ie anything
"""



































