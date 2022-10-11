#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  9 10:01:26 2019

@author: benjamin

isOscispine = True
    
    isDefaultOscispine = False
    
    startPerc = .30645 # % of chord
    
    endPerc = .87 # % of chord
    
    period = 0.0887 # % of chord
    
    spineWidth = 0.1 # % of thickness
    
    fillet_radius = 0.5 # a percentage of spine width
    
    firstUp = True
    
    isTlocalOsc = False

"""

from matplotlib import pyplot as plt
import gcode as gc



x,y = gc.NACAFileReader('/home/benjamin/Desktop/foilshapes/NACA 0016.txt')
#plt.axis('equal')
#plt.plot(x,y)
#plt.show()

def OsciSpine(x,y,chord,thickness,startPerc = .0,endPerc = .88,period = 0.09,spineWidth = 0.0,fillet_rad = 0.5,firstUp = True,isTlocal = False): # period is a % of the chordlength 
    from scipy.interpolate import interp1d as interp
    if(startPerc == 0):
        startPerc = x[y.index(max(y))] / chord
    if(spineWidth == 0):
        spineWidth = 0.2
    elif(spineWidth == 0 and isTlocal == True):
        spineWidth = 1
    sign = -1
    if(firstUp == True):
        sign = 1
    
    from numpy import cos, pi
    x_to_interpolate = x[:len(x)//2]
    y_to_interpolate = y[:len(x)//2]
    topInterp = interp(x_to_interpolate[::-1],y_to_interpolate[::-1])
    ox = []
    oy = []
    
    
    
    
    n = 300
    i = []
    ival = startPerc
    shift = (endPerc - startPerc) / n
    for w in range(n):
        i.append(ival + shift*w)
    
    k = 0
    frequency_over_chord = ((2*pi)/(period*chord))
    starting_point = startPerc*chord
    spine_thickness = spineWidth*thickness/2
    for val in i: # while loop that ends where x vals are greater than int(endPerc-startPerc / (period)) + (1/2) * endPerc-startPerc / (period) + startPerc  * chord
        val = i[k]
        ox.append(val*chord)
        interpolated_point = topInterp(ox[k])
        if(isTlocal):
            spine_thickness = spineWidth*interpolated_point
        if(ox[k] < (startPerc + period/2)*chord):
            oy.append(sign*((( interpolated_point/2 )*cos(frequency_over_chord*(ox[k] - starting_point) + pi)) + interpolated_point/2 - spine_thickness))
        elif((ox[k]/chord - startPerc)/period > int((endPerc - startPerc)/period)):
            if((ox[k]/chord - startPerc)/period < int((endPerc - startPerc)/period)+ 0.5):
                oy.append(sign*(interpolated_point - spine_thickness)*-cos(frequency_over_chord*(ox[k] - starting_point)))
            else:
                oy.append(sign*((( interpolated_point/2 )*cos(frequency_over_chord*(ox[k] - starting_point) + pi)) + interpolated_point/2 - spine_thickness))
        else:
            oy.append(sign*(interpolated_point - spine_thickness)*-cos(frequency_over_chord*(ox[k] - starting_point))) # create OsciSpine in such a way that the ending opposite
        k += 1
    
    ##
    jx = ox
    jy = oy
    ##
    
    
        #################################################################TURN THIS TEXT BLOCK INTO A SEPARATE FUNCTION WHICH ###################################
        ################################################ GIVEN A SHIFT VALUE MAKES TWO SETS OF DATA THAT ARE SHIFTED TO COVER THE LINE GIVEN####################
    
    topy = []
    topx = []
    botx = []
    boty = []
    for j in range(len(oy)-1):
        slope = (oy[j+1] - oy[j]) / (ox[j+1] - ox[j])
        inverse_slope = 1 / (-slope)
        xmid = ((ox[j+1] - ox[j]) / 2 ) + ox[j]
        if(isTlocal):
            spine_thickness = spineWidth*topInterp(xmid)
        ymid = ((oy[j+1] - oy[j]) / 2 ) + oy[j]
        xchange = ( spine_thickness**2 / (1 + (inverse_slope**2)) )**(1/2)
        ychange = ( spine_thickness**2 / (1 + (inverse_slope**-2)) )**(1/2)
        xtnew = xmid - (xchange if inverse_slope < 0 else -xchange)
        ytnew = ymid + ychange
        topy.append(ytnew)
        topx.append(xtnew)
        xbnew = xmid + (xchange if inverse_slope < 0 else -xchange)
        ybnew = ymid - ychange
        botx.append(xbnew)
        boty.append(ybnew)
        
    ###########################################################################################################################################################
    
    x_interpolate = x[len(x)//2:]
    y_interpolate = y[len(x)//2:]
    bottomInterp = interp(x_interpolate,y_interpolate)
    
    osx = []
    osy = []
    
    i = 0
    k = 0
    topxccw = topx[::-1]
    topyccw = topy[::-1]    
    while(i < len(x_to_interpolate)):
        if(topxccw[0] < x_to_interpolate[i] or (k != 0 and topxccw[len(topxccw)-1] > x_to_interpolate[i])):
            osx.append(x_to_interpolate[i])
            osy.append(y_to_interpolate[i])
        elif(topxccw[0] > x_to_interpolate[i] and k == 0):
            while(k < len(topxccw)):
                if(k == 0):
                    osx.append(topxccw[0])
                    osy.append(topInterp(topxccw[0]))
                osx.append(topxccw[k])
                osy.append(topyccw[k])
                if(k == len(topxccw)-1):
                    osx.append(topxccw[len(topxccw)-1])
                    osy.append(topInterp(topxccw[len(topxccw)-1]))
                k += 1
        i += 1
    
    i = 0
    k = 0
    
    # following 2 lines : add point at start of bottom squigley that matches the x location of the top last squigley point
    botx.insert(0,topxccw[len(topxccw)-1])
    boty.insert(0, oy[0] - ( osy[len(osy)-1] - oy[0]) )
    
    while(i < len(x_interpolate)):
        if(botx[0] > x_interpolate[i] or (k != 0 and botx[len(topxccw)-1] < x_interpolate[i])):            
            osx.append(x_interpolate[i])
            osy.append(y_interpolate[i])
        elif(botx[0] < x_interpolate[i] and k == 0):
            while(k < len(botx)):
                if(k == 0):
                    osx.append(botx[0])
                    osy.append(bottomInterp(botx[0]))
                osx.append(botx[k])
                osy.append(boty[k])
                if(k == len(botx)-1):
                    
                    # following 2 lines : add point at bottom of squigley that matches the x location of the top first squigley point
                    osx.append(topxccw[0])
                    osy.append(oy[len(oy)-1] - ( topyccw[0] - oy[len(oy)-1]))
                    
                    osx.append(osx[len(osx)-1])
                    osy.append(bottomInterp(osx[len(osx)-1]))
                k += 1            
        i += 1
    
    return osx,osy,jx,jy


c = 1
thck = 0.16
startPerc = 0.3
endPerc = .9
period = 0.09
spineWidth = 0.18
fillet_rad = 0.5
firstUp = True
isTlocal = False

ox,oy,jx,jy = OsciSpine(x,y,c,thck,startPerc,endPerc,period,spineWidth,fillet_rad,firstUp,isTlocal)

plt.axis('equal') #((0.85,1.0,-0.05,0.05))
plt.plot(ox,oy)
plt.plot(jx,jy)
plt.show()


















