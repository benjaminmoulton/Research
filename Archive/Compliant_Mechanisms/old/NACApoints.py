#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 11 15:16:31 2018

@author: benjamin
"""

def getNACAfour (foil): # used to find values for NACA four digit airfoil
    #               and returns the values to the user for computation
#    nacaFour = input ("What NACA 4-digit airfoil would you like to use? ")
#    nacaFour = int(nacaFour)
    #returns values for camber, location, and chord to length ratio
    _,_,nacaFour = foil.partition(' ')
    nacaFour = int(nacaFour)
    camber = int(nacaFour / 1000.0 )
    camber *= 100.0 #max camber returns decimal value of % camber  (max camber is _ % of chord)
    camber_location = 10.0 * int((nacaFour -(camber * 10.0)) / 100.0)
    # location of maximum camber (decimal value of % along wing)
    thickness_to_chord_ratio = int(nacaFour - (camber * 10.0) - (camber_location * 10.0))
    # ratio of chord to length ie, it is xx% thick as it is long
    m = camber * 0.0001
    p = camber_location / 100.0
    xx = thickness_to_chord_ratio / 100.0
    return m,p,xx

#given a certain number, this function calculates the locations
# for cosine clustering configuration from 0 to 1
# and returns a list of these numbers to the user (starting at 0 to 1)
def getCosineClustering(num):
    from math import cos, pi
    decimal_along_camber = []
    dtheta = 2.0* pi / (num - 1.0)
    for i in range(1,int(num/2)+1):
        decimal_along_camber.append(0.5 * (1. - cos((i - 0.5) * dtheta)) )
    return decimal_along_camber

#given the percent along the chord, the camber and the camber Location,
# this function finds the y values for the mean camber line and returns
# them to the user
def getMeanCamberLine(m,p,percentOfChord):
    camberLineYvalues = []
    for val in percentOfChord:
        yc = 0.0
        if val >= 0.0 and val <= p:
            yc = (m/(p**2.0)) * (2.0*p*val - (val**2.0))
        elif val > p and val <= 1:
            yc = (m/((1.0-p)**2.0)) * ((1.0 - 2.0*p) + 2.0*p*val - (val**2.0))
        camberLineYvalues.append(yc)
    return camberLineYvalues

#given the thickness to chord ratio, and the percent of chord
# this function returns a list of the y values (thickness function of
# the airfoil)
def getYht(t,x): #maxThickness, percentOfChord
    Yht = []
    for val in x:
        Yht.append(5*t*(0.2969*((val)**(1.0/2.0)) - 0.126*val - 0.3516*(val**2.0) + 
                        0.2843*(val**3.0) - 0.1036*(val**4.0)))
    return Yht

#given the camber, camber location, and a list of the percent of Chord
# values, this function returns a list of the dy/dx values for each value given
def getdydx(m,p,percentOfChord): #camber, camber location, percentofChord
    dydx = []
    for val in percentOfChord:
        dy = 0.0
        if val >= 0.0 and val <= p:
            dy = (2.0*m/(p**2.0)) * (p - val)
        elif val > p and val <= 1:
            dy = (2.0*m/((1.0-p)**2.0)) * (p - val)
        dydx.append(dy)
    return dydx

##this function takes the chord and the cosine clustering percentOfChord
## values and returns the xpoints
#def cosineToChord(x,chord):
#    xwithChord = []
#    for i in range(len(x)):                                        # not needed
#        xwithChord.append(x[i] * chord)
#    return xwithChord

#the function below takes all the information gathered thus far,
# and finds the x and y points of the foilshape. It begins at the
# trailing edge and follows the airfoilshape counter clockwise.
def getfoilShape(t,dydx,x,yc):
    from math import sqrt
    xpoints = []
    ypoints = []
    for j in range(len(x)-1,-1,-1):
        xpoints.append(x[j] - (t[j] * dydx[j]) / (2.0*sqrt(1.0 + (dydx[j]**2.0))) )
        ypoints.append(yc[j] + t[j] / (2.0*sqrt(1.0 + (dydx[j]**2.0))) )
    for k in range(0,len(x)):
        xpoints.append(x[k] + (t[k] * dydx[k]) / (2.0*sqrt(1.0 + (dydx[k]**2.0))) )
        ypoints.append(yc[k] - t[k] / (2.0*sqrt(1.0 + (dydx[k]**2.0))) )
    return xpoints,ypoints

#the function below takes a set of coordinates for a non-cambered wing and determines the 
# coordinates for a fishBAC design on the wing through cubic splines, and returns the 
# new coordinates to the user.
def FishboneMaker(x,y,xbone,ybone,x1,x2,y1,y2):
    fishx = []
    fishy = []
    nextBoneIndex = 0

    i = 0
    while(i < len(x1)):      
        if((xbone[nextBoneIndex] >= x1[i]) and nextBoneIndex < int(len(xbone)/2)):
            for h in range(4):
                fishx.append(xbone[nextBoneIndex])
                fishy.append(ybone[nextBoneIndex])
                nextBoneIndex += 1
        elif((nextBoneIndex != 0 and xbone[nextBoneIndex-1] > x1[i]) or nextBoneIndex == 0): # (xbone[nextBoneIndex] < x1[i]):
            fishx.append(x1[i])
            fishy.append(y1[i]) 
            i += 1
        else:
            i += 1
    turnpt = nextBoneIndex
    j = 0    
    while(j < len(x2)):
        if(nextBoneIndex < len(ybone) and (xbone[nextBoneIndex] <= x2[j])):
             for h in range(4):
                 fishx.append(xbone[nextBoneIndex])
                 fishy.append(ybone[nextBoneIndex])
                 nextBoneIndex += 1
        elif((xbone[nextBoneIndex-1] < x2[j]) or nextBoneIndex == turnpt): # (xbottombones[nextBoneIndex] > x2[j]):
            fishx.append(x2[j])
            fishy.append(y2[j])
            j += 1
        else:
            j += 1
    
    return fishx,fishy


#the function below gets the points to be used by the FishboneMaker
# as the locations for the fishbone points
# uses the FishboneMaker function to create the fish bone, and return the x and y values
# to the user
def FishBonePoints(x,y,thickness,chord,startperc = 0,endperc = 0.85,boneWidth = 0.0225,gapWidth = 0.0225,spineWidth = 0.25,isTlocal = False):
    from scipy.interpolate import interp1d as interp
    if(startperc == 0):
        startperc = x[y.index(max(y))] / chord
    xpost = endperc*chord
    topbonex = []
    while(xpost >= ((startperc + gapWidth + boneWidth)*chord)):
        topbonex.append(xpost)
        topbonex.append(xpost)
        xpost -= gapWidth*chord
        topbonex.append(xpost)
        topbonex.append(xpost)
        xpost -= boneWidth*chord
    
    bottombonex = list(reversed(topbonex))
    xone = x[:len(x)//2]
    xtwo = x[len(x)//2:]
    yone = y[:len(x)//2]
    ytwo = y[len(x)//2:]
    xonef = xone[::-1]
    yonef = yone[::-1]
    fone = interp(xonef,yonef)
    ftwo = interp(xtwo,ytwo)
    topboney = []
    bottomboney = []
    if(isTlocal == True):
        for g in range(0,len(topbonex),4):
            topboney.append(float(fone(topbonex[g])))
            topboney.append(spineWidth*float(fone(topbonex[g+1])) )#/ 2)
            topboney.append(spineWidth*float(fone(topbonex[g+2])) )#/ 2)
            topboney.append(float(fone(topbonex[g+3])))
            bottomboney.append(float(ftwo(bottombonex[g])))
            bottomboney.append(spineWidth*float(ftwo(bottombonex[g+1])) )#/ 2)
            bottomboney.append(spineWidth*float(ftwo(bottombonex[g+2])) )#/ 2)
            bottomboney.append(float(ftwo(bottombonex[g+3])))
    else:
        for g in range(0,len(topbonex),4):
            topboney.append(float(fone(topbonex[g])))
            topboney.append(spineWidth*thickness / 2.0)
            topboney.append(spineWidth*thickness / 2.0)
            topboney.append(float(fone(topbonex[g+3])))
            bottomboney.append(float(ftwo(bottombonex[g])))
            bottomboney.append(-spineWidth*thickness / 2.0)
            bottomboney.append(-spineWidth*thickness / 2.0)
            bottomboney.append(float(ftwo(bottombonex[g+3])))
    bonex = []
    bonex.extend(topbonex)
    bonex.extend(bottombonex)
    boney = []
    boney.extend(topboney)
    boney.extend(bottomboney)
    fx,fy = FishboneMaker(x,y,bonex,boney,xone,xtwo,yone,ytwo)
    
    return fx,fy

#the function below takes a set of x and y values for an airfoil and interpolates
# the mean line of camber based on the location of the y values assumed to associate
# with each x position. The function returns a new set of x values and yvalues to 
# the user that follow the line of mean camber for the airfoil.
def makeCamberLine(xval,yval):
    halflength = int(len(yval)/2.0)
    newyvals = []
    newxvals = []
    for l in range(halflength,-1,-1):
        newyvals.append((yval[l] + yval[len(yval)-l-1]) /2.0 )
        newxvals.append(xval[l])
    return newxvals,newyvals

#given the camber line and chord of a set of points, finds the camber of the foil
def getCambervalue(camx,camy,chord):
    m = max(camy) / chord
    p = camx[camy.index(m * chord)] / chord
    return m,p

#given the x values for a set of data points (supposing no twist) 
# returns the chord length of the airfoil
def getChord(x): 
    maxi = max(x)
    mini = min(x)
    chord = round(maxi - mini,2)
    return chord

#in order to create a bonefoil, this function creates a list of the values
# which are to be used for the thickness of the bone in the fishBone design
# this amount is calculated to be 3.2% of the chord.
def boneThickness(chord,x):
    thickness = []
    for val in x:
        thickness.append(chord * 0.032)
    return thickness

#make a function that enthickens an airfoil to a desired thickness....
#def enThicken(a,v):
#    
#
#
#

#this function creates the airfoil to the size of the chord given
# assuming that the wing currently has no twist
def enlarge(x,y,newChord):
    oldChord = getChord(x)
    ex = []
    ey = []
    for u in range(len(x)):
        ex.append(x[u] * (newChord/oldChord) )
        ey.append(y[u] * (newChord/oldChord) )
    return ex,ey

#this function finds the length of a line, given its x and y values
def lineLength(x,y):
    length = 0
    for g in range(len(x)-1):
        dx = x[g+1] - x[g]
        dy = y[g+1] - y[g]
        length += (((dx**2.0) + (dy**2.0))**(0.5))
    return length


#this function takes an airfoil list of coordinates (imported from text file)
# and determines the fishbone shape of the airfoil, given x,y
def txtFishboneMaker(x,y):
    chord = getChord(x)
    fx,fy,why,do = FishboneMaker(x,y,chord,None,None)
    camberx,cambery = makeCamberLine(x,y)
    m,p = getCambervalue(camberx,cambery,chord)
    dydx = getdydx(m,p,camberx)
    t = boneThickness(chord,dydx)
    
    newfx,newfy = getfoilShape(t,dydx,camberx,cambery)
    return fx,fy,camberx,cambery,newfx,newfy

#this function determines the outline of an airfoil shape (even fishbone) given 
# its data points and kerf value. shifts the foil shape to include this kerf value 
# for the whole geometry
def kerfMaker(x,y,kerf):  
    
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
    
    return kerfedx,kerfedy



def OsciSpine(x,y,chord,thickness,startPerc = .0,endPerc = .88,period = 0.09,spineWidth = 0.0,firstUp = True,isTlocal = False): # period is a % of the chordlength 
    from scipy.interpolate import interp1d as interp
    if(startPerc == 0):
        startPerc = x[y.index(max(y))] / chord
    if(spineWidth == 0):
        spineWidth = 0.2
    elif(spineWidth == 0 and isTlocal == True):
        spineWidth = 1
    sgn = -1
    if(firstUp == True):
        sgn = 1
    
    from numpy import cos, pi
    x_to_interpolate = x[:len(x)//2]
    y_to_interpolate = y[:len(x)//2]
    topInterp = interp(x_to_interpolate[::-1],y_to_interpolate[::-1])
    ox = []
    oy = []
    
    
    
    
    n = 150
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
            oy.append(sgn*((( interpolated_point/2 )*cos(frequency_over_chord*(ox[k] - starting_point) + pi)) + interpolated_point/2 - spine_thickness))
        elif((ox[k]/chord - startPerc)/period > int((endPerc - startPerc)/period)):
            if((ox[k]/chord - startPerc)/period < int((endPerc - startPerc)/period)+ 0.5):
                oy.append(sgn*(interpolated_point - spine_thickness)*-cos(frequency_over_chord*(ox[k] - starting_point)))
            else:
                oy.append(sgn*((( interpolated_point/2 )*cos(frequency_over_chord*(ox[k] - starting_point) + pi)) + interpolated_point/2 - spine_thickness))
        else:
            oy.append(sgn*(interpolated_point - spine_thickness)*-cos(frequency_over_chord*(ox[k] - starting_point))) # create OsciSpine in such a way that the ending opposite
        k += 1
    
        #################################################################TURN THIS TEXT BLOCK INTO A SEPARATE FUNCTION WHICH ###################################
        ################################################ GIVEN A SHIFT VALUE MAKES TWO SETS OF DATA THAT ARE SHIFTED TO COVER THE LINE GIVEN####################
    topy = []
    topx = []
    botx = []
    boty = []
    for j in range(len(oy)-1):
        slope = (oy[j+1] - oy[j]) / (ox[j+1] - ox[j])
        inslope = 1 / (-slope)
        xmid = ((ox[j+1] - ox[j]) / 2 ) + ox[j]
        if(isTlocal):
            spine_thickness = spineWidth*topInterp(xmid)
        ymid = ((oy[j+1] - oy[j]) / 2 ) + oy[j]
        xchange = ( spine_thickness**2 / (1 + (inslope**2)) )**(1/2)
        ychange = ( spine_thickness**2 / (1 + (inslope**-2)) )**(1/2)
        xtnew = xmid - (xchange if inslope < 0 else -xchange)
        ytnew = ymid + ychange
        topy.append(ytnew)
        topx.append(xtnew)
        xbnew = xmid + (xchange if inslope < 0 else -xchange)
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
                    osx.append(botx[len(botx)-1])
                    osy.append(bottomInterp(botx[len(botx)-1]))
                k += 1
        i += 1

    return osx,osy



def infoCompile(foil,enlargeFactor,isfish,fishbone_variables,oscispine_variables):
    chord = 1
    camber,camber_location,thickness_to_chord_ratio = getNACAfour(foil)
    cosine_clustering_x_values = getCosineClustering(100)
    camberliney = getMeanCamberLine(camber,camber_location,cosine_clustering_x_values)
    thickness = getYht(thickness_to_chord_ratio *chord * 2.0,cosine_clustering_x_values)
    derivatives = getdydx(camber,camber_location,cosine_clustering_x_values)
    xval,yval = getfoilShape(thickness,derivatives,cosine_clustering_x_values,camberliney)
    
    xval,yval = enlarge(xval,yval,enlargeFactor)
    from geometry import getThickness as getTh
    if(isfish == True):
        if(fishbone_variables[0] == True):
            fx,fy = FishBonePoints(xval,yval,getTh(xval,yval),enlargeFactor)
        else:
            start = fishbone_variables[1]
            end = fishbone_variables[2]
            bone = fishbone_variables[3]
            gap= fishbone_variables[4]
            spine = fishbone_variables[5]
            islocal = fishbone_variables[6]
            fx,fy = FishBonePoints(xval,yval,getTh(xval,yval),enlargeFactor,start,end,bone,gap,spine,islocal)
        ox = 0; oy = 0;
    
    elif(oscispine_variables[0] == True):
        if(oscispine_variables[1] == True):
            ox,oy = OsciSpine(xval,yval,enlargeFactor,getTh(xval,yval))
        else:
            start = oscispine_variables[2]
            end = oscispine_variables[3]
            period = oscispine_variables[4]
            spineWidth = oscispine_variables[5]
            firstUp = oscispine_variables[6]
            isTlocalOsc = oscispine_variables[7]
            ox,oy = OsciSpine(xval,yval,enlargeFactor,getTh(xval,yval),start,end,period,spineWidth,firstUp,isTlocalOsc)
        fx = 0; fy = 0;
    else:
        fx = 0; fy = 0;
        ox = 0; oy = 0;
    
    return xval,yval,fx,fy,ox,oy


#infoCompile('NACA 0018',5)

#def OsciSpine(x,y,startPerc,endPerc,frequency,spineWidth,firstUp = True):
#    from scipy.interpolate import interp1d as interp
#    from numpy import cos,exp
#    topInterp = interp(x[:len(x)//2],y[:len(x)//2])
#    bottomInterp = interp(x[len(x)//2:],y[len(x)//2:])
#    
#    
#    
#    return ox,oy

