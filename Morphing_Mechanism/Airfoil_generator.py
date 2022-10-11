#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 11 15:16:31 2018

@author: benjamin
"""
import point as pt

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
    
    camber_line = pt.points()
    
    dtheta = 2.0* pi / (num - 1.0)
    
    for i in range(1,int(num/2)+1):
        
        x_value = 0.5 * (1. - cos((i - 0.5) * dtheta))
        
        new_point = pt.point(x_value)
        
        camber_line.add(new_point)
#    camber_line.plot(display = ['bo'])
    return camber_line



#given the percent along the chord, the camber and the camber Location,
# this function finds the y values for the mean camber line and returns
# them to the user
def getNACAMeanCamberLine(m,p,camber_line):
    
    for val in camber_line.points:
        
        yc = 0.0
        
        if val.x >= 0.0 and val.x <= p:
            
            yc = (m/(p**2.0)) * (2.0*p*val.x - (val.x**2.0))
            
        elif val.x > p and val.x <= 1:
            
            yc = (m/((1.0-p)**2.0)) * ((1.0 - 2.0*p) + 2.0*p*val.x - (val.x**2.0))
        
        elif(m == 0 and p == 0):
            
            yc = 0.0
            
        val.y = yc
        
    return camber_line



##given the thickness to chord ratio, and the percent of chord
## this function returns a list of the y values (thickness function of
## the airfoil)
#def getYht(t,x): #maxThickness, percentOfChord
#    
#    # a values are set
#    
#    a0 = 0.2969
#    
#    a1 = -0.126
#    
#    a2 = -0.3516
#    
#    a3 = 0.2843
#    
#    a4 = -0.1036
#    
#    # if x is a list, run through for each point
#    
#    if(type(x) == list):
#        Yht = []
#        
#        for val in x:
#            
#            y = ( t / 0.2 ) * ( a0*val**0.5 + a1*val + a2*val**2 + a3*val**3 + a4*val**4 )
#            
#            Yht.append(y)
#    
#    elif(type(x) == float or type(x) == int):
#        
#        y = ( t / 0.2 ) * ( a0*x**0.5 + a1*x + a2*x**2 + a3*x**3 + a4*x**4 )
#        
#        Yht = y
##        Yht.append(0.5*t*( 2.98*(val**0.5) - 1.32*val - 3.286*(val**2) +
##                       2.441*(val**3) - 0.815*(val**4)))
#    return Yht

def getYht(t,x): #maxThickness, percentOfChord
    
    # a values are set
    
    a0 = 0.2969
    
    a1 = -0.126
    
    a2 = -0.3516
    
    a3 = 0.2843
    
    a4 = -0.1036
    
    
    Yht = 0
    
    # if x is a list, run through for each point
    
    if(type(x) == list):
        Yht = []
        
        for val in x:
            
            y = ( t / 0.2 ) * ( a0*val**0.5 + a1*val + a2*val**2 + a3*val**3 + a4*val**4 )
            
            Yht.append(y)
    
    elif(type(x) == float or type(x) == int):
        
        y = ( t / 0.2 ) * ( a0*x**0.5 + a1*x + a2*x**2 + a3*x**3 + a4*x**4 )
        
        Yht = y
    
    
    
    elif(type(x) == pt.point):
        
        xval = x.x
        
        x.y = ( t / 0.2 ) * ( a0*xval**0.5 + a1*xval + a2*xval**2 + a3*xval**3 + a4*xval**4 )
        
        Yht = x
    
    elif(type(x) == pt.points):
        
        Yht = pt.points()
        
        for val in x.points:
            
            new_point = pt.point()
            
            new_point.x = val.x
            
            new_point.y = ( t / 0.2 ) * ( a0*val.x**0.5 + a1*val.x + a2*val.x**2 + a3*val.x**3 + a4*val.x**4 )
            
            Yht.add(new_point)
    
    return Yht

#given the camber, camber location, and a list of the percent of Chord
# values, this function returns a list of the dy/dx values for each value given
def getdydx(m,p,camber_line): #camber, camber location, percentofChord in camber line
    
    dydx = pt.points()
    
    for val in camber_line.points:
        
        dy = 0.0
        
        if(val.x >= 0.0 and val.x <= p):
            
            dy = ( 2.0 * m / ( p**2.0 ) ) * ( p - val.x )
            
        elif(val.x > p and val.x <= 1):
            
            dy = ( 2.0 * m / ( ( 1.0 - p )**2.0 ) ) * ( p - val.x )
            
        dydx.add(pt.point(y = dy))
        
    return dydx

"""
the function directly below is not transferred over to pts object use
"""
###this function takes the chord and the cosine clustering percentOfChord
### values and returns the xpoints
##def cosineToChord(x,chord):
##    xwithChord = []
##    for i in range(len(x)):                                        # not needed
##        xwithChord.append(x[i] * chord)
##    return xwithChord



#the function below takes all the information gathered thus far,
# and finds the x and y points of the foilshape. It begins at the
# trailing edge and follows the airfoilshape counter clockwise.
def getfoilShape(t,dydx,camber_line):
    from math import sqrt
    airfoil = pt.points()
    
    for j in range(len(camber_line.points)-1,-1,-1):
        
        new_point = pt.point()
        
        new_point.x = camber_line.points[j].x - (t.points[j].y * dydx.points[j].y) / (2.0*sqrt(1.0 + (dydx.points[j].y**2.0)))
        
        new_point.y = camber_line.points[j].y + t.points[j].y / (2.0*sqrt(1.0 + (dydx.points[j].y**2.0)))
        
        airfoil.add(new_point)
        
    for k in range(0,len(camber_line.points)):
        
        new_point = pt.point()
        
        new_point.x = camber_line.points[k].x + (t.points[k].y * dydx.points[k].y) / (2.0*sqrt(1.0 + (dydx.points[k].y**2.0)))
        
        new_point.y = camber_line.points[k].y - t.points[k].y / (2.0*sqrt(1.0 + (dydx.points[k].y**2.0)))
        
        airfoil.add(new_point)
        
    return airfoil



def BoneMaker( airfoil, startperc = 0, endperc = 0.85, 
              boneWidth = 0.0225, gapWidth = 0.0225, spineWidth = 0.25,
              tribase = 0, isTlocal = False, ledgeWidth = 0.0, 
              ledgeDepth = 0.0, ledgeIsFlat = False,skipPoints = False,perc_r = 1.):
    import Fishbone as fb
    bonefoil = fb.BoneMaker(airfoil,startperc,endperc,boneWidth,
                         gapWidth,spineWidth,tribase,isTlocal,ledgeWidth,
                         ledgeDepth,ledgeIsFlat,skipPoints,perc_r)
    
    return bonefoil






#the function below takes a set of x and y values for an airfoil and interpolates
# the mean line of camber based on the location of the y values assumed to associate
# with each x position. The function returns a new set of x values and yvalues to 
# the user that follow the line of mean camber for the airfoil.
def makeCamberLine(airfoil):
    
    camber_line = pt.points()
    
    for i in range(airfoil.x_min_index):
        
        new_point = pt.point()
        
        top = airfoil[i].y
        
        bottom = airfoil.interp(airfoil[i].x,isTop = False)
        
        new_point.x = airfoil[i].x
        
        new_point.y = ( top + bottom ) / 2
        
        camber_line.add(new_point)
    
    airfoil.camber_line = camber_line
    
    return airfoil


"""
getCambervalue needs to be modified as it can not determine the exact location of max camber
in it's current setup
"""
#given the camber line and chord of a set of points, finds the camber of the foil
def getCambervalue(camber_line,chord):
    
    mx = max(camber_line.y()) / chord
    
    p = camber_line.x()[camber_line.y().index(mx * chord)] / chord
    return mx,p



#in order to create a bonefoil, this function creates a list of the values
# which are to be used for the thickness of the bone in the fishBone design
# this amount is calculated to be 3.2% of the chord.
def boneThickness(chord,airfoil):
    thickness = []
    for val in airfoil.x():
        thickness.append(chord * 0.032)
    return thickness


#this function creates the airfoil to the size of the chord given
# assuming that the wing currently has no twist
def change_chord(airfoil,newChord,oldChord = 0):
    if(oldChord == 0):
        oldChord = airfoil.chord()
    
    for u in range(len(airfoil.points)):
        
        airfoil.points[u].x = airfoil.points[u].x * ( newChord / oldChord )
        
        airfoil.points[u].y = airfoil.points[u].y * ( newChord / oldChord )
    
    return airfoil



#this function finds the length of a line, given its x and y values
def lineLength(line):
    
    length = 0
    
    for g in range(len(line.points)-1):
        
        dx = line.points[g+1].x - line.points[g].x
        
        dy = line.points[g+1].y - line.points[g].y
        
        length += (((dx**2.0) + (dy**2.0))**(0.5))
    
    return length


"""
the function directly below is not transferred over to pts object use
"""
##this function takes an airfoil list of coordinates (imported from text file)
## and determines the fishbone shape of the airfoil, given x,y
#def txtFishboneMaker(x,y):
#    chord = getChord(x)
#    fx,fy,why,do = FishboneMaker(x,y,chord,None,None)
#    camberx,cambery = makeCamberLine(x,y)
#    m,p = getCambervalue(camberx,cambery,chord)
#    dydx = getdydx(m,p,camberx)
#    t = boneThickness(chord,dydx)
#    
#    newfx,newfy = getfoilShape(t,dydx,camberx,cambery)
#    return fx,fy,camberx,cambery,newfx,newfy



#this function determines the outline of an airfoil shape (even fishbone) given 
# its data points and kerf value. shifts the foil shape to include this kerf value 
# for the whole geometry
def kerfMaker(airfoil,kerf):  
    
    import kerf as k
    
    kerfed_airfoil = k.Kerf_pt(airfoil,kerf)
    
    return kerfed_airfoil



"""
the function directly below is not transferred over to pts object use
"""
#def OsciSpine(x,y,chord,thickness,startPerc = .0,endPerc = .88,period = 0.09,
#              spineWidth = 0.0,fillet_rad = 0.5,firstUp = True,
#              isTlocal = False): # period is a % of the chordlength 
#    
#    import Oscispine as os
#    
#    osx, osy = os.OsciSpine(x,y,chord,thickness,startPerc,endPerc,period,
#                            spineWidth,fillet_rad,firstUp,isTlocal)
#    
#    return osx,osy






def Bezier_Spine(airfoil, start_percentage = 0.3, end_percentage = 0.9, 
                 P1_percentage = 0.75, P2_percentage = 0.25, T = 0.09, 
                 spine_width_LE = 0.12, spine_width_TE = 0.1, first_up = True, 
                 is_Tlocal = False):
    
    import Bezierspine as bz
    
    bezier_airfoil = bz.Bezier_Spine(airfoil, start_percentage, end_percentage, 
                            P1_percentage, P2_percentage, T, spine_width_LE, 
                            spine_width_TE, first_up, is_Tlocal)
    
    return bezier_airfoil


 # end need to change

def Airfoil_Creator(foil,enlargeFactor,p, cosine_clustering_value = 100): # change
    
    chord = 1
    
    camber,camber_location,thickness_to_chord_ratio = getNACAfour(foil)
    
    camber_line = getCosineClustering(cosine_clustering_value)
    
    camber_line = getNACAMeanCamberLine(camber,camber_location,camber_line)
    
    thickness = getYht(thickness_to_chord_ratio *chord * 2.0,camber_line)
    
    derivatives = getdydx(camber,camber_location,camber_line)
    
    airfoil = getfoilShape(thickness,derivatives,camber_line)
    
    
    airfoil = change_chord(airfoil,enlargeFactor)
    
    if(p['is fishbone'] and not p['is bezier']):
        
        if(p['default fishbone']):
            
            bonefoil = BoneMaker(airfoil)
            
        else:
            
            bonefoil = BoneMaker(airfoil,
                              p['fishbone start'],p['fishbone end'],
                              p['fishbone bone width'],p['fishbone gap width'],
                              p['fishbone spine width'],p['fishbone base'],
                              p['fishbone is t local'],p['fishbone ledge width'],
                              p['fishbone ledge depth'],p['fishbone ledge is flat'],
                              p['fishbone skip points'],p['fishbone percent ledge right'])
        
        bezierfoil = pt.points()
    
    elif(not p['is fishbone'] and p['is bezier']):
        
        if(p['default bezier']):
            
            bezierfoil = Bezier_Spine(airfoil) # change
        else:
            
            bezierfoil = Bezier_Spine(airfoil, 
                              p['bezier start'], p['bezier end'],
                              p['bezier P1 percentage'], p['bezier P2 percentage'],
                              p['bezier period'], p['bezier leading edge spine width'], 
                              p['bezier trailing edge spine width'], 
                              p['bezier first up'], p['bezier is t local'])
            
        bonefoil = pt.points()
        
    else:
        
        bonefoil = pt.points()
        bezierfoil = pt.points()
    
    return airfoil,bonefoil,bezierfoil

# g = {}

# g['is fishbone'] = False
# g['is bezier'] = True
# g["default bezier"] = True


# airfoil,bone,bez = Airfoil_Creator('NACA 0016',1.0,g)


# airfoil.plot( ) # display = 'bo')
