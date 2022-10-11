#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 10 16:24:26 2018

@author: benjamin
"""
from math import radians, cos, sin,tan
import point as pt
import numpy as n

""" function below not changed over to points """
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


# if a list is not whole, ( a closed loop ) make it so
def makeWhole(data, start_at_end = True): # verses ending at the start ie sectionfoil, start on bottom or top
    
    if (round(data.points[len(data.points)-1].x,5) != round(data.points[0].x,5) or 
                 round(data.points[len(data.points)-1].y,5) != round(data.points[0].y,5) or
                 round(data.points[len(data.points)-1].z,5) != round(data.points[0].z,5)):
        
        # copy the lists
        new_data = data.split()
        point = pt.point()
        
#        last_point.fix(data.points[len(data.points)-1].x,data.points[len(data.points)-1].y,data.points[len(data.points)-1].z)
        if(start_at_end):
            point.set_point(data.points[len(data.points)-1])
            
            new_data.add(point,True)
        else:# end at start
            point.set_point(data.points[0])
            new_data.add(point)
        
        return new_data
    
    else:
        return data



# given a list of points, find the area of the outlined shape
#  assume data is planar, ie. z values are all the same ( flat surface )
def Area(data):
    
    A = 0
    
    data.makeWhole()
    
    for i in range(len(data.points)-1):
        
        A += data.points[i].x * data.points[i+1].y
        
        A -= data.points[i+1].x * data.points[i].y
    
    A = 0.5 * abs(A)
    
    return A

# determine the centroid of a closed shape # assuming all at same z
def Centroid(data):
    
    #initialize centroid point
    C = pt.point()
    
    # make the shape whole
    whole_data = makeWhole(data)
    
    # find the area of the shape
    A = Area(whole_data)
    
    # iteratively find the centroid summations
    for i in range(len(whole_data.points)-1):
        
        C.x += ( whole_data.points[i].x + whole_data.points[i+1].x 
                ) * ( whole_data.points[i].x * whole_data.points[i+1].y - 
                whole_data.points[i+1].x * whole_data.points[i].y )
        
        C.y += ( whole_data.points[i].y + whole_data.points[i+1].y 
                ) * ( whole_data.points[i].x * whole_data.points[i+1].y - 
                whole_data.points[i+1].x * whole_data.points[i].y )
    
    # multiply each by 1/(6*A) to get the true centroid
    C.x /= 6 * A
    C.y /= 6 * A
    
    C.z = whole_data.points[0].z
    
    # return centroid
    return C

# given 3 points in cartesian coordinates, determine the normal 
#  to the plane defined by the 3 points
def Normal(point_zero,point_one,point_two):
    
    # set up vectors
    V1 = pt.point( point_one.x - point_zero.x, point_one.y - point_zero.y, point_one.z - point_zero.z )
    V2 = pt.point( point_two.x - point_zero.x, point_two.y - point_zero.y, point_two.z - point_zero.z )
    
    # initialize normal vector
    N = pt.point()
    
    # calculate V1 x V2
    N.x = V1.y * V2.z - V1.z * V2.y
    N.y = V1.z * V2.x - V1.x * V2.z
    N.z = V1.x * V2.y - V1.y * V2.x
    
    # return vector
    return N


#a function that has not yet been developed
def CenterOfMass():
    return pt.point()



#a function that is still to be developed, but is in use.
def quarter_chord(airfoil): # finds the cg of a given NACA 4 digit airfoil
    qc = pt.point()
    chord = airfoil.chord()
    qc.x = chord / 4
    qc.x += airfoil.points[airfoil.x_min_index].x
    #some code that finds the exact position of y of cg on wing using NACAfour number
    avg_y = 0;
    
    for i in range(len(airfoil.points)):
        avg_y += airfoil.points[i].y
    
    qc.y = avg_y / len(airfoil.points)
    return qc

""" function below not changed over to points """
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
def giveWashout(airfoil,twistangle = 10): 
    qc = quarter_chord(airfoil)
    
    # the angle is made negative and changed to radians for mathematical purposes
    twist = radians(-twistangle)
    
    washed_out = pt.points()
    
    # an iterative process to take each point in x,y and rotate it about the center of chord
    for i in range(len(airfoil.points)):
        
        wpt = pt.point()
        
        wpt.x = airfoil.points[i].x - qc.x # set about origin
        wpt.y = airfoil.points[i].y - qc.y
        
        wpt_new = pt.point()
        wpt_new.x = wpt.x*cos(twist) + wpt.y*sin(twist) # rotate by given angle
        wpt_new.y = wpt.y*cos(twist) - wpt.x*sin(twist)
        
        wpt.x = qc.x + wpt_new.x # return to position
        wpt.y = qc.y + wpt_new.y
        
        # append to list
        washed_out.add(wpt)
    
    # return airfoil
    return washed_out
 
    


#this function takes in a points object, and rounds each value in the list for Gcode manipulation
def roundValuesG(airfoil):
    
    rounded = pt.points()
    for point in airfoil.points:
        rnd_pt = pt.point()
        
        rnd_pt.x = n.round(point.x,3)
        rnd_pt.y = n.round(point.y,3)
        rnd_pt.z = n.round(point.z,3)
        
        rounded.add(rnd_pt)
    return rounded
    



def giveSweep(airfoil,sweepAngle=0,span=1): # if there is sweep, this function returns the foilshape of the UV side shifted for the sweep
    #assuming degree is from front of wing to parallel from old wing root ##### we dont need y vals, right?
   
    swept = pt.points()
    shiftValue = tan(radians(sweepAngle)) * span # the shift is calculated
    for i in range(len(airfoil.points)):
        swept_pt = pt.point()
        swept_pt.x = shiftValue + airfoil.points[i].x
        swept_pt.y = airfoil.points[i].y
        swept_pt.z = airfoil.points[i].z
        
        swept.add(swept_pt)
    
    # return
    return swept  # returns the airfoil that was swept



def giveDihedral(airfoil,dihedralAngle=0,span=1): # if there is dihedral, this function returns the foilshape coords of the uv side shifted
    
    dihedraled = pt.points()
    
    shiftValue = cos(radians(dihedralAngle)) * span
    for i in range(len(airfoil.points)):
        dih_pt = pt.point()
        dih_pt.x = airfoil.points[i].x
        dih_pt.y = shiftValue + airfoil.points[i].y
        dih_pt.z = airfoil.points[i].z
        
        dihedraled.add(dih_pt)
    return dihedraled



#aligns tip quarter chord with root quarter chord
def alignQuarterChord(xyzfoil,uvwfoil):
    rchord = xyzfoil.chord()
    tchord = uvwfoil.chord()
    shiftval = (rchord/4) - (tchord/4)
    for e in range(len(uvwfoil.points)):
        uvwfoil.points[e].x += shiftval
    return uvwfoil



"""

fix getThickness situation so that it solves for the thickness

"""


#this function takes a given set of points for the x,y and u,v shapes of
# a foam wing, and determines what points the CNC hotwire machine should follow
# in cutting out the wing given different shapes,offsets, and the wire length
def fixAxisPoints(xyzfoil,uvwfoil,span,offset,wireLength):
    #offset on the left side of a right wing
    rightSpace = wireLength - span - offset
    
    Axyzfoil = pt.points();  Auvwfoil = pt.points()
    
    # an iterative procedure is followed to determine where the points should be
    # based on the slope in the x and y directions between the two sets of points
    for r in range(len(xyzfoil.points)):
        slopex = (uvwfoil.points[r].x - xyzfoil.points[r].x) / span
        slopey = (uvwfoil.points[r].y - xyzfoil.points[r].y) / span
        Axyz_pt = pt.point();  Auvw_pt = pt.point()
        Axyz_pt.x = xyzfoil.points[r].x - slopex * offset
        Axyz_pt.y = xyzfoil.points[r].y - slopey * offset
        Auvw_pt.x = uvwfoil.points[r].x + slopex * rightSpace
        Auvw_pt.y = uvwfoil.points[r].y + slopey * rightSpace
        Axyzfoil.add(Axyz_pt); Auvwfoil.add(Auvw_pt);
    
    # return point objects
    return Axyzfoil, Auvwfoil

#this function shifts the y and v data up the axis so it will fit in the machine constraints
def shiftUp(xyzfoil,uvwfoil,height):
    maxheight = max(max(xyzfoil.y()),max(uvwfoil.y()))
    shift = height - maxheight - 0.125
    
    for i in range(len(xyzfoil.points)):
        xyzfoil.points[i].y += shift
        uvwfoil.points[i].y += shift
    return xyzfoil,uvwfoil

#this function shifts the x and u data right on axis so it fits in the machine
def shiftRight(xyzfoil,uvwfoil):
    shift = 1#inch
    
    for j in range(len(xyzfoil.points)):
        xyzfoil.points[j].x += shift
        uvwfoil.points[j].y += shift
    return xyzfoil,uvwfoil




# this function is for partioning. if the user only needs the trailing edge, or the leading edge, 
#   or some other section of the airfoil, this function cuts the points out and returns them to 
#   the user.  
def sectionMaker(airfoil,percentage_start = 9001,percentage_end = 1.,add_fillets = False):
    
    # a shift value is initialized in case the foil is shifted to the right.
    shift_value = airfoil.points[airfoil.x_min_index].x
    
    
    # if no start percentage is given( or 9001 is given) start at the location of max thickness
    if(percentage_start == 9001):
        percentage_start = airfoil.x()[airfoil.y().index(max(airfoil.y()))] / airfoil.chord()
    
    
    # the lists of x y coordinates to be returned are initialized
    section = pt.points()
    
    # a count factor is initialized to hold the place of the computation
    count = 0
    
    chord = airfoil.chord()
    
    # a while loop is run which skips the points at the trailing edge top of the foil
    #   which are not to be included
    while(airfoil.points[count].x > (percentage_end*chord + shift_value)):
        count += 1
    
    
    # the point at percentage_end*chord is interpolated
    new_pt = pt.point()
    if(percentage_end != 1.0):
        new_pt.x = percentage_end*chord + shift_value
        new_pt.y = airfoil.interp(new_pt.x)
        section.add(new_pt)
    
    # if the start perc is not at 0, do the following:
    if(percentage_start != 0):
        
        # a while loop adds all the points up to the percentage_start value
        while(airfoil.points[count].x > (percentage_start*chord + shift_value)):
            section.add(airfoil.points[count])
            
            count += 1
        
        
        # the end point of the top of the airfoil is interpolated
        new_pt = pt.point()
        new_pt.x = percentage_start*chord + shift_value
        new_pt.y = airfoil.interp(percentage_start*chord + shift_value)
        section.add(new_pt)
        
        
        # a while loop skips all the points between the top and bottom of the airfoil
        #    at the start percentage point. count is added to prevent errors
        count += 2
        while(airfoil.points[count].x < ( percentage_start*chord + shift_value )):
            count += 1
        
        
        # the start point of the bottom of the airfoil is interpolated
        new_pt = pt.point()
        new_pt.x = percentage_start*chord + shift_value
        new_pt.y = airfoil.interp(new_pt.x,False)
        section.add(new_pt)

    # else: continue 
    #    the code below adds all the points which should be added to the end of the section        
        
    # a while loop runs through and adds points which are between the start and
    #   end percentages on the bottom of the airfoil
    while(count < len(airfoil.points) and airfoil.points[count].x < (
            percentage_end*chord + shift_value)):
        section.add(airfoil.points[count])
        count += 1
    
    # the final interpolation points are added to the section coordinates
    # the point at percentage_end*chord is interpolated
    new_pt = pt.point()
    if(percentage_end != 1.0):
        new_pt.x = percentage_end*chord + shift_value
        new_pt.y = airfoil.interp(new_pt.x,False)
        section.add(new_pt)
    
#    section = makeWhole(section,start_at_end = False)
    section.makeWhole(start_at_end = False)
    # and the lists of coordinate values are returned  to the user
    return section



#import File_handler as fh
#
#
#ax,ay = fh.FileReader('/home/benjamin/Desktop/foilshapes/E335.txt')
#
#airfoil = pt.points(ax,ay)
#
#airfoil.plot()
#
#airfoil = makeWhole(airfoil)
#
#print(Area(airfoil))
#
#point = Centroid(airfoil)
#point.coord()
#
#Normal(airfoil.points[0],airfoil.points[1],airfoil.points[2]).coord()
#
#print(airfoil.x_max_index)
#airfoil.points[airfoil.x_max_index].coord()
#chord = getChord(airfoil)
#print(chord)
#
#qc = quarter_chord(airfoil,chord)
#qc.coord()
#
#
#wash = giveWashout(airfoil)
#wash.plot()
#
#rounded = roundValuesG(airfoil)
#rounded.plot()
#
#swept = giveSweep(airfoil,10)
#swept.plot()
#
#dihed = giveDihedral(airfoil,10)
#dihed.plot()
#
#swept,dihed = shiftRight(swept,swept)
#swept.plot()
#
#swept = alignQuarterChord(airfoil,swept)
#airfoil.plot([swept])
#
#print(getPercentThickness(airfoil,chord))
#
#section = sectionMaker(airfoil,chord,0.1,0.9)
#
#section.plot()
#
#one,two = fixAxisPoints(airfoil,airfoil,1,10,30)
    


