#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 20 10:27:56 2019

@author: benjamin
"""
import point as pt
from math import cos,sin,radians,acos,degrees

# given two points, find the distance between them
def get_distance(point_a,point_b):
    
    # find vertical separation
    rise = point_a.y - point_b.y
    
    # find horizontal separation
    run = point_a.x - point_b.x
    
    # find distance using pythagorean theorum
    distance = ( rise**2 + run**2 )**0.5
    
    # return value
    return distance

# given a point (x,y) and slope m, where is the y intercept?
def get_intercept(pt,m):
    
    k = pt.y - (m * pt.x)
    
    return k

# given two points, what is the slope between them?
def get_slope(point_a,point_b):
    # find rise
    rise = point_a.y - point_b.y
    
    # find run
    run = point_a.x - point_b.x
    
    # find slope
    if(not run == 0):
        m = rise / run
    else:
        m = 0
    
    # return slope
    return m
    

# given the slopes and intercepts (a and m   and   b and k respectively)
#  determine the point at which they intersect
def find_intersection(a,b,m,k):
    
    # if a == m, invalid operation
    if(a == m):
        if(b == k):
            print('these are the same line')
            x = 1
            y = m*x + k
            return x,y
        else:
            print('these lines are parallel')
            return 0,0
    
    # else, continue, lines will intersect
    
    # given y = m*x + k and y = a*x + b
    #  means m*x + k = a*x + b
    # thus (m - a)*x = b - k
    #      x = ( b - k ) / ( m - a )
    
    x = ( b - k ) / ( m - a )
    
    # this x value can be used with either equation to determine the y value
    
    y = m * x + k
    
    # return the point
    
    point = pt.point(x,y)
    
    return point



# function which, given either the slope and intercept (a,b) of a linear line, or
#  two points (a,b) and (u,v) which define a line, determine the point on the line
#  which lies closest to (where normal begins) the point (x,y)
def closest_point(off_point,point_1,point_2 = pt.point(),a = 0, b = 0):
    
    # if two points given
    if( a == 0 and b == 0):
        
        # find slope
        a = get_slope(point_1,point_2)
        
        # find intercept
        b = get_intercept(point_2,a)
    
    # else if not, continue with a as slope and b as intercept
        
    # find perpendicular slope
    _a = - 1 / a
    
    # find perpendicular y intercept
    
    _b = get_intercept(off_point,_a)
    
    # find intersection
    
    closest_point = find_intersection(a,b,_a,_b)
    
    # return point
    
    return closest_point


# determine the two quadratic roots of a quadratic equation
def quadratic_root(a,b,c):
    
    plus_x = ( - b + ( b**2 - 4 * a * c )**0.5 ) / ( 2 * a )
    
    minus_x = ( - b - ( b**2 - 4 * a * c )**0.5 ) / ( 2 * a )
    
    plus_pt = pt.point(plus_x,0,0)
    
    minus_pt = pt.point(minus_x,0,0)
    
    return plus_pt,minus_pt


# find the pt on a line y = m*x + b which is a set distance r from the point (x,y)
def pt_set_distance_on_line(mid_pt,m,b,r,right = True):
    
    # find intersection between line y = m*x + b
    #  and circle (x - xc)**2 + (y - yc)**2 = r**2
    # place y equation in circle equations,
    # yields:
    # (1 + m**2) * x**2 + (-2*xc + 2*m*b - 2*yc*m) * x + (xc**2 + b**2 -2*yc*b + yc**2 - r**2) = 0
    x_sqrd_terms = 1 + m**2
    x_terms = - 2 * mid_pt.x + 2 * m * b - 2 * mid_pt.y * m
    const_terms = mid_pt.x**2 + b**2 - 2 * mid_pt.y * b + mid_pt.y**2 - r**2
    
    # determine positive and negative x roots ( + and - due to how equations is solved)
    
    plus_pt,minus_pt = quadratic_root(x_sqrd_terms,x_terms,const_terms)
    
    # determine right and left x roots
    r_root = pt.point(); l_root = pt.point()
    
    if(plus_pt.x > mid_pt.x):
        r_root.set_point(plus_pt)
        l_root.set_point(minus_pt)
    else: # if(minus_x > x)
        r_root.set_point(minus_pt)
        l_root.set_point(plus_pt)
    
    # initialize solved points
    point_set_d = pt.point()
    
    # if right point desired
    if(right):
        point_set_d.set_point(r_root)
    else: # if left point desired
        point_set_d.set_point(l_root)
    
    point_set_d.y = m * point_set_d.x + b
    
    # return the points
    
    return point_set_d


#def pt_set_dist_on_line(x,y,m,b,r,right = True):
#    
#    # given that point is distance r along line
#    # if right
#    # then r**2 = (ptx - x)**2 + (pty - y)**2
#    # AND pty = m * ptx + b
#    # thus r**2 = ptx**2 - 2*x*ptx + x**2 + pty**2 - 2*y*pty + y**2
#    # inst r**2 = ptx**2 - 2*x*ptx + x**2 + ( m**2 * ptx**2 + 2*b*m*ptx + b**2 ) - 2*y*m*ptx - 2*y*b + y**2
#    # now 0 = ptx**2 - 2*x*ptx + x**2 + ( m**2 * ptx**2 + 2*b*m*ptx + b**2 ) - 2*y*m*ptx - 2*y*b + y**2 - r**2
#    # x_sqrd_t = 1 + m**2
#    # x_t = -2*x +2*b*m - 2*y*m
#    # c_t = x**2 + b**2 - 2*y*b + y**2 - r**2
#    
#    
#    
#    return 0,0

    
# find the angle from a point(x,y) to another point(u,v) from the +x relative axis
def get_theta(point_1,point_2,radius,right = True):
    
    m = get_slope(point_1,point_2)
    
    b = get_intercept(point_1,m)
    
    set_point = pt_set_distance_on_line(point_1,m,b,radius,right)
    
    run = 0
    
    if(right):
        run = set_point.x - point_1.x
    else:
        run = point_1.x - set_point.x
#    print(run / radius)
    theta_rad = acos( run / radius )
    
    theta = degrees(theta_rad)
    
    return theta

#this function creates an xy circle
def circle (center_point,diameter,percent_of_circle = 1.0, Theta = 0, n = 20, isccw = True): # defines circle who starts at rightmost of circle, moves ccw
    circle = pt.points()
    
    r = diameter / 2
    
    
    # a value is determined for the given set of points to be how far apart the points are interspersed, in degrees
    degStep = percent_of_circle * 360 / n
    
    # max cw if not ccw
    if(not isccw):
        degStep *= -1
    
#    # to create a closed hole
#    if(percent_of_circle == 1.0):
    n += 1
    
#    # to ensure that the control horn semicircle does not cut into TE of wing
#    if(percent_of_circle == 0.5):
#        n -= 0
    
    # an iterative process to find the x,y coordinates for each point along the circle
    for i in range(n):
        
        new_pt = pt.point()# center_point.x,center_point.y)
        new_pt.set_point(center_point)
#        center_point.coord()
        new_pt.move( r*cos(radians(Theta)), r*sin(radians(Theta)))
        
        Theta += degStep
        
        circle.add(new_pt)
#    circle.makeWhole()
    # the x,y points are returned to the user
    return circle



def control_horn_semicircle(airfoil,depth,height,diameter,add_gap):
    
    # set control horn values
    
    controlhorn = pt.points()
    
    # start at trailing edge of wing
    
    controlhorn.add(airfoil.points[0])
    
    
    # define hole
    
    center = pt.point()
    center.x = controlhorn.points[0].x - depth + diameter
    center.y = controlhorn.points[0].y + height;
    
    hole = circle(center,diameter)
    
    chord = airfoil.chord()
    
    # find two points that are near to the hole center x value (around 1/100 c away)
    
    pt_a = pt.point(); pt_b = pt.point();
    pt_a.x = center.x + ( chord / 100); pt_b.x = center.x - (chord / 100)
    pt_a.y = airfoil.interp(pt_a.x) ; pt_b.y = airfoil.interp(pt_b.x)
    
    # determine closest point
    
    closest_pt = closest_point(center, pt_a, pt_b )
    
    # determine distance to hole from closest point
    
    d = get_distance(closest_pt,center)
    
    circle_diameter = 2 * ( d + 1.5 * diameter )
    
    # determine theta value
    
    theta = get_theta(closest_pt,pt_a,circle_diameter / 2)
    
    half_circle = circle(closest_pt,circle_diameter,0.5,theta)
    
    # add circle points to chx and chy lists
    
    controlhorn.add(half_circle)
    
    
    
    # remaining airfoil points are added
    
    add_airfoil_points = True
    
    counter = 0
    
    while(add_airfoil_points):
        
        # if the control horn needs to be added, skip first few points
        
        if(airfoil.points[counter].x >= controlhorn.points[len(controlhorn.points)-1].x):
            
            counter += 1
        
        # else, add the airfoil points, turn off outer while loop
        
        else:
            controlhorn.add(airfoil.split(counter))
            #or   #   #  above is faster
#            while(counter < len(airfoil.x())):
#                controlhorn.add(airfoil.points[counter])
#                counter += 1
            
            add_airfoil_points = False
    
    
    
    controlhorn.add(hole)
    
    return controlhorn,hole


 ####################################################### stopped changing pts
def control_horn_sharp(airfoil,depth,height,diameter,add_gap):
    
    # set control horn vlaues
    
    control_horn = pt.points()
    
    # start at trailing edge of wing
    
    control_horn.add(airfoil.points[0])
    
    # move back depth and up height plus diam
    
    nextx = control_horn.points[0].x - depth + 2 * diameter
    
    nexty = control_horn.points[0].y + height + diameter
    
    control_horn.add(pt.point(nextx,nexty))
    
    # move to before hole
    
    nextx = nextx - 2 * diameter
    
    control_horn.add(pt.point(nextx,nexty))
    
    # move to below and before hole
    
    nexty = nexty - 2 * diameter
    
    control_horn.add(pt.point(nextx,nexty))
    
    # if gap is desired, add gap
    if(add_gap):
        
        # move to in front of and below hole
        
        nextx = nextx + 2 * diameter
        
        control_horn.add(pt.point(nextx,nexty))
    
    
    # move to on airfoil
    
    nexty = airfoil.interp(nextx)
    
    control_horn.add(pt.point(nextx,nexty))
    
    
    
    # remaining airfoil points are added
    
    add_airfoil_points = True
    
    counter = 0
    
    while(add_airfoil_points):
        
        # if the control horn needs to be added, skip first few points
        
        if(airfoil.points[counter].x >= control_horn.points[len(control_horn.points)-1].x):
            
            counter += 1
        
        # else, add the airfoil points, turn off outer while loop
        
        else:
            control_horn.add(airfoil.split(counter))
#            while(counter < len(x)):
#                chx.append(x[counter])
#                chy.append(y[counter])
#                counter += 1
            
            add_airfoil_points = False
    
    # define hole
    
    center = pt.point()
    center.x = control_horn.points[0].x - depth + diameter
    center.y = control_horn.points[0].y + height;
    
    hole = circle(center,diameter)
    
    control_horn.add(hole)
    
    return control_horn,hole
    
    


def control_horn(airfoil,depth,height,diameter,horn_type,add_gap):
    
    # if sharp horn picked
    if(horn_type == 'sharp'):
        control_horn,hole = control_horn_sharp(airfoil,depth,height,diameter,add_gap)
        
    # if half-circle horn picked
    elif(horn_type == 'semicircle'):
        control_horn,hole = control_horn_semicircle(airfoil,depth,height,diameter,add_gap)
    
    
    
    return control_horn,hole


#from matplotlib import pyplot as plt
#
#center = pt.point()
#circle = circle(center, 2)
#plt.axis('equal')
#plt.plot(circle.x(),circle.y())
#plt.show()





#from matplotlib import pyplot as plt
#import File_handler as fh
#import Airfoil_generator as ag
#
#ax,ay = fh.FileReader('/home/benjamin/Desktop/selig bezier NACA 0012.dat') # '/home/benjamin/Desktop/foilshapes/NACA 0015.txt'
#ax,ay = ag.change_chord(ax,ay,12)
#
#airfoil = pt.points(ax,ay)
#
#d = 1.2 / 2.54
#h = 0.5 / 2.54
#di = 0.2 / 2.54
#add_gap = True
#
#
#ch,h = control_horn(airfoil,d,h,di,'semicircle',add_gap)
#
#airfoil.plot([ch,h])
#
#
#axis = [0.9*max(ax),max(ax),-0.1 * max(ax),0.1 * max(ax)]
#
#airfoil.plot([ch,h],axis)






#from matplotlib import pyplot as plt
#import File_handler as fh
#import Airfoil_generator as ag
#
#ax,ay = fh.FileReader('/home/benjamin/Desktop/foilshapes/NACA 0015.txt')
###ax,ay = ag.change_chord(ax,ay,12)
#
#
#d = 1.2 / 2.54
#h = 0.5 / 2.54
#di = 0.2 / 2.54
#add_gap = True
#
#
#chx,chy,hx,hy = control_horn(ax,ay,d,h,di,add_gap)
#
#plt.axis('equal')
#plt.plot(ax,ay)
#plt.plot(chx,chy)
#plt.plot(hx,hy)
#plt.show()
#
#plt.axis([0.9*max(ax),max(ax),-0.1 * max(ax),0.1 * max(ax)])
#plt.plot(ax,ay)
#plt.plot(chx,chy)
#plt.plot(hx,hy)
#plt.show()
    

