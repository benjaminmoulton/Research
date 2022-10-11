#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 10 09:21:43 2019

@author: benjamin
"""

import point as pt

import geometry as geo


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

# first two points define a line, second two points define a line
def intersect(line1_pt1,line1_pt2,line2_pt1,line2_pt2):
    
    m0 = get_slope(line1_pt1,line1_pt2)
    m1 = get_slope(line2_pt1,line2_pt2)
    
    k0 = get_intercept(line1_pt1,m0)
    k1 = get_intercept(line2_pt1,m1)
    
    if(m0 == m1):
        return 'non-intersect'
    
    int_pt = find_intersection(m0,k0,m1,k1)
    
    in_x = int_pt.x > min(line1_pt1.x,line1_pt2.x) and int_pt.x < max(line1_pt1.x,
                         line1_pt2.x) and int_pt.x > min(line2_pt1.x,line2_pt2.x
                                   ) and int_pt.x < max(line2_pt1.x,line2_pt2.x)
    in_y = int_pt.y > min(line1_pt1.y,line1_pt2.y) and int_pt.y < max(line1_pt1.y,
                         line1_pt2.y) and int_pt.y > min(line2_pt1.y,line2_pt2.y
                                    ) and int_pt.y < max(line2_pt1.y,line2_pt2.y)
    
    if(in_x and in_y):
        return int_pt
    
    else:
        return 'non-intersect'
    
"""  VVVVVVVVVVVVVVVV This could be made much simpler.. collect all to remove points, 
        list(set()) them, and then remove them all in one fell swoop VVVVVVVVVVVVV"""

# given an xy list(s), remove any self intersections that may occur
def remove_intersections(airfoil, dont_keep_before_first = True):
    
    i_pts = pt.points(); i_vals = pt.points()
    
    duplicate_indices = pt.points()
    
    for i in range(len(airfoil.points)-1):
        for j in range(i+1,len(airfoil.points)-1):
            
            ip = airfoil.points[i]
            jp = airfoil.points[j]
            if(ip.x == jp.x and ip.y == jp.y and ip.z == jp.z):
                
                duplicate_indices.add(pt.point(j))
    
    values = pt.points(list(set(duplicate_indices.x())))
#    values.show()
    count = 0
    for i in range(len(values.points)):
        count += 1
        airfoil.remove(values.points[i].x - count)
    
    
    
            
    for i in range(len(airfoil.points)-1):
        for j in range(i,len(airfoil.points)-1):
            ip = intersect( airfoil.points[i], airfoil.points[i+1], 
                           airfoil.points[j], airfoil.points[j+1] )
            
            if(type(ip) != str):
                i_pts.add(ip); i_vals.add(pt.point(i,j))
#                """ """
#                ip.coord()
#                airfoil.plot([ip], axis = [ip.x -0.01,ip.x + 0.01, ip.y-0.01,ip.y+0.01] ,display = ['','ro'])
#                
#                
#                """ """
#    airfoil.show()
    
    pre_removal = pt.points()
    # check list to see if any double items will be removed
    for i in range(len(i_vals.points)):
        for j in range(len(i_vals.points)):
            if(not dont_keep_before_first or (not j == 0 and dont_keep_before_first)):
                if(i_vals.points[i].x > i_vals.points[j].x and
                   i_vals.points[i].y < i_vals.points[j].y):
                    pre_removal.add(pt.point(i))
    
    pre_removal = pt.points(list(set(pre_removal.x())))
    for k in range(len(pre_removal.points)):
#        print(len(i_vals.points),pre_removal.points[k].x)
        i_vals.remove(pre_removal.points[k].x)
    
    removed = 0
    
    for i in range(len(i_vals.points)-1,-1,-1):
        
        if(i == 0 and dont_keep_before_first):
            to_remove = pt.points()
            var = 0
            while(var <= i_vals.points[i].x):
                index = pt.point(var)
                to_remove.add(index)
                var += 1
            
            airfoil.points = [i for j, i in enumerate(airfoil.points) if j not in to_remove.x()]
            
            airfoil.add(i_pts.points[i],location = 0)
            
            removed += len(to_remove.points) - 1
            
            to_remove = pt.points()
            
            var = i_vals.points[i].y
            
            while(var - removed <= len(airfoil.points)):
                index = pt.point(var - removed)
                to_remove.add(index)
                var += 1
            
            airfoil.points = [i for j, i in enumerate(airfoil.points) if j not in to_remove.x()]
            
            airfoil.add(i_pts.points[i])
            
            removed += len(to_remove.points) - 1
            
        else:
            to_remove = pt.points()
            var = i_vals.points[i].x
            while(var <= i_vals.points[i].y):
                index = pt.point(var - removed)
                to_remove.add(index)
                var += 1
            airfoil.points = [i for j, i in enumerate(airfoil.points) if j not in to_remove.x()]
            
            airfoil.add(i_pts.points[i],location = i_vals.points[i].x)
            
            removed += len(to_remove.points) - 1
            
    return airfoil



def Kerf_pt(airfoil,kerf,dont_keep_b4_1st = True):
    
    x = airfoil.x(); y = airfoil.y()
    
    x0 = x[:len(x)]; y0 = y[:len(y)];
    x1 = x[1:]; y1 = y[1:];
    
    
    xm = [(one + two) / 2 for one,two in zip(x0,x1)]
    ym = [(one + two) / 2 for one,two in zip(y0,y1)]
    
    
    dx = [one - zero for one,zero in zip(x1,x0)]
    dy = [one - zero for one,zero in zip(y1,y0)]
    
    kx = [mid + ( kerf * y ) / ( x**2 + y**2 )**0.5 for mid,x,y in zip(xm,dx,dy)]
    ky = [mid - ( kerf * x ) / ( x**2 + y**2 )**0.5 for mid,x,y in zip(ym,dx,dy)]
    
    k = pt.points(kx,ky)
#    k.show()
    
    k = remove_intersections(k,dont_keep_b4_1st)
    
    if(airfoil.thickness() <= -kerf * 2):
        k = pt.points()
    
    return k



#import File_handler as fh
#import Airfoil_generator as ag
#af = fh.FileReader('/home/benjamin/Desktop/foilshapes/E335.txt')
#
#af = ag.change_chord(af,12)
#
#kerf = Kerf_pt(af,-0.3)
#
#kerf.plot()
#
#
#
#
#kerf = ag.change_chord(kerf,1,12)
#fh.SeligFileWriter('kerfed',kerf)
    



