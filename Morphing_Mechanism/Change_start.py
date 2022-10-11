#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 18 11:13:06 2019

@author: benjamin
"""

import point as pt

import kerf as k

from math import cos,sin,radians,acos,degrees

# find the angle from the positive x axis from main pt to sec_pt
def find_theta(main_pt,sec_pt):
    
    d = pt.point()
    d.x = sec_pt.x - main_pt.x; d.y = sec_pt.y - main_pt.y;
    
    h = ( d.x**2 + d.y**2 )**0.5
    
    rad = acos(d.x / h)
    
    theta = degrees(rad)
    
    return theta

def find_point(main_pt,hyp,theta):
    new_pt = pt.point()
    
    print(radians(theta))
    print(cos(radians(theta)))
    new_pt.x = main_pt.x + hyp * cos( radians(theta) )
    
    new_pt.y = main_pt.y + hyp * sin( radians(theta) )
    
    return new_pt

# change the starting value for the airfoil(running in the direction defined).
#  starting value us determined to be the first point encountered when rotating around the 
#  quarter chord point from the angle given in the direction given. this value is
#  set as the starting value
def change_start(airfoil,theta = 0,theta_true = False, ccw = True):
    
    airfoil.makeWhole()
    
    if(not airfoil.ccw):
        ccw = not ccw
    
    changed = pt.points()
    
    qc_pt = pt.point()
    
    qc_pt.x = airfoil.chord() / 4
    
    qc_pt.y = ( airfoil.interp(qc_pt.x) + airfoil.interp(qc_pt.x,False) ) / 2
    
    
    
    xmax_pt = pt.point()
    
    xmax_pt.set_point(airfoil.points[airfoil.x_max_index])
    
    if(theta_true):
        dtheta = find_theta(qc_pt,xmax_pt)
        
        theta -= dtheta
    
    
    
    
    
    angle_pt = find_point(qc_pt,airfoil.chord(),theta)
    
    start_index = 0
    
    for i in range(len(airfoil.points)-1):
        
        result = k.intersect(qc_pt,angle_pt,airfoil.points[i],airfoil.points[i+1])
        
        if(type(result) == pt.point):
            
            start_index = i + 1
            if(not ccw):
                start_index -= 1
    
    
    for i in range(len(airfoil.points) - start_index):
        changed.add(airfoil.points[i + start_index])
    
    for i in range( start_index ):
        changed.add(airfoil.points[i + 1])
    
    print(theta)
    return changed,qc_pt,xmax_pt




#import File_handler as fh
#
#af = fh.FileReader('/home/benjamin/Desktop/foilshapes/E335.txt')
#
##ch = change_start(af,60)
##ch.plot([ch.points[0]], display = ['','ro'])
##ch = change_start(af,90)
##ch.plot([ch.points[0]], display = ['','ro'])
##ch = change_start(af,30)
##ch.plot([ch.points[0]], display = ['','ro'])
##ch = change_start(af,270)
##ch.plot([ch.points[0]], display = ['','ro'])
##ch = change_start(af,-1)
##ch.plot([ch.points[0]], display = ['','ro'])
##ch = change_start(af,-90)
##ch.plot([ch.points[0]], display = ['','ro'])
#ch,qc,xm = change_start(af,180)
##fir = pt.points(); fir.add(qc); fir.add(xm);
##sec = pt.points(); sec.add(ch.points[0]); sec.add(xm)
##ch.plot([qc,ch.points[0],fir,sec], display = ['','bo','ro'])
#ch.plot([qc,ch.points[0]], display = ['','bo','ro'])
#
#






#import geometry as geo
##af = geo.makeWhole(af)
#af.makeWhole()
#
#af.show()
#print()
#af.add(pt.point(1,1,1))
#
##af = geo.makeWhole(af)
#af.makeWhole()
#
#af.show()