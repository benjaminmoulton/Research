#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 19 12:56:13 2019

@author: benjamin
"""

import point as pt
import geometry as geo
import File_handler as fh
import kerf as k

def length(data):
    
    data.makeWhole()
    
    length = 0
    
    for i in range(len(data.points)-1):
        
        dl = (  ( data.points[i+1].x - data.points[i].x )**2 +
              ( data.points[i+1].y - data.points[i].y )**2 +
              ( data.points[i+1].z - data.points[i].z )**2  )**0.5
        
        length += dl
        
    
    return length

skin_thickness = 0.4 / 25.4


# f = fh.FileReader('/home/benjamin/Desktop/selig fishbone NACA 2412.dat')
epp = fh.FileReader('/home/ben/Desktop/foilshapes/E335_200pts.dat')
print('eppler %%thickness is : %f' % epp.percentThickness())

# f.change_chord(12,1)
epp.change_chord(12,1)

# f_a = geo.Area(f)
epp_a = geo.Area(epp)

print(epp_a)

LE = geo.sectionMaker(f,percentage_start = 0.0,percentage_end = 0.5)
TE = geo.sectionMaker(f,percentage_start = 0.95,percentage_end = 1.0)
#LE_skin = k.Kerf_pt(LE,skin_thickness)
#LE.plot([LE_skin,TE])


l_LE = length(LE)
l_TE = length(TE)
print(l_LE)

#line = pt.points()
#
#line.add(pt.point(0,0,0))
#line.add(pt.point(1,0,0))
#line.add(pt.point(1,1,0))
#line.add(pt.point(0,1,0))
#line.add(pt.point(0,0,0))

border_LE = l_LE * skin_thickness
border_TE = l_TE * skin_thickness

netLE = geo.Area(LE) - border_LE
netTE = geo.Area(TE) - border_TE

net_f_a = f_a - netLE - netTE

print(net_f_a,epp_a)
print('%% material fish : %f' % (net_f_a / epp_a))

mid = geo.sectionMaker(epp,0.5,0.95)
l_mid = length(mid)
border_mid = l_mid * skin_thickness

net_f_scales = border_mid + border_LE + border_TE
print('%% material scales : %f' % (net_f_scales / epp_a))

border_epp = length(epp) * skin_thickness

print('%% material empty : %f' % (border_epp / epp_a))















