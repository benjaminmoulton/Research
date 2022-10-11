#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 21 13:27:17 2019

@author: benjamin
"""

import point as pt
import File_handler as fh

def length(data):
    
    length = 0
    
    for i in range(len(data.points)-1):
        
        dl = (  ( data.points[i+1].x - data.points[i].x )**2 +
              ( data.points[i+1].y - data.points[i].y )**2 +
              ( data.points[i+1].z - data.points[i].z )**2  )**0.5
        
        length += dl
        
    
    return length

af = fh.FileReader('/home/benjamin/Desktop/foilshapes/E335_200pts.dat'); af.change_chord(2,1)

line = pt.points()
line.add(pt.point(0.5*af.chord(),-0.2*af.chord()))
line.add(pt.point(0.5*af.chord(),0.2*af.chord()))


af.plot([line,af.points[50]])

af_50 = af.split(0,50)
#print(len(af_50.points))

lenth = length(af_50)
print(lenth)
l = 0.49207270376405715

print(lenth / l)

c = 12

size = c * l

print(size)
















