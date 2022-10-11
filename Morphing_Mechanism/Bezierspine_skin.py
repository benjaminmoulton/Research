#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun  3 09:00:25 2019

@author: benjamin

"""


import point as pt
import kerf as k
import Bezierspine as bz
 # adds top because the bezier shape being passed in has a wierd 0.1 at the start and end

#def distances(airfoil):
#    
#    
#    oth_distances = pt.points()
#    
#    for i in range(len(airfoil.points)-1):
#        
#        dist = ( (airfoil.points[i].x - airfoil.points[i+1].x)**2 + (airfoil.points[i].y -
#                airfoil.points[i+1].y)**2 + (airfoil.points[i].z - airfoil.points[i+1].z)**2 )**0.5
#        
#        oth_distances.add(pt.point(airfoil.points[i].x,dist))
#    
#    return oth_distances
#    
#
#def ydistances(airfoil):
#    
#    
#    oth_distances = pt.points()
#    
#    for i in range(len(airfoil.points)-1):
#        
#        dist = abs( airfoil.points[i].y - airfoil.points[i+1].y )
#        
#        oth_distances.add(pt.point(airfoil.points[i].x,dist))
#    
#    return oth_distances


def Bezier_skin(airfoil, start_percentage, end_percentage, P1_percentage, 
                 P2_percentage, T, spine_width_LE, spine_width_TE, first_up,
                 is_Tlocal,skin):
    
    airfoil.makeWhole()
    
    bezi = bz.Bezier_Spine(airfoil, start_percentage, end_percentage, P1_percentage, 
                 P2_percentage, T, spine_width_LE, spine_width_TE, first_up,
                 is_Tlocal)
    
    bs = pt.points()
    
    bs = k.Kerf_pt(airfoil,-skin )
    
    bs.points.pop(0)
    
    bs.makeWhole()
    
    other = bs.split()
    
    min_x = other.points[other.x_min_index].x
    
##    other.show()
#    for i in range(len(other.points)):
#        other.points[i].x -= min_x
##    other.shift(-min_x)
#    other.points[0].x += min_x
##    print();other.show()
    
    dechorded = airfoil.split()
    dechorded.unityChord()
    
    c = airfoil.chord()
    c_other = other.chord()
    
    start_percentage *= (c / c_other)
#    print('start perc %f' % start_percentage)
#    start_percentage -= min_x / c_other
    
    T *= c / c_other
#    print('T %f' % T)
    spine_width_LE *= c / c_other
    spine_width_TE *= c / c_other
    
    other = bz.Bezier_Spine(other, start_percentage, end_percentage, P1_percentage, 
                 P2_percentage, T, spine_width_LE, spine_width_TE, first_up,
                 is_Tlocal)
    
    
#    for i in range(len(other.points)):
#        other.points[i].x += min_x
    
    
    m,end_inner = other.ymax_in_range()
    
    
    m,start_inner = other.ymax_in_range(isTop = False)
    
    
    mb,start_outer = bezi.ymax_in_range()
    start_outer += 1
    
    
    mb,end_outer = bezi.ymax_in_range(isTop = False)
    end_outer -= 1
    
    
    skinned = pt.points()
    
    for i in range(end_inner + 1):
        skinned.add(other.points[i])
    
    for i in range(start_outer,end_outer+1):
        skinned.add(bezi.points[i])
    
    while(start_inner < len(other.points)):
        skinned.add(other.points[start_inner])
        start_inner += 1
    
    
    
    
    
    return skinned, start_percentage, T,other.chord()






#
#import File_handler as gc
#
#airfoil = gc.SeligFileReader('/home/benjamin/Desktop/foilshapes/E335.txt')
#
#airfoil.change_chord(12,1)
#
#
#start_perc = 0.5
#end_perc = 0.95
#
#T = 0.091
#
#P1_perc = 0.75   #1 # straight up and down # 2/3 # get closer to osci # 0.5 # original #
#P2_perc = 0.25   #0 # straight up and down # 1/3 # get closer to osci # 0.5 # original #
#
#spine_width_LE = 0.1
#spine_width_TE = 0.1  # spine_width_LE * 0.5
#
#first_up = True
#
#is_Tlocal = False
#
#is_skinned = True
#skin = 0.23
#
#ski,st,t,o = Bezier_skin(airfoil,start_perc,end_perc,P1_perc,P2_perc, T,spine_width_LE,spine_width_TE,first_up, is_Tlocal, skin)
#
#o.plot()
#
#ski.plot()
#
#
#
#ski.change_chord(1,12)
#gc.FileWriter('skinned',ski,True)
#
#gc.FileWriter('ofile',o,True)



