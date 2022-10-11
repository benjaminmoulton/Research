#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 12 10:01:57 2019

@author: benjamin
"""

import point as pt
import File_handler as fh




def doveTail(startpt,endpt,head_width,neck_width,tail_width):
    
    # initialize dovetail points
    dove_pts = pt.points()
    
    # set a floating point to position doveTail coordinates
    fl = pt.point()
    fl.set_point(startpt)
    
    # find distance from neck to head in x dir
    overhang =  - ( head_width - neck_width) / 2
    
    # while loop that adds dovetails sequentially till only a headwidth is left
    while(endpt.x - fl.x > head_width + neck_width):
        
        # add float point at top of bottom head
        np = pt.point(fl.x,fl.y,fl.z);dove_pts.add(np)
        
        # move across head
        fl.move(head_width,0)
        
        # add to dove_pts
        np = pt.point(fl.x,fl.y,fl.z); dove_pts.add(np)
        
        # move to left of top head
        fl.move(overhang,tail_width)
        
        # add to dove_pts
        np = pt.point(fl.x,fl.y,fl.z); dove_pts.add(np)
        
        # move across top of neck
        fl.move(head_width,0)
        
        # add to dove_pts
        np = pt.point(fl.x,fl.y,fl.z); dove_pts.add(np)
        
        # move to bottom left of bottom head
        fl.move(overhang,-tail_width)
    
    # add final fl pt
    dove_pts.add(fl)
    
    # continue to end
    dove_pts.add(endpt)
    
    
    
    
    return dove_pts


"""
switch y and z so that z is value that changes with dovetail in out tail width
make it so y values and z values  follow start point and end point y values and z values


"""
#start = pt.point()
#end = pt.point(12,0,0)
#
#h = 1
#n = 0.5
#t = 0.25
#
#dt = doveTail(start,end,h,n,t)
#
#dt.plot()
#
#fh.FileWriter('dovetail',dt,isSelig = True,skipFirst = False)


