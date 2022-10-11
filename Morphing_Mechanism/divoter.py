#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug  9 12:57:52 2019

@author: benjamin
"""

import point as pt
import Control_horn as ch
import kerf as k



def Halfcircle(center,diameter,num = 30,isTop = True):
    
    # use circle function from Control_horn
    
    # define as semicircle
    percent_of_circle = 1.0
    
    # set to begin at far right
    theta = 0
    
    # make clockwise so semicircle is concave into bezier spine
    isccw = False
    
    # run the function
    hc = ch.circle(center, diameter,percent_of_circle,theta,num,isccw)
    
    # return semicircle
    return hc

def Divoter(airfoil,start_percentage,end_percentage,period,every = 1,diameter = 2 / 25.4,isFirstup = True):
    
    # create points variable
    div = pt.points()
    
    # determine start and end, with period
    
    # find minimum x value
    min_x = airfoil.points[airfoil.x_min_index].x
    
    # get chord
    c = airfoil.chord()
    
    # find T in terms of chord given
    T = period * c
    
    # give tnum
    # value is set in Bezier function, so it is introduced here as final variable
    tnum = 30
    
    # find half T
    hT = T / 2
    
    # find start
    start = start_percentage * c + min_x + hT
    
    # find number half periods
    num_h = ( end_percentage - start_percentage ) // (period/2)
    
    # find end
    end = start + (num_h-2) * hT
    
    # determine where to start locating max's
    # set ind
    ind = 0
    
#    # while loop finds start of end
#    while(airfoil.points[ind].x >= end):
#        ind += 1
    
    # set list of points' for each max value
    circles = []
    
    # set list of points for each max point
    maxes = pt.points()
    
    # define number of circles to make on top and bottom
    num_c_top = int(num_h // 2)
    num_c_bot = int(num_h // 2)
    
    if(isFirstup and num_h % 2 == 0):
        num_c_bot -= 1
    
    elif(not isFirstup and num_h % 2 == 0):
        num_c_top -= 1
    
    
    # find max points, and make corresponding half circles
    # run through each period
    for i in range(num_c_top):
        
        xval = start + i * T
        yval = airfoil.interp(xval,True)
        
        max_pt = pt.point(xval,yval)
        
        maxes.add(max_pt)
        circles.append(Halfcircle(max_pt,diameter,tnum,True))
    
    # reverse the order of the circles so they go in ccw fashion around the airfoil
    circles.reverse()
    
    num_top = len(circles)
    
    for i in range(num_c_bot):
        
        xval = start + (i+0.5) * T
        yval = airfoil.interp(xval,False)
        
        max_pt = pt.point(xval,yval)
        
        maxes.add(max_pt)
        circles.append(Halfcircle(max_pt,diameter,tnum,False))
    
    
    a_count = 0; 
    
    
    # insert circles to airfoil
    for i in range(len(circles)):
        
        c_start = 900; c_end = 900;
        ipt_st = pt.point(); ipt_en = pt.point();
        
        for j in range(len(airfoil.points)-1):
            
            for l in range(len(circles[i].points)-1):
                
                itsc_pt = k.intersect(airfoil.points[j],airfoil.points[j+1],circles[i].points[l],circles[i].points[l+1])
                
                if(type(itsc_pt) != str and c_start == 900):
                    c_start = l+1
                    ipt_st = itsc_pt
                    
                elif(type(itsc_pt) != str and c_start != 900):
                    c_end = l
                    ipt_en = itsc_pt
        
        
        
        c_count = c_start;
#        print(ipt_st.x , airfoil.points[a_count].x , i , num_top ) 
        while( ipt_st.x <  airfoil.points[a_count].x and i < num_top):
            
            div.add(airfoil.points[a_count])
            a_count += 1
        
        while( ipt_st.x >  airfoil.points[a_count].x and i >= num_top):
            
            div.add(airfoil.points[a_count])
            a_count += 1
        
        # add interpolated point
        div.add(ipt_st)
        
        while( c_count <= c_end and i < num_top):
            div.add(circles[i].points[c_count])
            c_count += 1
        
        while( c_count <= c_end and i >= num_top):
            div.add(circles[i].points[c_count])
            c_count += 1
        
        # add interpolated point
        div.add(ipt_en)
        
        while( ipt_en.x < airfoil.points[a_count].x and i < num_top):
            a_count += 1
        
        while( ipt_en.x > airfoil.points[a_count].x and i >= num_top):
            a_count += 1
    
    while(a_count < len(airfoil.points)):
        div.add(airfoil.points[a_count])
        a_count += 1
            
        
    
    return div



#import File_handler as fh

#af = fh.FileReader('/home/benjamin/Desktop/foilshapes/selig bezier NACA 0012.dat')
#nb = fh.FileReader('/home/benjamin/Desktop/foilshapes/E335.txt')
#af.change_chord(12,1); nb.change_chord(12,1);
#d = Divoter(af,0.5,0.95,0.091,1,2/25.4,True)
##af.plot()
##d.plot(axiss = [9.75,9.9,0.175,0.25])
##d.plot(axiss = [7,7.2,-0.6,-0.5])
##d.plot(axiss = [8.1,8.3,-0.6,-0.4])
#d.plot()
##d.plot(circles) 
#
#d.makeWhole()
#d.change_chord(1,12)
#
#fh.FileWriter('divoted_foil',d,isSelig = True)




#af = fh.FileReader('/home/benjamin/Desktop/selig bezier NACA 2412.dat')
#af.change_chord(12,1);
#
#
#d = Divoter(af,0.5,0.8,0.091,1,2/25.4,True)
#
#d.plot()
##d.plot(circles) 
#
#d.makeWhole()
#d.change_chord(1,12)
#
#fh.FileWriter('divoted_foil',d,isSelig = True)




#af = fh.FileReader('/home/benjamin/Desktop/foilshapes/selig bezier NACA 0012.dat' ) # ('/home/benjamin/Desktop/skinned.dat') # 
#ag = af.split()
#af.change_chord(12,1);
#
#
#d = Divoter(af,0.5,0.95,0.091,1,2/25.4,True)
#ag.plot()
#d.plot([ag])
#
#d.makeWhole()
#d.change_chord(1,12)
#
#fh.FileWriter('divoted_foil',d,isSelig = True)


#for i in range(4,len(circles)):
#    af.plot(other_plotters = circles,display = [''], axiss = [ maxes.points[i].x - 0.1, maxes.points[i].x + 0.1, maxes.points[i].y - 0.1, maxes.points[i].y + 0.1])



#d = 2 / 25.4
#c = pt.point(0,1,0)
#
#hc = Halfcircle(c,d)
#
#hc.plot()
