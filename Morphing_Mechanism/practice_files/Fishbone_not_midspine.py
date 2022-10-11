#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 17 14:21:22 2019

@author: benjamin
"""

def BoneMaker( x, y, thickness, chord, startperc, endperc, boneWidth, gapWidth,
              spineWidth, tribase, isTlocal, ledgeWidth, ledgeDepth, ledgeIsFlat):
    
    # import scipy to be used for interpolating values
    from scipy.interpolate import interp1d as interp
    
    
    # if no start percentage is given( or 0 is given) start at the location of max thickness
    if(startperc == 0):
        startperc = x[y.index(max(y))] / chord
    
    
    # if no base percentage is given for the triBone values, assume 100%
    if(tribase == 0):
        tribase = 1.0
        
    # tribase is converted for mathematical purposes
    tribase = 1 - tribase
    
    # define the bone output arrays
    bx=[]; by=[]
    
    # The interpolation functions are defined
    top_interpolation = interp(x[:len(x)//2],y[:len(y)//2])
    bottom_interpolation = interp(x[len(x)//2:],y[len(y)//2:])
    
    # set up x points to be used to define the bone points
    x_points = []
    
    # set up thickness points for each x_points
    t = []
    
    #xposition value is defined to begin the while loop
    x_position = endperc * chord
    
    
    
    # while loop is set so that when the xposition is less than the start percentage and addt'l correction values, it will end
    while(x_position >= ((startperc + gapWidth + boneWidth)*chord)):
        
        # a point is placed in position for a bone section
        x_points.append(x_position)
        
        
        # the thickness value is added for the x_points based on whether the boneMade is Tlocal or not
        if(isTlocal):
            
            t.append(spineWidth * top_interpolation(x_position))
            
        else:
            
            t.append(spineWidth * thickness/2)
        
        
        # a gap is accounted for
        x_position -= gapWidth*chord
        
        
        # a second point is given
        x_points.append(x_position)
        
        
        # the thickness value is added for the x_points based on whether the boneMade is Tlocal or not
        if(isTlocal):
            t.append(spineWidth * top_interpolation(x_position))
        else:
            t.append(spineWidth * thickness/2)
        
        
        # a bone is accounted for
        x_position -= boneWidth*chord
    



    # with the x values found, a bone structure can begin
    index = 0
       
    
    
    # a for loop is developed which will add the bone points for the top of the airfoil
    for i in range(0,len(x_points),2):
        
        # a while loop is begun which adds the points till the next coordinate x_points
        while(x[index] > x_points[i]):
            # the values are added to the bone lists bx and by
            bx.append(x[index]); by.append(y[index])
            # the index is increased in order to approach endper*chord
            index += 1
        
        # the coordinate is passed to the bone lists for the first interpolation ridge
        bx.append(x_points[i]); by.append(top_interpolation(x_points[i]))
        
        # the coordinates for the first spine deep point are added
        bx.append(x_points[i] + tribase*boneWidth*chord/2); by.append(t[i])
        
        # secondly, the coordinates for the second spine deep point are added
        bx.append(x_points[i+1] - tribase*boneWidth*chord/2); by.append(t[i+1])
        
        # if no ledge is to be added, finish the bone
        
        if(ledgeWidth == 0.0 or ledgeDepth == 0.0):
            # third, the interpolation point is added which is for the top right of the bone
            bx.append(x_points[i+1]); by.append(top_interpolation(x_points[i+1]))
        
        
        # if a flat ledge is to be added
        elif(ledgeIsFlat):
            # third, the bottom left point of the ledge is added
            bx.append(x_points[i+1]); by.append(top_interpolation(x_points[i+1]) - ledgeDepth * thickness)
            
            # define the next point
            ledge_x = x_points[i+1] + ledgeWidth * gapWidth * chord
            ledge_y = by[len(by)-1]
            
            # add the point to the ledge
            bx.append(ledge_x); by.append(ledge_y)
            
            # add the top right point
            bx.append(ledge_x); by.append(top_interpolation(x_points[i+1]))
            
            # add the final top left point
            bx.append(x_points[i+1]); by.append(top_interpolation(x_points[i+1]))
         
        
        # if an airfoil continuous ledge is to be added
        else:
            # third, the bottom left point of the ledge is added
            bx.append(x_points[i+1]); by.append(top_interpolation(x_points[i+1]) - ledgeDepth * thickness)
            
            # define the next point
            ledge_x = x_points[i+1] + ledgeWidth * gapWidth * chord
            ledge_y = top_interpolation(ledge_x) - ledgeDepth * thickness
            
            # add the point to the ledge
            bx.append(ledge_x); by.append(ledge_y)
            
            # add the top right point
            bx.append(ledge_x); by.append(top_interpolation(ledge_x))
            
            # add the final top left point
            bx.append(x_points[i+1]); by.append(top_interpolation(x_points[i+1]))
        
        
        # a while loop runs through the x y coordinates which have now been superceded by the bone spine
        while(x[index] > x_points[i+1]):
            index += 1
           
    
    # a for loop is developed which will add the bone points for the bottom of the airfoil
    for i in range(len(x_points)-1,-1,-2):
        
        # a while loop is begun which adds the points till the next coordinate x_points
        while(x[index] < x_points[i]):
            # the values are added to the bone lists bx and by
            bx.append(x[index]); by.append(y[index])
            # the index is increased in order to approach endper*chord
            index += 1
        
        # the coordinate is passed to the bone lists for the first interpolation ridge
        bx.append(x_points[i]); by.append(bottom_interpolation(x_points[i]))
        
        
        if( ( ledgeWidth == 0.0 or ledgeDepth == 0.0 ) and ledgeIsFlat):
            
            # define bottom right point
            ledge_x = x_points[i] + ledgeWidth * gapWidth * chord
            ledge_y = by[len(by)-1]
            
            # add it
            bx.append(ledge_x); by.append(ledge_y)
            
            # add top right point
            bx.append(ledge_x); by.append(bottom_interpolation(x_points[i]) + ledgeDepth * thickness)
            
            # add top left point
            bx.append(x_points[i]); by.append(by[len(by)-1])
            
        else:
            
            # define bottome right point
            ledge_x = x_points[i] + ledgeWidth * gapWidth * chord
            ledge_y = bottom_interpolation(ledge_x)
            
            # add it
            bx.append(ledge_x); by.append(ledge_y)
            
            # add top right point
            bx.append(ledge_x); by.append(bottom_interpolation(ledge_x) + ledgeDepth * thickness)
            
            # add top left point
            bx.append(x_points[i]); by.append(bottom_interpolation(x_points[i]) + ledgeDepth * thickness)
        
        # the coordinates for the first spine deep point are added
        bx.append(x_points[i] - tribase*boneWidth*chord/2); by.append(-t[i])
        
        # secondly, the coordinates for the second spine deep point are added
        bx.append(x_points[i-1] + tribase*boneWidth*chord/2); by.append(-t[i-1])
        
        # third, the interpolation point is added which is for the bottom left of the bone
        bx.append(x_points[i-1]); by.append(bottom_interpolation(x_points[i-1]))
        
        # a while loop runs through the x y coordinates which have now been superceded by the bone spine
        while(x[index] < x_points[i-1]):
            index += 1
    
    # a while loop runs through the remaining x y coordinate points
    while(index < len(x)):
        bx.append(x[index]); by.append(y[index])
        index += 1
    
    # the values are returned, the function ends
    return bx,by






#import File_handler as fh
#from matplotlib import pyplot as plt
#
#ax,ay = fh.FileReader('/home/benjamin/Desktop/foilshapes/NACA 0015.txt')
#
#
#thickness = 0.15
#chord = 1
#startperc = 0.3
#endperc = 0.9
#boneWidth = 0.025
#gapWidth = 0.025
#spineWidth = 0.1
#tribase = 1.0
#isTlocal = False
## % of gapWidth
#ledgeWidth = 0.6 # ( gapWidth * 6 - 0.03 ) / ( gapWidth * 6 )
##print(ledgeWidth)
## % of thickness
#ledgeDepth = 0.025 / 6.0 / 0.15
#ledgeIsFlat = False
#
#fx,fy = BoneMaker( ax, ay, thickness, chord, startperc, endperc, boneWidth, gapWidth,
#              spineWidth, tribase, isTlocal, ledgeWidth, ledgeDepth, ledgeIsFlat)
#
#
#
#plt.axis('equal')
##plt.plot(ax,ay)
#plt.plot(fx,fy)
#plt.show()
#
#
#plt.axis([0.6,0.7,0.1,0])
##plt.plot(ax,ay)
#plt.plot(fx,fy)
#plt.show()








