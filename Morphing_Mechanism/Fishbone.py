#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 17 14:21:22 2019

@author: benjamin


if ledge depth is > 1, create over next ledge...
"""
import point as pt

def BoneMaker( airfoil, startperc, endperc, boneWidth, gapWidth,
              spineWidth, tribase, isTlocal, ledgeWidth, ledgeDepth, ledgeIsFlat, skipPoints, percent_right):
    
    # define chord
    chord = airfoil.chord()
    
    # define thickness
    
    thickness = airfoil.thickness()
    
    # if no start percentage is given( or 0 is given) start at the location of max thickness
    if(startperc == 0):
        startperc = airfoil.x()[airfoil.y().index(max(airfoil.y()))] / chord
    
    
    # if no base percentage is given for the triBone values, assume 100%
    if(tribase == 0):
        tribase = 1.0
        
    # tribase is converted for mathematical purposes
    tribase = 1 - tribase
    
    # define the bone output points objects
    bone = pt.points()
    
    # find index of min x value
    minxind = airfoil.x_min_index
    
    # define percent_left
    percent_left = 1 - percent_right
    
    # determine mid_spine points
    
    mid_spine = pt.points()
    
    for i in range(minxind):
        
        new_spinept = pt.point()
        
        ytop = airfoil.points[i].y
        
        ybot = airfoil.interp(airfoil.points[i].x,False)
        
        new_spinept.y = ( ytop + ybot ) / 2
        
        new_spinept.x = airfoil.points[i].x
        
        # save the midspine values to the lists
        
        mid_spine.add(new_spinept)
    
    # set up x and tpoints to be used to define the bone points
    # set up thickness points for each x_points
    xt_top = pt.points()
    
    xt_bot = pt.points()
    
    #xposition value is defined to begin the while loop
    x_position = endperc * chord
    
    
    # while loop is set so that when the xposition is less than the start percentage and addt'l correction values, it will end
    while(x_position >= ((startperc + gapWidth + boneWidth)*chord)):
        
        # a point is placed in position for a bone section
        top_bone_pt = pt.point(x_position)
        bot_bone_pt = pt.point(x_position)
        
        
        # the thickness value is added for the x_points based on whether the boneMade is Tlocal or not
        if(isTlocal):
            
            top_bone_pt.y = spineWidth * airfoil.interp(x_position) + mid_spine.interp(x_position)
            bot_bone_pt.y = mid_spine.interp(x_position) - spineWidth * airfoil.interp(x_position)
            
        else:
            top_bone_pt.y = spineWidth * thickness/2 + mid_spine.interp(x_position)
            bot_bone_pt.y = mid_spine.interp(x_position) - spineWidth * thickness/2
        
        
        # add point to xt lists
        xt_top.add(top_bone_pt)
        xt_bot.add(bot_bone_pt)
        
        # a gap is accounted for
        x_position -= gapWidth*chord
        
        
        # a second point is given
        top_bone_pt = pt.point(x_position)
        bot_bone_pt = pt.point(x_position)
        
        
        # the thickness value is added for the x_points based on whether the boneMade is Tlocal or not
        if(isTlocal):
            top_bone_pt.y = spineWidth * airfoil.interp(x_position) + mid_spine.interp(x_position)
            bot_bone_pt.y = mid_spine.interp(x_position) - spineWidth * airfoil.interp(x_position)
        else:
            top_bone_pt.y = spineWidth * thickness/2 + mid_spine.interp(x_position)
            bot_bone_pt.y = mid_spine.interp(x_position) - spineWidth * thickness/2
        
        
        # add point to xt lists
        xt_top.add(top_bone_pt)
        xt_bot.add(bot_bone_pt)
        
        
        # a bone is accounted for
        x_position -= boneWidth*chord
    
    # with the x values found, a bone structure can begin
    index = 0
    
    
    # a for loop is developed which will add the bone points for the top of the airfoil
    for i in range(0,len(xt_top.points),2):
        
        # a while loop is begun which adds the points till the next coordinate x_points
        while(airfoil.points[index].x > xt_top.points[i].x):
            
            if(skipPoints and i != 0 and i != len(xt_top.points)-2):
                # the index is increased in order to approach endper*chord
                index += 1
            
            else:
                # the values are added to the bone lists bx and by
                bone.add(airfoil.points[index])
                # the index is increased in order to approach endper*chord
                index += 1
        
        new_point = pt.point()
        
        # if no ledge is to be added, finish the bone
        
        if(ledgeWidth == 0.0 or ledgeDepth == 0.0):
            # third, the interpolation point is added which is for the top right of the bone
            # the coordinate is passed to the bone lists for the first interpolation ridge
            new_point = pt.point()
            new_point.x = xt_top.points[i].x; new_point.y = airfoil.interp(xt_top.points[i].x)
            bone.add(new_point)        
        
        # if a flat ledge is to be added
        elif(ledgeIsFlat):
            
            # third, the top right point of the ledge is added
            bone.add(pt.point(xt_top.points[i].x, airfoil.interp(xt_top.points[i].x)))
            
            # define the next point
            new_point = pt.point()
            new_point.x = xt_top.points[i].x - percent_left * ledgeWidth * gapWidth * chord
            new_point.y = airfoil.interp(xt_top.points[i].x) 
            
            # add the point to the ledge
            bone.add(new_point)
            
            # add the bottom left point
            bone.add(pt.point(new_point.x,airfoil.interp(xt_top.points[i].x) - ledgeDepth * thickness))
            
            # add the final bottom right point
            bone.add(pt.point(xt_top.points[i].x,airfoil.interp(xt_top.points[i].x) - ledgeDepth * thickness))
         
        
        # if an airfoil continuous ledge is to be added
        else:
            # third, the top right point of the ledge is added
            bone.add(pt.point(xt_top.points[i].x, airfoil.interp(xt_top.points[i].x)))
            
            # define the next point
            new_point = pt.point()
            new_point.x = xt_top.points[i].x - percent_left * ledgeWidth * gapWidth * chord
            new_point.y = airfoil.interp(new_point.x) 
            
            # add the point to the ledge
            bone.add(new_point)
            
            # add the bottom left point
            bone.add(pt.point(new_point.x,airfoil.interp(new_point.x) - ledgeDepth * thickness))
            
            # add the final bottom right point
            bone.add(pt.point(xt_top.points[i].x,airfoil.interp(xt_top.points[i].x) - ledgeDepth * thickness))
        
        
        # the coordinates for the first spine deep point are added
        new_point = pt.point()
        new_point.x = xt_top.points[i].x + tribase*boneWidth*chord/2; new_point.y = xt_top.points[i].y
        bone.add(new_point)
        
        # secondly, the coordinates for the second spine deep point are added
        new_point = pt.point()
        new_point.x = xt_top.points[i+1].x - tribase*boneWidth*chord/2; new_point.y = xt_top.points[i+1].y
        bone.add(new_point)
        
        # if no ledge is to be added, finish the bone
        
        if(ledgeWidth == 0.0 or ledgeDepth == 0.0):
            # third, the interpolation point is added which is for the top right of the bone
            new_point = pt.point()
            new_point.x = xt_top.points[i+1].x; new_point.y = airfoil.interp(xt_top.points[i+1].x)
            bone.add(new_point)        
        
        # if a flat ledge is to be added
        elif(ledgeIsFlat):
            # third, the bottom left point of the ledge is added
            new_point = pt.point()
            new_point.x = xt_top.points[i+1].x; new_point.y = airfoil.interp(xt_top.points[i+1].x) - ledgeDepth * thickness
            bone.add(new_point)
            
            # define the next point
            new_point = pt.point()
            new_point.x = xt_top.points[i+1].x + percent_right * ledgeWidth * gapWidth * chord
            new_point.y = bone.points[len(bone.points)-1].y
            
            # add the point to the ledge
            bone.add(new_point)
            
            # add the top right point
            bone.add(pt.point(new_point.x,airfoil.interp(xt_top.points[i+1].x)))
            
            # add the final top left point
            bone.add(pt.point(xt_top.points[i+1].x, airfoil.interp(xt_top.points[i+1].x)))
         
        
        # if an airfoil continuous ledge is to be added
        else:
            # third, the bottom left point of the ledge is added
            bone.add(pt.point(xt_top.points[i+1].x, airfoil.interp(xt_top.points[i+1].x) - ledgeDepth * thickness))
            
            # define the next point
            new_point = pt.point()
            new_point.x = xt_top.points[i+1].x + percent_right * ledgeWidth * gapWidth * chord
            new_point.y = airfoil.interp(new_point.x) - ledgeDepth * thickness
            
            # add the point to the ledge
            bone.add(new_point)
            
            # add the top right point
            bone.add(pt.point(new_point.x,airfoil.interp(new_point.x)))
            
            # add the final top left point
            bone.add(pt.point(xt_top.points[i+1].x,airfoil.interp(xt_top.points[i+1].x)))
        
        
        # a while loop runs through the x y coordinates which have now been superceded by the bone spine
        while(airfoil.points[index].x > xt_top.points[i+1].x):
            index += 1
    
    
    
    # a for loop is developed which will add the bone points for the bottom of the airfoil
    for i in range(len(xt_bot.points)-1,-1,-2):
        
        # a while loop is begun which adds the points till the next coordinate x_points
        while(airfoil.points[index].x < xt_bot.points[i].x):
            
            if(skipPoints and i != 1 and i != len(xt_top.points)-1):
                # the index is increased in order to approach endper*chord
                index += 1
            
            else:
                # the values are added to the bone lists bx and by
                bone.add(airfoil.points[index])
                # the index is increased in order to approach endper*chord
                index += 1
        
        new_point = pt.point()
        
        # the coordinate is passed to the bone lists for the first interpolation ridge
        new_point.x = xt_bot.points[i].x; new_point.y = airfoil.interp(xt_bot.points[i].x,False)
        bone.add(new_point)
        
        
        if( ( ledgeWidth == 0.0 or ledgeDepth == 0.0 ) and ledgeIsFlat):
            
            # define bottom right point
            new_point = pt.point()
            new_point.x = xt_bot.points[i].x + ledgeWidth * gapWidth * chord
            new_point.y = bone.points[len(bone.points)-1].y
            
            # add it
            bone.add(new_point)
            
            # add top right point
            bone.add(pt.point(new_point.x, airfoil.interp(xt_bot.points[i].x,False) + ledgeDepth * thickness))
            
            # add top left point
            bone.add(pt.point(xt_bot.points[i].x, bone.points[len(bone.points)-1].y))
            
        else:
            
            # define bottome right point
            new_point = pt.point()
            new_point.x = xt_bot.points[i].x + percent_right * ledgeWidth * gapWidth * chord; new_point.y = airfoil.interp(new_point.x,False)
            
            # add it
            bone.add(new_point)
            
            # add top right point
            new_point = pt.point()
            new_point.x = xt_bot.points[i].x + percent_right * ledgeWidth * gapWidth * chord; new_point.y = airfoil.interp(new_point.x,False) + ledgeDepth * thickness
            bone.add(new_point)
            
            # add top left point
            new_point = pt.point()
            new_point.x = xt_bot.points[i].x; new_point.y = airfoil.interp(xt_bot.points[i].x,False) + ledgeDepth * thickness
            bone.add(new_point)
        
        # the coordinates for the first spine deep point are added
        new_point = pt.point()
        new_point.x = xt_bot.points[i].x - tribase*boneWidth*chord/2; new_point.y = xt_bot.points[i].y
        bone.add(new_point)
        
        # secondly, the coordinates for the second spine deep point are added
        bone.add(pt.point(xt_bot.points[i-1].x + tribase*boneWidth*chord/2 , xt_bot.points[i-1].y))
        
        
        if( ledgeWidth == 0.0 or ledgeDepth == 0.0 ):
            
            # third, the interpolation point is added which is for the bottom left of the bone
            bone.add(pt.point(xt_bot.points[i-1].x , airfoil.interp(xt_bot.points[i-1].x,False)))
        
        elif(ledgeIsFlat):
            
            # define top right point
            new_point = pt.point()
            new_point.x = xt_bot.points[i-1].x #+ ledgeWidth * gapWidth * chord
            new_point.y = airfoil.interp(xt_bot.points[i-1].x,False) + ledgeDepth * thickness
            
            # add it
            bone.add(new_point)
            
            # add top left point
            new_point = pt.point()
            new_point.x = xt_bot.points[i-1].x - percent_left * ledgeWidth * gapWidth * chord
            new_point.y = airfoil.interp(xt_bot.points[i-1].x,False) + ledgeDepth * thickness
            
            # add it
            bone.add(new_point)
            
            # add bottom left point
            bone.add(pt.point(new_point.x, airfoil.interp(xt_bot.points[i-1].x,False) ))
            
            # finally add bottom right pt
            bone.add(pt.point(xt_bot.points[i-1].x , airfoil.interp(xt_bot.points[i-1].x,False)))
            
        else:
            
            # define top right point
            new_point = pt.point()
            new_point.x = xt_bot.points[i-1].x #+ ledgeWidth * gapWidth * chord
            new_point.y = airfoil.interp(xt_bot.points[i-1].x,False) + ledgeDepth * thickness
            
            # add it
            bone.add(new_point)
            
            # add top left point
            new_point = pt.point()
            new_point.x = xt_bot.points[i-1].x - percent_left * ledgeWidth * gapWidth * chord
            new_point.y = airfoil.interp(new_point.x,False) + ledgeDepth * thickness
            
            # add it
            bone.add(new_point)
            
            # add bottom left point
            bone.add(pt.point(new_point.x, airfoil.interp(new_point.x,False) ))
            
            # finally add bottom right pt
            bone.add(pt.point(xt_bot.points[i-1].x , airfoil.interp(xt_bot.points[i-1].x,False)))
            
        
        # a while loop runs through the x y coordinates which have now been superceded by the bone spine
        while(airfoil.points[index].x < xt_bot.points[i-1].x):
            index += 1
    
    # a while loop runs through the remaining x y coordinate points
    while(index < len(airfoil.points)):
        bone.add(airfoil.points[index])
        index += 1
    
    # check if airfoil shape is whole
    bone.makeWhole()
    
    # the values are returned, the function ends
    return bone






#import File_handler as fh
#
#airfoil = fh.FileReader('/home/benjamin/Desktop/foilshapes/E335.txt')
#
#thickness = 0.15
#chord = 1
#startperc = 0.5
#endperc = 0.95
#boneWidth = 0.025
#gapWidth = 0.025
#spineWidth = 0.1
#tribase = 1.0
#isTlocal = False
## % of gapWidth
#ledgeWidth =  0.9 # ( gapWidth * 6 - 0.03 ) / ( gapWidth * 6 )
##print(ledgeWidth)
## % of thickness
#ledgeDepth =  0.05 # 0.025 / 6.0 / 0.15
#ledgeIsFlat = False
#percent_r = .5
#
#fishbone = BoneMaker( airfoil, startperc, endperc, boneWidth, gapWidth,
#              spineWidth, tribase, isTlocal, ledgeWidth, ledgeDepth, ledgeIsFlat,False,percent_r)
#
#
#fishbone.plot()
#fishbone.plot(axis = [0.55,0.65,-0.1,0.1])
##fishbone.plot(axis = [0.95,1.05,-0.01,0.01])




