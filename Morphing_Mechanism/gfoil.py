#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 20 16:12:54 2018

@author: benjamin
This creates a gfoil

The user may edit inputs in the getParams() function such that the .txt file can be located 
to input point values (as from a seligdat file).
Else, parameters can be inserted for NACA 4 digit airfoils.
"""

# files are imported for use in the program ::::

import point as pt


# file_handler file consists of functions used in the manipulation of files
#  and file formats.

import File_handler as fh                   # seligdatfile.txt

# NACApoints file consists of functions used in the creation, modification,
#  and use of data points for naca four digit airfoils

import Airfoil_generator as ag

# geometry file consists of functions used in creating the geometry specified
# by the user in the getParams() function

import geometry as geo

# the pyplot class is imported for use in plotting airfoil data to show the user
#  how foilshapes are constructed after input.

from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D



###### add other files here

import Control_horn as ch
import dovetail as dt
from datetime import datetime as date



def getParams():#               
    
    # stm # what vars will modify which foilshape output: self(rootairfoil etc), .txt, .mu
    
    p = {} # 'p' for parameters
    
      #     # root foil shape # SHOULD be a naca four digit
    p['root airfoil'] = 'NACA 2412'
    
      #     # tip foil shape  # SHOULD be a naca four digit
    p['tip airfoil'] = 'NACA 2412' 
    
    p['NACA number of points'] = 100
    
      ###   # offset value in inches of the foam from the XY axis (z value) 
            # as well as how tall the piece of foam is (y height)( for use in the movements of the HotwireCNC) thickness
            # and the length of the foam peice (x height) ( from Trailing Edge to Leading Edge) chord 
            # foam must be positioned 1 inch out from the wire.
    p['offset'] = 2.
    p['depth'] = 2
    p['length'] = 10
    
      ###   # length of the wire from the attachment at the spring to the opposite eyescrew
    p['wire length'] = 10  #31.5 # is what it is like right now
         
      ##    # This is the chord length of the root foil
    p['root chord'] = 6.
           
      ##    # This is the chord length of the tip foil
    p['tip chord'] = 6. 
            
      ##    # This is the degree value of the rotation of the tipfoil downward, or CW on a right wing.
    p['washout'] = 0.  
    
      ##    # This is the degree value, where + is CW(right wing) toward the tail.             
    p['sweep'] = 0. 
                
      ##    # semispan of wing (ie. left or right only)
    p['span'] = 7  
               
      ###   # whether the values are to be in inches. the other option is milimeters
    p['is in inches'] = True  
       
      ##    # Whether the wing to cut is the left wing or not(right)
    p['is left wing'] = False  
      
      ###     # only to be modified if user wants the wire to cut differently than normal
    p['root wire kerf'] = -0.0#5
      #0.02
      ###     # only to be modified if user wants the wire to cut differently than normal
    p['tip wire kerf'] = -0.0#5       
    
    # offset for where the foilshape is? u axis is along the baseshape. determine how to find where the foam begins/ends
    
    
      ###    # if the air foil is a section, and if so, what are the start and end percentages of the airfoil
             #  In percent chord
    p['is section'] = False
    
    p['section start'] = 0.0
    
    p['section end'] = 0.3
    
    
      ###    # whether the airfoil desired will be a fishBac design !!!!!! only for symmetrical airfoils !!!!!!
    p['is fishbone'] = True
    
      ###   # whether the fishbone design will be a default setting, ie. (startperc = location of max thickness,endperc = 0.85,
              # boneWidth = 0.0225,gapWidth = 0.0225,spineWidth = 0.25,isTlocal = False) or not (params to be used are those following)
    p['default fishbone'] = False
    
      ###    # at what % chord should the fishbone start from the Leading Edge
    p['fishbone start'] = .4
    
      ###    # at what % chord should the fishbone end from the Leading Edge
    p['fishbone end'] = .9
    
    """ """
    
    ##    # This is the chord length of the root foil
    p['root chord'] = 1.0
      ##    # This is the chord length of the tip foil
    p['tip chord'] = p['root chord'] 
    """ """
    
    
    chord =  p['root chord']
    smallval = 0.017# / chord
      ###   # how wide should the bones be in % chord       
    p['fishbone bone width'] =  smallval / chord # 0.025 / chord # 0.4 / 0.0254 / chord
    
    # add 2 *  inside first parenth to get double width columns
      ###   # how wide should the gaps between bones be in % chord     
    p['fishbone gap width'] = 0.05 # ( ( (0.18 + 0.05) / 6 ) ) - p['fishbone bone width']
    
      ###   # how wide should the spine be in % thickness     
    p['fishbone spine width'] = smallval / chord / 0.2 # 0.03 / chord / 0.15 # 0.025 in with 6 in c and .15 th # .05 / 6 /0.15 # 0.05 in with 6 in c and .15 th # 
    
      ###   # what percentage of the bone top width will be the base boneWidth
    p['fishbone base'] = 1.0 ## was triboneperc, now fishbone base
    
      ###   # should the Thickness be based on the Max thickness, or the local thickness at each point    
    p['fishbone is t local'] = False
    
    # thinnest got was - 0.0025
    p['fishbone ledge width'] = 0.9 # ( p['fishbone gap width'] * chord - 0.02 ) / ( p['fishbone gap width'] * chord ) # needs to be less than .19 / .23 / 6
    
    p['fishbone ledge depth'] = smallval / chord / 0.3 # 0.03 / chord / 0.12
    
    p['fishbone ledge is flat'] = False
    
    p['fishbone percent ledge right'] = 0.5#1.0 # 0.5 # 1.0 # -1.0 # 
    
    p['fishbone skip points'] = True # if you want to skip interpolation points for ease of use
      
    
    
    p['is bezier'] = False
    
    p['default bezier'] = False
    
    p['bezier start'] = 0.4
    
    p['bezier end'] = 0.95
    
    p['bezier period'] = 0.091 #### there might be a discontinuity at the peaks when period is greater (ie 0.6 etc)
    
    p['bezier P1 percentage'] = 0.75       #1 # straight up and down # 2/3 # get closer to osci # 0.5 # original #
    
    p['bezier P2 percentage'] = 0.25       #0 # straight up and down # 1/3 # get closer to osci # 0.5 # original #
    
    thinnest = 0.1 #  0.03 / .12 / 4
    
    p['bezier leading edge spine width'] = thinnest
    # 0.13
    p['bezier trailing edge spine width'] = thinnest # p['leading edge spine width'] * 0.5
    
    p['bezier first up'] = True
    
    p['bezier is t local'] = False
    
    
    
    
    p['add control horn'] = False
    
    p['control horn depth'] = 1.2 / 2.54 # default # 1.2 / 2.54
    
    p['control horn height'] = 0.5 / 2.54 # default # 0.5 / 2.54 
    
    p['control horn diameter'] = 0.2 / 2.54 # default # 0.2 / 2.54
    
    p['control horn type'] = 'semicircle' # only other option - 'sharp'
    
    p['control horn add overhang'] = False
    
    
    
    p['make dovetail'] = False
    
    # p['dovetail start point'] = pt.point( 5.31863128 ,11.09910554 ,4 )
    
    # p['dovetail end point'] = pt.point( 11.00283128 ,11.09910554 ,4 )
    
    # p['dovetail head width'] = 0.5
    
    # p['dovetail neck width'] = 0.25
    
    # p['dovetail tail width'] = 0.125
    
    
    
    # selig datfile location if data points are to be used instead. for tip
    p['tip wingshape directory'] = 'None'  # '/home/ben/Desktop/input.mu' # 
    
    # selig datfile location if data points are to be used instead. for root
    p['root wingshape directory'] = 'None' # '/home/ben/Desktop/input.mu' #  
    
    # '/home/ben/Desktop/input.mu'
    # '/home/ben/Desktop/foilshapes/E335_200pts.dat'
    # '/home/ben/Desktop/foilshapes/E335_150pts.txt'
    # '/home/ben/Desktop/foilshapes/E335.txt'
    # '/home/ben/Desktop/foilshapes/E186.txt'
    # '/home/ben/Desktop/foilshapes/E591.txt'   
    # '/home/ben/Desktop/foilshapes/NACA 64A210.txt'
    # '/home/ben/Desktop/foilshapes/selig NACA 8406.dat'
    # '/home/ben/Desktop/foilshapes/grumman k2.txt'
    # '/home/ben/Desktop/foilshapes/kenmar.txt'
    # '/home/ben/Desktop/foilshapes/oaf095.txt'
    # '/home/ben/Desktop/foilshapes/ys-900.txt'
    # '/home/ben/Desktop/foilshapes/wing_one.mu'
    # '/home/ben/Desktop/foilshapes/test.mu'
    # '/home/ben/Desktop/foilshapes/test_two.mu'
    
    p['write lednicer'] = False
    
    p['write selig'] = False
    
    p['write Gcode'] = False
    
    p['write stl'] = False
    
    p['write z'] = True
    
    p['write skip first line'] = True
    
    # give file name. file must be .ngc to be a gcode file
    
    p['Gcode file name'] = 'gcode.ngc'
    
    
    
    return p



"""

figure out how to read a multiple wing machup file and create multiple gcode files...

make an option to have a bone ledge which direction it goes ie true if towards LE, false if towards TE

read in selig dat files from internet

make it so ie. only a root could be txt and the tip could be generated like a NACA 0012

why does bezier take so long????

make inner kerfer that creates a hole in the LE and or the TE of Fishbone and Bezier...

make it possible to make an stl multi wing like manta print...





"""



def gFoil(): 
    
    
    p = getParams()
    import json
    with open('into.json','w') as fp:
        print(json.dump(p,fp,indent=4))
    
    

    
    file_type = ''
    
    rt_d = p['root wingshape directory']; tp_d = p['tip wingshape directory']
    
    if((rt_d[-4:] == '.dat' and tp_d[-4:] == '.dat' )  or ( rt_d[-4:] == '.txt' and tp_d[-4:] == '.txt' ) ): # != 'none'
        
        file_type = '.txt file of data points'
        
        xyzfoil = fh.FileReader( p['root wingshape directory'] )
        
        uvwfoil = fh.FileReader( p['tip wingshape directory'] )
        
        xyzfoil = ag.change_chord( xyzfoil, p['root chord'] )
        
        uvwfoil = ag.change_chord( uvwfoil, p['tip chord'] )
        
        
        if(p['is fishbone'] and not p['is bezier']):
            
            if(p['default fishbone']):
                
                fishxyz = ag.BoneMaker(xyzfoil)
                
                fishuvw = ag.BoneMaker(uvwfoil)
                
            else:
                
                fishxyz = ag.BoneMaker(xyzfoil, p['fishbone start'], 
                                     p['fishbone end'], p['fishbone bone width'], p['fishbone gap width'], 
                                     p['fishbone spine width'], p['fishbone base'], p['fishbone is t local'],
                                     p['fishbone ledge width'],p['fishbone ledge depth'],p['fishbone ledge is flat'],
                                     p['fishbone skip points'],p['fishbone percent ledge right'])
                
                
                fishuvw = ag.BoneMaker(uvwfoil, p['fishbone start'], 
                                     p['fishbone end'], p['fishbone bone width'], p['fishbone gap width'], 
                                     p['fishbone spine width'], p['fishbone base'], p['fishbone is t local'],
                                     p['fishbone ledge width'],p['fishbone ledge depth'],p['fishbone ledge is flat'],
                                     p['fishbone skip points'],p['fishbone percent ledge right'])
            
            # assign the wing shape coordinates back to the used values
            xyzfoil = fishxyz; uvwfoil = fishuvw
                
        elif(not p['is fishbone'] and p['is bezier']):
            
            if(p['default bezier']):
                
                bezxyz = ag.Bezier_Spine(xyzfoil)
                
                bezuvw = ag.Bezier_Spine(uvwfoil)
                
            else:
                
                bezxyz = ag.Bezier_Spine(xyzfoil, p['bezier start'],
                                     p['bezier end'], p['bezier P1 percentage'], p['bezier P2 percentage'],
                                     p['bezier period'], p['bezier leading edge spine width'], 
                                     p['bezier trailing edge spine width'], 
                                     p['bezier first up'], p['bezier is t local'])
                
                bezuvw = ag.Bezier_Spine(uvwfoil, p['bezier start'],
                                     p['bezier end'], p['bezier P1 percentage'], p['bezier P2 percentage'],
                                     p['bezier period'], p['bezier leading edge spine width'], 
                                     p['bezier trailing edge spine width'], 
                                     p['bezier first up'], p['bezier is t local'])
            
            # assign the wing shape coordinates back to the used values
            xyzfoil = bezxyz; uvwfoil = bezuvw
    
    elif(p['root wingshape directory'][-3:] == '.mu' and p['root wingshape directory'][-3:] == '.mu'): 
        
        file_type = '.mu file'
        
        mu = fh.MachUpFileReader(p['root wingshape directory'])
        
        p['is left wing'] = mu['is left']; p['span'] = mu['span']
        
        p['sweep'] = mu['sweep']; p['dihedral'] = mu['dihedral']
        
        p['mounting angle'] = mu['mounting angle']; p['washout'] = mu['washout']
        
        p['root chord'] = mu['root chord']; p['tip chord'] = mu['tip chord']
        
        p['root airfoil'] = mu['root airfoil']; p['tip airfoil'] = mu['tip airfoil']

        
        # we won't use dihedral or mounting angle. Those parameters can be physically modified by the user after the wing is cut
        
        xyzfoil,fishxyz,bezxyz = ag.Airfoil_Creator(p['root airfoil'], p['root chord'],p ,p['NACA number of points'])
        
        uvwfoil,fishuvw,bezuvw = ag.Airfoil_Creator(p['tip airfoil'], p['tip chord'], p ,p['NACA number of points'])
        
        
    else: # reading in input values  # if(fileLocationRoot == 'none' and fileLocationTip == 'none'):
        
        file_type = '%s root airfoil and %s tip airfoil' % (p['root airfoil'], p['tip airfoil'])
        
        xyzfoil,fishxyz,bezxyz = ag.Airfoil_Creator(p['root airfoil'], p['root chord'], p ,p['NACA number of points'])
        
        uvwfoil,fishuvw,bezuvw = ag.Airfoil_Creator(p['tip airfoil'], p['tip chord'], p ,p['NACA number of points'])

    
    #aligns outer quarter chord with root quarter chord
    if(p['is fishbone'] and not p['is bezier']):
        xyzfoil = fishxyz; uvwfoil = fishuvw
    
    elif(not p['is fishbone'] and p['is bezier']):
        xyzfoil = bezxyz; uvwfoil = bezuvw
    
    # fix non wholeness in airfoil shapes
    
    xyzfoil = geo.makeWhole(xyzfoil)
    uvwfoil = geo.makeWhole(uvwfoil)
    
    
#    fh.NACACSVWriter('OSCI',x,y) # creates a file in the aerolab folder of the listed shape

    
    uvwfoil = geo.alignQuarterChord(xyzfoil,uvwfoil)
    
    if(p['washout']):
        
        uvwfoil = geo.giveWashout(uvwfoil,p['washout'])
    
    if(p['sweep']):
        
        uvwfoil = geo.giveSweep(uvwfoil,p['sweep'],p['span'])
        
    if(p['is left wing']): # this might need to be swapped becuase of the way the hotwireCutter positions 4 axes
        
        temp = xyzfoil
        
        xyzfoil = uvwfoil; uvwfoil = temp;
    
    if(p['root wire kerf'] and p['tip wire kerf']):
        
        kerfxyz = ag.kerfMaker(xyzfoil,p['root wire kerf'])
        
        kerfuvw = ag.kerfMaker(uvwfoil,p['tip wire kerf'])
        
    if(not p['root wire kerf'] and not p['tip wire kerf']):
        
        kerfxyz = xyzfoil
        
        kerfuvw = uvwfoil
    
    
    # if only a section of the foil is to be cut, the values are returned to the user.
    if(p['is section']):
        
        sectionxyz = geo.sectionMaker( kerfxyz, p['section start'], p['section end'])
        
        sectionuvw = geo.sectionMaker( kerfuvw, p['section start'], p['section end'])
    
    # otherwise pass along the values
    else:
        sectionxyz = kerfxyz; sectionuvw = kerfuvw;
    
    
    if(p['add control horn']):
        xyzfoil,hole = ch.control_horn(xyzfoil, p['control horn depth'], 
                                        p['control horn height'], 
                                        p['control horn diameter'],
                                        p['control horn type'],
                                        p['control horn add overhang'])
    
    
    
    w_xyz,w_uvw = geo.fixAxisPoints(sectionxyz,sectionuvw,p['span'],p['offset'],p['wire length'])
    
    # determine centroids
    cn_xyz = geo.Centroid(xyzfoil)
    cn_uvw = geo.Centroid(uvwfoil)
    
    
    # plot the airfoils
    other_plots = [[],[]]
    
    other_plots[0].append(cn_xyz)
    other_plots[1].append(cn_uvw)
    
    displays = ['k','','k','k','k'] # ["","ro","","",""]

    if(p['is section']):
        other_plots[0].append(sectionxyz)
        other_plots[1].append(sectionuvw)
        
    if(p['root wire kerf'] and p['tip wire kerf']):
        other_plots[0].append(kerfxyz)
        other_plots[1].append(kerfuvw)
    
    if(p['add control horn']):
        other_plots[0].append(hole)

    # ax = xyzfoil.x(); ay = xyzfoil.y()

    # # filename = "/home/ben/Desktop/airfoil.txt"
    # # filename = "/home/ben/Desktop/bezier.txt"
    # # filename = "/home/ben/Desktop/fishbone.txt"
    # filename = "/home/ben/Desktop/fishscale.txt"

    # with open(filename,"w") as f:

    #   for i in range(len(ax)):
    #     f.write("{:<16.14f} {:<16.14f} \n".format(ax[i],ay[i]))
      
    #   f.close()

    
    xyzfoil.plot(other_plots[0],display = displays, Subplot = 1)
    uvwfoil.plot(other_plots[1],display = displays, Subplot = 2)
    
    
    
    # plot a 3D image of the outer shapes traced by the wire cutter,
    #  and the inner shapes on the edges of the wing.
    
    fig = plt.figure()
    
    ax = fig.add_subplot(111, projection='3d')
    
    face = 'y'
    
    ax.plot(w_xyz.x(),w_xyz.y(),0,zdir= face)
    
    ax.plot(xyzfoil.x(),xyzfoil.y(),p['offset'],zdir= face)
    
    ax.plot(uvwfoil.x(),uvwfoil.y(),p['offset']+p['span'],zdir= face)
    
    ax.plot(w_uvw.x(),w_uvw.y(),p['wire length'],zdir= face)
    
    ax.view_init(20,-30)
    
    for i in range(0,len(w_xyz.points),4):
        
        ax.plot([w_xyz.points[i].x,w_uvw.points[i].x],[w_xyz.points[i].y,
                 w_uvw.points[i].y],[0,p['wire length']],zdir= face)
    
    w_xyz,w_uvw = geo.shiftUp(w_xyz,w_uvw,p['depth'])
    w_xyz,w_uvw = geo.shiftRight(w_xyz,w_uvw)
    
    
    w_xyz = geo.roundValuesG(w_xyz)
    w_uvw = geo.roundValuesG(w_uvw)
    
    ######
    # add Control horn if desired
    
    
    
    ##########
    
    
    
    new_chord = 1
    
    special = ''
    c1xyzfoil = ag.change_chord(xyzfoil, new_chord)
    
#    if(p['add control horn']):
#        
#        old_chord = 12
#        
#        hx_c1, hy_c1 = ag.change_chord(hx,hy, new_chord, old_chord)
    
    if(p['is fishbone'] and not p['is bezier']):
        
        special = 'fishbone '
        
    if(not p['is fishbone'] and p['is bezier']):
        
        special = 'bezier '
    
    if(p['add control horn']):
        
        special += 'ctrhorn '

    # get wingshape name
    airfoil = p['root airfoil']
    if p['root wingshape directory'] != 'None' :
      airfoil = p['root wingshape directory'].split('/')[-1].split('.')[0]
  
    # datetime object containing current date and time
    now = date.now()
    # dd/mm/YY H:M:S
    dt_string = now.strftime("%m-%d-%Y %H:%M:%S")

    if p['write selig']:
        
        fh.FileWriter('selig ' + special + airfoil + ' ' + dt_string ,c1xyzfoil,True,p['write z'],p['write skip first line']) 
        
        
    if(p['write lednicer']):
        
        fh.FileWriter('lednicer ' + special + airfoil + ' ' + dt_string ,c1xyzfoil,False,p['write z']) # + ' %.f' % p['root chord'] + 'c'    
    
    
    
    #insert into the first values repositioning so that the machine can cut the desired shape
    
    
    fishTEXT = ''
    if(p['is fishbone'] and not p['is bezier']):
        fishTEXT = 'Fishbone VCCW,'
    
    beziTEXT = ''
    if(not p['is fishbone'] and p['is bezier']):
        beziTEXT = 'Bezier Spine VCCW,'
    
    leftTEXT = 'right wing,'
    if(p['is left wing']):
        leftTEXT = 'left wing,'
    
    unit = 'in'
    if(p['is in inches']):
        unit = 'mm'
    
    file_explanation = """(This file is a G-Code generated wing shape for use in the Aerolab CNCHotwire)
(Wing shape is based off a %s)
(Wing Characteristics::)\n(%s %s %s root chord %.1f %s, tip chord %.1f %s,)
(%.1f degrees washout, %.1f degrees sweep,span %.1f %s)\n""" % (file_type,fishTEXT, beziTEXT,
    leftTEXT,p['root chord'],unit,p['tip chord'],unit,p['washout'],p['sweep'],p['span'],unit)  
    
    if(p['write Gcode']):
        
        # write the gcode file
        
        fh.GcodeFileWriter(w_xyz,w_uvw,p['is in inches'],file_explanation,p['depth'],p['length'],p['Gcode file name'])
    
    
#    if(p['make dovetail']):
#        
#        dove = dt.doveTail(p['dovetail start point'],p['dovetail end point'],p['dovetail head width'],
#                    p['dovetail neck width'], p['dovetail tail width'])
#        
#        dove.plot()
#        
#        fh.FileWriter('dovetail',dove,True,p['write z'],skipFirst = False)
        
        
    
    
    
    return 0


gFoil()



