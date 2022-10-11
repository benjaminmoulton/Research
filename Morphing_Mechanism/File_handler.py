#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 10 16:30:30 2018

@author: benjamin
"""
import point as pt
import json
import os


#this function takes two lists of data, and writes a Gcode file for CNC hotwiring
# HOWEVER, this is not translated to LinuxCNC, and is for non-Tapered wings
def ParallelFileWriter(airfoil):  # writes a given data file into Gcode for CNC hotwiring
    
    # file is opened for writing
    f = open("gcode.txt","w+")

    # a for loop loops through given data and creates a file with no taper (root airfoil = tip airfoil, same chords, no sweep, etc.)
    for i in range(len(airfoil.points)):
        
        # points are written to file
        f.write("G01 X" + str(airfoil.points[i].x) + " Y" + str(airfoil.points[i].y) 
                + " U" + str(airfoil.points[i].x) + " V" + str(airfoil.points[i].y)+ " F20\n")
    
    # Gcode requires file end with G2 to terminate the procedure
    f.write("G2")
    
    # file is closed
    f.close()

    return 0


#this function takes  a list of x and y data, and creates a simple cut out for 
# a two axis GCODE machine
def twoDimFileWriter(airfoil):
    
    # file is opened for writing
    f = open("gcode.txt","w+")

    for i in range(len(airfoil.points)):
        f.write("G01 X" + str(airfoil.points[i].x) + " Y" + str(airfoil.points[i].y) + " F20\n")
    f.write("G2")
    f.close()
    return 0


def DetermineFormat(file):
    
    is_selig = True
    
    f = open(file,'r')
    
    f.readline()
    lst = []
    for line in f:
        lst.append([ float(point) for point in line.split()])
    
    f.close()
    
    first_value = lst[0][0]
    
    if(first_value > 1):
        is_selig = False
    
    return is_selig



#this function takes a file location and reads the file to create a list of 
# the x and y values found for an airfoil shape
def SeligFileReader(file): # reads in a NACA dat file as on airfoiltools.com and converts the file to two lists of x and y values
    
    f = open(file,"r")
    
    lst = []
    
    f.readline()
    
    for line in f:
        
        lst.append([ float(point) for point in line.split()])
        
    x_values = [ point[0] for point in lst]
    
    y_values = [ point[1] for point in lst]
    
    airfoil = pt.points(x_values,y_values)
    
    return airfoil



def LednicerFileReader(file):
    
    f = open(file,"r")
    
    f.readline()
    
    lst = []
    
    for line in f:
        lst.append([ float(point) for point in line.split()])
    
    top_num = int( lst[0][0] )
    
    bottom_num = int( lst[0][1] )
    
    lst.pop(0)
    lst.pop(0)
    
    top_x = []; top_y = [];
    
    for i in range(top_num):
        
        top_x.append(lst[0][0])
        
        top_y.append(lst[0][1])
        
        lst.pop(0)
    
    lst.pop(0)
    
    bottom_x = []; bottom_y = [];
    
    for i in range(bottom_num):
        
        bottom_x.append(lst[0][0])
        
        bottom_y.append(lst[0][1])
        
        lst.pop(0)
    
    # reorder top lists
    
    top_x = top_x[::-1]; top_y = top_y[::-1];
    
    # remove duplicate value used by lednicer file
    
    bottom_x.pop(0); bottom_y.pop(0);
    
    x = top_x + bottom_x; y = top_y + bottom_y;
    
    airfoil = pt.points(x,y)
    
    return airfoil



def FileReader(file):
    
    is_selig = DetermineFormat(file)
    
    if(is_selig):
        
        airfoil = SeligFileReader(file)
    
    else:
        
        airfoil = LednicerFileReader(file)
    
    return airfoil


#this function MAKES a .txt file of a given set of points
def SeligFileWriter(fileName,airfoil,write_z = False,skipFirst = False):
    
    fullFileName = fileName + '.dat'
    
    file_description = fileName + ' AIRFOIL' + '\n'
    
    f = open(fullFileName,"w+")
    
    if(not skipFirst):
        f.write(file_description)
    
    if(not write_z):
        for xval,yval in zip(airfoil.x(),airfoil.y()):
            
            space = ' '
            
            if(yval < 0):
                
                space = ''
                
            f.write(' %.5f %s%.5f\n' % (xval,space,yval))
            
    else:
        for xval,yval,zval in zip(airfoil.x(),airfoil.y(),airfoil.z()):
            
            space = ' '
            
            if(yval < 0):
                
                space = ''
                
            f.write(' %.5f %s%.5f %.5f\n' % (xval,space,yval,zval))
    
    f.close()
    
    from path_maker import move_file_to_Desktop as mD
    
    mD(fullFileName)
    
    return



#this function MAKES a .txt file of a given set of points
def LednicerFileWriter(fileName,airfoil,write_z = False):
    
    fullFileName = fileName + '.dat'
    
    file_description = fileName + ' AIRFOIL' + '\n'
    
    f = open(fullFileName,"w+")
    
    
    f.write(file_description)
    
    min_index = airfoil.x_min_index
    
    top_num = min_index + 1
    
    bottom_num = len(airfoil.points) - min_index
    
    number_line = '%9d. %9d.\n' % (top_num, bottom_num)
    
    f.write(number_line)
    
    f.write('\n')
    
    top_x = airfoil.split(0,top_num).x(); top_y = airfoil.split(0,top_num).y();
    
    top_x = top_x[::-1]; top_y = top_y[::-1];
    
    bottom_x = airfoil.split(top_num,0).x(); bottom_y = airfoil.split(top_num,0).y();
    
    if(not write_z):
        
        for xval,yval in zip(top_x,top_y):
            
            y_coord = '%.7f' % yval
            
            if(yval < 0):
                
                y_coord = '-' + y_coord[2:]
                
            f.write(' %.7f %s\n' % (xval,y_coord))
        
        
        f.write('\n')
        
        
        for xval,yval in zip(bottom_x,bottom_y):
            
            y_coord = '%.7f' % yval
            
            if(yval < 0):
                
                y_coord = '-' + y_coord[2:]
                
            f.write(' %.7f %s\n' % (xval,y_coord))
            
            
    else:
        top_z = airfoil.split(0,top_num).z(); top_z = top_z[::-1];
        bottom_z = airfoil.split(top_num,0).z();
        
        
        for xval,yval,zval in zip(top_x,top_y,top_z):
            
            y_coord = '%.7f' % yval
            
            if(yval < 0):
                
                y_coord = '-' + y_coord[2:]
                
            f.write(' %.7f %s %.7f\n' % (xval,y_coord,zval))
        
        
        f.write('\n')
        
        
        for xval,yval,zval in zip(bottom_x,bottom_y,bottom_z):
            
            y_coord = '%.7f' % yval
            
            if(yval < 0):
                
                y_coord = '-' + y_coord[2:]
                
            f.write(' %.7f %s %.7f\n' % (xval,y_coord,zval))
            
    
    f.close()
    
    from path_maker import move_file_to_Desktop as mD
    
    mD(fullFileName)
    
    return


def FileWriter(fileName,airfoil,isSelig = True, write_z = False,skipFirst = False):
    
    if(isSelig):
        SeligFileWriter(fileName,airfoil,write_z,skipFirst)
    
    else:
        LednicerFileWriter(fileName,airfoil,write_z)
    
    return

#this function MAKES a .txt file of a given set of points
def NACACSVWriter(fileName,airfoil):
    fullFileName = fileName + '.csv'
    file_description = fileName + '\n'
    f = open(fullFileName,"w+")
    f.write(file_description)
    for xval,yval in zip(airfoil.x(),airfoil.y()):
        space = ' '
        if(yval < 0):
            space = ''
        f.write(' %.5f, %s%.5f,\n' % (xval,space,yval))
#    loc = '/home/benjamin/Desktop/PythonCode/aerolab/' + fullFileName
#    dest = '/home/benjamin/Desktop/foilshapes/' + fullFileName
#    copyRenaMove(loc,dest)
    return 0
    
    
#this function takes four lists and creates a tapered Gcode File for CNC hotwiring
def GcodeFileWriter(xy_airfoil,uv_airfoil,isInches,shapeExplained,depth,width,name):  # writes a given data file into Gcode for CNC hotwiring
    f = open(name,"w+")
    f.write(shapeExplained)
    f.write("G17\n")
    unit = 1
    if(isInches):
        f.write("G20\n")
    else:
        f.write("G21\n")
        unit = 25.4
    #enter code here for posthome, precut movements
    pt1x = str(0.5)
    pt1y = str(depth + 1.0)
    pt2x = str(1 + width + 1)
    pt2y = pt1y
    pt3x = pt2x
    pt3y = str(xy_airfoil.points[0].y)
    pt4x = str(round(xy_airfoil.points[0].x + 0.12*unit,3))
    pt4u = str(round(uv_airfoil.points[0].x + 0.12*unit,3))
    pt4y = pt3y
    
    f.write("G1 X" + pt1x + " Y" + pt1y + " U" + pt1x + " V" + pt1y + " F5\n")
    f.write("G4 P0.3\n")
    f.write("G1 X" + pt2x + " Y" + pt2y + " U" + pt2x + " V" + pt2y + "\n")
    f.write("G4 P0.3\n")
    f.write("G1 X" + pt3x + " Y" + pt3y + " U" + pt3x + " V" + pt3y + " F4\n")
    f.write("G4 P0.3\n")
    f.write("G1 X" + pt4x + " Y" + pt3y + " U" + pt4u + " V" + pt3y + " F1\n")
    f.write("G4 P0.3\n")
    
    Feedrate = " F5"
    for i in range(len(xy_airfoil.points)):
        if(i == 0):
            Feedrate = " F5"
        elif(i > 1):
            Feedrate = ""
        f.write("G1 X" + str(xy_airfoil.points[i].x) + " Y" + str(xy_airfoil.points[i].y) 
                + " U" + str(uv_airfoil.points[i].x) + " V" + str(uv_airfoil.points[i].y)+ Feedrate + "\n")
    
    
    # the point is added at the beginning of the airfoil in case the 
    #    piece to be cut out is an airfoil section
    f.write("G1 X" + str(xy_airfoil.points[0].x) + " Y" + str(xy_airfoil.points[0].y) 
                + " U" + str(uv_airfoil.points[0].x) + " V" + str(uv_airfoil.points[0].y)+ "\n")
    
    
    f.write("G1 X" + pt4x + " Y" + pt4y + " U" + pt4u + " V" + pt4y + " F5\n")
    f.write("G1 X" + pt3x + " Y" + pt3y + " U" + pt3x + " V" + pt3y + " F1\n")
    f.write("G1 X" + pt2x + " Y" + pt2y + " U" + pt2x + " V" + pt2y + " F5\n")
    f.write("G1 X" + pt1x + " Y" + pt1y + " U" + pt1x + " V" + pt1y + "\n")
    f.write("G1 X0 Y0 U0 V0\n")
    
    
    f.write("M2")
    f.close()
    
    import path_maker as pm
    
    pm.move_file_to_Desktop(name)
    
    return 0



#this function takes a MachUp file, and reads in information to determine 
# what foilshape is to be specified
def MachUpFileReader(file):

    # open json file
    json_string = open(file).read()
    json_vals = json.loads(json_string)

    side =          json_vals["geometries"][2]["side"]
    span =          json_vals["geometries"][2]["span"]
    sweep =         json_vals["geometries"][2]["sweep"]
    dihedral =      json_vals["geometries"][2]["dihedral"]
    mountingAngle = json_vals["geometries"][2]["mount"]
    washout =       json_vals["geometries"][2]["washout"]
    rchord =        json_vals["geometries"][2]["root_chord"]
    tchord =        json_vals["geometries"][2]["tip_chord"]
    rfoil = "NACA 2412" #json_vals["geometries"]["side"]
    tfoil = "NACA 2412" #json_vals["geometries"]["side"]
#    params = [side,span,sweep,dihedral,mountingAngle,washout,rchord,tchord,rfoil,tfoil]   
    isLeft = False
    
    if side == 'right':
        isLeft = False
    elif side == 'left':
        isLeft = True
    else: # side == 'both':
        printBoth = True
    
    mu = {}; mu['is left'] = isLeft; mu['span'] = span; mu['sweep'] = sweep
    
    mu['dihedral'] = dihedral; mu['mounting angle'] = mountingAngle
    
    mu['washout'] = washout; mu['root chord'] = rchord; mu['tip chord'] = tchord
    
    mu['root airfoil'] = rfoil; mu['tip airfoil'] = tfoil
    
    return mu



#this function takes a MachUp file, and reads in information to determine 
# what foilshape is to be specified
def MultiMachUpFileReader(file):
#    import os
#    for filename in os.listdir("."):
#        if filename.endswith("mu"):
#            os.rename(filename, filename[:-2] + 'txt')
    f = open(file,"r")
    for i in range(36):
        f.readline()
    side = f.readline().strip().replace(" ","")[8:-2]

    span = float(f.readline().strip().replace(" ","")[7:-1])
    sweep = float(f.readline().strip().replace(" ","")[8:-1])
    dihedral = float(f.readline().strip().replace(" ","")[11:-1])
    mountingAngle = float(f.readline().strip().replace(" ","")[8:-1])
    washout = float(f.readline().strip().replace(" ","")[10:-1])
    rchord = float(f.readline().strip().replace(" ","")[13:-1])
    tchord = float(f.readline().strip().replace(" ","")[12:-1])
    f.readline()
    rfoil = f.readline().strip()[1:-4]
    for j in range(14):
        f.readline()
    tfoil = f.readline().strip()[1:-4]
#    params = [side,span,sweep,dihedral,mountingAngle,washout,rchord,tchord,rfoil,tfoil]   
    isLeft = False
    
    if side == 'right':
        isLeft = False
    elif side == 'left':
        isLeft = True
    else: # side == 'both':
        printBoth = True
    
    mu = {}; mu['is left'] = isLeft; mu['span'] = span; mu['sweep'] = sweep
    
    mu['dihedral'] = dihedral; mu['mounting angle'] = mountingAngle
    
    mu['washout'] = washout; mu['root chord'] = rchord; mu['tip chord'] = tchord
    
    mu['root airfoil'] = rfoil; mu['tip airfoil'] = tfoil
    
    return mu



##this function moves a le from one location to another  
#def move(location):
#    import os
#    #import shutil
#    
#    os.rename(location, "/home/benjamin/Desktop/Hotwire/NACA 0009.txt")
#    #shutil.move("path/to/current/file.foo", "path/to/new/destination/for/file.foo")
#    return 0
#
#def copyRenaMove(location,destination):
#    from os import rename
#    from shutil import copy2
#    copy2(location,destination)
##    rename(destination,'/home/benjamin/Desktop/gcode.ngc')
#    
##    os.rename(location, "/home/benjamin/Desktop/Hotwire/NACA 0009.txt")
#    #shutil.move("path/to/current/file.foo", "path/to/new/destination/for/file.foo")
#    return 0


