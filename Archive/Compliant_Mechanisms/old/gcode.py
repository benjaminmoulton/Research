#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 10 16:30:30 2018

@author: benjamin
"""

#this function takes two lists of data, and writes a Gcode file for CNC hotwiring
# HOWEVER, this is not translated to LinuxCNC, and is for non-Tapered wings
def ParallelFileWriter(xdata,ydata):  # writes a given data file into Gcode for CNC hotwiring
    f = open("gcode.txt","w+")

    for i in range(len(xdata)):
        f.write("G01 X" + str(xdata[i]) + " Y" + str(ydata[i]) 
                + " U" + str(xdata[i]) + " V" + str(ydata[i])+ " F20\n")
    f.write("G2")
    f.close()

    return 0


#this function takes  a list of x and y data, and creates a simple cut out for 
# a two axis GCODE machine
def twoDimFileWriter(xdata,ydata):
    f = open("gcode.txt","w+")

    for i in range(len(xdata)):
        f.write("G01 X" + str(xdata[i]) + " Y" + str(ydata[i]) + " F20\n")
    f.write("G2")
    f.close()
    return 0


#this function takes a file location and reads the file to create a list of 
# the x and y values found for an airfoil shape
def NACAFileReader(file): # reads in a NACA dat file as on airfoiltools.com and converts the file to two lists of x and y values
    f = open(file,"r")
    lst = []
    f.readline()
    for line in f:
        lst.append([ float(x) for x in line.split()])
    x_values = [ x[0] for x in lst]
    y_values = [ x[1] for x in lst]
    return x_values,y_values



#this function MAKES a .txt file of a given set of points
def NACAFileWriter(fileName,x,y):
    fullFileName = fileName + '.txt'
    namestart = fileName + '\n'
    f = open(fullFileName,"w+")
    f.write(namestart)
    for xval,yval in zip(x,y):
        space = ' '
        if(yval < 0):
            space = ''
        f.write(' %.5f %s%.5f\n' % (xval,space,yval))
    loc = '/home/benjamin/Desktop/PythonCode/aerolab/' + fullFileName
    dest = '/home/benjamin/Desktop/foilshapes/' + fullFileName
    copyRenaMove(loc,dest)
    return 0
    
    
#this function takes four lists and creates a tapered Gcode File for CNC hotwiring
def taperedFileWriter(xdata,ydata,udata,vdata,isInches,shapeExplained,depth,width):  # writes a given data file into Gcode for CNC hotwiring
    f = open("gcode.txt","w+")
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
    pt3y = str(ydata[0])
    pt4x = str(round(xdata[0] + 0.12*unit,3))
    pt4u = str(round(udata[0] + 0.12*unit,3))
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
    for i in range(len(xdata)):
        if(i == 0):
            Feedrate = " F5"
        elif(i > 1):
            Feedrate = ""
        f.write("G1 X" + str(xdata[i]) + " Y" + str(ydata[i]) 
                + " U" + str(udata[i]) + " V" + str(vdata[i])+ Feedrate + "\n")
    
    f.write("G1 X" + pt4x + " Y" + pt3y + " U" + pt4u + " V" + pt3y + " F5\n")
    f.write("G1 X" + pt3x + " Y" + pt3y + " U" + pt3x + " V" + pt3y + " F1\n")
    f.write("G1 X" + pt2x + " Y" + pt2y + " U" + pt2x + " V" + pt2y + " F5\n")
    f.write("G1 X" + pt1x + " Y" + pt1y + " U" + pt1x + " V" + pt1y + "\n")
    f.write("G1 X0 Y0 U0 V0\n")
    
    
    f.write("M2")
    f.close()
    copyRenaMove("/home/benjamin/Desktop/PythonCode/aerolab/gcode.txt",'/home/benjamin/Desktop/gcode.ngc')
    return 0



#this function takes a MachUp file, and reads in information to determine 
# what foilshape is to be specified
def MachUpFileReader(file):
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
    
    if side == 'right':
        isLeft = False
    elif side == 'left':
        isLeft = True
    else: # side == 'both':
        printBoth = True

    return isLeft,span,sweep,dihedral,mountingAngle,washout,rchord,tchord,rfoil,tfoil


#this function moves a file from one location to another  
def move(location):
    import os
    #import shutil
    
    os.rename(location, "/home/benjamin/Desktop/Hotwire/NACA 0009.txt")
    #shutil.move("path/to/current/file.foo", "path/to/new/destination/for/file.foo")
    return 0

def copyRenaMove(location,destination):
    from os import rename
    from shutil import copy2
    copy2(location,destination)
#    rename(destination,'/home/benjamin/Desktop/gcode.ngc')
    
#    os.rename(location, "/home/benjamin/Desktop/Hotwire/NACA 0009.txt")
    #shutil.move("path/to/current/file.foo", "path/to/new/destination/for/file.foo")
    return 0
















