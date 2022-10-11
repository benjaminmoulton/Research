#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep  4 11:11:44 2018

@author: benjamin
"""

import NACAworks as ncw
from matplotlib import pyplot as plt

#def NACAFileReader(file): # reads in a NACA dat file as on airfoiltools.com and converts the file to two lists of x and y values
#    f = open(file,"r")
#    lst = []
#    f.readline()
#    for line in f:
#        lst.append([ float(x) for x in line.split()])
#    x_values = [ x[0] for x in lst]
#    y_values = [ x[1] for x in lst]
#    return x_values,y_values

#filetoimport = print("input file location: ")
#x,y = NACAFileReader(filetoimport)
#plt.plot(x,y,'bo')
#plt.show()
