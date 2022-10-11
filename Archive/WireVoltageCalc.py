#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  6 09:47:04 2018

@author: benjamin
"""

#Nichrome 80
wire_diam = 0.0799 # diam in mm
d = (wire_diam/1000.) * 39.3701
print("wire diameter: " + str(round(d,5)))

RperF = 67.63 # Resistance per foot of the wire
wireL = 31.5  # wire length
R = wireL * RperF # Resistance of the overall wire

V = 42. # voltage 
#V = I*R
I = V/R
print(str(round(I,5)) + " Amperes")



# thats a roll