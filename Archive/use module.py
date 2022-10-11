#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 10 16:28:31 2018

@author: benjamin
"""

from aerolab import gcode as gc
from aerolab import geometry as geo
from aerolab import NACApoints as nc

same = 'NACA 6410'

x,y,fx,fy = nc.infoCompile(same,1)
lets,eat,some,food = nc.infoCompile(same,0.85)
newx = geo.giveSweep(lets,44,1)





xaxisx = [0,1]
xaxisy = [0,0]
import matplotlib.pyplot as plt

plt.axis([-.01,max(p)*1.3,min(q)*1.2,max(q)])
plt.axis('equal')
plt.plot(x,y)
plt.plot(newx,eat)

plt.show()







#Set the program so that if a text file is input, then it calculates the chord, and if not, it uses the givine chord.









