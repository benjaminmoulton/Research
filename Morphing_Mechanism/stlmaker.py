#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul  2 15:18:33 2019

@author: benjamin
"""

#x = [3, 5,12, 9, 5]
#y = [4,11, 8, 5, 6]
#
import geometry as geo
#
#print(geo.Area(x,y))
#
import File_handler as fh

airfoil = fh.FileReader('/home/benjamin/Desktop/divoted_foil.dat') #  ('/home/benjamin/Desktop/foilshapes/E335.txt')  #
# 


#area = geo.Area(ax,ay)
#print(area)
#
#
#print(geo.Centroid(ax,ay))











import numpy as np
from stl import mesh


airfoil = geo.makeWhole(airfoil)

ax = airfoil.x(); ay = airfoil.y()


az = []; bz = []; bx = []; by = []
for i in range(len(ax)):
    az.append(0)
    bz.append(0.125 * 25.4)
    bx.append(ax[i] * 6 * 25.4)
    by.append(ay[i] * 6 * 25.4)
    ax[i] *= 6 * 25.4
    ay[i] *= 6 * 25.4

va = []
vb = []

for i in range(len(ax)):
    va.append( [ ax[i], ay[i], az[i] ] )
    vb.append( [ bx[i], by[i], bz[i] ] )


v = va + vb
#print(v)


lva = len(va)-1
f = []

# outer
for i in range(lva-1):
    tri = [     i,      i + 1,i + lva + 1 ]
    f.append(tri)
    tri = [ i + 1,i + lva + 2,i + lva + 1 ]
    f.append(tri)
# front side

topc = 0; botc = len(va)-2; top = True
while(topc < botc-1):
    if(top):
        tri = [topc,topc+1,botc]
        f.append(tri)
        topc += 1
        top = False
    else:
        tri = [topc,botc-1,botc]
        f.append(tri)
        botc -= 1
        top = True
#print(f)

#  back side
topc = lva + 1; botc = len(v)-2; top = True
while(topc < botc-1):
    if(top):
        tri = [topc,topc+1,botc]
        f.append(tri)
        topc += 1
        top = False
    else:
        tri = [topc,botc-1,botc]
        f.append(tri)
        botc -= 1
        top = True
#print(f)



# Define the 8 vertices of the cube
vertices = np.array(v)

# Define the 12 triangles composing the cube
faces = np.array(f)

# Create the mesh
obj = mesh.Mesh(np.zeros(faces.shape[0], dtype=mesh.Mesh.dtype))
for i, f in enumerate(faces):
    for j in range(3):
        obj.vectors[i][j] = vertices[f[j],:]

# Write the mesh to file "cube.stl"
obj.save('wing.stl')

import path_maker as pm

pm.move_file_to_Desktop('wing.stl')











#g = [[],[]]
#g[0].append('food')
#print(g[0][0])




