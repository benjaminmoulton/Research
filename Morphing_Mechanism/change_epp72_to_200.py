#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 12 15:09:03 2019

@author: benjamin
"""

import File_handler as fh
import point as pt

af = fh.FileReader('/home/benjamin/Desktop/foilshapes/E335.txt')
#ah = fh.FileReader('/home/benjamin/Desktop/foilshapes/NACA 2415.txt')
#ah.plot([af],display = ['o','o'])
af.makeWhole()
af.plot()
af.plot(display = ['o'])

from scipy.interpolate import interp1d as inte

epp_top = inte(af.x()[:af.x_min_index+1],af.y()[:af.x_min_index+1],kind = 'cubic')
epp_bot = inte(af.x()[af.x_min_index+1:],af.y()[af.x_min_index+1:],kind = 'cubic')

import Airfoil_generator as ag

cos_clust = ag.getCosineClustering(200)
print( af.points[af.x_min_index].x,  cos_clust.points[cos_clust.x_min_index].x )
#cos_clust.show()
HD = pt.points()

for i in range(len(cos_clust.points)-1,-1,-1):
#    print(i)
    newpt = pt.point(cos_clust.points[i].x,epp_top(cos_clust.points[i].x))
    
    HD.add(newpt)

for i in range(len(cos_clust.points)):
    newpt = pt.point(cos_clust.points[i].x,epp_bot(cos_clust.points[i].x))
    
    HD.add(newpt)

HD.plot(display = ['o'])
HD.plot([af])

fh.FileWriter('E335_200pts',HD)

#HD.show()