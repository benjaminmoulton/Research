#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 10 11:49:17 2019

@author: benjamin
"""
import point as pt
import geometry as geo
import Airfoil_generator as ag
import Change_start as cs
import kerf as k

def Hole_maker(airfoil,skin,isLE = True, isTE = True, p = {}):
    
    airfoil.makeWhole()
    holefoil = pt.points()
    LE = pt.points()
    TE = pt.points()
    
    if(isLE):
        
        start = 0; end = 0.5;
        
        if(p['is fishbone'] and not p['is bezier']):
            end = p['fishbone start']
        
        elif(not p['is fishbone'] and p['is bezier']):
            end = p['bezier start']
        
        LE = geo.sectionMaker(airfoil,start,end)
        LE.plot()
        
        KE = LE.split()
        LE = k.Kerf_pt(LE,-skin,dont_keep_b4_1st = False)
        
        LE.plot([KE])
            
    
    if(isTE):
        
        start = 0.5; end = 1;
        if(p['is fishbone'] and not p['is bezier']):
            start = p['fishbone end']
        
        elif(not p['is fishbone'] and p['is bezier']):
            start = p['bezier end']
        
        TE = geo.sectionMaker(airfoil,start,end)
        
        TE = k.Kerf_pt(TE,-skin)
        
        TE.plot()
    
    
    
    holefoil.add(airfoil); holefoil.add(TE); holefoil.add(LE);
    
    return holefoil,LE,TE




import File_handler as fh; import Airfoil_generator as ag;

af = fh.FileReader('/home/benjamin/Desktop/foilshapes/E335.txt')
#af = ag.change_chord(af,12,1)

#af.plot()


p = {}
p['is fishbone'] = True
p['is bezier'] = False
p['fishbone start'] = 0.3
p['bezier start'] = p['fishbone start']
p['fishbone end'] = 0.9
p['bezier end'] = p['fishbone end']

isLE = True
isTE = True
holed,le,te = Hole_maker(af,0.025/12,isLE,isTE,p)

holed.plot()










