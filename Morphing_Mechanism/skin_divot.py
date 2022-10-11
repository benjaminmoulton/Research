#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 14 15:42:21 2019

@author: benjamin
"""

import point as pt
import Bezierspine_skin as bzs
import divoter as div



def skinDivot(airfoil, start_percentage, end_percentage, P1_percentage, 
                 P2_percentage, T, spine_width_LE, spine_width_TE, first_up,
                 is_Tlocal,skin,every,diameter):
    
    skinned,new_st,new_T,small_bez_c = bzs.Bezier_skin(airfoil, start_percentage, end_percentage, P1_percentage, 
                              P2_percentage, T, spine_width_LE, spine_width_TE, first_up,
                              is_Tlocal,skin)
    
#    skinned.plot()
    
    skc = skinned.chord()
    sbc = small_bez_c
    
    new_st *= sbc / skc
    new_T *= sbc / skc
    
    
    divskin = div.Divoter(skinned, new_st, end_percentage, new_T, every,
                          diameter,first_up)
    
#    divoted.plot()
    
    
    
    
    return divskin





import File_handler as fh

af = fh.FileReader('/home/benjamin/Desktop/foilshapes/E335.txt')
#af.plot()

af.change_chord(12,1)

start = 0.5
end = 0.95
P1 = 0.75
P2 = 0.25
T = 0.091
LE = 0.1
TE = 0.1
firstup = True
isTloc = False
skin = 0.251
every = 1
diameter = 2.25 / 25.4

skd = skinDivot(af,start,end,P1,P2,T,LE,TE,firstup,isTloc,skin,every,diameter)


skd.plot()

skd.change_chord(1,12)


#fh.FileWriter('skin_divoted',skd,write_z = True, skipFirst = True)














