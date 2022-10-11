#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 23 13:03:42 2019

@author: benjamin
"""
import math

kc = 1
sut = 770
ka = 57.7*(sut**(-0.718))
kb = 0.8
n = 1.5
kb = 0.8
sep = 0.5*sut
num = 104
dnew = 0
se = ka*kb*kc*sep
fo = 2 / 1000
L = 0.6 

f = 0.83
a = ((f*sut)**2) / se
b = (-1/3)*math.log(f*sut/se)
sf = a*num**b

#first guess
bval = (  6*fo*L*n/sf  )**(1/3)
print(bval)

for i in range(15):
    bvalold = bval
    nkb = 0
    d = 0.808*bvalold*10
    if(d < 51):
        nkb = (d/7.62)**(-0.107)
    else:
        nkb = 1.51*(d**(-0.157))
    print("nkb is: %f" % nkb)
    se = ka*nkb*kc*sep
    
    a = ((f*sut)**2) / se
    b = (-1/3)*math.log(f*sut/se)
    sf = a*num**b
    
    bvalnew = (  6*fo*L*n/sf  )**(1/3)
    print(bvalnew)
    bval = bvalnew