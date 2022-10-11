#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun  7 10:06:10 2019

@author: benjamin
"""

di = {}

x_vals = [0,1,2,3,4]
y_vals = x_vals.copy()

di['x'] = x_vals

di['y'] = y_vals

di['spine_width'] = 0.01

print(di)

print(type(di))

print(type(di['x']))

x = di['x']

print(x)

if(type(di) == dict): print(True)