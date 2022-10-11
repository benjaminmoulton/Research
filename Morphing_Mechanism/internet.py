#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 17 14:47:40 2019

@author: benjamin
"""

import urllib as u
u.request

page = u.request.urlopen("http://airfoiltools.com/search/index")
#print (page.read())



e335page = u.request.urlopen('http://airfoiltools.com/airfoil/details?airfoil=e335-il')

selig_page = u.request.urlopen('http://airfoiltools.com/airfoil/seligdatfile?airfoil=naca0015-il')

print( selig_page.read())




import requests

def url_ok(url):
    r = requests.head(url)
    return r.status_code == 200



print(url_ok('http://aero.go.usu.edu')  )

#print(u.request.urlopen("http://airfoiltools.com/airfoil/seligdatfile?airfoil=e305-il").getcode())

#if(0.0000001):
#    print('yes!')

#import webbrowser as w
#
#w.open('http://aero.go.usu.edu')



















