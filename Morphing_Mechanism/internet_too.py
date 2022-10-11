#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 18 09:22:15 2019

@author: benjamin
"""

from selenium import webdriver

# Using Chrome to access web
driver = webdriver.Chrome('/home/benjamin/Documents/chromedriver')

# Open the website
driver.get("http://airfoiltools.com/search/index")

searchbox_name = 'MAirfoilSearchForm_textSearch'

find = 'naca 2415'


searchbox = driver.find_element_by_id(searchbox_name)

searchbox.send_keys(find)

# Find search button
search_button = driver.find_element_by_name('yt0')

# Click search
search_button.click()

if(find[:6] == 'eppler'):
    find = 'e' + find[7:]

elif(find[:5] == 'naca '):
    find = find[:4] + find[5:]

print(find)

content = driver.find_element_by_class_name('content')

content_id = content.find_element_by_id('content')


table = content_id.find_element_by_xpath('//table/tbody/tr/td')

print()
print(table.get_attribute)
print()
print(table.get_property)
print()
print(table.id)
print()
print(table.location)

#search_results = table.find_element_by_class_name('afSearchResult')





#title = table.find_element_by_class_name('cell12')

#click = table.find_element_by_link_text('Selig format dat file')

import time

time.sleep(2)

# close window
driver.close()