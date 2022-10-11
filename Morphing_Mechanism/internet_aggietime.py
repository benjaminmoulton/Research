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
driver.get("https://aggietime.usu.edu/login/auth")

searchbox_name = 'form-username'

find = 'a01849708'


searchbox = driver.find_element_by_id(searchbox_name)

searchbox.send_keys(find)

# Find search button
passbox_name = 'form-password'

passw = 'now take sir Francis Drake'


passbox = driver.find_element_by_id(passbox_name)

passbox.send_keys(passw)

# Find search button
search_button = driver.find_element_by_name('login-submit')

# Click search
search_button.click()

# find date box
datebox_name = 'entries[0].timeHolder'
datebox = driver.find_element_by_id(datebox_name)


import datetime

now = datetime.datetime.today()

date_today = now.strftime( "%a, %d %b %Y" )

# insert date
datebox.clear()
datebox.send_keys(date_today)



hrs = 0



hrsbox_name = 'entries[0].totalHours'
hrsbox = driver.find_element_by_id( hrsbox_name )

hrsbox.send_keys(hrs)

submit_name = 'bulk-submit'

submitbox = driver.find_element_by_id( submit_name )
submitbox.click()








#stay = True
#
#while(stay):
#    print('finished? (y/n)')
#    inp = input()
#    if(inp == 'y'):
#        stay = False




import time
# wait 0.2 seconds
time.sleep(0.2)

# close window
driver.close()