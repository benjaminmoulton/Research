#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 19 09:03:42 2019

@author: benjamin
"""



def get_Desktop_path():
    
    import os
    
    file_path = os.path.abspath('Desktop')
    
    split = file_path.split(os.path.sep)
    
    path_Desktop = ''
    
    for val in split:
        
        path_Desktop += val + '/'
        
        
        
        if(val == 'Desktop'):
            
            break
    
    return path_Desktop



def get_current_path():
    
    import os
    
    file_path = os.path.abspath('')
    
    split = file_path.split(os.path.sep)
    
    path_Desktop = ''
    
    for val in split:
        
        path_Desktop += val + '/'
    
    return path_Desktop



def move_file(current_location,destination):
    
    from os import rename
    
    rename(current_location,destination)
    
    return 0

#def rename_and_move_file(current_location,destination):
#    
#    from os import rename
#    
#    rename(current_location,destination)
#    
#    return 0

def move_file_to_Desktop(name):
    
    current_path = get_current_path() + name
    
    destination = get_Desktop_path() + name
    
    from os import rename
    
    rename(current_path,destination)
    
    return 0

#print(get_current_path() )
#
#
#current = get_current_path() + 'internet.py'
#
#destiny = get_Desktop_path() + 'interneter.py'
#
#move_file(current,destiny)



















