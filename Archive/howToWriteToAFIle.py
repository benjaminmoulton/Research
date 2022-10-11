#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug  8 14:15:47 2018

@author: benjamin
"""
def main():
    f = open("cap.txt","w+")

    for i in range(10):
        f.write("This is line %d\r\n" % (i+1))
    
    f.close()



    f = open("cap.txt","a+")
    for i in range(2):
        f.write("Appended line %d\r\n" % (i+1))
    
    f.close()
    
    
    f = open("cap.txt","r")
    if f.mode == 'r':
        contents = f.read() # or f.readlines() to read 1 line at a time, then for i in f: print(i)
        print(contents)
    
    
    return 0
main()
