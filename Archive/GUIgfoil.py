#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 13 11:22:43 2018

@author: benjamin
"""

from tkinter import *

base = Tk()
frame = Frame(base)
frame.pack()

bottomframe = Frame(base)
bottomframe.pack( side = BOTTOM )

redbutton = Button(frame, text="Red", fg="red")
redbutton.pack( side = LEFT)

brownbutton = Button(frame, text="Brown", fg="brown")
brownbutton.pack( side = LEFT )

bluebutton = Button(frame, text="Blue", fg="blue")
bluebutton.pack( side = LEFT )

blackbutton = Button(bottomframe, text="Black", fg="black")
blackbutton.pack( side = BOTTOM)

greenbutton = Button(bottomframe, text ="Green", fg="green")
greenbutton.pack( side = RIGHT)

purplebutton = Button(bottomframe,text="Purple",fg="purple")
purplebutton.pack( side = RIGHT)

# Code to add widgets will go here...
base.mainloop()