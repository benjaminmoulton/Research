#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun  3 09:00:25 2019

@author: benjamin

"""


import point as p

def Yht(t,x): #maxThickness, percentOfChord
    
    # a values are set
    
    a0 = 0.2969
    
    a1 = -0.126
    
    a2 = -0.3516
    
    a3 = 0.2843
    
    a4 = -0.1036
    
    
    Yht = 0
    
    # if x is a list, run through for each point
    
    if(type(x) == list):
        Yht = []
        
        for val in x:
            
            y = ( t / 0.2 ) * ( a0*val**0.5 + a1*val + a2*val**2 + a3*val**3 + a4*val**4 )
            
            Yht.append(y)
    
    elif(type(x) == float or type(x) == int):
        
        y = ( t / 0.2 ) * ( a0*x**0.5 + a1*x + a2*x**2 + a3*x**3 + a4*x**4 )
        
        Yht = y
    
    
    
    elif(type(x) == p.point):
        
        xval = x.x
        
        x.y = ( t / 0.2 ) * ( a0*xval**0.5 + a1*xval + a2*xval**2 + a3*xval**3 + a4*xval**4 )
        
        Yht = x
    
    return Yht



def Bezier_Cubic(points,n):
    
    # Determine the Bezier curve outlined by 4 points described in points, over
    #  n count locations along the curve.
    
    P0,P1,P2,P3 = points
    
    x = []; y = []
    
    # for loop to determine the x and y coordinates at each point
    
    for i in range(n):
        
        # determine t value as fraction of increments over n count values
        
        t = ( i / n )
        
        # determine the Bezier values in x and y coordinates for the given t value
        
        Bx = ( 1 - t )**3 * P0[0] + 3 * ( 1 - t)**2 * t * P1[0] + 3 * ( 1 - t ) * t**2 * P2[0] + t**3 * P3[0]
        
        By = ( 1 - t )**3 * P0[1] + 3 * ( 1 - t)**2 * t * P1[1] + 3 * ( 1 - t ) * t**2 * P2[1] + t**3 * P3[1]
        
        # append the x and y values to the x and y lists
        
        x.append(Bx)
        
        y.append(By)
    
    # return the lists to the user
    
    return x,y


def Bezier_Cubic_deriv(points,n):
    
    # determine the derivative of a bezier cubic function given the vertices and number of points
    
    # separate the points
    
    P0,P1,P2,P3 = points
    
    # initialize x and y lists
    
    x = []; y = []
    
    # using a for loop, cycle through each t value and determine the coorresponding 
    #  x and y values and add them to the x and y lists.
    
    for i in range(n):
        
        t = ( i / n )
        
        Bdx = 3 * ( 1 - t )**2 * (P1[0] - P0[0]) + 6 * ( 1 - t) * t * (P2[0] - P1[0]) + 3 * t**2 * (P3[0] - P2[0]) 
        
        Bdy = 3 * ( 1 - t )**2 * (P1[1] - P0[1]) + 6 * ( 1 - t) * t * (P2[1] - P1[1]) + 3 * t**2 * (P3[1] - P2[1])
        
        x.append(Bdx)
        
        y.append(Bdy)
    
    # return the x and y lists to the user
    
    return x,y



def addSpine(x,y,dxlist,dylist,high_d,low_d,t,isT):
    
    # initialize spine coordinate lists
    
    spinex = []
    
    spiney = []
    
    # find number of iterations to perform
    
    n = len(x)
    
    # d value saved for future use
    
    sw_high = high_d
    
    sw_low = low_d
    
    
    # run a for loop through each x y coordinated given, and determine
    #  the spined points on either side of the middle wave
    
    for i in range(n):
        
        # set individual points
        
        xval = x[i]
        
        yval = y[i]
        
        # if is T Local calculation, set d for each calculation
        
        if(isT):
            
            high_d = sw_high * Yht(t,xval)
            
            low_d  = sw_low  * Yht(t,xval)
        
        local_d =  high_d - ( i / n ) * (high_d - low_d)
        
        
        
        
        # set dx and dy values from list given at beginning
        
        dx = dxlist[i]
        
        dy = dylist[i]
        
        # find shift values
        
        xshift = local_d * dy / ( ( dx**2 + dy**2 )**0.5 )
        
        yshift = local_d * dx / ( ( dx**2 + dy**2 )**0.5 )
        
        # find d coordinates
        
        xd = xval + xshift
        
        yd = yval - yshift
        
        # add the coordinates to the spinex and spiney lists
        
        spinex.append(xd)
        
        spiney.append(yd)
    
    # return the lists to the user
    
    return spinex, spiney






def Bezier_Spine(x, y, c, t, start_percentage, end_percentage, P1_percentage, 
                 P2_percentage, T, spine_width_LE, spine_width_TE, first_up,
                 is_Tlocal):
    
    
    # insert default params
    
    
    
    
    # remove chord from data for ease of calculation
    
    dcx = x
    
    dcy = y
    
    if(not c == 1):
        for i in range(len(x)):
            dcx[i] /= c
            dcy[i] /= c
    
    # determine whether the first bend in the shape will be up or down,
    #  the sign will be used later in determining direction for Bezier curves
    
    sign = -1 # +1 going down, -1 going up
    
    if(first_up == False):
        
        sign = 1
    
    # import scipy.interpolate for interpolation
    
    from scipy.interpolate import interp1d as itp
    
    # determine interpolation points for the upper half of the air foil
    
    x_to_interpolate = dcx[:len(x)//2]
    
    y_to_interpolate = dcy[:len(x)//2]
    
    topitp = itp(x_to_interpolate[::-1],y_to_interpolate[::-1])
    
#    d = spine_width * t / 2
    
    nhp = int( ( end_percentage - start_percentage) // (T / 2) )
    
    bezx = []; bezy = []
    
    dx = []; dy = []
    
    # length of one quarter period
    
    qp = T / 4
    
    for i in range(nhp):
        
        # local spine width and local displacement
        
        l_sw = spine_width_LE - (i / nhp) * (spine_width_LE - spine_width_TE)
        
        l_d = l_sw * t / 2
        
        # local spine width increment and local displacement increment
        
        l_sw_i = (1 / nhp) * (spine_width_LE - spine_width_TE)
        
        l_d_i = l_sw_i * t / 2
        
        
        perc = start_percentage + i * (T / 2)
        
        x0 = perc; x1 = perc + P1_percentage * 2 * qp;
        x2 = perc + P2_percentage * 2 * qp; x3 = perc + 2 * qp;
        
        if(is_Tlocal):
            
            y0 = sign * topitp(perc) - sign * ( l_sw ) * topitp(perc)
            
            y1 = sign * topitp(perc + qp) - sign * ( l_sw + l_sw_i ) * topitp(perc + qp)
            
            y2 = -sign * topitp(perc + qp) + sign * ( l_sw + 2 * l_sw_i ) * topitp(perc + qp)
            
            y3 = -sign * topitp(perc + 2 * qp) + sign * ( l_sw + 3 * l_sw_i ) * topitp(perc + 2 * qp)
        
        else:
            
            y0 = sign * topitp(perc) - sign * ( l_d )
            
            y1 = sign * topitp(perc + qp) - sign * ( l_d + l_d_i )
            
            y2 = -sign * topitp(perc + qp) + sign * ( l_d + 2 * l_d_i )
            
            y3 = -sign * topitp(perc + 2 * qp) + sign * ( l_d + 3 * l_d_i )
        
        P0 = ( x0, y0 )
        
        P1 = ( x1, y1 )
        
        P2 = ( x2, y2 )
        
        P3 = ( x3, y3 )
        
        if(i == 0):
            
            P0 = P0[:1] + (0,)
            
            P1 = (perc + P1_percentage * 2 * qp,) + (0,) # P1[:1] + (0,)
            
            P2 = (perc + P2_percentage * 2 * qp,) + P2[1:]
            
        elif(i == nhp - 1):
            
            P1 = (perc + P1_percentage * 1 * qp,) + P1[1:]
            
            P2 = (perc + P2_percentage * 1 * qp,) + (0,)
            
            P3 = P3[:1] + (0,)
        
        points = [P0,P1,P2,P3]
        
        tnum = 30
        
        xp,yp = Bezier_Cubic(points,tnum)
        
        dxl,dyl = Bezier_Cubic_deriv(points,tnum)
        
        
        bezx = bezx + xp
        
        bezy = bezy + yp
        
        dx = dx + dxl
        
        dy = dy + dyl
        
        
        sign *= -1
    
    hi_d = spine_width_LE * t / 2
    
    lo_d = spine_width_TE * t / 2
    
    topx,topy = addSpine(bezx,bezy,dx,dy,-hi_d,-lo_d, t, is_Tlocal)
    
    bottomx,bottomy = addSpine(bezx,bezy,dx,dy,hi_d,lo_d, t, is_Tlocal)
    
    
    # add points to bx & by lists to form full bezier spine shape
    
    bx = []; by = []
    
    proceed = True
    
    top_added = False
    
    counter = 0
    
    while(proceed):
        
        # add the airfoil coordinates if adding top points of TE before bezier curve points
        
        if(dcx[counter] >= topx[len(topx)-1]):
            
            bx.append(dcx[counter])
            
            by.append(dcy[counter])
            
            counter += 1
        
        # add top bezier curve points
        
        elif(dcx[counter] < topx[len(topx)-1] and not top_added):
            
            itpTEx = topx[len(topx)-1] + 0.001
            
            itpTEy = topitp( itpTEx )
            
            bx.append(itpTEx); by.append(itpTEy)
            
            
            for i in range(len(topx)-1,-1,-1):
                
                bx.append(topx[i])
                
                by.append(topy[i])
            
            
            itpLEx = topx[0] - 0.001
            
            itpLEy = topitp( itpLEx )
            
            bx.append(itpLEx); by.append(itpLEy)
            
            
            top_added = True
            
            counter += 1
        
        # skip airfoil coordinates passed by squigley
        
        elif(dcx[counter] > bx[len(bx)-1] and top_added):
            
            counter += 1
          
        # add top LE points after adding bezier curve points
        
        else:
            while(dcx[counter] < dcx[counter-1]):
                
                bx.append(dcx[counter])
                
                by.append(dcy[counter])
                
                counter += 1
            
            proceed = False
    
    
    
    
    
    # bottom added
    
    proceed = True
    
    bottom_added = False
    
    while(proceed):
        
        if(dcx[counter] <= bottomx[0]):
            
            bx.append(dcx[counter])
            
            by.append(dcy[counter])
            
            counter += 1
        
        elif(dcx[counter] > bottomx[0] and not bottom_added):
            
            itpLEx = bottomx[0] - 0.001
            
            itpLEy = - topitp( itpLEx )
            
            bx.append(itpLEx); by.append(itpLEy)
            
            
            for i in range(len(bottomx)):
                
                bx.append(bottomx[i])
                
                by.append(bottomy[i])
            
            
            itpTEx = bottomx[len(bottomx)-1] + 0.001
            
            itpTEy = - topitp( itpTEx )
            
            bx.append(itpTEx); by.append(itpTEy)
            
            
            bottom_added = True
            
            counter += 1
        
        elif(dcx[counter] < bx[len(bx)-1] and bottom_added):
            
            counter += 1
            
        else:
            while(counter < len(dcx) and dcx[counter] > dcx[counter-1]):
                
                bx.append(dcx[counter])
                
                by.append(dcy[counter])
                
                counter += 1
            
            proceed = False
    
    
    # return chord to points
    
    for i in range(len(bx)):
        bx[i] *= c
        by[i] *= c
    
    # the Bezier curve coordinates are returned to the user
    
    return bx, by



#
#from matplotlib import pyplot as plt
#import File_handler as gc
#
#ax,ay = gc.SeligFileReader('/home/benjamin/Desktop/foilshapes/NACA 0015.txt')
#
#
#
#
#x, y = ax, ay
#c = 1
#t = 0.15
#start_perc = 0.6
#end_perc = 0.85
#
#T = 0.15
#
#P1_perc = 0.75   #1 # straight up and down # 2/3 # get closer to osci # 0.5 # original #
#P2_perc = 0.25   #0 # straight up and down # 1/3 # get closer to osci # 0.5 # original #
#
#spine_width_LE = 0.1
#spine_width_TE = 0.1  # spine_width_LE * 0.5
#
#first_up = True
#
#is_Tlocal = False
#
#bx,by = Bezier_Spine(ax,ay,c,t,start_perc,end_perc,P1_perc,P2_perc,T,spine_width_LE,spine_width_TE,first_up,is_Tlocal)
#"""
#
#Figure out how to input assymetric airfoils (ie E335) into Bezier Spine
#
#"""
##import NACApoints as nc
##
##bx,by = nc.change_chord(bx,by,6)
##ax,ay = nc.change_chord(ax,ay,6)
##
#plt.axis('equal')
#plt.plot(bx,by)
#plt.plot(ax,ay)
##plt.plot([0,1],[0,0])
#plt.show()

















#
#
#import gcode as gc
#
#gc.NACAFileWriter('Bezier_Spine',bx,by)





