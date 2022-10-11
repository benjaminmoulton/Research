#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun  3 09:00:25 2019

@author: benjamin

"""


import point as pt
import geometry as geo


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
    
    
    
    elif(type(x) == pt.point):
        
        xval = x.x
        
        x.y = ( t / 0.2 ) * ( a0*xval**0.5 + a1*xval + a2*xval**2 + a3*xval**3 + a4*xval**4 )
        
        Yht = x
    
    elif(type(x) == pt.points):
        
        Yht = pt.points()
        
        for val in x.points:
            
            new_point = pt.point()
            
            new_point.x = val.x
            
            new_point.y = ( t / 0.2 ) * ( a0*val.x**0.5 + a1*val.x + a2*val.x**2 + a3*val.x**3 + a4*val.x**4 )
            
            Yht.add(new_point)
    
    return Yht



def Bezier_Cubic(points,n = 30):
    
    # Determine the Bezier curve outlined by 4 points described in points, over
    #  n count locations along the curve.
    
    P0,P1,P2,P3 = points
    
    Bez_points = pt.points()
    
    # for loop to determine the x and y coordinates at each point
    
    for i in range(n):
        
        B = pt.point()
        
        # determine t value as fraction of increments over n count values
        
        t = ( i / n )
        
        # determine the Bezier values in x and y coordinates for the given t value
        
        B.x = ( 1 - t )**3 * P0.x + 3 * ( 1 - t)**2 * t * P1.x + 3 * ( 1 - t ) * t**2 * P2.x + t**3 * P3.x
        
        B.y = ( 1 - t )**3 * P0.y + 3 * ( 1 - t)**2 * t * P1.y + 3 * ( 1 - t ) * t**2 * P2.y + t**3 * P3.y
        
        # append the point to the points list
        
        Bez_points.add(B)
    
    # return the lists to the user
    
    return Bez_points


def Bezier_Cubic_deriv(points,n):
    
    # determine the derivative of a bezier cubic function given the vertices and number of points
    
    # separate the points
    
    P0,P1,P2,P3 = points
    
    # initialize pointslist
    
    Bez_derivs = pt.points()
    
    # using a for loop, cycle through each t value and determine the coorresponding 
    #  x and y values and add them to the x and y lists.
    
    for i in range(n):
        
        Bd = pt.point()
        
        t = ( i / n )
        
        Bd.x = 3 * ( 1 - t )**2 * (P1.x - P0.x) + 6 * ( 1 - t) * t * (P2.x - P1.x) + 3 * t**2 * (P3.x - P2.x) 
        
        Bd.y = 3 * ( 1 - t )**2 * (P1.y - P0.y) + 6 * ( 1 - t) * t * (P2.y - P1.y) + 3 * t**2 * (P3.y - P2.y)
        
        # add point to derivs list
        
        Bez_derivs.add(Bd)
    
    # return the x and y lists to the user
    
    return Bez_derivs



def addSpine(mid_spine,dxdy,high_d,low_d,t,isT):
    
    # initialize spine coordinate lists
    
    spine = pt.points()
    
    # find number of iterations to perform
    
    n = len(mid_spine.points)
    
    # d value saved for future use
    
    sw_high = high_d
    
    sw_low = low_d
    
    
    # run a for loop through each x y coordinated given, and determine
    #  the spined points on either side of the middle wave
    
    for i in range(n):
        
        new_point = pt.point()
        
        # set individual points
        
        new_point.x = mid_spine.points[i].x
        
        new_point.y = mid_spine.points[i].y
        
        # if is T Local calculation, set d for each calculation
        
        if(isT):
            
            high_d = sw_high * Yht(t,new_point.x)
            
            low_d  = sw_low  * Yht(t,new_point.x)
        
        local_d =  high_d - ( i / n ) * (high_d - low_d)
        
        
        
        
        # set dx and dy values from list given at beginning
        
        dx = dxdy.points[i].x
        
        dy = dxdy.points[i].y
        
        # find shift values
        
        xshift = local_d * dy / ( ( dx**2 + dy**2 )**0.5 )
        
        yshift = local_d * dx / ( ( dx**2 + dy**2 )**0.5 )
        
        # find d coordinates
        
        new_point.x += xshift
        
        new_point.y -= yshift
        
        # add the coordinates to the spine lists
        
        spine.add(new_point)
    
    # return the lists to the user
    
    return spine




def Bezier_Spine(airfoil, start_percentage, end_percentage, P1_percentage, 
                 P2_percentage, T, spine_width_LE, spine_width_TE, first_up,
                 is_Tlocal):
    
    # insert default params
    c = airfoil.chord()
    t = airfoil.percentThickness()
    
    
    # remove chord from data for ease of calculation
    
    dechorded_airfoil = airfoil.split()
    
#    dechorded_airfoil.show()
#    print('here are dechordeds')
#    print(dechorded_airfoil.points[dechorded_airfoil.x_max_index].x)
#    dechorded_airfoil.points[0].coord()
#    dechorded_airfoil.points[len(dechorded_airfoil.points)-1].coord()
    
    if(not c == 1):
        dechorded_airfoil.unityChord()
#        for i in range(len(airfoil.points)):
#            
#            dechorded_airfoil.points[i].x /= c
#            
#            dechorded_airfoil.points[i].y /= c
        
        
#    print( dechorded_airfoil.points[len(dechorded_airfoil.points)-1].x , (  airfoil.points[len(airfoil.points)-1].x / (2 * airfoil.chord()) ) )
#    if(dechorded_airfoil.points[len(dechorded_airfoil.points)-1].x == (  airfoil.points[len(airfoil.points)-1].x / (2 * airfoil.chord()) ) ):
#        print('allarum! allarum!')
#        dechorded_airfoil.points[len(dechorded_airfoil.points)-1].x *= airfoil.chord()
#        dechorded_airfoil.points[len(dechorded_airfoil.points)-1].y *= airfoil.chord()
        
    
#    dechorded_airfoil.points[0].coord()
#    dechorded_airfoil.points[len(dechorded_airfoil.points)-1].coord()
#    dechorded_airfoil.points[1].coord()
#    dechorded_airfoil.points[len(dechorded_airfoil.points)-2].coord()
#    print('here were dechordeds')
#    print();print();print();
    
    
    # determine whether the first bend in the shape will be up or down,
    #  the sign will be used later in determining direction for Bezier curves
    
    sign = 1 # +1 going down, -1 going up
    
    if(first_up == False):
        
        sign = -1
    
#    d = spine_width * t / 2
    
    nhp = int( ( end_percentage - start_percentage) // (T / 2) )
    
    bez = pt.points()
    
    dxdy = pt.points()
    
    # length of one quarter period
    
    qp = T / 4
    
    for i in range(nhp):
        
        # set up and down for assymetrical airfoils
        
        if(sign < 0):
            topi = True
            boti = False
            
        else:
            topi = False
            boti = True
                
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
            """
            This needs to be fixed for topi and boti ( assymetric airfoils )
            """
            y0 = sign * dechorded_airfoil.interp(perc,topi) - sign * ( l_sw ) * dechorded_airfoil.interp(perc)
            
            y1 = sign * dechorded_airfoil.interp(perc + qp,topi) - sign * ( l_sw + l_sw_i ) * dechorded_airfoil.interp(perc + qp)
            
            y2 = -sign * dechorded_airfoil.interp(perc + qp,boti) + sign * ( l_sw + 2 * l_sw_i ) * dechorded_airfoil.interp(perc + qp)
            
            y3 = -sign * dechorded_airfoil.interp(perc + 2 * qp,boti) + sign * ( l_sw + 3 * l_sw_i ) * dechorded_airfoil.interp(perc + 2 * qp)
        
        else:
            
##            y0 = sign * topitp(perc) - sign * ( l_d )
##            
##            y1 = sign * topitp(perc + qp) - sign * ( l_d + l_d_i )
##            
##            y2 = -sign * topitp(perc + qp) + sign * ( l_d + 2 * l_d_i )
##            
##            y3 = -sign * topitp(perc + 2 * qp) + sign * ( l_d + 3 * l_d_i )
            
            y0 = dechorded_airfoil.interp(perc,topi) + sign * ( l_d )
            
            y1 = dechorded_airfoil.interp(perc + qp,topi) + sign * ( l_d + l_d_i )
            
            y2 = dechorded_airfoil.interp(perc + qp,boti) - sign * ( l_d + 2 * l_d_i )
            
            y3 = dechorded_airfoil.interp(perc + 2 * qp,boti) - sign * ( l_d + 3 * l_d_i )
        
        P0 = pt.point(x0, y0 )
        
        P1 = pt.point( x1, y1 )
        
        P2 = pt.point( x2, y2 )
        
        P3 = pt.point( x3, y3 )
        
        if(i == 0):
            
            TE_loc = P0.x
            
            top_val = dechorded_airfoil.interp(TE_loc); bot_val = dechorded_airfoil.interp(TE_loc,False)
            
            mid_val = ( top_val + bot_val ) / 2
            
            P0.y = mid_val
            
            P1.x = perc + P1_percentage * 2 * qp; P1.y = mid_val
            
            P2.x = perc + P2_percentage * 2 * qp
            
        elif(i == nhp - 1):
            
            P1.x = perc + P1_percentage * 1 * qp
            
            TE_loc = P3.x
            
            top_val = dechorded_airfoil.interp(TE_loc); bot_val = dechorded_airfoil.interp(TE_loc,False)
            
            mid_val = ( top_val + bot_val ) / 2
            
            P2.x = perc + P2_percentage * 1 * qp; P2.y = mid_val
            
            P3.y = mid_val
        
        points = [P0,P1,P2,P3]
        
        tnum = 30
        
        bez_mini = Bezier_Cubic(points,tnum)
        
        dxdy_mini = Bezier_Cubic_deriv(points,tnum)
        
        
        bez.add(bez_mini)
        
        dxdy.add(dxdy_mini)
        
        
        sign *= -1
    
    hi_d = spine_width_LE * t / 2
    
    lo_d = spine_width_TE * t / 2
    
    top = addSpine(bez,dxdy,-hi_d,-lo_d, t, is_Tlocal)
    
    bottom = addSpine(bez,dxdy,hi_d,lo_d, t, is_Tlocal)
    
    
    # add points to bx & by lists to form full bezier spine shape
    
    bezier = pt.points()
    
    proceed = True
    
    top_added = False
    
    counter = 0
    
    while(proceed):
        
        # add the airfoil coordinates if adding top points of TE before bezier curve points
#        print(dechorded_airfoil.points[counter].x , top.points[len(top.points)-1].x)
        
        if(dechorded_airfoil.points[counter].x >= top.points[len(top.points)-1].x and not top_added):
            
            bezier.add(dechorded_airfoil.points[counter])
            
            counter += 1
        
        # add top bezier curve points
        
        elif(dechorded_airfoil.points[counter].x < top.points[len(top.points)-1].x and not top_added):
            
            interp_point = pt.point()
            
            interp_point.x = top.points[len(top.points)-1].x + 0.001
            
            interp_point.y = dechorded_airfoil.interp( interp_point.x )
            
            bezier.add(interp_point)
            
            
            top_reversed = top.split()
            
            top_reversed.points.reverse()
            
            bezier.add(top_reversed)
            
            
            interp_point = pt.point()
            
            interp_point.x = top.points[0].x - 0.001
            
            interp_point.y = dechorded_airfoil.interp( interp_point.x )
            
            bezier.add(interp_point)
            
            top_added = True
            
            counter += 1
        
        # skip airfoil coordinates passed by squigley
        
        elif(dechorded_airfoil.points[counter].x > bezier.points[len(bezier.points)-1].x and top_added):
            
            counter += 1
          
        # add top LE points after adding bezier curve points
        
        else:
            while(dechorded_airfoil.points[counter].x < bezier.points[len(bezier.points)-1].x):
                
                bezier.add(dechorded_airfoil.points[counter])
                
                counter += 1
            
            proceed = False
    
    
    
    
    
    # bottom added
    
    proceed = True
    
    bottom_added = False
    
    while(proceed):
        
        if(dechorded_airfoil.points[counter].x <= bottom.points[0].x):
            
            bezier.add(dechorded_airfoil.points[counter])
            
            counter += 1
        
        elif(dechorded_airfoil.points[counter].x > bottom.points[0].x and not bottom_added):
            
            interp_point = pt.point()
            
            interp_point.x = bottom.points[0].x - 0.001
            
            interp_point.y = dechorded_airfoil.interp( interp_point.x, False )
            
            bezier.add(interp_point)
            
            
            bezier.add(bottom)
            
            
            interp_point = pt.point()
            
            interp_point.x = bottom.points[len(bottom.points)-1].x + 0.001
            
            interp_point.y = dechorded_airfoil.interp( interp_point.x, False )
            
            bezier.add(interp_point)            
            
            bottom_added = True
            
            counter += 1
        
        elif(dechorded_airfoil.points[counter].x < bezier.points[len(bezier.points)-1].x and bottom_added):
            
            counter += 1
            
        else:
            while(counter < len(dechorded_airfoil.points) ): #and dechorded_airfoil.points[counter].x > dechorded_airfoil.points[counter-1].x):
                
                bezier.add(dechorded_airfoil.points[counter])
                
                counter += 1
            
            proceed = False
    
    
    # return chord to points
    
    for i in range(len(bezier.points)):
        bezier.points[i].x *= c; bezier.points[i].y *= c;
    
    # check if airfoil shape is whole
#    bezier.points[0].coord()
#    bezier.points[len(bezier.points)-1].coord()
    bezier.makeWhole()
#    bezier.points[0].coord()
#    bezier.points[len(bezier.points)-1].coord()
    
    
    # the Bezier curve coordinates are returned to the user
    
    return bezier


# import File_handler as gc

# ax,ay = gc.SeligFileReader('/home/ben/foilshapes/E335.txt')

# airfoil = pt.points(ax,ay)


# x, y = ax, ay
# c = 1
# t = 0.15
# start_perc = 0.3
# end_perc = 0.9

# T = 0.09

# P1_perc = 0.75   #1 # straight up and down # 2/3 # get closer to osci # 0.5 # original #
# P2_perc = 0.25   #0 # straight up and down # 1/3 # get closer to osci # 0.5 # original #

# spine_width_LE = 0.1
# spine_width_TE = 0.1  # spine_width_LE * 0.5

# first_up = True

# is_Tlocal = False

# bezier = Bezier_Spine(airfoil,c,t,start_perc,end_perc,P1_perc,P2_perc,T,spine_width_LE,spine_width_TE,first_up,is_Tlocal)

# airfoil.plot([bezier])
# bezier.plot(axis = [0.95,1.05,-0.01,0.01])
# #airfoil.plot([bezier],axis = [0.85,0.95,-0.05,0.05])





#
#
#p0 = pt.point(0,0)
#p1 = pt.point(0,0)
#p2 = pt.point(2,1)
#p3 = pt.point(2,3)
#points = [p0,p1,p2,p3]
#
#bz = Bezier_Cubic(points)
#
#bz.plot()

