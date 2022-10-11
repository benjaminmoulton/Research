#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 11 09:45:25 2019

@author: benjamin
"""

from matplotlib import pyplot as plt
import gcode as gc



x,y = gc.NACAFileReader('/home/benjamin/Desktop/foilshapes/NACA 0024.txt')
#plt.axis('equal')
#plt.plot(x,y)
#plt.show()

def shift_Spine(midx,midy,isTop,t,spine,isTlocal,interp, ppp):
    
    # if calculations are for the top, sign is set
    sign = -1
    if(isTop):
        sign = 1
        
    # initialize spine line to make
    
    shiftx = []
    shifty = []
    
    # start with the fillet radius 
    
    # set the fillet radius:
    fillet_radius = spine * t * 0.1
    fillet_points = 10
    
    for r in range(fillet_points):
        fillet_x = fillet_radius * r / fillet_points
        spine_width = t * spine / 2
        shiftx.append( fillet_x  + midx[0] )
        shifty.append( -sign*( fillet_radius**2 - ( fillet_x - fillet_radius )**2 )**0.5 + sign * fillet_radius + sign * spine_width )
    
    # spine thickness set
    
    spine_thickness = spine * t / 2
    
    # the points are set
    
    p = 0
    while( p < len(midx)-1 ):
        if(midx[p] < ( fillet_radius + midx[0] )):
            p += 1
        elif(midx[p] > ( fillet_radius + midx[0] ) and midx[p] < ( midx[len(midx)-1] - fillet_radius)):
            inverse_slope = 1 / - (  (midy[p+1] - midy[p]) / (midx[p+1] - midx[p])  )
            x_midpoint = ((midx[p+1] - midx[p]) / 2 ) + midx[p]
            if(isTlocal):
                spine_thickness = spineWidth*interp(x_midpoint)
            y_midpoint = ((midy[p+1] - midy[p]) / 2 ) + midy[p]
            xchange = ( spine_thickness**2 / (1 + (inverse_slope**2)) )**(1/2)
            ychange = ( spine_thickness**2 / (1 + (inverse_slope**-2)) )**(1/2)
            xtnew = x_midpoint - sign * (xchange if inverse_slope < 0 else -xchange)
            ytnew = y_midpoint + sign * ychange
            shiftx.append(xtnew)
            shifty.append(ytnew)
            p += 1
        else:
            p += 1
    
    dv = []
    for a in range(len(midx)):
        dv.append(spine_thickness)
    
    pphp = ppp // 2
    
    minus = 0
    if(isTop):
        minus = 1
    else:
        minus = -1
        
    shiftx,shifty = deintersect(midx,midy,shiftx,shifty,len(midx),pphp,spine_thickness,dv,minus)
    
    
    
    
    
    
    
    
    
    for r in range(fillet_points):
        fillet_x = fillet_radius * r / fillet_points
        spine_width = t * spine / 2
        shiftx.append( fillet_x  + midx[len(midx)-1] )
        shifty.append( -sign*( fillet_radius**2 - ( fillet_x )**2 )**0.5 + sign * fillet_radius + sign * spine_width )
    
    return shiftx,shifty



def deintersect(x,y,spinex,spiney,n,pphp,d,d_vals,minus):
    indices_to_remove = []
    
    # a nested for loop loops through each point of the spine line,
    #  and for every point the distance is calculated to the values in the previous
    #  half period, and the following half period.
    #  If the distance between any of these points is smaller than 99% of the 
    #  expected theoretical distance, the point index is added to a list
    #  for removal from the list of spine points
    
    for a in range(n):
        
        for b in range(-pphp,pphp):
            
            if( ( a + b >= 0 ) and ( a + b < n ) and minus * y[a+b] < 0): # ): #
                
                dist = ( ( spinex[a] - x[a+b] )**2 + ( spiney[a] - y[a+b] )**2 )**0.5
#                if(a < 5) : 
#                    print(dist < d) 
#                    if(dist < d): print(a+b)
                
                if(dist < d_vals[a]*0.95):
                    
                    indices_to_remove.append(a)
    
    # duplicate values are removed from the list
    
    indices = list(set(indices_to_remove))
    
    # and the list is sorted
    
    indices.sort()
    
#    print(indices)    
    
    # the values at the indices found are removed from the spine list
    #  the removals are done so that the order is unaffected of previous values
        
    for j in range(len(indices)):
        
        spinex.pop(indices[j] - j)
        
        spiney.pop(indices[j] - j)
#    
    return spinex,spiney







def deintersecter(ax,ay,bx,by,n,pphp,d):
    indices_to_remove = []
    
    # a nested for loop loops through each point of the spine line,
    #  and for every point the distance is calculated to the values in the previous
    #  half period, and the following half period.
    #  If the distance between any of these points is smaller than 99% of the 
    #  expected theoretical distance, the point index is added to a list
    #  for removal from the list of spine points
    
    for a in range(n):
        
        for b in range(-pphp,pphp):
            
            if( ( a + b >= 0 ) and ( a + b < n ) ): # ): #
                
                dist = ( ( bx[a] - ax[a+b] )**2 + ( by[a] - ay[a+b] )**2 )**0.5
#                if(a < 5) : 
#                    print(dist < d) 
#                    if(dist < d): print(a+b)
                
                if(dist < 2 * abs(d)*0.95):
                    
                    indices_to_remove.append(a)
    
    # duplicate values are removed from the list
    
    indices = list(set(indices_to_remove))
    
    # and the list is sorted
    
    indices.sort()
    
#    print(indices)    
    
    # the values at the indices found are removed from the spine list
    #  the removals are done so that the order is unaffected of previous values
        
    for j in range(len(indices)):
        
        if(ax[indices[j] - j] < 0):
            ax.pop(indices[j] - j)
            
            ay.pop(indices[j] - j)
            
        else:
            bx.pop(indices[j] - j)
            
            by.pop(indices[j] - j)
    
    return ax,ay,bx,by








def Yht(t,x): #maxThickness, percentOfChord
    
    # a values are set
    
    a0 = 0.2969
    
    a1 = -0.126
    
    a2 = -0.3516
    
    a3 = 0.2843
    
    a4 = -0.1036
    
    # if x is a list, run through for each point
    
    if(type(x) == list):
        Yht = []
        
        for val in x:
            
            y = ( t / 0.2 ) * ( a0*val**0.5 + a1*val + a2*val**2 + a3*val**3 + a4*val**4 )
            
            Yht.append(y)
    
    elif(type(x) == float or type(x) == int):
        
        y = ( t / 0.2 ) * ( a0*x**0.5 + a1*x + a2*x**2 + a3*x**3 + a4*x**4 )
        
        Yht = y
    
    return Yht


def dA(t,x):
    
    # a values are set
    
    a0 = 0.2969
    
    a1 = -0.126
    
    a2 = -0.3516
    
    a3 = 0.2843
    
    a4 = -0.1036
    
    dy = ( t / 0.2 ) * ( 0.5*a0*x**-0.5 + a1 + 2*a2*x + 3*a3*x**2 + 4*a4*x**3 )
    
    return dy


def addSpine(x,y,t,f,d,pphp,cos,sin,pi,sign,isT,addh):
    
    # the spine coordinate lists are initialized
    spinex = []
    
    spiney = []
    
    # number of points for iteration is determined
    
    n = len(x)
    
    # determine the period from the frequency
    
    T = 2 * pi / f
    
    # spinewidth is separated from d for use if is Tlocal shape
    
    sw = abs(d)
    
    minus = 0
    
    if(d < 0):
        minus = 1
    else:
        minus = -1
    
    # d values will be set to a list for use when determining self intersecting points
    
    d_vals = []
    
    # a for loop runs through each x y coordinated given, and determines
    #  the spined points on either side of the middle wave
    
    for i in range(n):
        
        # the individual points are set for calculation
        
        xval = x[i]
        
        yval = y[i]
    
        # a constant used when the half period pieces are used
        #  the factor is set to 1 for the middle, 0.5 on the edges
        
        # similarly, an additional value is created based on which section is being processed
        
        factor = 1
        
        addition = 0
        
        if(xval <= min(x) + T/2 or xval >= max(x) - T/2):
            
            factor = 0.5
            
            if(xval <= min(x) + T/2): addition =        sign * factor * dA(t,xval)
            if(xval >= max(x) - T/2): addition = addh * sign * factor * dA(t,xval)
        
        # if is T Local calculation, d is set for each calculation
        
        if(isT):
            
            d = sw * Yht(t,xval)
        
        # for use in later calculation, the d values are stored to a list
        
        d_vals.append(abs(d))
        
        # the dx and dy values are determined for future use
        
        dx = 1
        
        dy = factor * sign * (   dA(t,xval) * cos(f * xval) + Yht(t,xval) * f * sin(f * xval)   ) + addition
        
        # the shift values are found for use in determining the d coordinates
        
        xshift = d * dy / ( ( dx**2 + dy**2 )**0.5 )
        
        yshift = d * dx / ( ( dx**2 + dy**2 )**0.5 )
        
        # the d coordinates are found
        
        xd = xval + xshift
        
        yd = yval - yshift
        
        # the coordinates are added to the spinex and spiney lists
        
        spinex.append(xd)
        
        spiney.append(yd)
    
    # 2 lists are initialized to hold the indices where self intersection occurs
    #  on the spine wave
    
    spinex, spiney = deintersect(x,y,spinex,spiney,n,pphp,d,d_vals,minus)
    
    ##############################################################################3333
#    print(max(spinex))
    
    count = 0
    for i in range(len(spinex)):
        if(spinex[i - count] <= min(x)):
            spinex.pop(i - count)
            spiney.pop(i - count)
            count += 1
        elif(spinex[i - count] >= max(x)):
            spinex.pop(i - count)
            spiney.pop(i - count)
            count += 1
    ######################################################333333333333333333333333333
    # the point is added to the start and end of the spine so that the curve meets correctly

    
    x_begin = x[0]
    
    y_begin = y[0] - d
    
    spinex.insert(0,x_begin)
    
    spiney.insert(0,y_begin)
    
    x_end = x[len(x) - 1]
    
    y_end = y[len(y) - 1] - d
    
    spinex.append(x_end)
    
    spiney.append(y_end)
    
    # the spinex and spiney lists are then returned to the user
    
    return spinex, spiney


def filleter(x,y,r):
    
    return 0





def OsciSpine(x,y,chord,thickness,startPerc = .0,endPerc = .88,period = 0.09,spineWidth = 0.0,fillet = 0.0, firstUp = True,isTlocal = False): # period is a % of the chordlength 
    
    # the interpolation function to be used later in calculation
    
    from scipy.interpolate import interp1d as itp
    
    # if no start Percentage is specified, start at the location of max thickness on the airfoil
    
    if(startPerc == 0):
        
        startPerc = x[y.index(max(y))] / chord
    
    # if no spine width specified, use a spine width of 20% of the thickness
    
    if(spineWidth == 0):
        
        spineWidth = 0.2
    
    # if no spinewidth is given, and T-local thickness is desired, spinethickness
    # will have no effect on the shape
    
    elif(spineWidth == 0 and isTlocal == True):
        
        spineWidth = 1
    
    # if the spine oscillation is desired to begin by moving up, the sign factor 
    #  will modify the calculation
    
    sign = 1
    
    if(firstUp == True):
        
        sign = -1
    
    # cosine function and pi constant values are passed in to be used in calculation.
    
    from numpy import cos, sin, pi
    
    # the interpolation points for the upper half of the air foil are specified
    
    x_to_interpolate = x[:len(x)//2]
    
    y_to_interpolate = y[:len(x)//2]
    
    # concerns of thickness and chord are removed through dividing by chord for
    # all points
    
    for n in range(len(x_to_interpolate)):
        
        x_to_interpolate[n] /= chord
        
        y_to_interpolate[n] /= chord
    
    # an interpolation function is formed for use in future calculation
    
    topInterp = itp(x_to_interpolate[::-1],y_to_interpolate[::-1])
    
    # the number of periods within the oscilating spine is calculated
    
    n = (endPerc - startPerc - period) / period # of periods in range excluding start and end shifts
    
    # a boolean is made to determine whether an additional half period must be added to maximize the amount
    # of the oscilating spine within the given start and end perc's
    
    add_half = False
    
    if((n - (n//1)) >= 0.5):
        
        add_half = True
    
    # the addition is added to determine the new endPerc value
    
    addition = 0.5 if add_half else 0
    
    new_endPerc = startPerc + (( (n//1) + addition + 1) * period)
    
    # the sign for the last half period is determined to be used when calculating
    # the points of the sine curve.
    
    add_half_sign = -1 if add_half else 1
    
    # a value is set to determine how many points to use in calculation
    #    30 points for one period is determined a sufficient amount to give the accuracy needed
    
    points_per_period = 60
    
    # this value is multiplied by the number of periods to give the number of points to be found
    
    z = int( points_per_period * n )
    
    # the x points are determined for the sine curve
    
    i = []
    
    xp = []
    
    ival = startPerc
    
    shift = (endPerc - startPerc) / z
    
    for w in range(z):
        
        i.append(ival + shift*w)
        
        xp.append(shift*w)
    
    
    # the lists are initialized for the middle sine wave points
    
    sinex = []
    siney = []
    
    # a counter value is initialized
    
    j = 0
    
    # the frequency of the wave is found
    
    f = 2 * pi / period
    
    # the spine half thickness as a percentage of the chord is found
    
    spine_half_thick = spineWidth * thickness / 2
    
    # a while loop runs through till the coutner reaches the number of x points

    while(j < len(xp)):
        
        # if the oscispine is a tlocal oscispine, the iterate spine half thickness is
        #  determined based on the local point
        
        if(isTlocal):
            
            spine_half_thick = spineWidth * topInterp(xp[j] + startPerc) #Yht(thickness, xp[j] + startPerc)#
        
        # the "amplitude" for the calculated point is determined based on the 
        #  airfoil interpolation, subtracting the spine half thickness
        
        interpolation_point = topInterp(xp[j] + startPerc) - spine_half_thick #Yht(thickness, xp[j] + startPerc) - spine_half_thick#
        
        # if the point is in the first half period, the x and y coordinates are found
        
        if(xp[j] <= (0.5*period)):
            
            # the x point is assigned to the x list
            
            sinex.append(i[j])
            
            # the y point is assigned to the y list
            
            siney.append( sign * interpolation_point * 0.5 * cos(f * xp[j]) - sign * 0.5 * interpolation_point )
        
        # if the point is in middle section, the x and y coordinates are added to the sine lists
        
        elif(xp[j] > (0.5*period) and xp[j] <= (new_endPerc - startPerc - 0.5*period)):
            
            # the x point is assigned to the x list
            
            sinex.append(i[j])
            
            # the y point is assigned to the y list
            
            siney.append( sign * interpolation_point * cos(f * xp[j]) )
        
        # if the point is in the last portion of the wave, the x and y coords are found
        
        elif(xp[j] > (new_endPerc - startPerc - 0.5*period) and xp[j] <= (new_endPerc - startPerc)):
            
            # the x point is assigned to the x list
            
            sinex.append(i[j])
            
            # the y point is assigned to the y list
            
            siney.append( sign * interpolation_point * 0.5 * cos(f * xp[j]) - sign * add_half_sign * 0.5 * interpolation_point )
        
        # the counter is increased for the next iteration
        
        j += 1
    
    # set the half thickness displacer value
    
    d = spineWidth * thickness / 2
    
    # set d to spineWidth if isTlocal for use later on
    
    if(isTlocal):
        
        d = spineWidth
    
    # set the number of points in a half period
    
    num = points_per_period //2 # * 10  #
    
    # the above and below sine spine lists are found
    
    ax,ay = addSpine(sinex,siney,thickness,f, d,num,cos,sin,pi,sign,isTlocal,add_half_sign)
    
    bx,by = addSpine(sinex,siney,thickness,f,-d,num,cos,sin,pi,sign,isTlocal,add_half_sign)
    
#    ax,ay,bx,by = deintersecter(ax,ay,bx,by,len(sinex),num,d)
    
    
    abovex,abovey = shift_Spine(sinex,siney,True,thickness,spineWidth,isTlocal,topInterp  , points_per_period)
    belowx,belowy = shift_Spine(sinex,siney,False,thickness,spineWidth,isTlocal,topInterp  , points_per_period)
    
    
    
    
    
    
    
    
    
    
    
    
#        if(xp[j] < r_fillet):
#            for l in range(10):
#                sinex.append(i[j] + l * r_fillet / 10)
#                siney.append(r_fillet - (r_fillet**2 - (sinex[len(sinex)-1] - r_fillet)**2 )**0.5 )
#            while(xp[j] < r_fillet):
#                j += 1
#        
    
    
    
    
    ox = []
    oy = []
    
    k = 0
    frequency_over_chord = ((2*pi)/(period*chord))
    starting_point = startPerc*chord
    spine_thickness = spineWidth*thickness/2
    for val in i: # while loop that ends where x vals are greater than int(endPerc-startPerc / (period)) + (1/2) * endPerc-startPerc / (period) + startPerc  * chord
        val = i[k]
        ox.append(val)
        interpolated_point = topInterp(ox[k])
        if(isTlocal):
            spine_thickness = spineWidth*interpolated_point
        if(ox[k] < (startPerc + period/2)*chord):
            oy.append(sign*((( interpolated_point/2 )*cos(frequency_over_chord*(ox[k] - starting_point) + pi)) + interpolated_point/2 - spine_thickness))
        elif((ox[k]/chord - startPerc)/period > int((endPerc - startPerc)/period)):
            if((ox[k]/chord - startPerc)/period < int((endPerc - startPerc)/period)+ 0.5):
                oy.append(sign*(interpolated_point - spine_thickness)*-cos(frequency_over_chord*(ox[k] - starting_point)))
            else:
                oy.append(sign*((( interpolated_point/2 )*cos(frequency_over_chord*(ox[k] - starting_point) + pi)) + interpolated_point/2 - spine_thickness))
        else:
            oy.append(sign*(interpolated_point - spine_thickness)*-cos(frequency_over_chord*(ox[k] - starting_point))) # create OsciSpine in such a way that the ending opposite
        k += 1
    
    
    osx = ox
    osy = oy
    return osx,osy,sinex,siney,abovex,abovey,belowx,belowy,ax,ay,bx,by  


#c = 1
#thck = 0.24
#startPerc = 0.30645
#endPerc = .87
#period = 0.0887
#spineWidth = 0.1
#fillet = 0.5
#firstUp = True
#isTlocal = False

c = 1
thck = 0.24
startPerc = 0.3
endPerc = .9
period = 0.09
spineWidth = 0.1
fillet = 0.5
firstUp = True
isTlocal = False

ox,oy,sinex,siney,ax,ay,bx,by,jx,jy,kx,ky  = OsciSpine(x,y,c,thck,startPerc,endPerc,period,spineWidth,fillet,firstUp,isTlocal) ## 


gc.NACAFileWriter('sharp_curve_points',jx,jy)
gc.NACAFileWriter('sharp_curve_points_a',kx,ky)
gc.NACAFileWriter('middle',sinex,siney)




plt.axis('equal') #((0.298,startPerc + spineWidth * thck * 0.1,0.01,0.025)) #
#plt.plot(ox,oy)
plt.plot(sinex,siney)
plt.plot([startPerc,endPerc],[0,0])
plt.plot(ax,ay)
plt.plot(bx,by)
plt.show()


plt.axis('equal')
plt.plot(sinex,siney)
plt.plot([startPerc,endPerc],[0,0])
plt.plot(kx,ky)
plt.plot(jx,jy)#,'ro')
#plt.plot(kx[:4],ky[:4])
plt.show()
#
#for index in range(len(ai)):
#    plt.axis('equal')
#    plt.plot(sinex,siney)
#    plt.plot(jx,jy)
#    plt.plot(jx[ai[index]],jy[ai[index]],'bo')
#plt.show()



plt.show()


#import gcode as gc
#gc.NACAFileWriter('withA',jx,jy)





# #     This portion of code shows what needs to be fixed with the fillets
#plt.axis((0.29,0.33,-0.02,0.02)) #((0.298,startPerc + spineWidth * thck * 0.1,0.01,0.025)) #
##plt.plot(ox,oy)
#plt.plot(sinex,siney)
#plt.plot([startPerc,endPerc],[0,0])
#plt.plot(ax,ay)
#plt.plot(bx,by)
#plt.show()
#plt.axis((0.838,0.842,-0.02,0.02)) #((0.298,startPerc + spineWidth * thck * 0.1,0.01,0.025)) #
##plt.plot(ox,oy)
#plt.plot(sinex,siney)
#plt.plot([startPerc,endPerc],[0,0])
#plt.plot(ax,ay)
#plt.plot(bx,by)
#plt.show()


#plt.plot(x,y)
#plt.plot([startPerc + spineWidth*thck * 0.1,startPerc + spineWidth*thck * 0.1],[-0.1,0.1])
#plt.show()
#
#plt.axis('equal') # ((0.845,0.85,-0.02,0.02))#
#plt.plot(ax,ay)
#plt.plot(bx,by)
#plt.show()
#
#plt.axis('equal')
#plt.plot(ax[-10:],ay[-10:])
#plt.plot(bx[-10:],by[-10:])
#plt.show()
#
#plt.axis('equal')
#plt.plot(ax[:10],ay[:10])
#plt.plot(bx[:10],by[:10])
#plt.show()















### the below is used to determine math for fillets
#sharp_fillet = []
#round_fillet = []
#r_sharp = spineWidth*thck * 0.02
#r_round = spineWidth*thck * 0.1
##print(r_round)
#xvals = [a * r_round / 100 for a in range(100)]
#ysharp = []
#yround = []
#for val in xvals:
#    ysharp.append( -(r_sharp**2 - (val - r_sharp)**2)**0.5 + r_sharp)
#    yround.append( -(r_round**2 - (val - r_round)**2)**0.5  + r_round)
#    
##plt.axis((0,r_round + 0.001,0,r_round))
##plt.plot(xvals,ysharp)
##plt.plot(xvals,yround)
##plt.show()
