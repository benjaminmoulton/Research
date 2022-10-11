#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 29 13:48:17 2019

@author: benjamin
"""
from matplotlib import pyplot as p
from math import cos,pi
import gcode as gc

ax,ay = gc.NACAFileReader('/home/benjamin/Desktop/foilshapes/NACA 0015.txt')
#p.axis('equal')
#p.plot(ax,ay)
#p.show()



def Bezier_Cubic(points,n):
    P0,P1,P2,P3 = points
    x = []
    
    y = []
    
    ts = []
    
    for i in range(n):
        t = ( i / n )
        Bx = ( 1 - t )**3 * P0[0] + 3 * ( 1 - t)**2 * t * P1[0] + 3 * ( 1 - t ) * t**2 * P2[0] + t**3 * P3[0]
        By = ( 1 - t )**3 * P0[1] + 3 * ( 1 - t)**2 * t * P1[1] + 3 * ( 1 - t ) * t**2 * P2[1] + t**3 * P3[1]
        
        ts.append(t)    # do we need t's?
        x.append(Bx)
        y.append(By)
    
    return x,y


def Bezier_Cubic_deriv(points,n):
    
    P0,P1,P2,P3 = points
    x = []
    
    y = []
    
    ts = []
    
    for i in range(n):
        t = ( i / n )
        Bx = 3 * ( 1 - t )**2 * (P1[0] - P0[0]) + 6 * ( 1 - t) * t * (P2[0] - P1[0]) + 3 * t**2 * (P3[0] - P2[0]) 
        By = 3 * ( 1 - t )**2 * (P1[1] - P0[1]) + 6 * ( 1 - t) * t * (P2[1] - P1[1]) + 3 * t**2 * (P3[1] - P2[1]) 
        
        ts.append(t)    # do we need t's?
        x.append(Bx)
        y.append(By)
    
    return x,y




#ts = []
#cx = []
#cy = []
#
#tmax = 2
#tmin = 0
#n = 100
#for i in range(n):
#    tval = (i / n) * ( tmax - tmin )
#    
#    ts.append(tval)
#
#P0 = ( 0 ,-2 )
#P1 = ( 1 ,-2 )
#P2 = ( 1 , 2 )
#P3 = ( 2 , 2 )
#pts = [P0,P1,P2,P3]
#
#x,y = Bezier_Cubic(pts,n)
#
#T = 4
#f = 2 * pi / T
#for i in range(n):
#    t = ts[i]
#    
#    cxval = t
#    cyval = -2*cos(f * t)
#    
#    cx.append(cxval)
#    cy.append(cyval)
##    Bx = ( 1 - t )**3 * P0[0] + 3 * ( 1 - t)**2 * t * P1[0] + 3 * ( 1 - t ) * t**2 * P2[0] + t**3 * P3[0]
##    By = ( 1 - t )**3 * P0[1] + 3 * ( 1 - t)**2 * t * P1[1] + 3 * ( 1 - t ) * t**2 * P2[1] + t**3 * P3[1]
##    
##    x.append(Bx)
##    y.append(By)

#p.axis('equal')
#p.plot(x,y)
#p.plot(cx,cy)
#p.legend(('Bezier','Sine'))
#p.show()


from scipy.interpolate import interp1d as itp
# the interpolation points for the upper half of the air foil are specified

x_to_interpolate = ax[:len(ax)//2]

y_to_interpolate = ay[:len(ax)//2]

topitp = itp(x_to_interpolate[::-1],y_to_interpolate[::-1])

c = 1
th = 0.15

start = 0.3
end = 0.87
T = 0.0887
spineWidth = 0.1
d = spineWidth * th / 2
nhp = int( ( end - start) // (T / 2) )
#print(nhp)

bezx = []
bezy = []

abezx = []
abezy = []

bbezx = []
bbezy = []

dx = []
dy = []

sign = -1 # +1 going down, -1 going up
for i in range(nhp):
    perc = start + i * (T / 2)
    qp = T / 4
    
    j = 2 # straight up and down # 4/3 # get closer to osci # 1 # original #
    k = 0 # straight up and down # 2/3 # get closer to osci # 1 # original #
    
    P0 = ( perc          ,           sign * topitp(perc) - sign * d )
    P1 = ( perc + j*qp     ,      sign * topitp(perc + qp) - sign * d )
    P2 = ( perc + k*qp     ,     -sign * topitp(perc + qp) + sign * d )
    P3 = ( perc + 2 * qp , -sign * topitp(perc + 2 * qp) + sign * d )
    
    if(i == 0):
        P0 = P0[:1] + (0,)
        P1 = P1[:1] + (0,)
    elif(i == nhp - 1):
        P2 = P2[:1] + (0,)
        P3 = P3[:1] + (0,)
    
    points = [P0,P1,P2,P3]
    tnum = 30
    xp,yp = Bezier_Cubic(points,tnum)
    dxl,dyl = Bezier_Cubic_deriv(points,tnum)
    
#    bezx.append(xp)
#    bezy.append(yp)
    bezx = bezx + xp
    bezy = bezy + yp
    
    dx = dx + dxl
    dy = dy + dyl
    
#    if(i == 1):
#        print(points)
    
    sign *= -1

#print(bezx)

import bez_funs as bf

abovebezx,abovebezy = bf.addSpine(bezx,bezy,dx,dy,-d)
belowbezx,belowbezy = bf.addSpine(bezx,bezy,dx,dy,d)


sinex,siney = gc.NACAFileReader('/home/benjamin/Desktop/PythonCode/aerolab/middle.txt')

#p.axis('equal')
#p.plot(bezx,bezy)
#p.plot(sinex,siney)
#p.legend(('Bezier','Sine'))
#p.show()

p.axis('equal')
p.plot(bezx,bezy)
p.plot(abovebezx,abovebezy)
p.plot(belowbezx,belowbezy)
#p.plot(ax,ay)
p.show()
#
#p.axis([0.4,0.6,0,0.1])
#p.plot(bezx,bezy)
#p.plot(abovebezx,abovebezy)
#p.plot(belowbezx,belowbezy)
#p.plot(ax,ay)
#p.show()  
    



abezx = []
abezy = []

bbezx = []
bbezy = []

sign = -1 # +1 going down, -1 going up
for i in range(nhp):
    perc = start + i * (T / 2)
    qp = T / 4
    
    P0 = ( perc          ,           sign * topitp(perc) - sign * d + d )
    P1 = ( perc + qp    ,      sign * topitp(perc + qp) - sign * d + d )
    P2 = ( perc + qp     ,     -sign * topitp(perc + qp) + sign * d + d )
    P3 = ( perc + 2 * qp , -sign * topitp(perc + 2 * qp) + sign * d + d )
    
    if(i == 0):
        P0 = P0[:1] + (d,)
        P1 = P1[:1] + (d,)
    elif(i == nhp - 1):
        P2 = P2[:1] + (d,)
        P3 = P3[:1] + (d,)
    
    points = [P0,P1,P2,P3]
    tnum = 30
    xp,yp = Bezier_Cubic(points,tnum)
    
    abezx = abezx + xp
    abezy = abezy + yp
    
    sign *= -1

sign = -1 # +1 going down, -1 going up
for i in range(nhp):
    perc = start + i * (T / 2)
    qp = T / 4
    
    P0 = ( perc          ,           sign * topitp(perc) - sign * d - d )
    P1 = ( perc + qp     ,      sign * topitp(perc + qp) - sign * d - d )
    P2 = ( perc + qp     ,     -sign * topitp(perc + qp) + sign * d - d )
    P3 = ( perc + 2 * qp , -sign * topitp(perc + 2 * qp) + sign * d - d )
    
    if(i == 0):
        P0 = P0[:1] + (0,)
        P1 = P1[:1] + (0,)
    elif(i == nhp - 1):
        P2 = P2[:1] + (0,)
        P3 = P3[:1] + (0,)
    
    points = [P0,P1,P2,P3]
    tnum = 30
    xp,yp = Bezier_Cubic(points,tnum)
    
    bbezx = bbezx + xp
    bbezy = bbezy + yp
    
    sign *= -1




#p.axis('equal')
#p.plot(ax,ay)
#p.plot(abezx,abezy)
#p.plot(bbezx,bbezy)
#p.show()
#
#p.axis([0.35,0.4,-0.1,0])
#p.plot(bezx,bezy)
#p.plot(abezx,abezy)
#p.plot(bbezx,bbezy)
#p.plot(ax,ay)
#p.show() 
















P0 = (0.345, 0.0669147111257841)
P1 = (0.3675, 0.06631048638132295)
P2 = (0.3675, -0.06631048638132295) 
P3 = (0.38999999999999996, -0.06541311284046691)
n = 30
lin = []
points = [P0,P1,P2,P3]
thex,they = Bezier_Cubic(points,n)
d = -d

slp01 = ( P1[1] - P0[1] ) / ( P1[0] - P0[0] )
slp12 = 'eat' #( P2[1] - P1[1] ) / ( P2[0] - P1[0] )
slp23 = ( P3[1] - P2[1] ) / ( P3[0] - P2[0] )

lin.append([slp01,slp12,slp23])


int01 = P0[1] - slp01 * P0[0]
int12 = ( P2[0] + P1[0] ) / 2
int23 = P2[1] - slp23 * P2[0]

lin.append([int01,int12,int23])

#print(lin)

for i in range(len(lin[1])):
    lin[1][i] += d


#print(lin)

px = lin[1][1]
py1 = lin[0][0] * px + lin[1][0]
py2 = lin[0][2] * px + lin[1][2]

pP0 = (P0[0], P0[1] + d)
pP1 = ( px , py1 )
pP2 = ( px , py2 )
pP3 = (P3[0], P3[1] + d)

ppoints = [pP0,pP1,pP2,pP3]

pax,pay = Bezier_Cubic(ppoints,n)





indices = []

#print(d)
for j in range(len(thex)):
    
    nd = ( ( thex[j] - pax[j] )**2 + ( they[j] - pay[j] )**2 )**0.5
    
    if(d - nd < d/100):
#        print(nd)
        indices.append(j)

#print(len(indices) - n)



apax = [0.345, 0.3479016666666666, 0.3506133333333333, 0.35314500000000004, 0.3555066666666667, 0.3577083333333334, 
        0.3597600000000001, 0.36167166666666656, 0.36345333333333335, 0.3651149999999999, 0.3666666666666667, 
        0.36811833333333327, 0.36948000000000003, 0.3707616666666666, 0.37197333333333327, 0.373125, 
        0.3742266666666666, 0.3752883333333332, 0.37632000000000004, 0.3773316666666667, 0.3783333333333333, 0.379335, 
        0.3803466666666666, 0.3813783333333333, 0.38244, 0.38354166666666667, 0.38469333333333333, 0.385905, 
        0.3871866666666667, 0.3885483333333333]

apay = [0.0744147111257841, 0.07390621082153348, 0.07258095258079143, 0.07049761051378826, 0.0677148587307543, 
        0.0642913713419198, 0.06028582245751512, 0.05575688618777051, 0.05076323664291635, 0.045363547933182854, 
        0.03961649416880042, 0.03358074945999927, 0.027314987917009766, 0.020877883650062164, 0.01432811076938682, 
        0.007724343385214011, 0.0011252556077740328, -0.0054104784527027935, -0.011824184685986177, 
        -0.018057188981845804, -0.02405081723005137, -0.029746395320372574, -0.035085249142579114, 
        -0.040008704586440694, -0.04445808754172699, -0.048374723898207704, -0.05169993954565254, -0.05437506037383116, 
        -0.05634141227251332, -0.05754032113146866]








#p.axis('equal')
##p.plot(thex,they)
#p.plot(pax,pay)
#p.plot(apax,apay)
#
#npax = []
#npay = []
#napax = []
#napay = []
#
#for i in range(len(pax)-1,-1,-1):
#    xpax = ( pax[i] - apax[len(apax)-1] )* -1 + apax[len(apax)-1]
#    ypay = pay[i]
#    
#    xapax = ( apax[i] - apax[len(apax)-1] )* -1 + apax[len(apax)-1]
#    yapay = apay[i]
#    
#    npax.append(xpax)
#    npay.append(ypay)
#    napax.append(xapax)
#    napay.append(yapay)
#
#p.plot(npax,npay)
#p.plot(napax,napay)
#
#
##for i in range(len(indices)):
##    p.plot(pax[i],pay[i],'bo')
#p.show()
#
##print(pax,pay)

















#### the following shows that a list can be easily taken apart
#Apples = [0,2,1,0]
#a,b,c,d = Apples
#print(a,b,c,d)












a = (0,4)
b = (4,4)
c = (0,0)
d = (4,0)
z = [a,b,c,d]
n = 30
xx,yy = Bezier_Cubic(z,n)

#p.axis('equal')
#p.plot(xx,yy)
#p.show()

































