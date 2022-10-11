# Program to create the shapes and guide curves necessary
#    for the internal chamber of the ducted fan area


import numpy as np 
import math

# import transition_cubic_ben as transcub 



def behgeo(gci,xgc,ygc,zgc,cg,r,besh,zpos,ax):
    ####---------------End Contour------------------------------####
    behgeox = []
    behgeoy = []
    behgeoz = []
    for i in range(len(ygc[gci])):
        if ygc[gci][i] > (cg[1]-(r+besh)):
            start_i = i
            break
    for i in range(len(ygc[gci])):
        if ygc[gci][i] > (cg[1]+(r+besh)):
            end_i = (i-1)
            break
    #print(start_i,end_i)
    n = ((xgc[gci][end_i]-xgc[gci][start_i]) / (ygc[gci][end_i]-ygc[gci][start_i]))
    a = xgc[gci][start_i] - (n*(ygc[gci][start_i]))
    alpha = math.atan(n)

    for i in range(start_i+1,end_i,1):
        behgeox.append((n*ygc[gci][i]) +a)
        behgeoy.append(ygc[gci][i])
        behgeoz.append(zgc[gci][i]) 
    for i in range(end_i-2,start_i+3,-1):
        behgeox.append((n*ygc[gci][i]) +a)
        behgeoy.append(ygc[gci][i])
        behgeoz.append(zgc[gci+1][i]+0.125)
    behgeox.append((n*ygc[gci][start_i+1]) +a)
    behgeoy.append(ygc[gci][start_i+1])
    behgeoz.append(zgc[gci][start_i+1])
    ax.plot(behgeox,behgeoy,behgeoz, label = 'Chamber',color = 'red')
    bx = np.zeros((len(behgeox)))
    by = np.zeros((len(behgeoy)))
    bz = np.zeros((len(behgeoz)))
    bx[:] = behgeox
    by[:] = behgeoy
    bz[:] = behgeoz
    ####--------------------------------------------------------####
    ####--------------------------------------------------------####
    # Create points for Cubic Guide Curve Endpoints
    behptx = []
    behpty = []
    behptz = []
    # Order of input into pt lists is Clock-Wise from 6:00.
    oldone = 10
    oldtwo = 10
    oldthr = 10
    oldfou = 10
    i1 = 0
    i2 = 0
    i3 = 0
    i4 = 0
    for i in range(len(behgeoz)):
        one = abs(-0.7-behgeoz[i])
        two = abs(6.5-behgeoy[i])
        thr = abs(-0.5-behgeoz[i])
        fou = abs(6.5-behgeoy[i])
        if ((one < oldone) and (behgeoy[i] < cg[1])):
            i1 = i
            oldone = one
        if ((two < oldtwo) and (behgeoz[i] < min(zpos))):
            i2 = i
            oldtwo = two
        if ((thr < oldthr) and (behgeoy[i] > cg[1])):
            i3 = i
            oldthr = thr
        if ((fou < oldfou) and (behgeoz[i] > min(zpos))):
            i4 = i
            oldfou = fou
    behptx.append(behgeox[i4])
    behpty.append(behgeoy[i4])
    behptz.append(behgeoz[i4])
    behptx.append(behgeox[i3])
    behpty.append(behgeoy[i3])
    behptz.append(behgeoz[i3])
    behptx.append(behgeox[i2])
    behpty.append(behgeoy[i2])
    behptz.append(behgeoz[i2])
    behptx.append(behgeox[i1])
    behpty.append(behgeoy[i1])
    behptz.append(behgeoz[i1])
    #ax.plot(behptx,behpty,behptz,color = 'blue')
    #print('The Behind Points\n',behptx,'\n',behpty,'\n',behptz)

    # combine i1, i2, i3, i4
    indices = (i1,i2,i3,i4)

    return bx,by,bz,behptx,behpty,behptz,alpha,indices

def radcheck(r,cg,y,z):
    ydist = cg[1]-y
    zdist = cg[2]-z
    dist = math.sqrt((ydist**2) + (zdist**2))
    if dist <= r:
        return True
    else:
        return False

def frogeo(gci,xgc,ygc,zgc,cg,r,ax):
    ####---------------End Contour------------------------------####
    frogeox = []
    frogeoy = []
    frogeoz = []
    ind = []
    for i in range(len(ygc[gci])):
        if radcheck(r,cg,ygc[gci][i],zgc[gci][i]) == True:
            ind.append(i)
    start_i = min(ind)
    end_i = max(ind)
    n = ((xgc[gci][end_i]-xgc[gci][start_i]) / (ygc[gci][end_i]-ygc[gci][start_i]))
    a = xgc[gci][start_i] - (n*(ygc[gci][start_i]))
    alpha = math.atan(n)

    for i in ind:
        frogeox.append((n*ygc[gci][i]) +a)
        frogeoy.append(ygc[gci][i])
        frogeoz.append(zgc[gci][i]) 
    # Calculate the starting and ending angles from which to complete the circle
    s_yd = cg[1] - ygc[gci][start_i]
    s_zd = cg[2] - zgc[gci][start_i]
    e_yd = cg[1] - ygc[gci][end_i]
    e_zd = cg[2] - zgc[gci][end_i]
    #print('\nStarts and ends: ',s_yd,s_zd,e_yd,e_zd)
    s_theta = round((math.atan(s_zd/s_yd)) * (180/math.pi) )
    e_theta = round(((math.atan(e_zd/e_yd)) * (180/math.pi) ) + 180 )
    #print('\nthetas: ',s_theta,e_theta)
    if s_theta%6 != 0:
        s_theta += (6-s_theta%6)
    if e_theta%6 != 0:
        e_theta -= (e_theta%6)
    #print('\nthetas: ',s_theta,e_theta)
    thets = np.linspace(e_theta,s_theta,round(((e_theta-s_theta)+1)/6))
    #print('thets linspace is: ',thets)
    # Fill in the rest of the circle
    for j in range(thets.shape[0]):
        #print(j,'\t',thets[j])
        rad = thets[j] * (math.pi/180)
        frogeoy.append(cg[1] - (r*math.cos(rad)) )
        frogeoz.append(cg[2] - (r*math.sin(rad)) )
        frogeox.append((n*   (cg[1] - (r*math.cos(rad)))   ) +a)
    frogeox.append(xgc[gci][start_i])
    frogeoy.append(ygc[gci][start_i])
    frogeoz.append(zgc[gci][start_i])
    ax.plot(frogeox,frogeoy,frogeoz, label = 'Chamber',color = 'red')
    fx = np.zeros((len(frogeox)))
    fy = np.zeros((len(frogeoy)))
    fz = np.zeros((len(frogeoz)))
    fx[:] = frogeox
    fy[:] = frogeoy
    fz[:] = frogeoz
    ####--------------------------------------------------------####
    ####--------------------------------------------------------####
    # Create points for Cubic Guide Curve Endpoints
    froptx = []
    fropty = []
    froptz = []
    # Order of input into pt lists is Clock-Wise from 6:00.
    for i in range(len(frogeoz)):
        if frogeoz[i] == min(frogeoz):
            froptx.append(frogeox[i])
            fropty.append(frogeoy[i])
            froptz.append(frogeoz[i])
            i1 = i
        if frogeoy[i] == min(frogeoy):
            froptx.append(frogeox[i])
            fropty.append(frogeoy[i])
            froptz.append(frogeoz[i])
            i2 = i
        if frogeoz[i] == max(frogeoz):
            froptx.append(frogeox[i])
            fropty.append(frogeoy[i])
            froptz.append(frogeoz[i])
            i3 = i
        if frogeoy[i] == max(frogeoy):
            froptx.append(frogeox[i])
            fropty.append(frogeoy[i])
            froptz.append(frogeoz[i])
            i4 = i

    # combine i1, i2, i3, i4
    indices = (i1,i2,i3,i4)

    #print('The Front Points:\n',froptx,'\n',fropty,'\n',froptz)
    #ax.scatter(froptx,fropty,froptz, color = 'black',marker = 'o')
    return fx,fy,fz,froptx,fropty,froptz,alpha,indices

# The matrix filler subroutine for the cubic solver
def mat_fill(x,z,m,mat,b,adder):
    mat[0+adder][0] = (x**3)
    mat[0+adder][1] = (x**2)
    mat[0+adder][2] = x
    mat[0+adder][3] = 1
    b[0+adder][0] = z

    mat[1+adder][0] = 3*(x**2)
    mat[1+adder][1] = 2*x
    mat[1+adder][2] = 1
    mat[1+adder][3] = 0
    b[1+adder][0] = m
    return (mat,b)

# Recieves the x,y coordinate of the 2 points and
#   the 2 endpoints slopes. Outputs the cubic coefficients
def cubic_solver(x1,y1,x2,y2,m1,m2):

    mat = np.zeros((4,4))
    b = np.zeros((4,1))

    (mat,b) = mat_fill(x1,y1,m1,mat,b,0)
    (mat,b) = mat_fill(x2,y2,m2,mat,b,2)

    inv = np.linalg.inv(mat)
    cub = np.dot(inv,b)
    return cub

# Solves the cubic at a given x-position
def val_finder(exam,cub):
    out0 = cub[0]*(exam**3) + cub[1]*(exam**2) + cub[2]*exam + cub[3]
    out = float(out0[0])
    return out 

def gc_cubs(x1,y1,z1,x2,y2,z2,ax,cg,r,shift):
    if cg == 0:
        x = np.linspace(x2,x1,100)
        y = np.zeros((x.shape))
        z = np.zeros((x.shape))
        cub_y = cubic_solver(x1,y1,x2,y2,0,0)
        cub_z = cubic_solver(x1,z1,x2,z2,0,0)
        for i in range(x.shape[0]):
            y[i] = val_finder(x[i],cub_y)
            z[i] = val_finder(x[i],cub_z)
        #print('\nThe x y and z: ',x,'\ny is:',y,'\nz is:',z)
        #print(type(x),type(y),type(z))
        #print(x.shape,y.shape,z.shape)
    else:
        x = np.zeros((102))
        y = np.zeros((x.shape))
        z = np.zeros((x.shape))
        x[0] = x2
        y[0] = y2
        z[0] = z2
        x[1] = cg[0]+shift
        y[1] = cg[1]
        z[1] = cg[2]+r
        x[2:] = np.linspace(cg[0],x1,100)
        cub_y = cubic_solver(x1,y1,cg[0],cg[1],0,0)
        cub_z = cubic_solver(x1,z1,cg[0],cg[2]+r,0,0)
        for i in range(x.shape[0]-2):
            y[i+2] = val_finder(x[i+2],cub_y)
            z[i+2] = val_finder(x[i+2],cub_z)
        #print(type(x),type(y),type(z))
        #print(x.shape,y.shape,z.shape)

    #ax.plot(x,y,z, label = 'Chamber_GC',color = 'orange')
    return x,y,z

def main(xgc,ygc,zgc,over,cg,r,shift,zpos,ax,vals):
    # Create the Chamber end geometries
    bx,by,bz,bptx,bpty,bptz,alpha_b,behind_i = behgeo(min(over)-1,xgc,ygc,zgc,cg,r,vals["z bump"]["span range Behind [in]"][0],zpos,ax)
    fx,fy,fz,fptx,fpty,fptz,alpha_f,front_i = frogeo(max(over)+1,xgc,ygc,zgc,cg,r,ax)
    # Create the Chamber Guide Curves
    gc1x,gc1y,gc1z = gc_cubs(bptx[0],bpty[0],bptz[0],fptx[0],fpty[0],fptz[0],ax,cg,r,shift)
    gc2x,gc2y,gc2z = gc_cubs(bptx[1],bpty[1],bptz[1],fptx[1],fpty[1],fptz[1],ax,0,0,0)
    gc3x,gc3y,gc3z = gc_cubs(bptx[2],bpty[2],bptz[2],fptx[2],fpty[2],fptz[2],ax,cg,-r,shift)
    gc4x,gc4y,gc4z = gc_cubs(bptx[3],bpty[3],bptz[3],fptx[3],fpty[3],fptz[3],ax,0,0,0)

    gc1 = (gc1x,gc1y,gc1z)
    gc2 = (gc2x,gc2y,gc2z)
    gc3 = (gc3x,gc3y,gc3z)
    gc4 = (gc4x,gc4y,gc4z)

    #print('\n',bx.shape,by.shape,bz.shape)
    #print(type(bx),type(by),type(bz))
    #print(fx.shape,fy.shape,fz.shape)
    #print(type(fx),type(fy),type(fz))
    #print(gc1x.shape,gc1y.shape,gc1z.shape)
    #print(type(gc1x),type(gc1y),type(gc1z))
    #print(gc2x.shape,gc2y.shape,gc2z.shape)
    #print(type(gc2x),type(gc2y),type(gc2z))
    #print(gc3x.shape,gc3y.shape,gc3z.shape)
    #print(type(gc3x),type(gc3y),type(gc3z))
    #print(gc4x.shape,gc4y.shape,gc4z.shape)
    #print(type(gc4x),type(gc4y),type(gc4z))

    #print('The behind angle is: ',alpha_b,alpha_b*(180/math.pi))
    #print('The front angle is: ',alpha_f,alpha_f*(180/math.pi))

    return bx,by,bz,alpha_b,fx,fy,fz,alpha_f,gc1,gc2,gc3,gc4,behind_i,front_i