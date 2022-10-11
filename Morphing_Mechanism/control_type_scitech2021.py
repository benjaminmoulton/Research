import numpy as np
import json
import dxf
import os
import C14_geoPRECISE as C14
from scipy import interpolate as interp
import airfoil_db as adb
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# function that reads in a json file as a dictionary
def read_json(filename):
    # import json file
    json_string=open(filename).read()
    vals = json.loads(json_string)
    return vals

# read in MachUpX distributions file
def read_dist(filename):
    # initialize data dictionary
    info = []

    # open file
    with open(filename,"r") as fn:
        col_names = fn.readline()

        # read in data
        for line in fn:
            info.append([float(point) for point in line.split()[2:]])
        
        # make a numpy 2d array
        info = np.matrix(info)

        # close file
        fn.close()
    
    return info

# read in text file
def read_text(filename):
    # initialize data dictionary
    info = []

    # open file
    with open(filename,"r") as fn:
        # read column names
        col_names = fn.readline()

        # read in data
        for line in fn:
            info.append([float(point) for point in line.split()])
        
        # make a numpy 2d array
        info = np.array(info)

        # split off x and y
        x = np.array([info[:,0]])
        y = np.array([info[:,1]])
        z = np.array([info[:,2]])

        # close file
        fn.close()
    
    return x,y,z

# read in text file
def read_c14(filename):
    # initialize arrays
    dx = []; dy = []; dz = []

    # open file
    with open(filename,"r") as fn:
        col_names = fn.readline()

        # read in data
        # get number of lines
        num_lines = int(fn.readline().split()[0])
        for i in range(num_lines):
            # initialize data dictionary
            info = []

            # get number of points
            num_points = int(fn.readline().split()[0])
            for j in range(num_points):
                line = fn.readline()
                info.append([float(point) for point in line.split()])
        
            # make a numpy 2d array
            info = np.array(info)

            # split off x and y
            x = np.array(info[:,0])
            y = np.array(info[:,1])
            z = np.array(info[:,2])

            # save to d arrays
            dx.append(x);dy.append(y);dz.append(z)
        
        # turn into numpy arrays
        dx = np.array(dx)
        dy = np.array(dy)
        dz = np.array(dz)

        # close file
        fn.close()
    
    return dx,dy,dz

# create dist dictionary
def get_dist(vals):
    # read in distributions file
    dist = read_dist(vals["distributions file"])

    # remove all but control point, chord, and dihedral
    dist = dist[:,1:7]
    dist = np.delete(dist,4,axis=1)

    # where the y value is negative (ie, left wing)
    # make dihedral negative
    for i in range(dist.shape[0]):
        if dist[i,1] < 0.0:
            dist[i,4] *= -1.
    
    # put units for cpx,cpy,cpz,chord in terms of inches and scale them
    dist[:,:4] *= 12. * vals["scale"]

    # scale specific values in the vals dictionary
    for i in range(len(vals["ducted fan cg [in]"])):
        vals["ducted fan cg [in]"][i] *= vals["scale"]
    for i in range(len(vals["ducted fan diameter, width, and length [in]"])):
        vals["ducted fan diameter, width, and length [in]"][i] *= vals["scale"]
    vals["ducted fan cover thickness [in]"] *= vals["scale"]
    vals["z bump"]["span range Behind [in]"] *= vals["scale"]
    vals["z bump"]["span range Over [in]"] *= vals["scale"]
    vals["z bump"]["span range Front [in]"] *= vals["scale"]
    vals["wing tip y start [in]"] *= vals["scale"]
    
    
    # create dxf file of lines that can be used to define planes
    # determine where the "mid point" lies -- index of last negative y value
    index = np.max(np.argwhere(dist[:,1] < 0.0)[:,0])
    
    # input a distribution for 0 if not there
    if not dist[index+1,1] == 0.0:
        # create cubic interpolations of each unit in dist as function of y (span)
        cubic = interp.interp1d
        x = []; y = []; z = []; c = []; d = []
        for i in range(dist.shape[0]):
            x.append(dist[i,0])
            y.append(dist[i,1])
            z.append(dist[i,2])
            c.append(dist[i,3])
            d.append(dist[i,4])
        x_y = cubic(y,x,kind="cubic")
        z_y = cubic(y,z,kind="cubic")
        c_y = cubic(y,c,kind="cubic")
        d_y = cubic(y,d,kind="cubic")
        dist = np.insert(dist,index+1,np.array([x_y(0.0),0.0,z_y(0.0),\
            c_y(0.0),0.0]),axis=0)

    # set start from index
    start = index + 1

    return dist,start

# translate guides
def translate_guides(guides,dist,i):
    # translate c14 guides
    for q in range(guides.shape[0]):
        # mirror across x axis
        guides[q][:,0] *= -1

        # flip iy and iz
        ia = guides[q][:,1] * -1.; guides[q][:,1] = guides[q][:,2] * 1.
        guides[q][:,2] = ia * 1.

        # rotate by the dihedral
        ry = guides[q][:,1]*np.cos(dist[i,4]) + \
            guides[q][:,2]*np.sin(dist[i,4])
        rz = guides[q][:,2]*np.cos(dist[i,4]) - \
            guides[q][:,1]*np.sin(dist[i,4])
        guides[q][:,2] = rz; guides[q][:,1] = ry

        # shift to the c/4 point
        guides[q][:,0] += dist[i,0]
        guides[q][:,1] += dist[i,1]
        guides[q][:,2] += dist[i,2]

    return guides

# determine x shift value
def get_x_shift(i,dist):
    # determine slope of plane
    m = -np.cos(dist[i,4]) / np.sin(dist[i,4])

    # read in y and z values
    yval = dist[i,1]
    zval = -dist[i,2]

    # determine intercept and diagonal values
    b = zval - m * yval
    d = b * np.cos(dist[i,4])

    # calculate normal point values ( of origin)
    vval = d * np.sin(dist[i,4])
    wval = b - d * np.cos(dist[i,4])

    # calculate x shift
    xshift = ((vval-yval)**2. + (wval-zval)**2.)**0.5

    return xshift

# transform the 2D x y arrays for the 2D dxf file
def transform_2D(x,y,dist,i,index,out_arrays=False):
    ### determine 2D values of the shape
    # mirror airfoil across y axis
    x2d = -1. * x

    # rotate by 90 degrees around origin
    rx2d = y*-1.
    ry2d = x2d*1.
    x2d = rx2d; y2d = ry2d

    # determine how much to shift down due to dihedral
    # shift back x and y
    if not dist[i,4] == 0.0:
        # solve for shift value
        xshift  = get_x_shift(i,dist)
        xsshift = get_x_shift(index,dist)

        # add x shift and y shift to proper arrays
        if out_arrays:
            x2d += xsshift
            y2d += dist[index,0]
        else:
            x2d  += xshift
            y2d  += dist[i,0]

    return x2d,y2d

# transform the 3D x y z arrays for the 3D dxf file and plotting
def transform_3D(x,y,z,dist,i,index,out_arrays=False):
    ### determine 3D values of the shape
    # mirror across x axis
    x3d = -1. * x
    
    # flip iy and iz
    a3d = y * -1.; y3d = z * 1.; z3d = a3d * 1.

    # rotate by the dihedral
    ry3d = y3d*np.cos(dist[i,4]) + z3d*np.sin(dist[i,4])
    rz3d = z3d*np.cos(dist[i,4]) - y3d*np.sin(dist[i,4])
    z3d = rz3d; y3d = ry3d

    # shift to the c/4 point
    if out_arrays:
        x3d += dist[index,0]
        y3d += dist[index,1]
        z3d += dist[index,2]
    else:
        x3d += dist[i,0]
        y3d += dist[i,1]
        z3d += dist[i,2]
    
    return x3d,y3d,z3d

# make cubic interpolations for later use
def make_interpolations(dist,vals):
    # create cubic interpolations of each unit in dist as function of y (span)
    cubic = interp.interp1d
    x = []; y = []; z = []; c = []; d = []
    for i in range(dist.shape[0]):
        x.append(dist[i,0])
        y.append(dist[i,1])
        z.append(dist[i,2])
        c.append(dist[i,3])
        d.append(dist[i,4])
    vals["x(y)"] = cubic(y,x,kind="cubic")
    vals["z(y)"] = cubic(y,z,kind="cubic")
    vals["c(y)"] = cubic(y,c,kind="cubic")
    vals["d(y)"] = cubic(y,d,kind="cubic")

    # initialize lengths var list, ismorphing bool list, and section anme list
    b = [0.0]
    c = [0.0]
    ismorphing = [True]
    name = ["bay"]
    spec = ["bay"]

    # determine lengths
    # determine where the fan portion starts
    b.append(vals["ducted fan cg [in]"][1]-\
        vals["ducted fan diameter, width, and length [in]"][1]/2.)
    c.append(vals["ducted fan cg [in]"][1]-\
        vals["ducted fan diameter, width, and length [in]"][1]/2.)
    ismorphing.append(False); name.append("fan"); spec.append("fan")
    # determine where the fan portion ends
    b.append(vals["ducted fan cg [in]"][1]+\
        vals["ducted fan diameter, width, and length [in]"][1]/2.)
    c.append(vals["ducted fan cg [in]"][1]+\
        vals["ducted fan diameter, width, and length [in]"][1]/2.)
    ismorphing.append(True); name.append("servo00"); spec.append("sv0")

    # determine the servo lengths
    n = vals["actuators per semispan"]
    for i in range(1,n):
        b.append((vals["wing tip y start [in]"] - b[2])/n*i + b[2])
        c.append((vals["wing tip y start [in]"] - c[2])/n*i + c[2])
        ismorphing.append(True)
        name.append("servo{}".format(str(i).zfill(2)))
        spec.append("sv{}".format(str(i).zfill(1)))
    
    # add the last two points, the start and end of the wing tips
    b.append(vals["wing tip y start [in]"])
    c.append(vals["wing tip y start [in]"])
    ismorphing.append(False); name.append("wingtip"); spec.append("wtp")
    b.append(dist[-1,1])
    c.append(dist[-1,1])
    # save to vals
    b = np.array(b)
    c = np.array(b)
    
    # initialize important things
    vals["b"] = b
    vals["c"] = c
    vals["is morphing"] = ismorphing
    vals["name"] = name
    vals["spec"] = spec

    return

def get_bist(vals):
    # initialize important things
    b = vals["b"]
    ismorphing = vals["is morphing"]
    name = vals["name"]
    spec = vals["spec"]

    # create new dist array which can have all the b values interpolated
    bist = np.zeros((b.shape[0],5))
    bist[:,0] = vals["x(y)"](b)
    bist[:,1] = b
    bist[:,2] = vals["z(y)"](b)
    bist[:,3] = vals["c(y)"](b)
    bist[:,4] = vals["d(y)"](b)

    return bist

def make_af(vals):
    # create airfoil shape
    add = {}
    add["geometry"] = {}
    add["geometry"]["type"] = "outline_points"
    add["geometry"]["outline_points"] = \
        vals["C14"]["airfoil dat file location"]

    # initialize airfoil database
    foil = adb.Airfoil("foil",add)
    coords = foil.get_outline_points()
    x = coords[:,0]; y = coords[:,1]
    # ensure the airfoil is closed
    x0 = x[0]; x1 = x[-1]
    if not x0 == x1:
        if x0 > x1:
            x = np.append(x,x[0])
            y = np.append(y,y[0])
        else:
            x = np.insert(x,0,x[-1])
            y = np.insert(y,0,y[-1])
    # determine camberline values
    x_cam = 0.0; y_cam = 0.0#foil.get_camber(x_cam)
    vals["x_cam"] = 0.25

    # Dallin airfoil shape
    vals["x zb-af"] = x - x_cam
    vals["y zb-af"] = y - y_cam

    # save to airfoil shape
    x = np.array([x[:100],x[99:]]) - x_cam
    y = np.array([y[:100],y[99:]]) - y_cam
    # save to vals
    vals["x af"] = x
    vals["y af"] = y

    return

def deflect_af(bist,vals):
    # intialize x y z   in d and c arrays
    xd = []; yd = []; zd = []
    xc = []; yc = []; zc = []

    # set deflection values for each
    dd = [10.0,0.001,10.0,-15.0,0.001,-20.0,25.0,0.001]
    dc = []
    for i in range(len(dd)-1):
        if i == 1:
            dc.append(dd[i])
            dc.append(dd[i])
        dc.append(dd[i])
        dc.append(dd[i+1])
    dc.append(dd[-1])
    dc.append(dd[-1])
    
    # run through each segment section
    for i in range(bist.shape[0]-1):
        # set deflection
        vals["da"] = 15.0 # deflection in deg

        # resize airfoil and add to dict
        xdi = {}; ydi = {}
        xdi["afo"] = np.array([vals["x zb-af"],vals["x zb-af"]*0.])# * bist[0,3]
        ydi["afo"] = np.array([vals["y zb-af"],vals["y zb-af"]*0.])# * bist[0,3]

        # create c14 camber dict
        O = C14.interpolate(vals["x zb-af"],vals["y zb-af"],make_camberline=True)

        # create root and tip discrete and continuous x and y values
        xrd = {}; yrd = {}; xrd["afo"] = xdi["afo"] * 1.; yrd["afo"] = ydi["afo"] * 1.
        xtd = {}; ytd = {}; xtd["afo"] = xdi["afo"] * 1.; ytd["afo"] = ydi["afo"] * 1.
        xrc = {}; yrc = {}; xrc["afo"] = xdi["afo"] * 1.; yrc["afo"] = ydi["afo"] * 1.
        xtc = {}; ytc = {}; xtc["afo"] = xdi["afo"] * 1.; ytc["afo"] = ydi["afo"] * 1.

        # deflect each
        vals["da"] = dd[i]
        C14.deflection(xrd,yrd,O,vals)
        C14.deflection(xtd,ytd,O,vals)
        vals["da"] = dc[2*i]
        C14.deflection(xrc,yrc,O,vals)
        vals["da"] = dc[2*i+1]
        C14.deflection(xtc,ytc,O,vals)

        # replace the middle few
        xrd["afo"][0,90:110] = vals["x zb-af"][90:110]
        yrd["afo"][0,90:110] = vals["y zb-af"][90:110]
        xtd["afo"][0,90:110] = vals["x zb-af"][90:110]
        ytd["afo"][0,90:110] = vals["y zb-af"][90:110]
        xrc["afo"][0,90:110] = vals["x zb-af"][90:110]
        yrc["afo"][0,90:110] = vals["y zb-af"][90:110]
        xtc["afo"][0,90:110] = vals["x zb-af"][90:110]
        ytc["afo"][0,90:110] = vals["y zb-af"][90:110]

        # shift by camber

        # retrieve x y z
        xrd = xrd["afo"][0] * bist[i,3]
        yrd = yrd["afo"][0] * bist[i,3]
        zrd = yrd * 0.
        xtd = xtd["afo"][0] * bist[i+1,3]
        ytd = ytd["afo"][0] * bist[i+1,3]
        ztd = ytd * 0.
        xrc = xrc["afo"][0] * bist[i,3]
        yrc = yrc["afo"][0] * bist[i,3]
        zrc = yrc * 0.
        xtc = xtc["afo"][0] * bist[i+1,3]
        ytc = ytc["afo"][0] * bist[i+1,3]
        ztc = ytc * 0.

        # transform
        xrd3,yrd3,zrd3 = transform_3D(xrd,yrd,zrd,bist,i,i)
        xtd3,ytd3,ztd3 = transform_3D(xtd,ytd,ztd,bist,i+1,i+1)
        xrc3,yrc3,zrc3 = transform_3D(xrc,yrc,zrc,bist,i,i)
        xtc3,ytc3,ztc3 = transform_3D(xtc,ytc,ztc,bist,i+1,i+1)

        # add to array
        xd.append(xrd3); xd.append(xtd3)
        yd.append(yrd3); yd.append(ytd3)
        zd.append(zrd3); zd.append(ztd3)
        xc.append(xrc3); xc.append(xtc3)
        yc.append(yrc3); yc.append(ytc3)
        zc.append(zrc3); zc.append(ztc3)
    
    # turn into numpy arrays
    xd = np.array(xd); yd = np.array(yd); zd = np.array(zd)
    xc = np.array(xc); yc = np.array(yc); zc = np.array(zc)

    # combine for return
    dis = [xd,yd,zd]
    con = [xc,yc,zc]

    return dis,con

def plottype(dis):

    xd,yd,zd = dis

    # plot a 3d plot of the wing
    fig = plt.figure()
    axe = fig.add_subplot(111, projection='3d')
    face = 'z'

    
    # initialize min and max vars
    xmin = xd[0][0]; xmax = xd[0][0]
    ymin = yd[0][0]; ymax = xd[0][0]
    zmin = zd[0][0]; zmax = zd[0][0]

    # plot discrete
    for i in range(xd.shape[0]):
        axe.plot(xd[i],yd[i],zd[i],c="k",zdir=face)

        # solve min and maxes
        axmin = np.min(xd[i]); axmax = np.max(xd[i])
        aymin = np.min(yd[i]); aymax = np.max(yd[i])
        azmin = np.min(zd[i]); azmax = np.max(zd[i])

        if axmin < xmin: xmin = axmin
        if axmax > xmax: xmax = axmax
        if aymin < ymin: ymin = aymin
        if aymax > ymax: ymax = aymax
        if azmin < zmin: zmin = azmin
        if azmax > zmax: zmax = azmax
        
    for i in range(0,xd.shape[0],2):
        ind = [0,50,100,150]
        for j in range(len(ind)):
            axe.plot([xd[i][ind[j]],xd[i+1][ind[j]]],\
                [yd[i][ind[j]],yd[i+1][ind[j]]],\
                    [zd[i][ind[j]],zd[i+1][ind[j]]],c="k",zdir=face)
    
    # solve for center
    xcent = (xmax + xmin) / 2.
    ycent = (ymax + ymin) / 2.
    zcent = (zmax + zmin) / 2.

    # solve for differences
    xdiff = np.abs(xmax - xmin)
    ydiff = np.abs(ymax - ymin)
    zdiff = np.abs(zmax - zmin)

    # solve for max difference
    max_diff = max([xdiff,ydiff,zdiff])

    # define limits
    x_lims = [xcent + 0.5*max_diff,xcent - 0.5*max_diff]
    y_lims = [ycent + 0.5*max_diff,ycent - 0.5*max_diff]
    z_lims = [zcent + 0.5*max_diff,zcent - 0.5*max_diff]

    # set limits
    axe.set_xlim3d(x_lims[1], x_lims[0])
    axe.set_ylim3d(y_lims[0], y_lims[1])
    axe.set_zlim3d(z_lims[1], z_lims[0])

    # Hide grid lines and axes
    plt.axis('off')
    plt.grid(b=None)
    
    axe.view_init(20,-30)
    plt.show()

    return

def plotter(dis,con):

    # plots
    plottype(dis)
    plottype(con)    

    return

def main(input_file):
    # read in json file
    if type(input_file) == str:
        vals = read_json(input_file)
    elif type(input_file) == dict:
        vals = input_file

    # get dist dictionary
    dist,start = get_dist(vals)
    # 0 - cpx, 1 - cpy, 2 - cpz, 3 - chord, 4 - dihedral
    # body fixed coordinates

    # make interpolation functions
    make_interpolations(dist,vals)

    # make bist array
    bist = get_bist(vals)

    # create airfoil
    make_af(vals)

    # transfer some values for deflection
    vals["f"] = vals["C14"]["flex start [x/c]"]
    vals["th"] = False
    vals["dty"] = "parabolic"

    # create a 3D airfoil at each begin and end segment
    dis,con = deflect_af(bist,vals)

    # plot
    plotter(dis,con)

    return


inpud = {}
inpud["distributions file"] = "horizonV0_distributions_Aug17_2020.txt"
inpud["scale"] = 1.0
inpud["ducted fan cg [in]"] = [0.00,6.5,0.00]
inpud["ducted fan diameter, width, and length [in]"] = [0.00,4.0,0.00]
inpud["ducted fan cover thickness [in]"] = 0.00
inpud["z bump"] = {}
inpud["z bump"]["span range Behind [in]"] = 0.0
inpud["z bump"]["span range Over [in]"] = 0.0
inpud["z bump"]["span range Front [in]"] = 0.0
inpud["wing tip y start [in]"] = 43.2
inpud["actuators per semispan"] = 5
inpud["C14"] = {}
inpud["C14"]["flex start [x/c]"] = 0.5
inpud["C14"]["airfoil dat file location"] ="/home/ben/foilshapes/NACA 2412.dat"

main(inpud)