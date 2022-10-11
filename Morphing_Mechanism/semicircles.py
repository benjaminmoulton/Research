"""
This code piece creates the Concept 14 geometric shape. The input json
file should have the following shape:
{
    "C14" : {
        "airfoil dat file location" : "/home/ben/foilshapes/E335.txt",
        "chord [in]" : 12,
        "num arcs" : 5,
        "num arc points" : 30,
        "arc type" : "exponential",
        "arc type notes" : "'exponential', '90%t', or 'hinge_point'",
        "shell thickness [mm]" : 0.8,
        "flex start [x/c]" : 0.5,
        "tongue start [in]" : 1.0,
        "tongue start [in]_notes" : "distance before flex start",
        "mouth start [in]" : 2.0,
        "mouth start [in]_notes" : "distance before flex start",
        "TE wall start [x/c]" : 0.95,
        "mouth clearance [mm]" :  1.0,
        "fillet side length [mm]" : 1.0,
        "Notes" : " ",
        "show plot" : false,
        "write dxf" : false,
        "show legend" : true,
        "all black" : false,
        "dxf file path" : "/home/ben/Desktop/",
        "dxf file name" : "c14"
    }
}
"""
import airfoil_db as adb
import numpy as np
import json
import dxf
import deflect
import os
from scipy import interpolate as itp
from scipy import optimize as opt
from shapely.geometry import polygon as poly

from matplotlib import pyplot as plt

# function that reads in a json file as a dictionary
def read_json(filename):
    # import json file
    json_string=open(filename).read()
    vals = json.loads(json_string)
    return vals

# simplify a json dictionary
def simplify(vals):
    # initialize dictionary
    dicti = {}

    # transfer vals
    in2mm = 25.4
    dicti["c"]  = vals["C14"]["chord [in]"]
    dicti["n"]  = vals["C14"]["num semicircles"]
    dicti["np"]  = vals["C14"]["num semicircle points"]
    dicti["smr"]  = vals["C14"]["semicircles radius [mm]"] / dicti["c"] / in2mm
    dicti["at"]  = vals["C14"]["arc type"]
    dicti["ac"] = vals["C14"]["arc concave"]
    dicti["t"]  = vals["C14"]["shell thickness [mm]"] / dicti["c"] / in2mm
    dicti["f"] = vals["C14"]["flex start [x/c]"]
    dicti["to"] = dicti["f"] - np.abs(vals["C14"]["tongue start [in]"]) \
        / dicti["c"]
    dicti["ts"] = dicti["f"] + np.abs(vals["C14"]["tongue end [in]"]) \
        / dicti["c"]
    dicti["mo"] = dicti["f"] - np.abs(vals["C14"]["mouth start [in]"]) \
        / dicti["c"]
    if dicti["to"] < dicti["mo"]:
        raise ValueError("Tongue cannot extend beyond mouth start")
    dicti["tw"] = vals["C14"]["TE wall start [x/c]"]
    dicti["mt"] = vals["C14"]["mouth clearance [mm]"] / dicti["c"] / in2mm
    dicti["fi"] = vals["C14"]["fillet side length [mm]"] / dicti["c"] / in2mm
    dicti["da"] = vals["C14"]["deflection angle [deg]"]
    dicti["dty"] = vals["C14"]["deflection type"]
    dicti["th"] = vals["C14"]["hinge at top"]
    dicti["sk"] = vals["C14"]["skip wing box"]

    return dicti

# function that determines the midpoint between points
def midpt(x,y):

    # create x[i] and x[i+1] arrays
    x0 = x[:-1]
    x1 = x[1:]

    # create y[i] and y[i+1] arrays
    y0 = y[:-1]
    y1 = y[1:]

    # find the midpoint between each point
    midx = (x0+x1)/2.0
    midy = (y0+y1)/2.0
    return midx,midy

# function that creates an (+) inward offset of a given line
def offset(x,y,d):

    # create x[i] and x[i+1] arrays
    x0 = x[:-1]
    x1 = x[1:]

    # create y[i] and y[i+1] arrays
    y0 = y[:-1]
    y1 = y[1:]

    # get the midpoints of the given line
    midx,midy = midpt(x,y)

    # determine the offset x and y values for each midpoint
    offx = midx - d*(y1-y0)/( (x0-x1)**2.0 + (y0-y1)**2.0 )**0.5
    offy = midy + d*(x1-x0)/( (x0-x1)**2.0 + (y0-y1)**2.0 )**0.5
    
    return offx,offy

# determine where a set of x y coordinates intersect
def intersections(x,y):
    # initialize intersection info
    insct = []

    # run through each panel and determine if the panel intersects any other
    for i in range(x.shape[0]-1): # 40,41): # 
        for j in range(x.shape[0]-1): # 157,158): # 
            # ignore if same panel
            if not i == j and not (i == j+1 or i==j-1):
                # determine the slope and intercept of both lines
                mi = (y[i+1]-y[i])/(x[i+1]-x[i])
                mj = (y[j+1]-y[j])/(x[j+1]-x[j])
                bi = y[i] - mi * x[i]
                bj = y[j] - mj * x[j]

                # determine where these lines intersect
                x_int =  (bj - bi) / (mi - mj)
                y_int = mi * x_int + bi

                # find which direction is forward
                if x[i] < x[i+1]:
                    i0 = i; i1 = i+1
                else:
                    i0 = i+1; i1 = i
                if x[j] < x[j+1]:
                    j0 = j; j1 = j+1
                else:
                    j0 = j+1; j1 = j

                # determine if the intersection occurs between the two lines
                if y[i] < y[i+1]:
                    iy0 = i; iy1 = i+1
                else:
                    iy0 = i+1; iy1 = i
                if y[j] < y[j+1]:
                    jy0 = i; jy1 = j+1
                else:
                    jy0 = j+1; jy1 = j
                
                # define boolean checks
                in_ix = x_int >= x[i0] and x_int <= x[i1]
                in_iy = y_int >= y[iy0] and y_int <= y[iy1]
                in_jx = x_int >= x[j0] and x_int <= x[j1]
                in_jy = y_int >= y[jy0] and y_int <= y[jy1]

                if in_ix and in_iy and in_jx and in_jy:
                    insct.append(np.array([x_int,y_int,i,j]))
    
    # make insct a numpy array
    insct = np.array(insct)

    # remove duplicate x values
    k = 0
    while k < insct.shape[0]:
        repeats = np.argwhere(insct[k,0]==insct[:,0])
        for index in repeats:
            if not index[0] == k:
                insct = np.delete(insct,index[0],axis=0)
        k += 1
    
    return insct

# remove intersections
def de_intersect(x,y,it):
    # check if only one intersection, or more
    if it.shape[0] == 1:
        # initialize de-intersected lists
        num = int(np.abs(it[0][2]-it[0][3]))
        ix = np.zeros((num+2,))
        iy = np.zeros((num+2,))

        # add the intersection point
        ix[0] = it[0][0]
        iy[0] = it[0][1]

        # add other points
        ind_max = int(np.max(np.array([it[0][2],it[0][3]])))
        ind_min = int(np.min(np.array([it[0][2],it[0][3]])))
        ix[1:-1] = x[ind_min+1:ind_max+1]
        iy[1:-1] = y[ind_min+1:ind_max+1]

        # add the intersection point
        ix[-1] = it[0][0]
        iy[-1] = it[0][1]

    else:
        # determine where to start and end
        if it[0][0] < it[1][0]:
            end_max = int(np.max(np.array([it[0][2],it[0][3]])))
            end_min = int(np.min(np.array([it[0][2],it[0][3]])))
            str_max = int(np.max(np.array([it[1][2],it[1][3]])))
            str_min = int(np.min(np.array([it[1][2],it[1][3]])))
        else:
            end_max = int(np.max(np.array([it[1][2],it[1][3]])))
            end_min = int(np.min(np.array([it[1][2],it[1][3]])))
            str_max = int(np.max(np.array([it[0][2],it[0][3]])))
            str_min = int(np.min(np.array([it[0][2],it[0][3]])))

        # initialize de-intersected lists
        top = end_min - str_min; bot = str_max - end_max
        ix = np.zeros((top + bot + 3,))
        iy = np.zeros((top + bot + 3,))

        # add top points
        ix[1:top+1] = x[str_min+1:end_min+1]
        iy[1:top+1] = y[str_min+1:end_min+1]

        # add top points
        ix[top+2:-1] = x[end_max+1:str_max+1]
        iy[top+2:-1] = y[end_max+1:str_max+1]

        if it[0][0] < it[1][0]:
            # add the intersection point
            ix[0] = it[1][0]
            iy[0] = it[1][1]

            # add the intersection point
            ix[top+1] = it[0][0]
            iy[top+1] = it[0][1]

            # add the intersection point
            ix[-1] = it[1][0]
            iy[-1] = it[1][1]
        else:
            # add the intersection point
            ix[0] = it[0][0]
            iy[0] = it[0][1]

            # add the intersection point
            ix[top+1] = it[1][0]
            iy[top+1] = it[1][1]

            # add the intersection point
            ix[-1] = it[0][0]
            iy[-1] = it[0][1]

    return ix,iy

# function that creates an offset with the kinks removed
# assume closed shape
def inner(x,y,d):

    # find offset curve
    offx,offy = offset(x,y,d)

    # plt.axis("equal")
    # plt.plot(x,y)
    # plt.plot(offx,offy)
    # plt.show()

    # find where the intersections occur
    itrsc = intersections(offx,offy)

    # remove them
    ix,iy = de_intersect(offx,offy,itrsc)

    return ix,iy

# create camberline and interpolation thereoof
def camberline(xy):
    # initialize camberline arrays
    xy["cx"]    = np.zeros((xy["tx"].shape[0],))
    xy["cy"]    = np.zeros((xy["ty"].shape[0],))
    xy["cl"]    = np.zeros((xy["tx"].shape[0],))

    summation = 0
    # cycle through the points and determine where the camberline is situated
    for i in range(xy["cx"].shape[0]):
        # determine x and y values
        xmid = (xy["tx"][i] + xy["bx"][i]) /2.
        ymid = (xy["ty"][i] + xy["by"][i]) /2.

        # save to cx and cy arrays
        xy["cx"][i] = xmid
        xy["cy"][i] = ymid

        # add length
        if not i == 0:
            summation += ((xy["cx"][i] - xy["cx"][i-1])**2. + \
                (xy["cy"][i] - xy["cy"][i-1])**2.)**0.5
            xy["cl"][i] = summation
    
    # create interpolation of camberline
    xy["c"] = itp.interp1d(xy["cx"],xy["cy"],kind="cubic")

    # create derivative of camberline
    xy["dc"] = (xy["cy"][1:] - xy["cy"][:-1]) / (xy["cx"][1:] - xy["cx"][:-1])

    return

# create top and bottom interpolation of a set of xy points ( assume airfoil)
def interpolate(x,y,make_camberline=False):
    # find middle point
    i_mid = int(x.shape[0]/2)

    # split arrays accordingly
    if x.shape[0] % 2 == 0:
        top_x = np.flip(x[:i_mid])
        top_y = np.flip(y[:i_mid])
        bot_x = x[i_mid:]
        bot_y = y[i_mid:]
    else:
        top_x = np.flip(x[:i_mid+1])
        top_y = np.flip(y[:i_mid+1])
        bot_x = x[i_mid:]
        bot_y = y[i_mid:]
    if top_x[0] == top_x[1]:
        top_x[1] -= 1e-10
    
    # points = np.zeros((x.shape[0],2))
    # points[:,0] = x
    # points[:,1] = y  
    # # create airfoil shape
    # add = {}
    # add["geometry"] = {}
    # add["geometry"]["type"] = "outline_points"
    # add["geometry"]["outline_points"] = points

    # # initialize airfoil database
    # foil = adb.Airfoil("foil",add)
    # def ytop(x):
    #     return foil.get_camber(x) + foil.get_thickness(x)
    # def ybot(x):
    #     return foil.get_camber(x) - foil.get_thickness(x)

    # create interpolations
    T = itp.interp1d(top_x,top_y,kind="cubic")
    B = itp.interp1d(bot_x,bot_y,kind="cubic")

    # save to dictionary
    I = {}
    I["t"] = T # ytop # 
    I["b"] = B # ybot # 
    I["x"] = x
    I["y"] = y
    I["tx"] = top_x
    I["ty"] = top_y
    I["bx"] = bot_x
    I["by"] = bot_y
    I["imid"] = i_mid

    # determine camber line of airfoil
    if make_camberline:
        camberline(I)

    return I

# determine where te line is 
def te_triangle(xy,info):
    # initialize triangle arrays
    tri_x = np.zeros((4,),dtype=np.ndarray)
    tri_y = np.zeros((4,),dtype=np.ndarray)

    # find linear interpolated point where wall meets skin
    yw_top = xy["t"](info["tw"])
    yw_bot = xy["b"](info["tw"])

    # find linear interpolated point where skinned wall meets skin
    yskw_top = xy["t"](info["tw"]-info["t"])
    yskw_bot = xy["b"](info["tw"]-info["t"])

    # save points to tri arrays
    for i in range(tri_x.shape[0]):
        if i == 0:
            # add top part of triangle
            x = np.linspace(np.max(xy["tx"]),info["tw"],num=10)
            y = xy["t"](x)           
            
            # add to tri
            tri_x[i] = x
            tri_y[i] = y
            
        elif i == 1:
            # add side part of triangle
            tri_x[i] = np.array([info["tw"],info["tw"]])
            tri_y[i] = np.array([yw_top,yw_bot])

        elif i == 2:
            # add bottom part of triangle
            x = np.linspace(np.max(xy["bx"]),info["tw"],num=10)
            y = xy["b"](x)           
            
            # add to tri
            tri_x[i] = x
            tri_y[i] = y
        else:
            # add skin thickness wall to left of triangle
            xvalue = info["tw"] - info["t"]
            tri_x[i] = np.array([xvalue,xvalue])
            tri_y[i] = np.array([yskw_top,yskw_bot])

    return tri_x,tri_y

# determine where a y location is on a line
def get_y_loc(x_loc,x,y):
    # determine which indices block in the arc
    it = 0
    while x[it] < x_loc:
        it += 1

    # determine percent along between camberline indices
    perc = (x_loc - x[it-1]) / (x[it] - x[it-1])
    y_loc = perc*(y[it] - y[it-1]) + y[it-1]

    return y_loc,perc,it

# a function that implements the bisection method to find a root
def Bisect(Fun,xl,xu, Err = 1e-5, maxiter = 1000):
    
    # initialize a large number so the err is less than
    E = 10000000
    # initialize the icounter
    i = 0
    # if the given values do not encompass the root, throw error
    if not Fun(xl)*Fun(xu) < 0:
        raise ValueError("Upper and lower limits do not encompass a root")

    # make xrold
    xrold = 100

    # till threshold error is reached, continue
    while(E > Err and i < maxiter):
        
        # estimate xr
        xr = (xl + xu) / 2

        # determine which subinterval has the root
        funmult = Fun(xl) * Fun(xr)
        if funmult < 0:
            xu = xr
        elif funmult > 0:
            xl = xr
        else: # funmult = 0
            return xr
        
        # calculate Err
        if xrold != 0:
            E = abs( (xr - xrold) / xrold )
        
        # return xrold
        xrold = xr

        # add counter value
        i += 1
    # return the root
    return xr

# function that given a point on the camberline, finds where the normal
#  to this "line" (with dc/dx) will pass the top or bottom of the innerskin
def intersect(xy,it,xp,yp,is_top=True):
    # determine inverse slope
    a = -1. / xy["dc"][it]

    # determine y intercept of this line
    b = yp - a*xp

    # solve for if top or bottom
    if is_top:
        loc = "t"
    else:
        loc = "b"
    
    # define function to use with Bisect method
    def func(x):
        return (a*x+b) - xy[loc](x)
    
    # determine the x value at which this intersection occurs
    low = np.min(xy["cx"]+0.05)
    hi  = np.max(xy["cx"])
    xval = Bisect(func,low,hi)
    
    return xval

# function that determines where a circle center is given a slope and point
# that passes through the center, as well as two points on the circle
def center(xt,yt,xb,yb,r,get_concave):
    # determine center using two points and a radius
    x3 = (xt + xb) / 2.
    y3 = (yt + yb) / 2.
    dx = xb - xt; dy = yb - yt
    q = (dx**2. + dy**2.)**0.5
    d = ((r**2.)-(q/2)**2.)**0.5

    x1 = x3 - d*dy/q
    y1 = y3 + d*dx/q
    x2 = x3 + d*dy/q
    y2 = y3 - d*dx/q

    # initialize boolean to check for where the circle center should be
    x_before = x1 < xt

    if (x_before and get_concave) or (not x_before and not get_concave):
        x_cent = x1; y_cent = y1
    else:
        x_cent = x2; y_cent = y2
    
    return x_cent, y_cent

# function that determines the angle on a circle
# from np.pi to -np.pi
def get_theta(xc,yc,xp,yp):
    # determine the x and y lengths to be used to find the theta value
    xlen = xp - xc
    ylen = yp - yc

    # determine theta value
    theta = np.arctan2(ylen,xlen)

    return theta

# function that determines where a point is on an interpolation
# range given a point and radius from that point to the interpolated point
def get_point(xy,xc,yc,r,is_concave,is_top=True):
    # solve for if top or bottom
    if is_top:
        loc = "t"
    else:
        loc = "b"
    
    # get minimum and maximum of the search range
    if is_concave:
        minx = xc
        maxx = np.max(xy[loc+"x"])
    else:
        minx = np.min(xy[loc+"x"])
        maxx = xc

    # create function to be used to bracket the root
    def func(x):
        return (xc-x)**2. + (yc-xy[loc](x))**2. - r**2.
    
    # determine x location
    x_root = Bisect(func,minx,maxx)
    y_root = xy[loc](x_root)

    return x_root,y_root

# create arcs
def arcs(xy,info,shift=0.0):
    # initialize arcs arrays first index is 0 arc 1 skin
    x_arc = np.zeros((2,info["n"],info["np"]))
    y_arc = np.zeros((2,info["n"],info["np"]))
    # initialize arcs top and bottom coordinates arrays
    # [i][ ][ ] 0 - top, 1 - bottom
    # [ ][i][ ] 0 -> 1 from LE to TE
    # [ ][ ][i] 0 - arc on left (LE), 1 - arc on right (TE)
    info["xpt"] = np.zeros((2,info["n"],2))
    info["ypt"] = np.zeros((2,info["n"],2))

    # determine where flex starts on camberline
    x_flex = info["f"]
    y_flex,p_flex,ifx = get_y_loc(x_flex,xy["cx"],xy["cy"])

    # determine "length" value at flex start
    l_flex = p_flex*(xy["cl"][ifx] - xy["cl"][ifx-1]) + xy["cl"][ifx-1]

    # determine where te wall starts on camberline
    x_end = info["tw"] - info["t"]
    y_end,p_end,ie = get_y_loc(x_end,xy["cx"],xy["cy"])

    # determine "length" value at te wall
    l_end = p_end*(xy["cl"][ie] - xy["cl"][ie-1]) + xy["cl"][ie-1]

    # cycle through each arc
    for i in range(info["n"]): # 2): # 
        # set up x location
        x_loc = info["f"] + xy["arcdx"] * (i+1) + shift
        
        # determine y location
        y_loc,perc,it = get_y_loc(x_loc,xy["cx"],xy["cy"])

        # get x and y value where the intersection occurs
        x_top = intersect(xy,it,x_loc,y_loc)
        x_bot = intersect(xy,it,x_loc,y_loc,is_top=False)
        y_top = xy["t"](x_top)
        y_bot = xy["b"](x_bot)

        # various types of arcs:
        if info["at"][-2:] == "%t":
            # determine local arc radius
            r_loc = ((x_top-x_bot)**2. + (y_top-y_bot)**2.)**0.5
            r_loc /= 2.

            size_arc = float(info["at"].split("%")[0])/100.
            
            # determine resized arc radius as local is 90% of the diameter
            r_res = r_loc / size_arc
            
            # determine where the circles center is 
            x_cent,y_cent = center(x_top,y_top,x_bot,y_bot,r_res,info["ac"])
        elif info["at"] == "hinge_point":
            # get flex point of camberline
            x_cam_flex = info["f"]
            # if convex desired, put hinge point at te wall... fix this
            if not info["ac"]:
                x_cam_flex = info["tw"]
            y_cam_flex = xy["c"](x_cam_flex)

            # determine arc radius as avg of top and bottom to hingepoint
            r_res_t = ((x_top-x_cam_flex)**2. + (y_top-y_cam_flex)**2.)**0.5
            r_res_b = ((x_bot-x_cam_flex)**2. + (y_bot-y_cam_flex)**2.)**0.5
            r_res = (r_res_t + r_res_b) / 2.
            
            # determine where the circles center is
            x_cent,y_cent = center(x_top,y_top,x_bot,y_bot,r_res,info["ac"])
        elif info["at"] == "reversed_exp":
            # determine local arc radius
            r_loc = ((x_top-x_bot)**2. + (y_top-y_bot)**2.)**0.5
            r_loc /= 2.

            # determine "length" value at locale
            l_loc = perc*(xy["cl"][it] - xy["cl"][it-1]) + xy["cl"][it-1]
            
            # determine resized arc radius
            r_res = r_loc * (l_end-l_flex) / (l_loc-l_flex)
            
            # determine where the circles center is
            x_cent,y_cent = center(x_top,y_top,x_bot,y_bot,r_res,info["ac"])
        else:
            # determine local arc radius
            r_loc = ((x_top-x_bot)**2. + (y_top-y_bot)**2.)**0.5
            r_loc /= 2.

            # determine "length" value at locale
            l_loc = perc*(xy["cl"][it] - xy["cl"][it-1]) + xy["cl"][it-1]
            
            # determine resized arc radius
            r_res = r_loc * (l_end-l_flex) / (l_end-l_loc)
            
            # determine where the circles center is
            x_cent,y_cent = center(x_top,y_top,x_bot,y_bot,r_res,info["ac"])
        
        # determine where circle starts
        theta_t = get_theta(x_cent,y_cent,x_top,y_top)
        theta_b = get_theta(x_cent,y_cent,x_bot,y_bot)

        # if convex desired, shift theta values
        if not info["ac"]:
            theta_b += 2.*np.pi

        # create array of angles and radii
        theta = np.linspace(theta_b,theta_t,info["np"])
        r = np.full(theta.shape,r_res)

        # calculate arc x y coordinates
        # if convex desired, switch theta values
        mult = 1.
        if not info["ac"]: mult *= -1.
        x_arc[0][i] = r * np.cos(theta) + x_cent
        y_arc[0][i] = r * np.sin(theta) + y_cent

        # find where the "thickened" arcs lie
        r_thick = r_res + mult * info["t"]

        # determine where top and bottom of skinned arcs are located
        xs_top,ys_top = get_point(xy,x_cent,y_cent,r_thick,info["ac"])
        xs_bot,ys_bot = get_point(xy,x_cent,y_cent,r_thick,info["ac"],\
            is_top=False)

        # determine where skinned arc starts
        theta_sk_t = get_theta(x_cent,y_cent,xs_top,ys_top)
        theta_sk_b = get_theta(x_cent,y_cent,xs_bot,ys_bot)

        # # if convex desired, shift theta values
        if not info["ac"]:
            theta_sk_b += 2.*np.pi

        # create array of angles and radii
        theta_sk = np.linspace(theta_sk_b,theta_sk_t,info["np"])
        r_thick_arr = np.full(theta_sk.shape,r_thick)

        # calculate arc skin x y coordinates
        x_arc[1][i] = r_thick_arr * np.cos(theta_sk) + x_cent
        y_arc[1][i] = r_thick_arr * np.sin(theta_sk) + y_cent

        # # save top and bottom x y values to the points arrays in info
        # [i][ ][ ] 0 - top, 1 - bottom
        # [ ][i][ ] 0 -> 1 from LE to TE
        # [ ][ ][i] 0 - arc on left (LE), 1 - arc on right (TE)
        info["xpt"][0][i][0] = x_arc[0][i][-1]
        info["xpt"][1][i][0] = x_arc[0][i][0]
        info["xpt"][0][i][1] = x_arc[1][i][-1]
        info["xpt"][1][i][1] = x_arc[1][i][0]
        info["ypt"][0][i][0] = y_arc[0][i][-1]
        info["ypt"][1][i][0] = y_arc[0][i][0]
        info["ypt"][0][i][1] = y_arc[1][i][-1]
        info["ypt"][1][i][1] = y_arc[1][i][0]

    return x_arc,y_arc

# create a cubic bezier curve from two xy arrays over 4 points
def Bezier_Cubic(x,y,n = 30):
    # initialize bezier arrays
    xbz = np.zeros((n,))
    ybz = np.zeros((n,))
    
    # for loop to determine the x and y coordinates at each point
    for i in range(n):
        # determine t value (fraction of number of points)
        t = ( i / (n-1) )
        
        # determine the Bezier values in x and y coordinates
        xbz[i] = (1-t)**3*x[0] + 3*(1-t)**2*t*x[1] + 3*(1-t)*t**2*x[2] + t**3*x[3]
        ybz[i] = (1-t)**3*y[0] + 3*(1-t)**2*t*y[1] + 3*(1-t)*t**2*y[2] + t**3*y[3]
        
    return xbz,ybz

# create outer mouth structure
def make_mouth(O,I,om,info):
    # initialize mouth array
    x_mouth = np.zeros((6,),dtype=np.ndarray)
    y_mouth = np.zeros((6,),dtype=np.ndarray)
    
    # create lower lip line
    # find where the flex starts as an interpolated point outer
    x_fx_o = info["f"]
    y_fx_o,p_fx_o,i_fx_o = get_y_loc(x_fx_o,O["bx"],O["by"])
    # save value so info to be used in a later function
    info["outer"] = {}; info["outer"]["start"] = [x_fx_o,y_fx_o,i_fx_o]

    # find where the flex starts as an interpolated point inner
    x_fx_i = info["f"]
    y_fx_i,p_fx_i,i_fx_i = get_y_loc(x_fx_i,I["bx"],I["by"])

    # create lower lip
    x_mouth[0] = np.array([x_fx_o,x_fx_i])
    y_mouth[0] = np.array([y_fx_o,y_fx_i])

    # create bottom of inner mouth
    x_mouth[1] = np.linspace(info["mo"],info["f"],num=15)
    y_mouth[1] = I["b"](x_mouth[1])
    y_mouth[1][-1] = y_fx_i

    # create back of inner mouth
    # find where the mouth starts as an interpolated point outer mouth
    x_ms_om = info["mo"]
    y_ms_om,p_ms_om,i_ms_om = get_y_loc(x_ms_om,om["bx"],om["by"])
    
    # find where the mouth starts as an interpolated point inner
    x_ms_i = info["mo"]
    y_ms_i,p_ms_i,i_ms_i = get_y_loc(x_ms_i,I["bx"],I["by"])
    y_ms_i = y_mouth[1][0]

    # initialize back outer mouth
    x_mouth[2] = np.array([x_ms_om,x_ms_i])
    y_mouth[2] = np.array([y_ms_om,y_ms_i])

    # create top of inner mouth
    x_mouth[3] = np.linspace(info["mo"],info["f"],num=15)
    y_mouth[3] = om["b"](x_mouth[3])
    # find where flex starts as an interpolated point outer mouth
    x_fx_om = info["f"]
    y_fx_om,p_fx_om,i_fx_om = get_y_loc(x_fx_om,om["bx"],om["by"])
    y_mouth[3][-1] = y_fx_om
    y_mouth[3][0] = y_ms_om
    
    # create "upper lip" line
    # find where flex starts as an interpolated point outer mouth
    x_fx_om = info["f"]
    y_fx_om,p_fx_om,i_fx_om = get_y_loc(x_fx_om,om["bx"],om["by"])
    
    # find where the flex starts as an interpolated point top inner
    x_fx_ti = info["f"]
    y_fx_ti,p_fx_ti,i_fx_ti = get_y_loc(x_fx_ti,I["tx"],I["ty"])

    # create upper lip
    x_mouth[4] = np.array([x_fx_ti,x_fx_om])
    y_mouth[4] = np.array([y_fx_ti - info["fi"],y_fx_om])

    # create fillet from upper lip to inner airfoil
    # find where the fillet ends as an interpolated point top inner
    x_fi_ti = info["f"] + info["fi"]*1.1
    y_fi_ti,p_fi_ti,i_fi_ti = get_y_loc(x_fi_ti,I["tx"],I["ty"])

    # initialze point arrays to get the bezier curve
    x_points = np.array([x_fx_ti,x_fx_ti,(x_fx_ti+x_fi_ti)/2.,x_fi_ti])
    y_points = np.array([y_fx_ti - info["fi"],(y_fx_ti - info["fi"]+y_fi_ti)\
        /2.,(y_fx_ti+y_fi_ti)/2.,y_fi_ti])

    # create fillet using a bezier curve maker
    x_mouth[5],y_mouth[5] = Bezier_Cubic(x_points,y_points)

    return x_mouth,y_mouth

# create tongue tip geometry
def tongue_tip(O,I,ot,it,info):
    # initialize mouth array
    x_tongue = np.zeros((6,),dtype=np.ndarray)
    y_tongue = np.zeros((6,),dtype=np.ndarray)

    # create a bezier curve from tongue tip to the bottom of first arc
    # find interpolation point where flex starts on inner tongue
    x_fx_it = info["f"]
    y_fx_it,p_fx_it,i_fx_it = get_y_loc(x_fx_it,it["bx"],it["by"])

    # create points for the bezier curve
    y_start,_,_ = get_y_loc(             x_fx_it+.03,it["bx"],it["by"])
    y_end,__,__ = get_y_loc(info["ts"]-.03,I["bx"],I["by"])

    # create bezier curve
    xpts = np.array([x_fx_it,x_fx_it+.03,info["ts"]-.03,\
        info["ts"]])
    ypts = np.array([y_fx_it,y_start,y_end,\
        I["b"](info["ts"])])
    x_tongue[0],y_tongue[0] = Bezier_Cubic(xpts,ypts,n=60)

    # create top of tongue tip
    # find where the tongue starts inner tongue
    x_ts_it = info["to"]
    y_ts_it,p_ts_it,i_ts_it = get_y_loc(x_ts_it,it["bx"],it["by"])
    x_tongue[1] = np.linspace(info["to"],info["f"],num=15)
    y_tongue[1] = it["b"](x_tongue[1])
    y_tongue[1][0] = y_ts_it
    y_tongue[1][-1] = y_fx_it

    # create very tip of tongue (flat face towards LE)
    # find where the tongue starts outer tongue
    x_ts_ot = info["to"]
    y_ts_ot,p_ts_ot,i_ts_ot = get_y_loc(x_ts_ot,ot["bx"],ot["by"])

    # save to tongue tip
    x_tongue[2] = np.array([x_ts_it,x_ts_ot])
    y_tongue[2] = np.array([y_ts_it,y_ts_ot])

    # create curve of tongue bottom
    # create an offset curve
    x_tongue[3],y_tongue[3] = offset(x_tongue[0],y_tongue[0],-info["t"])

    # create outer tongue bottom
    x_tongue[4] = np.linspace(info["to"],x_tongue[3][0],num=15)
    y_tongue[4] = ot["b"](x_tongue[4])
    y_tongue[4][ 0] = y_ts_ot
    y_tongue[4][-1] = y_tongue[3][0]

    # create line to fix discontinuity between bezier tongue bottom curve
    # and outer airfoil shape
    # determine the outer airfoil interpolation point from first arc bottom pt
    x_j_o = info["ts"]
    y_j_o,p_j_o,i_j_o = get_y_loc(x_j_o,O["bx"],O["by"])

    # save to tongue tip
    x_tongue[5] = np.array([x_tongue[3][-1],x_j_o])
    y_tongue[5] = np.array([y_tongue[3][-1],y_j_o])
    # save value so info to be used in a later function
    info["outer"]["end"] = [x_j_o,y_j_o,i_j_o]

    return x_tongue,y_tongue

# create a length array and interpolation functions
def length(I,side):
    # create length array
    l_sum = 0
    l = [0]
    # run through all values
    for i in range(1,I[side+"x"].shape[0]):
        l_sum += ((I[side+"x"][i]-I[side+"x"][i-1])**2. + \
            (I[side+"y"][i]-I[side+"y"][i-1])**2.)**0.5
        l.append(l_sum*1.)
    
    # turn l array into numpy
    l = np.array(l)
    # make max on l 1.
    l /= l_sum

    # create interpolation functions
    l_from_x = itp.interp1d(I[side+"x"],l,kind="cubic")
    x_from_l = itp.interp1d(l,I[side+"x"],kind="cubic")
    y_from_l = itp.interp1d(l,I[side+"y"],kind="cubic")

    # create two functions, one to find l from x, and one to find x,y from l
    def get_l(xval):
        return l_from_x(xval)
    def get_xy(lval):
        return x_from_l(lval),y_from_l(lval)
    
    # # test
    # print(get_l(0.5))
    # x1,y1 = get_xy(get_l(0.5))
    # x2,y2 = get_xy(0.5)
    # print(x2,y2,x1,y1)

    return get_l,get_xy

# create a semicircles arrays side
def make_half_semi(I,info,side="t"):
    # determine where to start the arcs and how many semicircles to create
    x_start = info[side+" start x"]

    # determine length from start to end
    lfun,xyfun = length(I,side)

    # determine start value l and end value l
    l_start = lfun(x_start)
    l_end   = lfun(info[side+" end x"])
    if side == "t":
        n = info["n"]
        l_step  = (l_end - l_start) / (n+1)
        info["top l step"] = l_step
    else:
        l_step = info["top l step"]
        n = int((l_end - l_start) / l_step  - 1) 
        l_step  = (l_end - l_start) / (n+1)

    # initialize arrays to return
    x_sem = np.zeros((2*n,),dtype=np.ndarray)
    y_sem = np.zeros((2*n,),dtype=np.ndarray)
    x_ski = np.zeros((2*n+1,),dtype=np.ndarray)
    y_ski = np.zeros((2*n+1,),dtype=np.ndarray)

    # initialize step center value
    lc = l_start * 1.

    # run through each semicircle to be made
    for i in range(n):
        # initialize l value
        lc += l_step

        # determine circle center
        xc,yc = xyfun(lc)

        # determine start skin intersection points, and thus angle
        is_top = side == "t"
        xs,ys = get_point(I,xc,yc,info["smr"],is_concave=False,is_top=is_top)
        theta_s = get_theta(xc,yc,xs,ys)
        if is_top:
            theta_s -= 2*np.pi
        else:
            theta_s += 2*np.pi

        # determine end skin intersection points, and thus angle
        xe,ye = get_point(I,xc,yc,info["smr"],is_concave=True,is_top=is_top)
        theta_e = get_theta(xc,yc,xe,ye)

        # create linspaced theta array and radius array
        theta = np.linspace(theta_s,theta_e,info["np"])
        r     = np.full(theta.shape,info["smr"])
        
        # create x y arrays
        k = int(2*i)
        x_sem[k] = r * np.cos(theta) + xc
        y_sem[k] = r * np.sin(theta) + yc

        # determine indices of start and end locations to create semicircle
        # flat sides
        y_s,_,i_s = get_y_loc(xs,I[side+"x"],I[side+"y"])

        # save index for next arc
        y_e,_,i_e = get_y_loc(xe,I[side+"x"],I[side+"y"])
        
        # create run through
        k += 1
        x_sem[k] = np.zeros((i_e - i_s + 2,))
        y_sem[k] = np.zeros((i_e - i_s + 2,))

        # start with first point
        x_sem[k][0] = x_sem[k-1][-1]
        y_sem[k][0] = y_sem[k-1][-1]

        # run through each additional point
        m = 1
        for j in range(i_s,i_e):
            x_sem[k][m] = I[side+"x"][j]
            y_sem[k][m] = I[side+"y"][j]
            m += 1
        
        # add final point
        x_sem[k][m] = x_sem[k-1][0]
        y_sem[k][m] = y_sem[k-1][0]










        
        # create outer run around
        # determine start skin intersection points, and thus angle
        xks,yks = get_point(I,xc,yc,info["smr"]+info["t"],is_concave=False,\
            is_top=is_top)
        theta_ks = get_theta(xc,yc,xks,yks)
        if is_top:
            theta_ks -= 2*np.pi
        else:
            theta_ks += 2*np.pi

        # determine end skin intersection points, and thus angle
        xke,yke = get_point(I,xc,yc,info["smr"]+info["t"],is_concave=True,\
            is_top=is_top)
        theta_ke = get_theta(xc,yc,xke,yke)

        # create linspaced theta array and radius array
        theta = np.linspace(theta_ks,theta_ke,info["np"])
        r     = np.full(theta.shape,info["smr"]+info["t"])
        
        # create x y arrays
        x_ski[k] = r * np.cos(theta) + xc
        y_ski[k] = r * np.sin(theta) + yc

        # determine indices of start and end locations to run to skin semicircl
        y_ks,_,i_ks = get_y_loc(xks,I[side+"x"],I[side+"y"])

        # initialize end index if first semicircle
        if i == 0:
            x_ke = info[side+" start x"]
            y_ke,_,i_ke = get_y_loc(x_ke,I[side+"x"],I[side+"y"])
            y_ke = info[side+" start y"]
        
        # create run through
        k -= 1
        x_ski[k] = np.zeros((i_ks - i_ke + 2,))
        y_ski[k] = np.zeros((i_ks - i_ke + 2,))

        # start with first point
        x_ski[k][0] = x_ke
        y_ski[k][0] = y_ke

        # run through each additional point
        m = 1
        for j in range(i_ke,i_ks):
            x_ski[k][m] = I[side+"x"][j]
            y_ski[k][m] = I[side+"y"][j]
            m += 1
        
        # add final point
        x_ski[k][m] = xks
        y_ski[k][m] = yks

        # determine end info for next semicircle
        x_ke = xke
        y_ke,_,i_ke = get_y_loc(x_ke,I[side+"x"],I[side+"y"])
        y_ke = yke
        
    # add in final line
    # determine indices of start and end locations to run to skin semicircl
    xks = info[side+" end x"]
    y_ks,_,i_ks = get_y_loc(xks,I[side+"x"],I[side+"y"])
    yks = info[side+" end y"]
    
    # create run through
    k += 2
    x_ski[k] = np.zeros((i_ks - i_ke + 2,))
    y_ski[k] = np.zeros((i_ks - i_ke + 2,))

    # start with first point
    x_ski[k][0] = x_ke
    y_ski[k][0] = y_ke

    # run through each additional point
    m = 1
    for j in range(i_ke,i_ks):
        x_ski[k][m] = I[side+"x"][j]
        y_ski[k][m] = I[side+"y"][j]
        m += 1
    
    # add final point
    x_ski[k][m] = xks
    y_ski[k][m] = yks

    # plot semicircles
    for j in range(2*n):
        plt.plot(x_sem[j],y_sem[j])
    plt.gca().set_prop_cycle(None)
    for j in range(2*n+1):
        plt.plot(x_ski[j],y_ski[j])


    return x_ski,y_ski,x_sem,y_sem

# create semicircles
def make_semicircles(I,info):
    # find semicircles on top surface
    xkt,ykt,xst,yst = make_half_semi(I,info)

    # find semicricles on bottom surface
    xkb,ykb,xsb,ysb = make_half_semi(I,info,"b")

    # initialize new arrays and add the results thereto
    xk = np.zeros((xkt.shape[0]+xkb.shape[0],),dtype=np.ndarray)
    yk = np.zeros((xkt.shape[0]+xkb.shape[0],),dtype=np.ndarray)
    xs = np.zeros((xst.shape[0]+xsb.shape[0],),dtype=np.ndarray)
    ys = np.zeros((xst.shape[0]+xsb.shape[0],),dtype=np.ndarray)

    # add in skin arrays
    xk[:xkt.shape[0]] = xkt
    xk[xkt.shape[0]:] = xkb
    yk[:xkt.shape[0]] = ykt
    yk[xkt.shape[0]:] = ykb

    # add in semicircles arrays
    xs[:xst.shape[0]] = xst
    xs[xst.shape[0]:] = xsb
    ys[:xst.shape[0]] = yst
    ys[xst.shape[0]:] = ysb

    return xk,yk,xs,ys

# create roofs and floors
def roofs_n_floors(I,info,x_fillet,y_tri_top,y_tri_bot):
    # initialize mouth array
    # [i][ ] 0 - roof, 1 - floor
    # [ ][i] 0 -> 1, LE -> TE
    # NOTE [1][0] will remain empty
    x_rofl = np.zeros((2,info["n"]+1),dtype=np.ndarray)
    y_rofl = np.zeros((2,info["n"]+1),dtype=np.ndarray)

    # cycle between roof and floor
    for i in range(x_rofl.shape[0]):
        # cycle through each arc and determine the roof and or floor
        for j in range(x_rofl.shape[1]):
            # if solving the roof
            if i == 0:
                # determine start
                if j == 0:
                    # determine the xy location of the fillet end
                    x_f_i = x_fillet
                    y_f_i,p_f_i,i_f_i = get_y_loc(x_f_i,I["tx"],I["ty"])
                    xstart = x_fillet
                    ystart = y_f_i
                    istart = i_f_i
                else:
                    # determin xy location and index of arc begin
                    x_f_i = info["xpt"][0][j-1][1]
                    y_f_i,p_f_i,i_f_i = get_y_loc(x_f_i,I["tx"],I["ty"])
                    xstart = x_f_i
                    ystart = info["ypt"][0][j-1][1]
                    istart = i_f_i
                
                # determine end
                if j == info["n"]:
                    # determine the xy location of the fillet end
                    x_f_i = info["tw"] - info["t"]
                    y_f_i,p_f_i,i_f_i = get_y_loc(x_f_i,I["tx"],I["ty"])
                    xend = x_f_i
                    yend = y_tri_top
                    iend = i_f_i
                else:
                    # determin xy location and index of arc begin
                    x_f_i = info["xpt"][0][j][0]
                    y_f_i,p_f_i,i_f_i = get_y_loc(x_f_i,I["tx"],I["ty"])
                    xend = x_f_i
                    yend = info["ypt"][0][j][0]
                    iend = i_f_i
                
                # solve for points
                x_rofl[i][j] = np.linspace(xstart,xend,num=10)
                y_rofl[i][j] = I["t"](x_rofl[i][j])
                
                # add the start y
                y_rofl[i][j][0] = ystart
                # add end y
                y_rofl[i][j][-1] = yend

            # if solving for the floors (excluding first)
            elif j >= 1:
                # determine start
                # determin xy location and index of arc begin
                x_f_i = info["xpt"][1][j-1][1]
                y_f_i,p_f_i,i_f_i = get_y_loc(x_f_i,I["bx"],I["by"])
                xstart = x_f_i
                ystart = info["ypt"][1][j-1][1]
                istart = i_f_i
                
                # determine end
                if j == info["n"]:
                    # determine the xy location of the fillet end
                    x_f_i = info["tw"] - info["t"]
                    y_f_i,p_f_i,i_f_i = get_y_loc(x_f_i,I["bx"],I["by"])
                    xend = x_f_i
                    yend = y_tri_bot
                    iend = i_f_i
                else:
                    # determin xy location and index of arc begin
                    x_f_i = info["xpt"][1][j][0]
                    y_f_i,p_f_i,i_f_i = get_y_loc(x_f_i,I["bx"],I["by"])
                    xend = x_f_i
                    yend = info["ypt"][1][j][0]
                    iend = i_f_i
                
                # solve for points
                x_rofl[i][j] = np.linspace(xstart,xend,num=10)
                y_rofl[i][j] = I["b"](x_rofl[i][j])
                
                # add the start y
                y_rofl[i][j][0] = ystart
                # add end y
                y_rofl[i][j][-1] = yend

    return x_rofl,y_rofl

# determine where a linear line crosses an interpolation line
def crossing(xy,m,b,is_top=True):
    # solve for if top or bottom
    if is_top:
        loc = "t"
    else:
        loc = "b"
    
    # create function to find where the crossing occurs
    def func(x):
        return m*x+b - xy[loc](x)

    # set min and max values for bisect
    minx = np.min(xy[loc+"x"])
    maxx = np.max(xy[loc+"x"])
    
    # run bisect function to find point
    x_cross = Bisect(func,minx,maxx)

    # solve for the y value and the index of the crossing
    y_cross = xy[loc](x_cross)
    _,_,i_cross = get_y_loc(x_cross,xy[loc+"x"],xy[loc+"y"])

    return x_cross,y_cross,i_cross

# create wing box
def wing_box(I,im,info):
    # initialize wing box info
    x_box = np.zeros((26,),dtype=np.ndarray)
    y_box = np.zeros((26,),dtype=np.ndarray)

    ### create TE triangle of wingbox
    # create diagonal
    # find top right of triangle
    x_trTE = info["f"] - info["t"]/2.
    y_trTE,p_trTE,i_trTE = get_y_loc(x_trTE,I["tx"],I["ty"])
    # find bottom left of triangle
    x_blTE = info["mo"] - info["t"]/2.
    y_blTE,p_blTE,i_blTE = get_y_loc(x_blTE,im["bx"],im["by"])

    # find offset of this line
    x_tri = np.array([x_trTE,x_blTE]); y_tri = np.array([y_trTE,y_blTE])
    # slope of line
    m = (y_tri[0] - y_tri[1]) / (x_tri[0] - x_tri[1])
    # x value of shift
    xvalue = ((info["t"]/2.)**2./(1+(-1./m)**2.))**0.5
    # offset line
    offxTE = x_tri + xvalue
    offyTE = y_tri + xvalue * -1./m
    # determine offset line intercept
    bTE = offyTE[0] - offxTE[0]*m

    # fix end point to meet with flex - skin
    offxTE[0] = info["f"] - info["t"]
    offyTE[0] = offxTE[0] * m + bTE

    # fix start point to meet with inner mouth interpolation, save
    x_cross_im,y_cross_im,i_cross_im = crossing(im,m,bTE,is_top=False)
    offxTE[1] = x_cross_im; offyTE[1] = y_cross_im
    x_box[0] = offxTE; y_box[0] = offyTE

    # create bottom piece
    x_box[1] = np.linspace(x_cross_im,info["f"]-info["t"],num=10)
    y_box[1] = im["b"](x_box[1])
    # determine index where meet up with flex - skin
    x_TE_end = info["f"] - info["t"]
    y_TE_end,p_TE_end,i_TE_end = get_y_loc(x_TE_end,im["bx"],im["by"])
    y_box[1][-1] = y_TE_end
    
    # create side piece, save to array, but first:
    x_box[2] = np.array([offxTE[0],x_TE_end])
    y_box[2] = np.array([offyTE[0],y_TE_end])

    ### create LE triangle of wingbox
    # offset lines
    offxLE = x_tri - xvalue
    offyLE = y_tri - xvalue * -1./m
    # determine offset line intercept
    bLE = offyLE[0] - offxLE[0]*m

    # fix start point to meet with flex - skin
    offxLE[1] = info["mo"]
    offyLE[1] = offxLE[1] * m + bLE

    # fix end point to meet with inner mouth interpolation, save
    x_cross_i,y_cross_i,i_cross_i = crossing(I,m,bLE)
    offxLE[0] = x_cross_i; offyLE[0] = y_cross_i
    x_box[3] = offxLE; y_box[3] = offyLE

    # create top piece
    x_box[4] = np.linspace(info["mo"],x_cross_i,num=10)
    y_box[4] = I["t"](x_box[4])
    # determine index where meet up with flex - skin
    x_LE_str = info["mo"]
    y_LE_str,p_LE_str,i_LE_str = get_y_loc(x_LE_str,I["tx"],I["ty"])
    y_box[4][0] = y_LE_str

    # create side piece, save to array, but first:
    x_box[5] = np.array([x_LE_str,offxLE[1]])
    y_box[5] = np.array([y_LE_str,offyLE[1]])

    ### create LE mini triangle
    # find start point meeting bottom inner
    x_bot_i,y_bot_i,i_bot_i = crossing(I,m,bTE,is_top=False)

    # find end point, save
    x_out = info["mo"] - info["t"]
    y_out = m * x_out + bTE
    x_box[6] = np.array([x_bot_i,x_out])
    y_box[6] = np.array([y_bot_i,y_out])

    # create bottom piece
    x_box[7] = np.linspace(x_bot_i,info["mo"]-info["t"],num=4)
    y_box[7] = I["b"](x_box[7])

    # create side piece, save to array, but first:
    # determine index where meet up with flex - skin
    x_TE_mini = info["mo"] - info["t"]
    y_TE_mini,p_TE_mini,i_TE_mini = get_y_loc(x_TE_mini,I["bx"],I["by"])
    x_box[8] = np.array([x_out,x_TE_mini])
    y_box[8] = np.array([y_out,y_TE_mini])

    ### create outer two lines
    # create vertical line
    # find top point
    x_top_vt = info["mo"] - info["t"]
    y_top_vt,p_top_vt,i_top_vt = get_y_loc(x_top_vt,I["tx"],I["ty"])

    # find bottom point, save
    x_bot_vt = info["mo"] - info["t"]
    y_bot_vt = m*x_bot_vt + bLE
    x_box[9] = np.array([x_top_vt,x_bot_vt])
    y_box[9] = np.array([y_top_vt,y_bot_vt])

    # create diagonal line
    x_bot_dia,y_bot_dia,i_bot_dia = crossing(I,m,bLE,is_top=False)

    # save
    x_box[10] = np.array([x_bot_vt,x_bot_dia])
    y_box[10] = np.array([y_bot_vt,y_bot_dia])

    ### create inner airfoil runaround
    # set min points to prevent running outside interpolation ranges
    t_min = I["tx"][10]
    b_min = I["bx"][10]

    # create top array that runs to t_min value, cosine clustered
    theta = np.linspace(0.0,np.pi/2.,30)
    txvals = np.flip((x_top_vt-t_min)*-np.cos(theta) + (x_top_vt))
    tyvals = I["t"](txvals)

    # create an array of the "skipped" values due to interpolation range
    bfxvals = np.append(np.flip(I["tx"][:11]),I["bx"][:11])
    bfyvals = np.append(np.flip(I["ty"][:11]),I["by"][:11])

    # create bottom array that runs from b_min value, cosine clustered
    theta = np.linspace(0.0,np.pi/2.,15)
    bxvals = (x_bot_dia-b_min)*(-np.cos(theta)) + (x_bot_dia)
    byvals = I["b"](bxvals)

    # create a full appended array
    runx = np.append(np.append(txvals,bfxvals),bxvals)
    runy = np.append(np.append(tyvals,bfyvals),byvals)
    runy[0] = y_top_vt

    # save these values to the box arrays
    x_box[11] = runx[:5];      y_box[11] = runy[:5]
    x_box[12] = runx[4:10];    y_box[12] = runy[4:10]
    x_box[13] = runx[9:12];    y_box[13] = runy[9:12]
    x_box[14] = runx[11:14];   y_box[14] = runy[11:14]
    x_box[15] = runx[13:15];   y_box[15] = runy[13:15]
    x_box[16] = runx[14:20];   y_box[16] = runy[14:20]
    x_box[17] = runx[19:25];   y_box[17] = runy[19:25]
    x_box[18] = runx[24:30];   y_box[18] = runy[24:30]
    x_box[19] = runx[29:35];   y_box[19] = runy[29:35]
    x_box[20] = runx[34:40];   y_box[20] = runy[34:40]
    x_box[21] = runx[39:45];   y_box[21] = runy[39:45]
    x_box[22] = runx[44:50];   y_box[22] = runy[44:50]
    x_box[23] = runx[49:55];   y_box[23] = runy[49:55]
    x_box[24] = runx[54:60];   y_box[24] = runy[54:60]
    x_box[25] = runx[59:];     y_box[25] = runy[59:]

    return x_box,y_box

# create outer airfoil shape
def outer(O,info):
    # initialize outer arrays
    x_outer = np.zeros((14,),dtype=np.ndarray)
    y_outer = np.zeros((14,),dtype=np.ndarray)

    ### create outer airfoil top runaround
    x_outer[0]  = O["x af"][:55];     y_outer[0]  = O["y af"][:55]
    x_outer[1]  = O["x af"][54:65];   y_outer[1]  = O["y af"][54:65]
    x_outer[2]  = O["x af"][64:70];   y_outer[2]  = O["y af"][64:70]
    x_outer[3]  = O["x af"][69:75];   y_outer[3]  = O["y af"][69:75]
    x_outer[4]  = O["x af"][74:80];   y_outer[4]  = O["y af"][74:80]
    x_outer[5]  = O["x af"][79:85];   y_outer[5]  = O["y af"][79:85]
    x_outer[6]  = O["x af"][84:90];   y_outer[6]  = O["y af"][84:90]
    x_outer[7]  = O["x af"][89:95];   y_outer[7]  = O["y af"][89:95]
    x_outer[8]  = O["x af"][94:100];  y_outer[8]  = O["y af"][94:100]
    x_outer[9]  = O["x af"][99:105];  y_outer[9]  = O["y af"][99:105]
    x_outer[10] = O["x af"][104:125]; y_outer[10] = O["y af"][104:125]

    # set min points to prevent running outside interpolation ranges
    b_min = O["x af"][124]; b_max = info["outer"]["start"][0]

    # create bottom array that runs from b_min value
    x_outer[-3] = np.linspace(b_min,b_max,25)
    y_outer[-3] = O["b"](x_outer[-3])
    y_outer[-3][-1] = info["outer"]["start"][1]

    # initialize second run
    # set min points to prevent running outside interpolation ranges
    b_min = info["outer"]["end"][0]; b_max = np.max(O["bx"])

    # create bottom array that runs from b_min value, cosine clustered
    theta = np.linspace(0.0,np.pi/2.,50)
    bxvals = (b_max-b_min)*(-np.cos(theta)) + (b_max)
    byvals = O["b"](bxvals)

    # create a full appended array
    x_outer[-2] = bxvals[0:25]
    y_outer[-2] = byvals[0:25]
    y_outer[-2][0] = info["outer"]["end"][1]
    x_outer[-1] = bxvals[24:]
    y_outer[-1] = byvals[24:]
    
    # connect TE if necessary
    top_max = x_outer[0][0]; bot_max = x_outer[-1][-1]
    if not (top_max == bot_max):
        if top_max > bot_max:
            x_outer[-1] = np.append(x_outer[-1],x_outer[0][0])
            y_outer[-1] = np.append(y_outer[-1],y_outer[0][0])
        else:
            x_outer[0] = np.insert(x_outer[0],0,x_outer[-1][-1])
            y_outer[0] = np.insert(y_outer[0],0,y_outer[-1][-1])

    return x_outer,y_outer

# create dxf file
def make_dxf(x,y,filename,write_to_file=True):
    # determine number of arrays spots to make
    size = 0
    for group in x:
        shape = x[group].shape
        if group == "arc" or group == "rfl":
            for i in range(shape[0]):
                for j in range(shape[1]):
                    if type(x[group][i][j]) == np.ndarray:
                        size += 1
        else:
            for i in range(shape[0]):
                if type(x[group][i]) == np.ndarray:
                    size += 1

    # initialize dxf inputs
    x_dxf = np.zeros((size,),dtype=np.ndarray)
    y_dxf = np.zeros((size,),dtype=np.ndarray)
    index = 0

    # input dxf points
    for group in x:
        shape = x[group].shape
        if group == "arc" or group == "rfl":
            for i in range(shape[0]):
                for j in range(shape[1]):
                    if type(x[group][i][j]) == np.ndarray:
                        x_dxf[index] = x[group][i][j]
                        y_dxf[index] = y[group][i][j]
                        index += 1
        else:
            for i in range(shape[0]):
                if type(x[group][i]) == np.ndarray:
                    x_dxf[index] = x[group][i]
                    y_dxf[index] = y[group][i]
                    index += 1
    
    # create z arrays
    z_dxf = y_dxf * 0.

    # write dxf file
    if write_to_file:
        dxf.dxf(filename,x_dxf,y_dxf,z_dxf,geometry="spline")

    return x_dxf,y_dxf,z_dxf

# deflect airfoil function
def deflection(x,y,O,info):
    # run through each group
    for group in x:
        if group == "arc" or group == "rfl":
            for i in range(x["arc"].shape[0]):
                x[group][i],y[group][i] = deflect.main(O,x[group][i],\
                    y[group][i],info)
        else:
            x[group],y[group] = deflect.main(O,x[group],y[group],info)
    
    return

# shift to quarter chord
def shift(x,y,O):
    # determine shift values
    x_shift = 0.25
    y_shift = O["c"](x_shift)

    # shift each group/index
    if type(x) == dict:
        for group in x:
            x[group] -= x_shift
            y[group] -= y_shift
    elif type(x) == np.ndarray:
        x -= x_shift
        y -= y_shift
    else:
        x -= x_shift
        y -= y_shift
    
    return

# resize by a chord value
def resize(x,y,c):
    # run through each group and resize
    if type(x) == dict:
        for group in x:
            x[group] *= c
            y[group] *= c
    elif type(x) == np.ndarray:
        x *= c
        y *= c
    else:
        x *= c
        y *= c
    
    return

# split off outer shape for ease in modeling
def split(x,y):
    # determine how many spots to place in sketch
    num = x["ton"].shape[0] + x["mou"].shape[0] + x["out"].shape[0] + 2

    # initialize x and y arrays
    xsplit = np.zeros((num,),dtype=np.ndarray)
    ysplit = np.zeros((num,),dtype=np.ndarray)
    j = 0

    # pass in tongue
    for i in range(x["ton"].shape[0]):
        # pass in values
        xsplit[j] = x["ton"][i]
        ysplit[j] = y["ton"][i]
        # increase counter
        j += 1

    # pass in mouth
    for i in range(x["mou"].shape[0]):
        # pass in values
        xsplit[j] = x["mou"][i]
        ysplit[j] = y["mou"][i]
        # increase counter
        j += 1

    # pass in outer
    for i in range(x["out"].shape[0]):
        # pass in values
        xsplit[j] = x["out"][i]
        ysplit[j] = y["out"][i]
        # increase counter
        j += 1

    # pass in roof and arc to complete shape
    xsplit[j] = x["rfl"][0][0]
    ysplit[j] = y["rfl"][0][0]
    xsplit[j+1] = x["arc"][0][0]
    ysplit[j+1] = y["arc"][0][0]

    # set z array as zero
    zsplit = ysplit * 0.

    return xsplit,ysplit,zsplit

# determine the max and min of each "hole"
def guide_curves(x,y,I,vals,info):
    # initialize points matrix.
    # i 0 --> holes, with 0 = outer skin
    # j guide curve number
    # k 0 - x, 1 - y, 2 - z
    number_to_make = 6 + info["n"] - (info["n"] + 4) * info["sk"]
    guides = np.zeros((number_to_make,),dtype=np.ndarray)
    xpts = np.zeros((number_to_make,),dtype=np.ndarray)
    ypts = np.zeros((number_to_make,),dtype=np.ndarray)

    # initialize group and point tuples lists
    grp = [["mou","mou","mou","mou","mou","mou","arc","arc","ton","ton","out",\
        "out"]]
    point = [
        [(0,-1),(1,0),(2,0),(3,-1),\
            (5,0),(5,-1),((0,0),-1),((0,0),0),(2,0),(2,-1),(-2,-1),(-1,-1)]]

    if not info["sk"]:
        grp = grp + [["box","box"],
            ["box","box","box"],
            ["box","box","box"],
            ["box","box","box"]
        ]
        point = point + [[(9,-1),(9,0)],
            [(7,0),(7,-1),(8,0)],
            [(5,-1),(4,-1),(4,0)],
            [(1,0),(1,-1),(2,0)]
        ]

    # initialize number of inner box points to add for hole 0
    num_outer = x["out"].shape[0] - 1

    # run through each number and add that many to the to be added lists
    grp0 = []; point0 = []
    for i in range(num_outer):
        grp0.append("out")
        point0.append((i,-1))
    
    # add to 2d lists
    grp[0] = grp0 + grp[0]; point[0] = point0 + point[0]

    if not info["sk"]:
        # initialize number of outer points to add for hole 1
        num_box = x["box"].shape[0] - 11

        # run through each number and add that many to the to be added lists
        grp1 = []; point1 = []
        for i in range(11,11+num_box):
            grp1.append("box")
            point1.append((i,-1))
        
        # add to 2d lists
        grp[1] = grp1 + grp[1]; point[1] = point1 + point[1]

        # add a new index to the grp and point lists for each arc (minus 1)
        j = 0
        for i in range(info["n"]-1):
            grp.append(["arc","arc","arc","arc"])
            point.append([((1,i),0),((0,i+1),0),((0,i+1),-1),((1,i),-1)])
            j += 1
        
        # add for final arc and te triangle
        grp.append(["arc","tri","tri","arc"])
        point.append([((1,j),0),(3,-1),(3,0),((1,j),-1)])
        grp.append(["tri","tri","tri"])
        point.append([(1,-1),(0,0),(0,-1)])

    # run through each hole
    for i in range(len(grp)):
        # create x y arrays
        x_hole = np.zeros((len(grp[i]),))
        y_hole = np.zeros((len(grp[i]),))

        # run through each guide curve
        for j in range(len(grp[i])):
            # set spline value
            spline = point[i][j][0]
            # set point value
            index = point[i][j][1]

            # save to x and y arrays
            x_hole[j] = x[grp[i][j]][spline][index]
            y_hole[j] = y[grp[i][j]][spline][index]

        # and subsequently save to arrays
        xpts[i] = x_hole
        ypts[i] = y_hole

    # if shift origin to quarter chord
    if "shift to c/4" in vals["C14"] and vals["C14"]["shift to c/4"]:
        # shift
        shift(xpts,ypts,I["O"])

    # resize chord of each segment
    resize(xpts,ypts,info["c"])

    # add points back to guides array
    q = 0
    for i in range(len(grp)):
        # initialize guides array
        guides[i] = np.zeros((len(grp[i]),3))

        # save x and y arrays to guides array
        guides[i][:,0] = xpts[i]
        guides[i][:,1] = ypts[i]
        q += 1
    
    # remove later
    for i in range(q,guides.shape[0]):
        guides[i] = np.zeros((1,3))
    
    return guides

# determine the split guide curves and dxf files 
def split_part(x,y,I,vals,info):
    # create result dictionary
    results = {}

    # create dxf file of inner hole for loft cut
    # initialize inner dxf file
    num_inner = 4*info["n"] + 2
    x_inner = np.zeros((num_inner,),dtype=np.ndarray)
    y_inner = np.zeros((num_inner,),dtype=np.ndarray)

    # run through roofs and floors and add the required lines
    k = 0 # counter for inner arrays
    for i in range(x["rfl"].shape[0]):
        # run through each roof or floor
        for j in range(1,x["rfl"].shape[1]):
            # determine arc index
            if i == 0:
                ai = -1; an = j
            else: ai = 0; an = 6 - j

            # add through arc line or roof/floor
            if i == 0:
                x_inner[k] = np.array([x["arc"][0][an-1][ai],\
                    x["arc"][1][an-1][ai]])
                y_inner[k] = np.array([y["arc"][0][an-1][ai],\
                    y["arc"][1][an-1][ai]])
            else:
                x_inner[k] = np.flip(x["rfl"][i][an])
                y_inner[k] = np.flip(y["rfl"][i][an])
            k += 1

            # add roof/floor or arc
            if i == 0:
                x_inner[k] = x["rfl"][i][an]
                y_inner[k] = y["rfl"][i][an]
            else:
                x_inner[k] = np.array([x["arc"][1][an-1][ai],\
                    x["arc"][0][an-1][ai]])
                y_inner[k] = np.array([y["arc"][1][an-1][ai],\
                    y["arc"][0][an-1][ai]])
            k += 1
        
        # add te wall
        if i == 0:
            x_inner[k] = x["tri"][3]
            y_inner[k] = y["tri"][3]
            k += 1
        else:
            x_inner[k] = np.array([x["arc"][0][0][0],x["arc"][0][0][-1]])
            y_inner[k] = np.array([y["arc"][0][0][0],y["arc"][0][0][-1]])
            k += 1  

    # save to results dictionary
    results["x inner"] = x_inner
    results["y inner"] = y_inner
    
    # create dxf file of te hole for loft cut
    # initialize te dxf file
    x_te = np.zeros((3,),dtype=np.ndarray)
    y_te = np.zeros((3,),dtype=np.ndarray)

    # add side of triangle
    x_te[0] = np.array([x["tri"][1][0],x["tri"][1][1]])
    y_te[0] = np.array([y["tri"][1][0]+3*info["t"],y["tri"][1][1]-3*info["t"]])

    # add top part
    x_te[1] = np.array([x["tri"][1][0],x["out"][0][0]+3*info["t"]])
    y_te[1] = np.array([y["tri"][1][0]+3*info["t"],y["out"][0][0]])

    # add bottom part
    x_te[2] = np.array([x["tri"][1][1],x["out"][0][0]+3*info["t"]])
    y_te[2] = np.array([y["tri"][1][1]-3*info["t"],y["out"][0][0]])

    # save to results dictionary
    results["x te"] = x_te
    results["y te"] = y_te

    # create guide curves for the inner, [0] index for all 
    # initialize guide curves for inner hole
    x_gc_in = np.zeros((num_inner,))
    y_gc_in = np.zeros((num_inner,))

    # run through the inner shape and find the points at each [0] index
    j = 0
    for i in range(results["x inner"].shape[0]):
        x_gc_in[j] = results["x inner"][i][0]
        y_gc_in[j] = results["y inner"][i][0]
        j += 1
    
    # save to results dictionary
    results["x inn gc"] = x_gc_in
    results["y inn gc"] = y_gc_in

    # create dxf file of arcs for lofting thereof
    # initialize dxf arrays
    x_arcs = np.zeros((4*info["n"],),dtype=np.ndarray)
    y_arcs = np.zeros((4*info["n"],),dtype=np.ndarray)

    # run through each arc and add it to the shape
    k = 0
    for i in range(info["n"]):
        # add front arc piece
        x_arcs[k] = x["arc"][0][i]
        y_arcs[k] = y["arc"][0][i]
        k += 1

        # add top arc piece
        x_arcs[k] = np.array([x["arc"][0][i][-1],x["arc"][1][i][-1]])
        y_arcs[k] = np.array([y["arc"][0][i][-1],y["arc"][1][i][-1]])
        k += 1

        # add back arc piece
        x_arcs[k] = np.flip(x["arc"][1][i])
        y_arcs[k] = np.flip(y["arc"][1][i])
        k += 1

        # add bottom arc piece
        x_arcs[k] = np.array([x["arc"][1][i][0],x["arc"][0][i][0]])
        y_arcs[k] = np.array([y["arc"][1][i][0],y["arc"][0][i][0]])
        k += 1
    
    # save to results dictionary
    results["x arcs"] = x_arcs
    results["y arcs"] = y_arcs

    # create guide curves for the arcs, [0] index for all 
    # initialize guide curves for inner hole
    x_gc_arc = np.zeros((info["n"],4))
    y_gc_arc = np.zeros((info["n"],4))

    # run through the inner shape and find the points at each [0] index
    k = 0
    for i in range(info["n"]):
        for j in range(4):
            x_gc_arc[i,j] = results["x arcs"][k][0]
            y_gc_arc[i,j] = results["y arcs"][k][0]
            k += 1
    
    # save to results dictionary
    results["x arcs gc"] = x_gc_arc
    results["y arcs gc"] = y_gc_arc

    # save tongue tip midpoint
    xton = (x["ton"][2][0]+x["ton"][2][1])/2.
    yton = (y["ton"][2][0]+y["ton"][2][1])/2.

    # if shift origin to quarter chord
    if "shift to c/4" in vals["C14"] and vals["C14"]["shift to c/4"]:
        # shift
        shift(results["x inner"],results["y inner"],I["O"])
        shift(results["x te"],results["y te"],I["O"])
        shift(results["x inn gc"],results["y inn gc"],I["O"])
        shift(results["x arcs"],results["y arcs"],I["O"])
        shift(results["x arcs gc"],results["y arcs gc"],I["O"])
        shift(xton,yton,I["O"])

    # resize chord of each segment
    resize(results["x inner"],results["y inner"],info["c"])
    resize(results["x te"],results["y te"],info["c"])
    resize(results["x inn gc"],results["y inn gc"],info["c"])
    resize(results["x arcs"],results["y arcs"],info["c"])
    resize(results["x arcs gc"],results["y arcs gc"],info["c"])
    resize(xton,yton,info["c"])

    # create inn guide curve dictionary element
    results["inn gc"] = np.zeros((num_inner,3))
    results["inn gc"][:,0] = results["x inn gc"]
    results["inn gc"][:,1] = results["y inn gc"]

    # create arcs guide curve dictionary element
    results["arcs gc"] = np.zeros((info["n"],4,3))
    results["arcs gc"][:,:,0] = results["x arcs gc"]
    results["arcs gc"][:,:,1] = results["y arcs gc"]

    results["ton tip"] = [xton,yton,0.0]

    return results

# create actuation hole
def actuate_hole(x,y,I,vals,info,results):
    # create actuation hole # initialize array
    num = 10
    x_hole = np.zeros((num,),dtype=np.ndarray)
    y_hole = np.zeros((num,),dtype=np.ndarray)
    j = 0

    # create empty arrays for each 
    for i in range(x_hole.shape[0]):
        x_hole[i] = np.array([])
        y_hole[i] = np.array([])

    # start at beginning of fillet at top of mouth
    x_hole[j] = np.append(x_hole[j],x["mou"][-1][0])
    y_hole[j] = np.append(y_hole[j],y["mou"][-1][0])

    # add top inner interp point
    x_start = x["mou"][-1][0]
    y_start,_,i_start = get_y_loc(x_start,I["I"]["tx"],I["I"]["ty"])
    x_hole[j] = np.append(x_hole[j],x_start)
    y_hole[j] = np.append(y_hole[j],y_start)
    j += 1

    # create start
    x_hole[j] = np.append(x_hole[j],x_start)
    y_hole[j] = np.append(y_hole[j],y_start)

    # add top front inner portion
    for i in range(i_start-1,-1,-1):
        x_hole[j] = np.append(x_hole[j],I["I"]["tx"][i])
        y_hole[j] = np.append(y_hole[j],I["I"]["ty"][i])

    # split between two more arrays
    x_hole[j+1] = x_hole[j][int(i_start/3)-1:int(2*i_start/3)]
    y_hole[j+1] = y_hole[j][int(i_start/3)-1:int(2*i_start/3)]
    x_hole[j+2] = x_hole[j][int(2*i_start/3)-1:]
    y_hole[j+2] = y_hole[j][int(2*i_start/3)-1:]
    x_hole[j] = x_hole[j][:int(i_start/3)]
    y_hole[j] = y_hole[j][:int(i_start/3)]
    j += 3

    # determine end
    x_end = x["mou"][2][0]
    y_end,_,i_end = get_y_loc(x_end,I["I"]["bx"],I["I"]["by"])    
    # add bottom front portion
    for i in range(i_end):
        x_hole[j] = np.append(x_hole[j],I["I"]["bx"][i])
        y_hole[j] = np.append(y_hole[j],I["I"]["by"][i])
    
    # add final point
    x_hole[j] = np.append(x_hole[j],x["mou"][1][0])
    y_hole[j] = np.append(y_hole[j],y["mou"][1][0])

    # split between two more arrays
    x_hole[j+1] = x_hole[j][int(i_end/3)-1:int(2*i_end/3)]
    y_hole[j+1] = y_hole[j][int(i_end/3)-1:int(2*i_end/3)]
    x_hole[j+2] = x_hole[j][int(2*i_end/3)-1:]
    y_hole[j+2] = y_hole[j][int(2*i_end/3)-1:]
    x_hole[j] = x_hole[j][:int(i_end/3)]
    y_hole[j] = y_hole[j][:int(i_end/3)]
    j += 3

    # append to array
    x_hole[j] = np.append(x_hole[j],x_end)
    y_hole[j] = np.append(y_hole[j],y_end)
    j += 1

    # create x array to be used later
    x_toptongue = np.linspace(info["to"],info["f"]+.01)
    y_toptongue = (I["om"]["b"](x_toptongue) + I["it"]["b"](x_toptongue)) / 2.

    # create array
    x_hole[j] = np.append(x_hole[j],x_toptongue)
    y_hole[j] = np.append(y_hole[j],y_toptongue)

    # append to array
    x_hole[j-1] = np.append(x_hole[j-1],x_toptongue[0])
    y_hole[j-1] = np.append(y_hole[j-1],y_toptongue[0])
    j += 1

    # create array
    x_hole[j] = np.append(x_hole[j],x_hole[j-1][-1])
    y_hole[j] = np.append(y_hole[j],y_hole[j-1][-1])
    x_hole[j] = np.append(x_hole[j],x_hole[0][0])
    y_hole[j] = np.append(y_hole[j],y_hole[0][0])

    # add to results
    results["x hole"] = x_hole; results["y hole"] = y_hole

    # create guide curves for the hole, [0] index for all 
    # initialize guide curves for inner hole
    x_gc_hl = np.zeros((num,))
    y_gc_hl = np.zeros((num,))

    # run through the inner shape and find the points at each [0] index
    for i in range(results["x hole"].shape[0]):
        x_gc_hl[i] = results["x hole"][i][0]
        y_gc_hl[i] = results["y hole"][i][0]
    
    # save to results dictionary
    results["x hole gc"] = x_gc_hl
    results["y hole gc"] = y_gc_hl

    # if shift origin to quarter chord
    if "shift to c/4" in vals["C14"] and vals["C14"]["shift to c/4"]:
        # shift
        shift(results["x hole"],results["y hole"],I["O"])
        shift(results["x hole gc"],results["y hole gc"],I["O"])

    # resize chord of each segment
    resize(results["x hole"],results["y hole"],info["c"])
    resize(results["x hole gc"],results["y hole gc"],info["c"])

    # create inn guide curve dictionary element
    results["hole gc"] = np.zeros((num,3))
    results["hole gc"][:,0] = results["x hole gc"]
    results["hole gc"][:,1] = results["y hole gc"]

    return

# main function
def main(jsonfile):
    """
    This is an explanation of what this file does
    """
    # import json file if a dictionary or json file, OR raise error
    if type(jsonfile) == dict:
        vals = jsonfile
    elif os.path.isfile(jsonfile):
        vals = read_json(jsonfile)
    else:
        raise ValueError("C14 input must be dictionary or file path")

    # simplify dictionary
    info = simplify(vals)

    # create airfoil shape
    add = {}
    add["geometry"] = {}
    add["geometry"]["type"] = "outline_points"
    afname = vals["C14"]["airfoil dat file location"]
    if afname[:4] == "NACA":
        add["geometry"]["NACA"] = afname.split(" ")[-1]
        isNACA = True
    else:
        add["geometry"]["outline_points"] = afname
        isNACA = False

    # initialize airfoil database
    foil = adb.Airfoil("foil",add)
    coords = foil.get_outline_points(close_te=not(isNACA),top_first=True,N=200)
    x_af = coords[:,0]; y_af = coords[:,1]

    # # close te
    # if x_af[-1] >= x_af[0]:
    #     x_af[0] = x_af[-1] * 1.
    # else:
    #     x_af[-1] = x_af[0] * 1.

    # enthicken the airfoil
    y_af *= vals["C14"]["enthicken"]

    # create inner airfoil
    polyx,polyy = inner(x_af,y_af,info["t"])

    # initialize interpolation dictionary
    I = {}

    # create interpolation of inner surface
    I["I"] = interpolate(polyx,polyy,make_camberline=True)

    # initialize x and y dictionaries
    x = {}; y = {}

    # find te triangle
    x["tri"],y["tri"] = te_triangle(I["I"],info)

    # # determine where the arcs should start-ish
    # I["I"]["arcdx"] = (info["tw"] - info["t"] - info["f"]) / (info["n"] + 1)

    # # find where each arc is
    # x["arc"],y["arc"] = arcs(I["I"],info)

    # create interpolation of outer surface
    I["O"] = interpolate(x_af,y_af,make_camberline=True)
    # save airfoil outline for future use
    I["O"]["x af"] = x_af; I["O"]["y af"] = y_af

    # create further thicknesses inner airfoils, as well as interpolations
    # outer tongue
    out_tx,out_ty = inner(x_af,y_af,   info["t"]+   info["mt"])
    I["ot"] = interpolate(out_tx,out_ty)
    # inner tongue
    inn_tx,inn_ty = inner(x_af,y_af,2.*info["t"]+   info["mt"])
    I["it"] = interpolate(inn_tx,inn_ty)
    # outer mouth
    out_mx,out_my = inner(x_af,y_af,2.*info["t"]+2.*info["mt"])
    I["om"] = interpolate(out_mx,out_my)
    # inner mouth
    inn_mx,inn_my = inner(x_af,y_af,3.*info["t"]+2.*info["mt"])
    I["im"] = interpolate(inn_mx,inn_my)
    # print(); print("here"); print()

    # create outer mouth structure
    x["mou"],y["mou"] = make_mouth(I["O"],I["I"],I["om"],info)

    # create tongue tip
    x["ton"],y["ton"] = tongue_tip(I["O"],I["I"],I["ot"],I["it"],info)

    # save certain values for future use
    info["t start x"] = x["mou"][5][-1]
    info["t start y"] = y["mou"][5][-1]
    info["b start x"] = x["ton"][0][-1]
    info["b start y"] = y["ton"][0][-1]
    info["t end x"]   = x["tri"][3][0]
    info["t end y"]   = y["tri"][3][0]
    info["b end x"]   = x["tri"][3][-1]
    info["b end y"]   = y["tri"][3][-1]

    x["ski"],y["ski"],x["sem"],y["sem"] = make_semicircles(I["I"],info)

    # # create roofs and floors
    # x["rfl"],y["rfl"] = roofs_n_floors(I["I"],info,x["mou"][5][-1],\
    #     y["tri"][3][0],y["tri"][3][1])

    # # create wing box
    # if not info["sk"]:
    #     x["box"],y["box"] = wing_box(I["I"],I["im"],info)

    # create outer run. start and end values are where break starts and ends
    x["out"],y["out"] = outer(I["O"],info)
    






    # # if desired to deflect, do so
    # if not info["dty"] == "none":
    #     # run through each set and deflect
    #     deflection(u,v,I["O"],info)

    # # if shift origin to quarter chord
    # if "shift to c/4" in vals["C14"] and vals["C14"]["shift to c/4"]:
    #     # shift
    #     shift(u,v,I["O"])

    # # resize chord of each segment
    # resize(u,v,info["c"])




    # initialize colors and labels
    g = {}
    # g["tri c"] = "#6a0dad"; g["tri l"] = "Triangle"
    # g["sem c"] = "#ff8c00"; g["sem l"] = "Semicircles"
    # g["mou c"] = "#0073ff"; g["mou l"] = "Mouth"
    # g["ton c"] = "#ff0000"; g["ton l"] = "Tongue"
    # g["ski c"] = "#3cb371"; g["ski l"] = "Skins"
    # g["box c"] = "#aa6c39"; g["box l"] = "Wing Box"
    # g["out c"] =       "k"; g["out l"] = "Outer"
    g["tri c"] = "b"; g["tri l"] = "PLA/TPU"
    g["sem c"] = "b"; g["sem l"] = ""
    g["mou c"] = "k"; g["mou l"] = "PLA"
    g["ton c"] = "b"; g["ton l"] = ""
    g["ski c"] = "b"; g["ski l"] = ""
    g["box c"] = "k"; g["box l"] = ""
    g["out c"] = "k"; g["out l"] = ""

    # plot to be saved items
    for group in x:
        # color
        if vals["C14"]["all black"]:
            c = "k"
        else:
            c = g[group+" c"]
        
        # if roofs and floors or arcs, do i and j
        if group == "arc": 
            for i in range(info["n"]): # 2): # 
                # 0 -> 4 from x = 0 -> 1.0
                for j in range(2):
                    # 0 - LE facing arc, 1 - TE facing arc
                    if i==0 and j==0:
                        l = g[group+" l"]
                    else:
                        l = ""
                    plt.plot(x[group][j][i],y[group][j][i],c=c,label=l)
        elif group == "rfl":
            for i in range(x["rfl"].shape[0]):
                for j in range(u["rfl"].shape[1]):
                    if i==0 and j==0:
                        plt.plot(x[group][i][j],y[group][i][j],c,\
                            label=g[group+" l"])
                    elif not(i==1 and j==0):
                        plt.plot(x[group][i][j],y[group][i][j],c)
        else:
            for i in range(x[group].shape[0]):
                if i == 0:
                    l = g[group+" l"]
                else:
                    l = ""

                if group == "out":
                    if i in [0,x[group].shape[0]-1,x[group].shape[0]-2]:
                        c = "b"
                    else:
                        c = "k"
                    
                    if i != 0:
                        plt.plot(x[group][i],y[group][i],c,label=l)
                    else:
                        d = 51
                        plt.plot(x[group][i][:d],y[group][i][:d],"b",label=l)
                        plt.plot(x[group][i][d-1:],y[group][i][d-1:],"k",label=l)

                
                else:
                    plt.plot(x[group][i],y[group][i],c,label=l)

    # if vals["C14"]["show legend"]:
    #     plt.legend()
    plt.xlabel("x/c")
    plt.ylabel("y/c")
    plt.axis("equal")
    plt.show()

    return #xvals, yvals, zvals

# run file
jsonfile = 'input_semicircles.json'
main(jsonfile)