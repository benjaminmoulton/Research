import airfoil_db as adb
import numpy as np
import json
import dxf
import svg
import deflect
import os
from scipy import interpolate as itp
from shapely.geometry import polygon as poly

from matplotlib import pyplot as plt

# function that reads in a json file as a dictionary
def read_json(filename):
    # import json file
    json_string=open(filename).read()
    vals = json.loads(json_string)
    return vals

# simplify a json dictionary
def simplify(vals,run):
    # initialize dictionary
    dicti = {}

    # transfer vals
    in2mm = 25.4
    dicti["c"]  = vals[run]["chord [in]"]

    dicti["nk"]  = vals[run]["num kinks"]
    dicti["knkns"] = vals[run]["kinkiness forward and aft"]
    dicti["kd"] = vals[run]["kink degree"]
    dicti["fkf"]  = vals[run]["first kink forward"]
    dicti["ktr"] = vals[run]["kinks to remove"]
    dicti["chkks"] = vals[run]["fillet kinks"]
    dicti["chkksh"] = vals[run]["fillet kinks sharps"]
    dicti["nkp"]  = vals[run]["num kink points"]
    dicti["kth"]  = [a / dicti["c"] / in2mm for a in \
        vals[run]["kink thickness [mm]"]]

    dicti["t"]  = vals[run]["shell thickness [mm]"] / dicti["c"] / in2mm
    dicti["f"]  = vals[run]["flex start [x/c]"]
    dicti["to"] = dicti["f"] - np.abs(vals[run]["tongue start [in]"]) \
        / dicti["c"]
    dicti["mo"] = dicti["f"] - np.abs(vals[run]["mouth start [in]"]) \
        / dicti["c"]
    if dicti["to"] < dicti["mo"]:
        raise ValueError("Tongue cannot extend beyond mouth start")
    dicti["tw"] = vals[run]["TE wall start [x/c]"]
    dicti["ts"] = vals[run]["Shift TE wall for hole"]
    dicti["mt"] = vals[run]["mouth clearance [mm]"] / dicti["c"] / in2mm
    dicti["fi"] = vals[run]["fillet side length [mm]"] / dicti["c"] / in2mm
    dicti["td"] = vals[run]["TE hole diameter [mm]"] / dicti["c"] / in2mm
    dicti["np"] = vals[run]["number of hole points"]

    dicti["srv"] = vals[run]["servo rotation pt [x,y] [x/c]"]
    dicti["x+s"] = vals[run]["servo x +shift [in]"] / dicti["c"]
    dicti["x-s"] = vals[run]["servo x -shift [in]"] / dicti["c"]
    dicti["sys"] = vals[run]["servo y shift [in]"] / dicti["c"]
    dicti["msh"] = vals[run]["make servo hole"]
    dicti["mch"] = vals[run]["make cord hole"]
    dicti["sp"]  = vals[run]["spar diameter [in]"] / dicti["c"]
    dicti["ch"]  = vals[run]["make control hole"]
    dicti["chd"] = vals[run]["control hole diameter [mm]"] / dicti["c"] / in2mm
    dicti["cwd"] = vals[run]["control washer diameter [mm]"] /dicti["c"] /in2mm

    dicti["da"] = vals[run]["deflection angle [deg]"]
    dicti["dty"] = vals[run]["deflection type"]
    dicti["th"] = vals[run]["hinge at top"]

    dicti["txt"] = vals[run]["svg text"]

    # hole shift vars
    if "TE hole x shift [in]" in vals[run]:
        dicti["tehxs"] = vals[run]["TE hole x shift [in]"] / dicti["c"]
    else:
        dicti["tehxs"] = 0.0
    if "TE hole y shift [in]" in vals[run]:
        dicti["tehys"] = vals[run]["TE hole y shift [in]"] / dicti["c"]
    else:
        dicti["tehys"] = 0.0
    if "tongue hole shift [in]" in vals[run]:
        dicti["tohs"] = vals[run]["tongue hole shift [in]"] / dicti["c"]
    else:
        dicti["tohs"] = 0.0

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
                    ind_max = int(np.max(np.array([i,j])))
                    ind_min = int(np.min(np.array([i,j])))
                    insct.append(np.array([x_int,y_int,ind_min,ind_max]))
    
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
        # if is_kink:
        #     num = x.shape[0] - num - 1
        ix = np.zeros((num+2,))
        iy = np.zeros((num+2,))

        # # add other points
        # if is_kink:
        #     # find indices max and min
        #     ind_max = int(np.max(np.array([it[0][2],it[0][3]])))
        #     ind_min = int(np.min(np.array([it[0][2],it[0][3]])))

        #     # add up to intersection
        #     ix[:ind_min] = x[:ind_min]
        #     iy[:ind_min] = y[:ind_min]

        #     # add intersection
        #     ix[ind_min] = it[0][0]
        #     iy[ind_min] = it[0][1]

        #     # add after intersection
        #     ix[ind_min+1:] = x[ind_max:]
        #     iy[ind_min+1:] = y[ind_max:]            
        # else:
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

# deintersect function
def deintersect_kink(x,y,it,fillet=False):
    # initialize x and y deintersected lists
    ix = []; iy = []

    # intialize counter for intersections
    j = 0
    
    # add extra intersections array to prevent failure
    itnew = np.zeros((it.shape[0]+1,it.shape[1]))
    itnew[:-1,:] = it * 1.
    itnew[-1,:] = np.array([0.0,0.0,x.shape[0],x.shape[0]])
    it = itnew * 1.

    if fillet:
        start = 0#3
    else:
        start = 0

    # run through each point
    for i in range(x.shape[0]):
        if it[j,2] <= i:
            ix.append(it[j,0]); iy.append(it[j,1])
            j += 1
        if j != 0 and j != it.shape[0]-1 and i >= it[j-1,3]+1+start:
            ix.append(x[i]); iy.append(y[i])
            if i == it[j,2] - 1:
                ix.append(x[i+1]); iy.append(y[i+1])


    # turn into numpy arrays
    ix = np.array(ix); iy = np.array(iy)
    
    return ix,iy

# deintersect function
def de_intersect_kink(x,y,it):
    # run through each interection and remove it
    for k in range(it.shape[0]):
        # initialize indices
        i = int(it[k,2]); j = int(it[k,3]) + 1

        # create a list of indices to remove
        inds = [a for a in range(i+1,j)]

        # delete indices, insert intersection pt
        x = np.delete(x,inds); y = np.delete(y,inds)
        x = np.insert(x,i+1,it[k,0]); y = np.insert(y,i+1,it[k,1])

        # decrease indices overall to 
        it[k:,2:] -= len(inds) - 1
    
    return x,y

# function that creates an offset with the kinks removed
# assume closed shape
def inner(x,y,d):

    # find offset curve
    offx,offy = offset(x,y,d)

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

# create top and bottom interpolation of xy points ( assume airfoil)
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
    # if x.shape[0] % 2 == 0:
    #     bot_x = x[i_mid-1:]
    #     bot_y = y[i_mid-1:]
    if top_x[0] == top_x[1]:
        top_x[1] -= 1e-10

    # create interpolations
    T = itp.interp1d(top_x,top_y,kind="cubic")
    B = itp.interp1d(bot_x,bot_y,kind="cubic")

    # save to dictionary
    I = {}
    I["t"] = T
    I["b"] = B
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
def get_point(xy,xc,yc,r,is_concave,is_top=True,bounds=[]):
    # solve for if top or bottom
    if is_top:
        loc = "t"
    else:
        loc = "b"
    
    # get minimum and maximum of the search range
    if not(len(bounds) == 0):
        minx = bounds[0]
        maxx = bounds[1]
    elif is_concave:
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

# define factorial function
def factorial(var):
    if var == 0 or var == 1:
        return 1
    else:
        return var * factorial(var-1)

# define binomial coefficient function
def binomial(n,k):
    # error prevention
    if not type(n) == int or not type(k) == int:
        raise ValueError("n and k variables must be integers")

    # calculate
    return factorial(n) / (factorial(k) * factorial(n-k))

# create bezier curve from two xy arrays over n points
def Bezier(x,y,k=3,n=30):
    # ensure method is used properly
    if x.shape[0] != y.shape[0]:
        raise ValueError("X and Y points arrays MUST be the same shape")
    if x.shape[0] != k+1:
        raise ValueError("X array must be 1 greater than degree in length")
    # initialize bezier arrays
    xbz = np.zeros((n,))
    ybz = np.zeros((n,))

    # determine t value (fraction of number of points)
    t   = np.linspace(0,1,num=n)

    # run through each degree and add the appropriate amount to the bezier arrs
    for i in range(k+1):
        factor = binomial(k,i) * (1-t)**(k-i) * t**i
        xbz += factor * x[i]
        ybz += factor * y[i]

    return xbz,ybz

# create arcs
def arcs(xy,info,shift=0.0):
    # initialize arcs arrays first index is 0 arc 1 skin
    x_arc = np.zeros((2,info["n"]),dtype=np.ndarray)
    y_arc = np.zeros((2,info["n"]),dtype=np.ndarray)
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
    x_end = info["tw"]# - info["t"]
    y_end,p_end,ie = get_y_loc(x_end,xy["cx"],xy["cy"])

    # determine "length" value at te wall
    l_end = p_end*(xy["cl"][ie] - xy["cl"][ie-1]) + xy["cl"][ie-1]

    # cycle through each arc
    for i in range(info["n"]): # 2): # 
        # initialize arrays
        x_arc[0][i] = np.zeros((info["np"],))
        y_arc[0][i] = np.zeros((info["np"],))
        x_arc[1][i] = np.zeros((info["np"],))
        y_arc[1][i] = np.zeros((info["np"],))

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
        r_thick = r_res + mult * info["ath"]

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
        x_arc[1][i] = np.flip(r_thick_arr * np.cos(theta_sk) + x_cent)
        y_arc[1][i] = np.flip(r_thick_arr * np.sin(theta_sk) + y_cent)

        if info["charc"]:
            # create bezier curve to remove arc stress concentration
            startit = 5; onei = 3; twoi = 4
            x_bez = np.array([x_arc[1][i][startit],
                            x_arc[1][i][onei],
                            x_arc[1][i][twoi],
                            x_arc[1][i][startit]+info["fi"]])
            y_bez = np.array([y_arc[1][i][startit],
                            y_arc[1][i][onei],
                            xy["t"](x_arc[1][i][twoi]),
                            xy["t"](x_arc[1][i][startit]+info["fi"])])
            xbzt,ybzt = Bezier(x_bez,y_bez)
            
            # flip to go in correct direction
            xbzt = np.flip(xbzt); ybzt = np.flip(ybzt)

            # bottom
            startib = -5; onei = -3; twoi = -4
            x_bez = np.array([x_arc[1][i][startib],
                            x_arc[1][i][onei],
                            x_arc[1][i][twoi],
                            x_arc[1][i][startib]+info["fi"]])
            y_bez = np.array([y_arc[1][i][startib],
                            y_arc[1][i][onei],
                            xy["b"](x_arc[1][i][twoi]),
                            xy["b"](x_arc[1][i][startib]+info["fi"])])
            xbzb,ybzb = Bezier(x_bez,y_bez)

            # create final array
            x_arc[1][i] = np.append(np.insert(x_arc[1][i][startit+1:startib],\
                0,xbzt),xbzb)
            y_arc[1][i] = np.append(np.insert(y_arc[1][i][startit+1:startib],\
                0,ybzt),ybzb)

        # # save top and bottom x y values to the points arrays in info
        # [i][ ][ ] 0 - top, 1 - bottom
        # [ ][i][ ] 0 -> 1 from LE to TE
        # [ ][ ][i] 0 - arc on left (LE), 1 - arc on right (TE)
        info["xpt"][0][i][0] = x_arc[0][i][-1]
        info["xpt"][1][i][0] = x_arc[0][i][0]
        info["xpt"][0][i][1] = x_arc[1][i][0]
        info["xpt"][1][i][1] = x_arc[1][i][-1]
        info["ypt"][0][i][0] = y_arc[0][i][-1]
        info["ypt"][1][i][0] = y_arc[0][i][0]
        info["ypt"][0][i][1] = y_arc[1][i][0]
        info["ypt"][1][i][1] = y_arc[1][i][-1]

    return x_arc,y_arc

# create arcs
def kincs(xy,info,shift=0.0):
    # initialize arcs arrays first index is 0 arc 1 skin
    x_spr = np.zeros((2,info["nk"]),dtype=np.ndarray)
    y_spr = np.zeros((2,info["nk"]),dtype=np.ndarray)
    # initialize arcs top and bottom coordinates arrays
    # [i][ ][ ] 0 - top, 1 - bottom
    # [ ][i][ ] 0 -> 1 from LE to TE
    # [ ][ ][i] 0 - arc on left (LE), 1 - arc on right (TE)
    
    info["xpt"] = np.zeros((2,info["nk"],2))
    info["ypt"] = np.zeros((2,info["nk"],2))

    # determine where flex starts on camberline
    x_flex = info["f"]
    y_flex,p_flex,ifx = get_y_loc(x_flex,xy["cx"],xy["cy"])

    # determine "length" value at flex start
    l_flex = p_flex*(xy["cl"][ifx] - xy["cl"][ifx-1]) + xy["cl"][ifx-1]

    # determine where te wall starts on camberline
    x_end = info["tw"]# - info["t"]
    y_end,p_end,ie = get_y_loc(x_end,xy["cx"],xy["cy"])

    # determine "length" value at te wall
    l_end = p_end*(xy["cl"][ie] - xy["cl"][ie-1]) + xy["cl"][ie-1]

    # cycle through each arc
    for i in range(info["nk"]): 
        # initialize arrays
        x_spr[0][i] = np.zeros((info["nkp"],))
        y_spr[0][i] = np.zeros((info["nkp"],))
        x_spr[1][i] = np.zeros((info["nkp"],))
        y_spr[1][i] = np.zeros((info["nkp"],))

        # set up x location
        x_loc = info["f"] + xy["arcdx"] * (i+1) + shift
        
        # determine y location
        y_loc,perc,it = get_y_loc(x_loc,xy["cx"],xy["cy"])

        # get x and y value where the intersection occurs
        x_top = intersect(xy,it,x_loc,y_loc)
        x_bot = intersect(xy,it,x_loc,y_loc,is_top=False)
        y_top = xy["t"](x_top)
        y_bot = xy["b"](x_bot)

        # determine "length" value at locale
        l_loc = perc*(xy["cl"][it] - xy["cl"][it-1]) + xy["cl"][it-1]

        # set lperc
        lperc = (l_end-l_loc) / (l_end-l_flex)

        # calculate kink thickness
        if len(info["kth"]) == 1:
            kink_thickness = info["kth"][0]
        else:
            # caclulate fraction length
            frac_l = (l_loc-l_flex) / (l_end-l_flex)

            # calculate kink_thickness
            kink_thickness = info["kth"][0] + frac_l * \
                (info["kth"][1] - info["kth"][0])
            

        # determine start and end points for kinks
        xpts = np.linspace(x_top,x_bot,num=2+info["kd"])
        ypts = np.linspace(y_top,y_bot,num=2+info["kd"])

        # if start forward
        if info["fkf"]:
            neg = 1
        else:
            neg = -1

        # shift to create kink
        for j in range(1,1+info["kd"]):
            if j % 2 == 1:
                xpts[j] -= xy["arcdx"] * info["knkns"][0] * lperc * neg
            else:
                xpts[j] += xy["arcdx"] * info["knkns"][1] * lperc * neg

        # run through each degree kink
        midx = np.array([]); midy = np.array([])
        for j in range(info["kd"]+1):
            # initialize bezier points
            xbpt = np.linspace(xpts[j],xpts[j+1],num=4)
            ybpt = np.linspace(ypts[j],ypts[j+1],num=4)

            # shift
            xbpt[1] = xbpt[0] * 1.;  ybpt[1] = (ybpt[0] + 3.*ybpt[-1]) / 4.
            xbpt[2] = xbpt[-1] * 1.; ybpt[2] = (3.*ybpt[0] + ybpt[-1]) / 4.

            # create bezier curve
            xb,yb = Bezier(xbpt,ybpt,k=3,n=30)
            # add to array
            midx = np.append(midx,xb[:-1])
            midy = np.append(midy,yb[:-1])
        # add final point
        midx = np.append(midx,xb[-1])
        midy = np.append(midy,yb[-1])

        # create offset kink
        x_spr[0][i],y_spr[0][i] = offset(midx,midy,-kink_thickness/2.)
        x_spr[1][i],y_spr[1][i] = offset(midx,midy,kink_thickness/2.)
        
        # find where the intersections occur
        it0 = intersections(x_spr[0][i],y_spr[0][i])
        it1 = intersections(x_spr[1][i],y_spr[1][i])

        # remove them
        x_spr[0][i],y_spr[0][i] = de_intersect_kink(x_spr[0][i],y_spr[0][i],\
            it0)
        x_spr[1][i],y_spr[1][i] = de_intersect_kink(x_spr[1][i],y_spr[1][i],\
            it1)

        if info["chkks"]:
            # create bezier curve to remove arc stress concentration
            if info["fkf"]:
                a = 0; b = 1
                shifty = -info["fi"]
                ai = [3,4,5]; bi = [0,0,1]
            else:
                a = 1; b = 0
                shifty = info["fi"]
                ai = [3,4,5]; bi = [0,0,1]
            xa = np.array([x_spr[a][i][ai[-1]]+shifty, x_spr[a][i][ai[2]],
                            x_spr[a][i][ai[1]], x_spr[a][i][ai[-1]]])
            ya = np.array([xy["t"](x_spr[a][i][ai[-1]]+shifty),
                            xy["t"](x_spr[a][i][ai[2]]),
                            y_spr[a][i][ai[1]], y_spr[a][i][ai[-1]]])
            xb = np.array([x_spr[b][i][bi[-1]]-shifty, x_spr[b][i][bi[2]],
                            x_spr[b][i][bi[1]], x_spr[b][i][bi[-1]]])
            yb = np.array([xy["t"](x_spr[b][i][bi[-1]]-shifty),
                            xy["t"](x_spr[b][i][bi[2]]),
                            y_spr[b][i][bi[1]], y_spr[b][i][bi[-1]]])
            xazt,yazt = Bezier(xa,ya,n=10)
            xbzt,ybzt = Bezier(xb,yb,n=10)

            # create final array
            x_spr[a][i] = np.insert(x_spr[a][i][ai[-1]+1:],0,xazt)
            y_spr[a][i] = np.insert(y_spr[a][i][ai[-1]+1:],0,yazt) 
            x_spr[b][i] = np.insert(x_spr[b][i][bi[-1]+1:],0,xbzt)
            y_spr[b][i] = np.insert(y_spr[b][i][bi[-1]+1:],0,ybzt)

            # bottom
            is_odd = info["kd"] % 2 == 1
            if is_odd and info["fkf"] or not is_odd and not info["fkf"]:
                a = 0; b = 1
                shifty = -info["fi"]
                ai = [-4,-5,-6]; bi = [-1,-1,-2]
            else:
                a = 1; b = 0
                shifty = info["fi"]
                ai = [-4,-5,-6]; bi = [-1,-1,-2]
            xa = np.array([x_spr[a][i][ai[-1]], x_spr[a][i][ai[1]],
                            x_spr[a][i][ai[2]], x_spr[a][i][ai[-1]]+shifty])
            ya = np.array([y_spr[a][i][ai[-1]], y_spr[a][i][ai[1]],
                            xy["b"](x_spr[a][i][ai[2]]),
                            xy["b"](x_spr[a][i][ai[-1]]+shifty)])
            xb = np.array([x_spr[b][i][bi[-1]], x_spr[b][i][bi[1]],
                            x_spr[b][i][bi[2]], x_spr[b][i][bi[-1]]-shifty])
            yb = np.array([y_spr[b][i][bi[-1]], y_spr[b][i][bi[1]],
                            xy["b"](x_spr[b][i][bi[2]]),
                            xy["b"](x_spr[b][i][bi[-1]]-shifty)])
                            
            xazt,yazt = Bezier(xa,ya,n=10)
            xbzt,ybzt = Bezier(xb,yb,n=10)

            # create final array
            x_spr[a][i] = np.append(x_spr[a][i][:ai[-1]],xazt)
            y_spr[a][i] = np.append(y_spr[a][i][:ai[-1]],yazt) 
            x_spr[b][i] = np.append(x_spr[b][i][:bi[-1]],xbzt)
            y_spr[b][i] = np.append(y_spr[b][i][:bi[-1]],ybzt)
        
        # flip
        x_spr[0][i] = np.flip(x_spr[0][i]); y_spr[0][i] = np.flip(y_spr[0][i])

        # # save top and bottom x y values to the points arrays in info
        # [i][ ][ ] 0 - top, 1 - bottom
        # [ ][i][ ] 0 -> 1 from LE to TE
        # [ ][ ][i] 0 - arc on left (LE), 1 - arc on right (TE)
        top = -1; bot = 0
        info["xpt"][0][i][0] = x_spr[0][i][top]
        info["xpt"][1][i][0] = x_spr[0][i][bot]
        info["xpt"][0][i][1] = x_spr[1][i][bot]
        info["xpt"][1][i][1] = x_spr[1][i][top]
        info["ypt"][0][i][0] = y_spr[0][i][top]
        info["ypt"][1][i][0] = y_spr[0][i][bot]
        info["ypt"][0][i][1] = y_spr[1][i][bot]
        info["ypt"][1][i][1] = y_spr[1][i][top]




    return x_spr,y_spr

# reorder a set of x y arrays as defined
def reorder(x,y,order,flip=[]):
    # initialize new arrays
    new_x = np.zeros((x.shape[0],),dtype=np.ndarray)
    new_y = np.zeros((y.shape[0],),dtype=np.ndarray)

    if len(order) == 0:
        n = new_x.shape[0]
        order = np.linspace(0,n-1,n).astype(int)

    # reorder the arrays
    for i in range(new_x.shape[0]):
        new_x[i] = x[order[i]] * 1.
        new_y[i] = y[order[i]] * 1.

    # flip if necessary
    for j in range(len(flip)):
        new_x[flip[j]] = np.flip(new_x[flip[j]])
        new_y[flip[j]] = np.flip(new_y[flip[j]])

    return new_x,new_y

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
    # find where the mouth starts as an interpolated point inner
    x_ms_i = info["mo"]
    y_ms_i,p_ms_i,i_ms_i = get_y_loc(x_ms_i,I["bx"],I["by"])

    # initialize bottom of inner mouth
    x_mouth[1] = np.zeros((i_fx_i-i_ms_i+2,))
    y_mouth[1] = np.zeros((i_fx_i-i_ms_i+2,))

    # set first point
    x_mouth[1][0] = x_ms_i
    y_mouth[1][0] = y_ms_i

    # set in between points
    j = 1
    for i in range(i_ms_i,i_fx_i):
        x_mouth[1][j] = I["bx"][i]
        y_mouth[1][j] = I["by"][i]
        j += 1
    
    # set last point
    x_mouth[1][-1] = x_fx_i
    y_mouth[1][-1] = y_fx_i

    # create back of inner mouth
    # find where the mouth starts as an interpolated point outer mouth
    x_ms_om = info["mo"]
    y_ms_om,p_ms_om,i_ms_om = get_y_loc(x_ms_om,om["bx"],om["by"])

    # initialize back outer mouth
    x_mouth[2] = np.array([x_ms_om,x_ms_i])
    y_mouth[2] = np.array([y_ms_om,y_ms_i])

    # create top of inner mouth
    # find where flex starts as an interpolated point outer mouth
    x_fx_om = info["f"]
    y_fx_om,p_fx_om,i_fx_om = get_y_loc(x_fx_om,om["bx"],om["by"])

    # initialize bottom of inner mouth
    x_mouth[3] = np.zeros((i_fx_om-i_ms_om+2,))
    y_mouth[3] = np.zeros((i_fx_om-i_ms_om+2,))

    # set first point
    x_mouth[3][0] = x_ms_om
    y_mouth[3][0] = y_ms_om

    # set in between points
    j = 1
    for i in range(i_ms_om,i_fx_om):
        x_mouth[3][j] = om["bx"][i]
        y_mouth[3][j] = om["by"][i]
        j += 1
    
    # set last point
    x_mouth[3][-1] = x_fx_om
    y_mouth[3][-1] = y_fx_om

    # create "upper lip" line
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

    # save fillet end info
    info["fnd"] = [x_fi_ti,y_fi_ti,p_fi_ti,i_fi_ti]

    # initialze point arrays to get the bezier curve
    x_points = np.array([x_fx_ti,x_fx_ti,(x_fx_ti+x_fi_ti)/2.,x_fi_ti])
    y_points = np.array([y_fx_ti - info["fi"],(y_fx_ti - info["fi"]+y_fi_ti)\
        /2.,(y_fx_ti+y_fi_ti)/2.,y_fi_ti])

    # create fillet using a bezier curve maker
    x_mouth[5],y_mouth[5] = Bezier(x_points,y_points)

    reorder
    x_mouth,y_mouth = reorder(x_mouth,y_mouth,[],[1,2,4])

    return x_mouth,y_mouth

# create tongue tip geometry
def tongue_tip(O,I,ot,it,ttt,info):
    # initialize mouth array
    x_tongue = np.zeros((9,),dtype=np.ndarray)
    y_tongue = np.zeros((9,),dtype=np.ndarray)

    # create a bezier curve from tongue tip to the bottom of first arc
    # find interpolation point where flex starts on inner tongue
    x_fx_it = info["f"]
    y_fx_it,p_fx_it,i_fx_it = get_y_loc(x_fx_it,it["bx"],it["by"])

    # create points for the bezier curve
    y_start,_,_ = get_y_loc(             x_fx_it+.03,it["bx"],it["by"])
    y_end,__,__ = get_y_loc(info["xpt"][1][0][0]-.03,I["bx"],I["by"])

    # create bezier curve
    xpts = np.array([x_fx_it,x_fx_it+.03,info["xpt"][1][0][0]-.03,\
        info["xpt"][1][0][0]])
    ypts = np.array([y_fx_it,y_start,y_end,\
        info["ypt"][1][0][0]])
    x_tongue[0],y_tongue[0] = Bezier(xpts,ypts,n=60)

    # create top of tongue tip
    # find where the tongue starts inner tongue
    x_ts_it = info["to"]
    y_ts_it,p_ts_it,i_ts_it = get_y_loc(x_ts_it,ttt["bx"],ttt["by"])
    if info["ch"]:
        x_pos = info["to"] + 3*info["t"] + info["td"] + info["cwd"]/2. +\
            info["chd"]
    else:
        x_pos = info["to"] + 2*info["t"] + info["td"]

    y_pos,p_pos,i_pos = get_y_loc(x_pos,ttt["bx"],ttt["by"])

    # initialize bottom of inner mouth
    x_tongue[1] = np.zeros((i_pos-i_ts_it+2,))
    y_tongue[1] = np.zeros((i_pos-i_ts_it+2,))

    # set first point
    x_tongue[1][0] = x_ts_it
    y_tongue[1][0] = y_ts_it

    # set in between points
    j = 1
    for i in range(i_ts_it,i_pos):
        x_tongue[1][j] = ttt["bx"][i]
        y_tongue[1][j] = ttt["by"][i]
        j += 1
    
    # set last point
    x_tongue[1][-1] = x_pos
    y_tongue[1][-1] = y_pos

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
    # find where the tongue joints outer tongue
    x_j_ot = x_tongue[3][0]
    y_j_ot,p_j_ot,i_j_ot = get_y_loc(x_j_ot,ot["bx"],ot["by"])

    # initialize bottom of inner mouth
    x_tongue[4] = np.zeros((i_j_ot-i_ts_ot+2,))
    y_tongue[4] = np.zeros((i_j_ot-i_ts_ot+2,))

    # set first point
    x_tongue[4][0] = x_ts_ot
    y_tongue[4][0] = y_ts_ot

    # set in between points
    j = 1
    for i in range(i_ts_ot,i_j_ot):
        x_tongue[4][j] = ot["bx"][i]
        y_tongue[4][j] = ot["by"][i]
        j += 1
    
    # set last point
    x_tongue[4][-1] = x_j_ot
    y_tongue[4][-1] = y_j_ot

    # create line to fix discontinuity between bezier tongue bottom curve
    # and bottom tongue imported points
    # save to tongue tip
    x_tongue[5] = np.array([x_j_ot,x_j_ot])
    y_tongue[5] = np.array([y_j_ot,y_tongue[3][0]])

    # create line to fix discontinuity between bezier tongue bottom curve
    # and outer airfoil shape
    # determine the outer airfoil interpolation point from first arc bottom pt
    x_j_o = info["xpt"][1][0][0]
    y_j_o,p_j_o,i_j_o = get_y_loc(x_j_o,O["bx"],O["by"])

    # save to tongue tip
    x_tongue[6] = np.array([x_tongue[3][-1],x_j_o])
    y_tongue[6] = np.array([y_tongue[3][-1],y_j_o])
    # save value so info to be used in a later function
    info["outer"]["end"] = [x_j_o,y_j_o,i_j_o]

    # create line up on top of tongue over control hole
    y_bos,p_bos,i_bos = get_y_loc(x_pos,it["bx"],it["by"])
    x_tongue[7] = np.array([x_pos,x_pos])
    y_tongue[7] = np.array([y_bos,y_pos])

    # initialize bottom of inner mouth
    x_tongue[8] = np.zeros((i_fx_it-i_bos+2,))
    y_tongue[8] = np.zeros((i_fx_it-i_bos+2,))

    # set first point
    x_tongue[8][0] = x_pos
    y_tongue[8][0] = y_bos

    # set in between points
    j = 1
    for i in range(i_bos,i_fx_it):
        x_tongue[8][j] = it["bx"][i]
        y_tongue[8][j] = it["by"][i]
        j += 1

    # set last point
    x_tongue[8][-1] = x_fx_it
    y_tongue[8][-1] = y_fx_it

    # reorder
    x_tongue,y_tongue = reorder(x_tongue,y_tongue,[0,8,7,1,2,4,5,3,6],[0,1,3])

    return x_tongue,y_tongue

# create roofs and floors
def roofs_n_floors(I,info,x_fillet):
    # initialize mouth array
    # [i][ ] 0 - roof, 1 - floor
    # [ ][i] 0 -> 1, LE -> TE
    # NOTE [1][0] will remain empty
    x_rofl = np.zeros((2,info["nk"]+1),dtype=np.ndarray)
    y_rofl = np.zeros((2,info["nk"]+1),dtype=np.ndarray)

    # initialize te wall top and bottom
    y_tri_top = I["t"](info["tw"])
    y_tri_bot = I["b"](info["tw"])

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
                if j == info["nk"]:
                    # determine the xy location of the fillet end
                    x_f_i = info["tw"]
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
                
                # initialize array for points
                x_rofl[i][j] = np.zeros((iend-istart+2,))
                y_rofl[i][j] = np.zeros((iend-istart+2,))
                
                # add the start to the points
                x_rofl[i][j][0] = xend
                y_rofl[i][j][0] = yend

                # add inbetween points
                k = 1
                for l in range(iend-1,istart-1,-1):
                    x_rofl[i][j][k] = I["tx"][l]
                    y_rofl[i][j][k] = I["ty"][l]
                    k += 1
                
                # add end points
                x_rofl[i][j][-1] = xstart
                y_rofl[i][j][-1] = ystart
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
                if j == info["nk"]:
                    # determine the xy location of the fillet end
                    x_f_i = info["tw"]
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
                
                # initialize array for points
                x_rofl[i][j] = np.zeros((iend-istart+2,))
                y_rofl[i][j] = np.zeros((iend-istart+2,))
                
                # add the start to the points
                x_rofl[i][j][0] = xstart
                y_rofl[i][j][0] = ystart

                # add inbetween points
                k = 1
                for l in range(istart,iend):
                    x_rofl[i][j][k] = I["bx"][l]
                    y_rofl[i][j][k] = I["by"][l]
                    k += 1
                
                # add end points
                x_rofl[i][j][-1] = xend
                y_rofl[i][j][-1] = yend
            # else do the te wall
            else:
                # initialize and save arrays
                x_rofl[i][j] = np.array([info["tw"],info["tw"]])
                y_rofl[i][j] = np.array([y_tri_bot,y_tri_top])

    return x_rofl,y_rofl

# remove any arcs desired
def remove_arcs(x,y,info):
    # check to make sure all arcs to remove are available
    for i in range(len(info["atr"])):
        if info["atr"][i] < 1:
            raise ValueError("ARC to remove value {}".format(info["atr"][i]) +\
                " below allowable values (1)")
        elif info["atr"][i] > info["n"]:
            raise ValueError("ARC to remove value {}".format(info["atr"][i]) +\
                " above allowable value {}".format(info["n"]))
        elif type(info["atr"][i]) != int:
            raise ValueError("ARC to remove value {}".format(info["atr"][i])\
                + " is not an integer")
    
    # remove arcs
    info["atr"] = np.array(info["atr"])
    
    # set to python indexing
    info["atr"] -= 1

    # run through each arc and remove it
    for i in range(len(info["atr"])):
        # set index
        j = info["atr"][i]

        # create line at top to remove arcs
        x["arc"][0][j] = np.array([x["rfl"][0][j+1][-1],x["rfl"][0][j][0]])
        y["arc"][0][j] = np.array([y["rfl"][0][j+1][-1],y["rfl"][0][j][0]])
        # create line at bottom to remove arcs
        if j == 0:
            x["arc"][1][j] = np.array([x["ton"][0][0],x["rfl"][1][j+1][0]])
            y["arc"][1][j] = np.array([y["ton"][0][0],y["rfl"][1][j+1][0]])
        else:
            x["arc"][1][j] = np.array([x["rfl"][1][j][-1],x["rfl"][1][j+1][0]])
            y["arc"][1][j] = np.array([y["rfl"][1][j][-1],y["rfl"][1][j+1][0]])

    return

# remove any arcs desired
def remove_springs(x,y,info):
    # check to make sure all arcs to remove are available
    for i in range(len(info["ktr"])):
        if info["ktr"][i] < 1:
            raise ValueError("Kink to remove value {}".format(info["ktr"][i])+\
                " below allowable values (1)")
        elif info["ktr"][i] > info["n"]:
            raise ValueError("Kink to remove value {}".format(info["ktr"][i])+\
                " above allowable value {}".format(info["n"]))
        elif type(info["ktr"][i]) != int:
            raise ValueError("Kink to remove value {}".format(info["ktr"][i])\
                + " is not an integer")
    
    # remove arcs
    info["ktr"] = np.array(info["ktr"])
    
    # set to python indexing
    info["ktr"] -= 1

    # run through each arc and remove it
    for i in range(len(info["ktr"])):
        # set index
        j = info["ktr"][i]

        # create line at top to remove arcs
        x["spr"][0][j] = np.array([x["rfl"][0][j+1][-1],x["rfl"][0][j][0]])
        y["spr"][0][j] = np.array([y["rfl"][0][j+1][-1],y["rfl"][0][j][0]])
        # create line at bottom to remove arcs
        if j == 0:
            x["spr"][1][j] = np.array([x["ton"][0][0],x["rfl"][1][j+1][0]])
            y["spr"][1][j] = np.array([y["ton"][0][0],y["rfl"][1][j+1][0]])
        else:
            x["spr"][1][j] = np.array([x["rfl"][1][j][-1],x["rfl"][1][j+1][0]])
            y["spr"][1][j] = np.array([y["rfl"][1][j][-1],y["rfl"][1][j+1][0]])

    return

# determine a point which is a set thickness away from both skins
def find_thickness_x(foil,thickness,xl,xu):
    # define function to be used
    def fun(x):
        return foil.get_thickness(x) - thickness
    
    # find where thickness is equal to that desired
    thick_x = Bisect(fun,xl,xu, Err = 1e-10, maxiter = 10000)

    return thick_x

# determine where the hole can be placed inside the inner airfoil
def find_wall(O,foil,info,report = False):
    # find where wall should start
    xl = 0.5; xu = np.max(O["tx"]); thick = info["td"]/2. + info["t"]
    info["tw"] = find_thickness_x(foil,thick,xl,xu) - thick - info["tehxs"]
    if report:
        print("New TE wall value : {:<7.5} x/c".format(info["tw"]))
    
    return info["tw"]

# create te hole
def te_hole(O,I,it,foil,info):
    # initialize outer arrays
    number = 4
    if info["ch"]:
        number += 1
    x_hole = np.zeros((number,),dtype=np.ndarray)
    y_hole = np.zeros((number,),dtype=np.ndarray)

    # determine front spar hole center
    thickness = info["t"] + info["sp"]*0.8#/2.
    x_center = find_thickness_x(foil,thickness,np.min(O["bx"]),0.5)
    y_center = foil.get_camber(x_center)

    # save circle center
    info["front spar"] = [x_center,y_center]

    # initialize theta, radius values
    theta = np.linspace(0,2.*np.pi,num=info["np"])
    radius = np.full(theta.shape,info["sp"]/2.)

    # input array values
    x_hole[0] = radius * np.cos(theta) + x_center
    y_hole[0] = radius * np.sin(theta) + y_center


    # determine aft spar hole center
    x_center = info["f"] - info["t"] - info["sp"]/2.
    y_center = O["t"](x_center) - info["t"] - info["sp"]/2.
    # save circle center
    info["aft spar"] = [x_center,y_center]

    # initialize theta, radius values
    radius = np.full(theta.shape,info["sp"]/2.)

    # input array values
    x_hole[1] = radius * np.cos(theta) + x_center
    y_hole[1] = radius * np.sin(theta) + y_center

    # determine hole center
    x_center = info["tw"] + info["t"] + info["td"]/2. + info["tehxs"]
    y_center = O["c"](x_center)
    print("y center",y_center,O["c"](0.5))

    # initialize theta, radius values
    radius = np.full(theta.shape,info["td"]/2.)

    # create shift array to shift hole center
    gr_pi  = theta > np.pi / 2.
    le_2pi = theta < 3*np.pi / 2.
    shift = np.logical_and(gr_pi,le_2pi)

    # input array values
    x_hole[2] = radius * np.cos(theta) + x_center - shift * info["tehxs"]
    y_hole[2] = radius * np.sin(theta) + y_center - shift * info["tehys"]

    # determine tongue connection hole center
    x_center = info["to"] + info["t"] + info["td"]/2.
    y_center = it["b"](x_center) + info["td"]/2.

    # input array values
    x_hole[3] = radius * np.cos(theta) + x_center
    y_hole[3] = radius * np.sin(theta) + y_center
    

    # if control hole desired
    if info["ch"]:
        # determine tongue control hole center
        x_center += info["td"]/2. + info["t"] + info["cwd"]/2.
        y_center  = it["b"](x_center) + info["cwd"]/2.

        # initialize radius values
        radius = np.full(theta.shape,info["chd"]/2.)

        # input array values
        x_hole[-1] = radius * np.cos(theta) + x_center
        y_hole[-1] = radius * np.sin(theta) + y_center

        #######################################################
        info["chxy"] = [x_center,y_center]
        #######################################################

    return x_hole,y_hole

# find x given y
def find_x_from_y(xy,yvalue,bounds,is_top=True):
    # determine which interp to use
    if is_top:
        interp = xy["t"]
    else:
        interp = xy["b"]
    
    # initialize function to find x value
    def fun(x):
        return interp(x) - yvalue
    
    # run bisect function
    xvalue = Bisect(fun,bounds[0],bounds[1])

    return xvalue

# create front inner runaround
def front_inner(I,im,info):
    # initialize arrays # add more if servo, remove some if servo and cord
    arrays_num = 7 + 4*info["msh"] - 2*(info["mch"] and info["msh"])
    x_inner = np.zeros((arrays_num,),dtype=np.ndarray)
    y_inner = np.zeros((arrays_num,),dtype=np.ndarray)

    # determine angle locations around aft spar
    x0,y0 = get_point(I,info["aft spar"][0],info["aft spar"][1],\
        info["t"]+info["sp"]/2.,is_concave=False)
    theta0 =  np.arctan2(y0-info["aft spar"][1],x0-info["aft spar"][0])
    theta0 -= 2.*np.pi
    y1 = info["aft spar"][1] - ((info["t"]+info["sp"]/2.)**2. - \
        (info["aft spar"][0]-(info["f"]-info["t"]))**2.)**0.5
    theta1 = np.arctan2(y1-info["aft spar"][1],\
        (info["f"]-info["t"])-info["aft spar"][0])    
    
    # calculate runaround circle
    theta  = np.linspace(theta0,theta1,num=info["np"])
    radius = np.full(theta.shape,info["t"]+info["sp"]/2.)
    x_inner[0] = np.flip(radius * np.cos(theta) + info["aft spar"][0])
    y_inner[0] = np.flip(radius * np.sin(theta) + info["aft spar"][1])

    # determine straight run down flex wall
    x_inner[1] = np.array([info["f"]-info["t"],x_inner[0][0]])
    y_inner[1] = np.array([im["b"](info["f"]-info["t"]),y_inner[0][0]])

    # determine run along inner mouth to back of mouth + skin
    # determine point aft
    x_im_ft = info["f"] - info["t"]
    y_im_ft,p_im_ft,i_im_ft = get_y_loc(x_im_ft,im["bx"],im["by"])
    # determine forward point
    x_im_mt = info["mo"] - info["t"]
    # if text is to be added to this place, add it here
    x_im_mt -= (len(info["txt"])*.0131*int(info["c"]))/info["c"] + info["t"]
    # save text start location to info
    info["txt start"] = [(x_im_mt+info["t"])*info["c"],
                        (I["b"](x_im_mt+info["t"])+info["t"])*info["c"]]
    y_im_mt,p_im_mt,i_im_mt = get_y_loc(x_im_mt,im["bx"],im["by"])

    # initialize array
    x_inner[2] = np.zeros((i_im_ft-i_im_mt+2,))
    y_inner[2] = np.zeros((i_im_ft-i_im_mt+2,))

    # add first point
    x_inner[2][0] = x_im_mt
    y_inner[2][0] = y_im_mt

    # run through up to final point
    j = 1
    for i in range(i_im_mt,i_im_ft):
        x_inner[2][j] = im["bx"][i]
        y_inner[2][j] = im["by"][i]
        j += 1
    
    # add final point
    x_inner[2][-1] = x_inner[1][0]
    y_inner[2][-1] = y_inner[1][0]

    # determine back of mouth line
    # determine bottom point
    x_i_mt = x_im_mt
    y_i_mt,p_i_mt,i_i_mt = get_y_loc(x_i_mt,I["bx"],I["by"])

    # save to array
    x_inner[3] = np.array([x_i_mt,x_im_mt])
    y_inner[3] = np.array([y_i_mt,y_im_mt])

    # determine theta values around front spar
    x0,y0 = get_point(I,info["front spar"][0],info["front spar"][1],\
        info["t"]+info["sp"]/2.,is_concave=True,is_top=False)
    theta0 =  np.arctan2(y0-info["front spar"][1],x0-info["front spar"][0])
    x1,y1 = get_point(I,info["front spar"][0],info["front spar"][1],\
        info["t"]+info["sp"]/2.,is_concave=True)
    theta1 =  np.arctan2(y1-info["front spar"][1],x1-info["front spar"][0])

    # create theta array, and consequently x and y arrays
    theta = np.linspace(theta0,theta1,info["np"])
    rad = np.full(theta.shape,info["t"]+info["sp"]/2.)
    x_inner[4] = np.flip(rad * np.cos(theta) + info["front spar"][0])
    y_inner[4] = np.flip(rad * np.sin(theta) + info["front spar"][1])

    # run from run around to front spar along bottom to back of mouth
    # find index of bottom part
    x_i_fs = x_inner[4][-1]
    y_i_fs,p_i_fs,i_i_fs = get_y_loc(x_i_fs,I["bx"],I["by"])
    y_i_fs = y_inner[4][-1]

    # initialize array
    x_inner[5] = np.zeros((i_i_mt-i_i_fs+2,))
    y_inner[5] = np.zeros((i_i_mt-i_i_fs+2,))

    # add first point
    x_inner[5][0] = x_i_fs
    y_inner[5][0] = y_i_fs

    # run through up to final point
    j = 1
    for i in range(i_i_fs,i_i_mt):
        x_inner[5][j] = I["bx"][i]
        y_inner[5][j] = I["by"][i]
        j += 1
    
    # add final point
    x_inner[5][-1] = x_i_mt
    y_inner[5][-1] = y_i_mt

    # determine points for top of spar circle runarounds
    # find index of front spar
    x_i_fst = x_inner[4][0]
    y_i_fst,p_i_fst,i_i_fst = get_y_loc(x_i_fst,I["tx"],I["ty"])
    y_i_fst = y_inner[4][0]
    # find index of aft spar
    x_i_ast = x_inner[0][-1]
    y_i_ast,p_i_ast,i_i_ast = get_y_loc(x_i_ast,I["tx"],I["ty"])
    y_i_ast = y_inner[0][-1]

    # if a servo hole is not to be added, run all around
    if not info["msh"]:
        x_a = x_i_fst; y_a = y_i_fst; i_a = i_i_fst
    else: # go to end of servo
        x_a = info["srv"][0] + info["x+s"] + info["t"]
        y_a, p_a, i_a = get_y_loc(x_a,I["tx"],I["ty"])
    
    # initialize array
    x_inner[6] = np.zeros((i_i_ast-i_a+2,))
    y_inner[6] = np.zeros((i_i_ast-i_a+2,))

    # add first point
    x_inner[6][0] = x_i_ast
    y_inner[6][0] = y_i_ast

    # run through up to final point
    j = 1
    for i in range(i_i_ast-1,i_a-1,-1):
        x_inner[6][j] = I["tx"][i]
        y_inner[6][j] = I["ty"][i]
        j += 1
    
    # add final point
    x_inner[6][-1] = x_a
    y_inner[6][-1] = y_a

    # add last arrays if servo hole to be added
    if info["msh"]:
        x_inner[7] = np.array([x_a,x_a])
        y_inner[7] = np.array([y_a,info["srv"][1] - info["sys"] - info["t"]])

        # add arrays if cord hole to be added or not
        if not info["mch"]:
            # create run along bottom of servo
            x_inner[8] = np.array([info["srv"][0]+info["x+s"]+info["t"],\
                info["srv"][0]-info["x-s"]-info["t"]])
            y_inner[8] = np.array([info["srv"][1]-info["sys"]-info["t"],\
                info["srv"][1]-info["sys"]-info["t"]])
            
            # find point where servo front side intersects inner skin
            x_i_srv = info["srv"][0]-info["x-s"]-info["t"]
            y_i_srv,p_i_srv,i_i_srv = get_y_loc(x_i_srv,I["tx"],I["ty"])

            # create servo front side array
            x_inner[9] = np.array([x_i_srv,x_i_srv])
            y_inner[9] = np.array([info["srv"][1]-info["sys"]-info["t"],\
                y_i_srv])
    
            # initialize array
            x_inner[10] = np.zeros((i_i_srv-i_i_fst+2,))
            y_inner[10] = np.zeros((i_i_srv-i_i_fst+2,))

            # add first point
            x_inner[10][0] = x_i_srv
            y_inner[10][0] = y_i_srv

            # run through up to final point
            j = 1
            for i in range(i_i_srv-1,i_i_fst-1,-1):
                x_inner[10][j] = I["tx"][i]
                y_inner[10][j] = I["ty"][i]
                j += 1
            
            # add final point
            x_inner[10][-1] = x_i_fst
            y_inner[10][-1] = y_i_fst
        else: # add cord hole
            y_int = info["srv"][1]-info["sys"]-info["t"]
            if np.max(y_inner[4]) <= y_int:
                # fix number of arrays
                x_innew = np.zeros((arrays_num+1,),dtype=np.ndarray)
                y_innew = np.zeros((arrays_num+1,),dtype=np.ndarray)
                x_innew[:-1] = x_inner *1.; y_innew[:-1] = y_inner *1.
                x_inner = x_innew; y_inner = y_innew

                # find where the line will intersect the inner airfoil
                x_int = find_x_from_y(I,y_int,[np.min(I["tx"]),0.5])

                # create a line to that point
                # create run along bottom of servo
                x_inner[8] = np.array([info["srv"][0]+info["x+s"]+info["t"],\
                    x_int])
                y_inner[8] = np.array([y_int,y_int])

                # run around airfoil for this short stretch yur
                # determine index of this point yur
                _,p_int,i_int = get_y_loc(x_int,I["tx"],I["ty"])
    
                # initialize array
                x_inner[9] = np.zeros((i_int-i_i_fst+2,))
                y_inner[9] = np.zeros((i_int-i_i_fst+2,))

                # add first point
                x_inner[9][0] = x_i_fst
                y_inner[9][0] = y_i_fst

                # run through up to final point
                j = 1
                for i in range(i_i_fst,i_int):
                    x_inner[9][j] = I["tx"][i]
                    y_inner[9][j] = I["ty"][i]
                    j += 1
                
                # add final point
                x_inner[9][-1] = x_int
                y_inner[9][-1] = y_int
            elif np.min(y_inner[4]) <= y_int:
                # determine the x value of the y height on the front spar arc
                x_arc = info["front spar"][0] + ((info["sp"]+info["t"])**2. -\
                    (info["front spar"][1] - y_int)**2.)**0.5
                theta1 =  np.arctan2(y_int-info["front spar"][1],\
                    x_arc-info["front spar"][0])

                # create theta array, and consequently x and y arrays
                theta = np.linspace(theta0,theta1,info["np"])
                rad = np.full(theta.shape,info["t"]+info["sp"]/2.)
                x_inner[4] = np.flip(rad *np.cos(theta) +info["front spar"][0])
                y_inner[4] = np.flip(rad *np.sin(theta) +info["front spar"][1])

                # create bottom of servo line to arc
                x_inner[8] = np.array([info["srv"][0]+info["x+s"]+info["t"],\
                    x_inner[4][0]])
                y_inner[8] = np.array([info["srv"][1]-info["sys"]-info["t"],\
                    y_inner[4][0]])
            else:
                raise ValueError("Servo hole bottom is lower than front spar")
    
    # reorder the arrays
    if not info["msh"]:
        x_inner,y_inner = reorder(x_inner,y_inner,[6,4,5,3,2,1,0])
    else:
        if not info["mch"]:
            x_inner,y_inner = reorder(x_inner,y_inner,[6,7,8,9,10,4,5,3,2,1,0])
        else:
            if np.max(y_inner[4]) <= y_int:
                x_inner,y_inner = reorder(x_inner,y_inner,\
                    [6,7,8,9,4,5,3,2,1,0])
            else:
                x_inner,y_inner = reorder(x_inner,y_inner,[6,7,8,4,5,3,2,1,0])

    return x_inner,y_inner

# create outer airfoil shape
def outer(O,info):
    # initialize outer arrays
    num_arrays = 2 + 4*info["msh"] + 1*(info["msh"] and info["mch"])
    x_outer = np.zeros((num_arrays,),dtype=np.ndarray)
    y_outer = np.zeros((num_arrays,),dtype=np.ndarray)

    # initialize bottom of aft run
    x_outer[0] = np.zeros((O["bx"].shape[0]-info["outer"]["end"][2]+1,))
    y_outer[0] = np.zeros((O["bx"].shape[0]-info["outer"]["end"][2]+1,))

    # add first point
    x_outer[0][0] = info["outer"]["end"][0]
    y_outer[0][0] = info["outer"]["end"][1]

    # run through top and add points
    j = 1
    for i in range(info["outer"]["end"][2],O["bx"].shape[0]):
        x_outer[0][j] = O["bx"][i]
        y_outer[0][j] = O["by"][i]
        j += 1

    # split top run around depending of if making servo and if cord hole
    if not info["msh"]:
        # initialize around run
        x_outer[1] = np.zeros((O["tx"].shape[0]+info["outer"]["start"][2]))
        y_outer[1] = np.zeros((O["ty"].shape[0]+info["outer"]["start"][2]))

        # run through top and add points
        j = 0
        for i in range(O["tx"].shape[0]-1,-1,-1):
            x_outer[1][j] = O["tx"][i]
            y_outer[1][j] = O["ty"][i]
            j += 1
        
        # add bottom points
        for i in range(1,info["outer"]["start"][2]):
            x_outer[1][j] = O["bx"][i]
            y_outer[1][j] = O["by"][i]
            j += 1
        
        # add final point
        x_outer[1][-1] = info["outer"]["start"][0]
        y_outer[1][-1] = info["outer"]["start"][1]
    else:
        # determine where servo point is location
        x_sv_te = info["srv"][0] + info["x+s"]
        y_sv_te,p_sv_te,i_sv_te = get_y_loc(x_sv_te,O["tx"],O["ty"])

        # initialize around run
        x_outer[1] = np.zeros((O["tx"].shape[0]-i_sv_te+1))
        y_outer[1] = np.zeros((O["ty"].shape[0]-i_sv_te+1))

        # run through top and add points
        j = 0
        for i in range(O["tx"].shape[0]-1,i_sv_te-1,-1):
            x_outer[1][j] = O["tx"][i]
            y_outer[1][j] = O["ty"][i]
            j += 1
        
        # add final point
        x_outer[1][-1] = x_sv_te
        y_outer[1][-1] = y_sv_te

        # add drop servo side wall
        x_outer[2] = np.array([x_sv_te,x_sv_te])
        y_outer[2] = np.array([y_sv_te,info["srv"][1] - info["sys"]])

        # split results between cord hole and ont
        if not info["mch"]:
            # add servo bottom
            x_outer[3] = np.array([x_sv_te,info["srv"][0] - info["x-s"]])
            y_outer[3] = np.array([info["srv"][1]-info["sys"],\
                info["srv"][1]-info["sys"]])

            # find where servo LE meets outer skin
            x_sv_le = info["srv"][0] - info["x-s"]
            y_sv_le,p_sv_le,i_sv_le = get_y_loc(x_sv_le,O["tx"],O["ty"])

            # add servo side front wall
            x_outer[4] = np.array([x_sv_le,x_sv_le])
            y_outer[4] = np.array([info["srv"][1]-info["sys"],y_sv_le])

            # initialize around run
            x_outer[5] = np.zeros((1+i_sv_le+info["outer"]["start"][2]))
            y_outer[5] = np.zeros((1+i_sv_le+info["outer"]["start"][2]))

            # add first point
            x_outer[5][0] = x_sv_le
            y_outer[5][0] = y_sv_le

            # run through top and add points
            j = 1
            for i in range(i_sv_le-1,-1,-1):
                x_outer[5][j] = O["tx"][i]
                y_outer[5][j] = O["ty"][i]
                j += 1
            
            # add bottom points
            for i in range(1,info["outer"]["start"][2]):
                x_outer[5][j] = O["bx"][i]
                y_outer[5][j] = O["by"][i]
                j += 1
            
            # add final point
            x_outer[5][-1] = info["outer"]["start"][0]
            y_outer[5][-1] = info["outer"]["start"][1]
        else: # if adding cord hole
            # find the y value of the x location at the front spare center
            x_o_sp = info["front spar"][0]
            y_o_sp,p_o_sp,i_o_sp = get_y_loc(x_o_sp,O["tx"],O["ty"])

            # set x and y points a radius and skin thickness away
            x_sp_rad = x_o_sp + info["sp"] + info["t"]
            y_sp_rad = info["srv"][1]-info["sys"]

            # set array for bottom of servo
            x_outer[3] = np.array([x_sv_te,x_sp_rad])
            y_outer[3] = np.array([y_sp_rad,y_sp_rad])

            # set wall up
            x_outer[4] = np.array([x_sp_rad,x_sp_rad])
            y_outer[4] = np.array([y_sp_rad,y_o_sp])

            # set flat to le top
            x_outer[5] = np.array([x_sp_rad,x_o_sp])
            y_outer[5] = np.array([y_o_sp,y_o_sp])

            # initialize around run
            x_outer[6] = np.zeros((1+i_o_sp+info["outer"]["start"][2]))
            y_outer[6] = np.zeros((1+i_o_sp+info["outer"]["start"][2]))

            # add first point
            x_outer[6][0] = x_o_sp
            y_outer[6][0] = y_o_sp

            # run through top and add points
            j = 1
            for i in range(i_o_sp-1,-1,-1):
                x_outer[6][j] = O["tx"][i]
                y_outer[6][j] = O["ty"][i]
                j += 1
            
            # add bottom points
            for i in range(1,info["outer"]["start"][2]):
                x_outer[6][j] = O["bx"][i]
                y_outer[6][j] = O["by"][i]
                j += 1
            
            # add final point
            x_outer[6][-1] = info["outer"]["start"][0]
            y_outer[6][-1] = info["outer"]["start"][1]
    

    # reorder the arrays
    if not info["msh"]:
        x_outer,y_outer = reorder(x_outer,y_outer,[1,0])
    else:
        if not info["mch"]:
            x_outer,y_outer = reorder(x_outer,y_outer,[1,2,3,4,5,0])
        else:
            x_outer,y_outer = reorder(x_outer,y_outer,[1,2,3,4,5,6,0])

    # plt.plot(info["srv"][0]*info["c"],info["srv"][1]*info["c"],"ro")
    # x = info["srv"][0]; y = info["srv"][1]
    # xl= [x+info["x+s"],x+info["x+s"],x-info["x-s"],x-info["x-s"],x+info["x+s"]]
    # yl= [y-info["sys"],y+info["sys"],y+info["sys"],y-info["sys"],y-info["sys"]]
    # xl = [a*info["c"] for a in xl]; yl = [a*info["c"] for a in yl]
    # plt.plot(xl,yl,"b")

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

#  re-organize c17 shape for polygon writing if desired
def make_polygon(x,y,vals):
    # determine how many polygons to create
    n = x["hol"].shape[0] + 2 +(vals["num kinks"]-len(vals["kinks to remove"]))
    # initialize arrays
    x_svg = np.zeros((n,),dtype=np.ndarray)
    y_svg = np.zeros((n,),dtype=np.ndarray)

    # add in hol arrays
    j = 0
    for i in range(x["hol"].shape[0]):
        x_svg[j] = x["hol"][i][:-1]
        y_svg[j] = y["hol"][i][:-1]

        j += 1
    
    # add inner hole array
    x_svg[j] = np.array([]); y_svg[j] = np.array([])
    for i in range(x["inn"].shape[0]):
        x_svg[j] = np.append(x_svg[j],x["inn"][i][:-1])
        y_svg[j] = np.append(y_svg[j],y["inn"][i][:-1])
    j += 1

    # initialize arcs numbers to remove and total
    arcs = np.linspace(1,vals["num kinks"],vals["num kinks"]).astype(int) - 1
    arem = np.array(vals["kinks to remove"]) - 1
    arkp = np.array([a for a in arcs if a not in arem])
    
    # run through each part to be added and add it!
    rfli = vals["num kinks"]
    for i in range(n-j-1):
        # initialize arrays
        x_svg[j] = np.array([]); y_svg[j] = np.array([])

        # initailize arc keep number
        arkpi = arkp[-1-i]

        # run through the removed ones ones
        for k in range(rfli,arkpi,-1):
            x_svg[j] = np.append(x_svg[j],x["rfl"][0][k][:-1])
            y_svg[j] = np.append(y_svg[j],y["rfl"][0][k][:-1])
            if k > arkpi+1:
                x_svg[j] = np.append(x_svg[j],x["spr"][0][k-1][:-1])
                y_svg[j] = np.append(y_svg[j],y["spr"][0][k-1][:-1])
        
        # add side arc
        x_svg[j] = np.append(x_svg[j],x["spr"][1][arkpi][:-1])
        y_svg[j] = np.append(y_svg[j],y["spr"][1][arkpi][:-1])

        # initialize return towards te line
        for k in range(arkpi+1,rfli+1):
            x_svg[j] = np.append(x_svg[j],x["rfl"][1][k][:-1])
            y_svg[j] = np.append(y_svg[j],y["rfl"][1][k][:-1])
            if k < rfli:
                x_svg[j] = np.append(x_svg[j],x["spr"][1][k][:-1])
                y_svg[j] = np.append(y_svg[j],y["spr"][1][k][:-1])

        # add side arc, or end cap depending on if it is the last
        if rfli == vals["num kinks"]:
            x_svg[j] = np.append(x_svg[j],x["rfl"][1][0][:-1])
            y_svg[j] = np.append(y_svg[j],y["rfl"][1][0][:-1])
        else:
            x_svg[j] = np.append(x_svg[j],x["spr"][0][rfli][:-1])
            y_svg[j] = np.append(y_svg[j],y["spr"][0][rfli][:-1])

        # move to next index
        rfli = arkpi
        j += 1
    
    # add final full runaround
    # initialize arrays
    x_svg[j] = np.array([]); y_svg[j] = np.array([])
    # add outer
    for i in range(x["out"].shape[0]-1):
        x_svg[j] = np.append(x_svg[j],x["out"][i][:-1])
        y_svg[j] = np.append(y_svg[j],y["out"][i][:-1])

    # add mouth
    for i in range(x["mou"].shape[0]):
        x_svg[j] = np.append(x_svg[j],x["mou"][i][:-1])
        y_svg[j] = np.append(y_svg[j],y["mou"][i][:-1])

    # add arcs and rfl part
    if len(arkp) == 0:
        rfli = vals["num kinks"]
    else:
        rfli = arkp[0]

    # run through the removed ones ones
    for k in range(rfli+1):
        x_svg[j] = np.append(x_svg[j],np.flip(x["rfl"][0][k])[:-1])
        y_svg[j] = np.append(y_svg[j],np.flip(y["rfl"][0][k])[:-1])
        if k < rfli-1:
            x_svg[j] = np.append(x_svg[j],np.flip(x["spr"][0][k])[:-1])
            y_svg[j] = np.append(y_svg[j],np.flip(y["spr"][0][k])[:-1])
    
    # add side part
    if len(arkp) == 0:
        # add end
        x_svg[j] = np.append(x_svg[j],np.flip(x["rfl"][1][0])[:-1])
        y_svg[j] = np.append(y_svg[j],np.flip(y["rfl"][1][0])[:-1])
    else:
        # add side arc
        x_svg[j] = np.append(x_svg[j],np.flip(x["spr"][0][arkpi])[:-1])
        y_svg[j] = np.append(y_svg[j],np.flip(y["spr"][0][arkpi])[:-1])

    # run through the removed ones ones
    for k in range(rfli,0,-1):
        x_svg[j] = np.append(x_svg[j],np.flip(x["rfl"][1][k])[:-1])
        y_svg[j] = np.append(y_svg[j],np.flip(y["rfl"][1][k])[:-1])
        
        # add bottom arc flat piece
        x_svg[j] = np.append(x_svg[j],np.flip(x["spr"][1][k-1])[:-1])
        y_svg[j] = np.append(y_svg[j],np.flip(y["spr"][1][k-1])[:-1])
    
    # add tongue
    for i in range(x["ton"].shape[0]):
        x_svg[j] = np.append(x_svg[j],x["ton"][i][:-1])
        y_svg[j] = np.append(y_svg[j],y["ton"][i][:-1])
    
    # add last of outer
    x_svg[j] = np.append(x_svg[j],x["out"][-1][:-1])
    y_svg[j] = np.append(y_svg[j],y["out"][-1][:-1])

    return x_svg,y_svg

# make polyline arrays for svg file
def make_polyline(x,y):
    # determine number of arrays spots to make
    size = 0
    for group in x:
        shape = x[group].shape
        if group == "arc" or group == "rfl" or group == "spr":
            for i in range(shape[0]):
                for j in range(shape[1]):
                    if type(x[group][i][j]) == np.ndarray:
                        size += 1
        else:
            for i in range(shape[0]):
                if type(x[group][i]) == np.ndarray:
                    size += 1

    # initialize svg inputs
    x_svg = np.zeros((size,),dtype=np.ndarray)
    y_svg = np.zeros((size,),dtype=np.ndarray)
    index = 0

    # input svg points
    for group in x:
        shape = x[group].shape
        if group == "arc" or group == "rfl" or group == "spr":
            for i in range(shape[0]):
                for j in range(shape[1]):
                    if type(x[group][i][j]) == np.ndarray:
                        x_svg[index] = x[group][i][j]
                        y_svg[index] = y[group][i][j]
                        index += 1
        else:
            for i in range(shape[0]):
                if type(x[group][i]) == np.ndarray:
                    x_svg[index] = x[group][i]
                    y_svg[index] = y[group][i]
                    index += 1

    return x_svg,y_svg

# create svg file
def make_svg(x,y,vals,info,filename):
    # initialize new vals
    txt = vals["svg text"]
    svg_type = vals["svg type"]
    write_to_file = vals["write svg"]

    # make polygon shape if desired
    if svg_type == "polygon":
        x_svg,y_svg = make_polygon(x,y,vals)
        # or make polyline
    elif svg_type == "polyline":
        x_svg,y_svg = make_polyline(x,y)
    else:
        raise ValueError("SVG file type must be 'polygon' or 'polyline'")

    # create text dictionary
    text_dict = {}
    text_dict["text"] = txt

    # pass in x y of start location
    text_dict["start"] = info["txt start"]

    # write svg file
    if write_to_file:
        svg.svg(filename,x_svg,y_svg,geometry=svg_type,text=text_dict)

    return x_svg,y_svg,text_dict

# deflect airfoil function
def deflection(x,y,O,info):
    # run through each group
    for group in x:
        x[group],y[group] = deflect.main(O,x[group],y[group],info)
    
    return

# shift to quarter chord
def shift(x,y,O):
    # determine shift values
    x_shift = 0.25
    y_shift = O["c"](x_shift)

    # shift each group
    for group in x:
        x[group] -= x_shift
        y[group] -= y_shift
    
    return

# resize by a chord value
def resize(x,y,c):
    # run through each group and resize
    for group in x:
        x[group] *= c
        y[group] *= c
    
    return

# points maker
def points(x,y,filename):
    # initialize points infos
    t = {}
    t["mou"] = [[0,True],[1,False],[2,False],[3,True],[4,False],[5,True]]
    t["ton"] = [[0,False],[1,False],[2,True],[4,True],[5,True],[3,True],\
        [6,True]]
    t["spr"] = [[0,True],[1,True],[2,True],[3,True]]
    t["hol"] = [[0,True]]
    t["inn"] = [[2,False],[1,False],[0,False]]
    t["out"] = [[1,True],[0,True]]
    t["order"] = ["mou","inn","ton","out"]
    t["new order"] = ["spr","hol"]

    # run through each and add to x y
    xp = []; yp = []
    for i in range(len(t["order"])):
        group = t["order"][i]
        for j in range(len(t[group])):
            sub = t[group][j]
            if sub[1]:
                for k in range(0,x[group][sub[0]].shape[0]-1):
                    xp.append(x[group][sub[0]][k])
                    yp.append(y[group][sub[0]][k])
            else:
                for k in range(x[group][sub[0]].shape[0]-1,0,-1):
                    xp.append(x[group][sub[0]][k])
                    yp.append(y[group][sub[0]][k])
    
    # add final point of outer
    xp.append(x["out"][0][-1])
    yp.append(y["out"][0][-1])

    # add spar and hole
    for i in range(len(t["new order"])):
        group = t["new order"][i]
        for j in range(len(t[group])):
            sub = t[group][j]
            for k in range(0,x[group][sub[0]].shape[0]):
                xp.append(x[group][sub[0]][k])
                yp.append(y[group][sub[0]][k])
    
    # make numpy arrays
    xp = np.array(xp); yp = np.array(yp)

    # create a csv file
    with open(filename+".txt","w") as f:
        for i in range(len(xp)):
            f.write("{:<10.8f}\t{:<10.8f}\n".format(xp[i],yp[i]))
        f.close()
    
    return xp,yp

def main(jsonfile,run="C17",foil="hallelujah"):
    # import json file if a dictionary or json file, OR raise error
    if type(jsonfile) == dict:
        vals = jsonfile
    elif os.path.isfile(jsonfile):
        vals = read_json(jsonfile)
    else:
        raise ValueError("C17 input must be dictionary or file path")

    # simplify dictionary
    info = simplify(vals,run)

    # create airfoil shape
    if type(foil) != adb.airfoil.Airfoil:
        add = {}
        add["geometry"] = {}
        add["geometry"]["type"] = "outline_points"
        add["geometry"]["outline_points"] = \
            vals[run]["airfoil dat file location"]

        # initialize airfoil database
        foil = adb.Airfoil("foil",add)
    
    coords = foil.get_outline_points()
    x_af = coords[:,0]; y_af = coords[:,1]
    # print(x_af)

    # if enthicken, do now
    if "enthicken" in vals[run]:
        y_af *= vals[run]["enthicken"]

    # create inner airfoil
    polyx,polyy = inner(x_af,y_af,info["t"])

    # initialize interpolation dictionary
    I = {}

    # for further use, if control hole, determine max diameter
    if info["ch"]:
        max_diam = max([info["td"],info["cwd"]])
    else: # use td diameter
        max_diam = info["cwd"] - info["t"]

    # create further thicknesses inner airfoils, as well as interpolations
    # create interpolation of inner surface
    I["I"] = interpolate(polyx,polyy,make_camberline=True)
    # create interpolation of outer surface
    I["O"] = interpolate(x_af,y_af,make_camberline=True)
    # outer tongue
    out_tx,out_ty = inner(x_af,y_af,   info["t"] + info["mt"])
    I["ot"] = interpolate(out_tx,out_ty)
    # inner tongue
    inn_tx,inn_ty = inner(x_af,y_af,2.*info["t"] + info["mt"])
    I["it"] = interpolate(inn_tx,inn_ty)
    # outer mouth
    out_mx,out_my = inner(x_af,y_af,3.*info["t"] + info["mt"] + max_diam)
    I["om"] = interpolate(out_mx,out_my)
    # outer mouth
    inn_mx,inn_my = inner(x_af,y_af,4.*info["t"] + info["mt"] + max_diam)
    I["im"] = interpolate(inn_mx,inn_my)
    # tip tongue top
    ttt_tx,ttt_ty = inner(x_af,y_af,3.*info["t"] + max_diam)
    I["ttt"] = interpolate(ttt_tx,ttt_ty)

    # if shift desired, shift hole
    if info["ts"]:
        find_wall(I["O"],foil,info)

    # plt.axis("equal")
    # plt.plot(x_af,y_af)
    # plt.plot(polyx,polyy)
    # plt.plot(out_tx,out_ty)
    # plt.plot(inn_tx,inn_ty)
    # plt.plot(out_mx,out_my)
    # plt.plot(inn_mx,inn_my)
    # plt.plot(ttt_tx,ttt_ty)
    # plt.show()

    # initialize x and y dictionaries
    x = {}; y = {}

    # determine where the arcs should start-ish
    I["I"]["arcdx"] = (info["tw"] - info["t"] - info["f"]) / (info["nk"] + 1)

    # create outer mouth structure
    x["mou"],y["mou"] = make_mouth(I["O"],I["I"],I["om"],info)

    # # create arcs
    # x["arc"],y["arc"] = arcs(I["I"],info)

    # create springs
    x["spr"],y["spr"] = kincs(I["I"],info)

    # create roofs and floors
    x["rfl"],y["rfl"] = roofs_n_floors(I["I"],info,x["mou"][5][-1])

    # create tongue tip
    x["ton"],y["ton"] = tongue_tip(I["O"],I["I"],I["ot"],I["it"],I["ttt"],info)

    # # remove arcs desired
    # remove_arcs(x,y,info)

    # remove kinks desired
    remove_springs(x,y,info)

    # create TE hole
    x["hol"],y["hol"] = te_hole(I["O"],I["I"],I["it"],foil,info)

    # create inner run around spar holes in wing nose
    x["inn"],y["inn"] = front_inner(I["I"],I["im"],info)

    # create outer run. start and end values are where break starts and ends
    x["out"],y["out"] = outer(I["O"],info)


    # initialize df array
    if not info["dty"] == "none" and type(info["da"]) == list:
        df = info["da"]
    else:
        df = [1.0]
        info["dty"] = "none"

    # run through each deflection angle
    for k in range(len(df)):
        if len(df) != 1:
            print("da={:>5.1f}, {:>5.1f}% complete".format(df[k],\
                (k+1)/len(df)*100.))
        # set deflection angle
        info["da"] = df[k]
        # reset x y set
        u = {}; v = {}
        for group in x:
            u[group] = x[group] * 1.; v[group] = y[group] * 1.
        
        # if desired to deflect, do so
        if not info["dty"] == "none":
            # run through each set and deflect
            deflection(u,v,I["O"],info)

        # if shift origin to quarter chord
        if "shift to c/4" in vals[run] and vals[run]["shift to c/4"]:
            # shift
            shift(u,v,I["O"])

        # resize chord of each segment
        resize(u,v,info["c"])

        # # create dxf file
        # filename = vals[run]["dxf file path"]+vals[run]["dxf file name"]
        # xvals,yvals,zvals = make_dxf(u,v,filename,vals[run]["write dxf"])

        # create svg file
        filename = vals[run]["svg file path"]+vals[run]["svg file name"]
        xsvg,ysvg,text = make_svg(u,v,vals[run],info,filename)
        
        # # if points file desired, do so
        # fname = vals[run]["points file path"]+vals[run]["points file name"]
        # if vals[run]["write points"]:
        #     points(u,v,fname)

        if vals[run]["show plot"]:
            # initialize colors and labels
            g = {}
            g["spr c"] = "#ff8c00"; g["spr l"] = "Spring"
            g["rfl c"] = "#3cb371"; g["rfl l"] = "roofsNfloors"
            g["mou c"] = "#0073ff"; g["mou l"] = "Mouth"
            g["ton c"] = "#ff0000"; g["ton l"] = "Tongue"
            g["inn c"] = "#aa6c39"; g["inn l"] = "Spar"
            g["hol c"] = "#6a0dad"; g["hol l"] = "Holes"
            g["out c"] =       "k"; g["out l"] = "Outer"

            cols = ["#FF0000","#FFA500","#FFFF00","#008000","#ADD8E6",
                    "#0000FF","#800080","#FFC0CB","0.6","#A52A2A","k",
                    "#008080"]

            # plot to be saved items
            for group in u:
                # color
                if vals[run]["all black"]:
                    c = "k"
                else:
                    c = g[group+" c"]
                
                # if roofs and floors or arcs, do i and j
                if group == "arc" or group == "spr": 
                    for i in range(info["nk"]): # 2): # 
                        # 0 -> 4 from x = 0 -> 1.0
                        for j in range(2):
                            # 0 - LE facing arc, 1 - TE facing arc
                            if i==0 and j==0:
                                l = g[group+" l"]
                            else:
                                l = ""
                            plt.plot(u[group][j][i],v[group][j][i],c=c,label=l)
                elif group == "rfl":
                    for i in range(u["rfl"].shape[0]):
                        for j in range(u["rfl"].shape[1]):
                            if i==0 and j==0:
                                plt.plot(u[group][i][j],v[group][i][j],c,\
                                    label=g[group+" l"])
                            else:
                                plt.plot(u[group][i][j],v[group][i][j],c)
                else:
                    for i in range(u[group].shape[0]):
                        if i == 0:
                            l = g[group+" l"]
                        else:
                            l = ""
                        # c = cols[i]
                        plt.plot(u[group][i],v[group][i],c,label=l)
                        # plt.plot(u[group][i][-2],v[group][i][-2],"k+")

    # for i in range(x_spr.shape[0]):
    #     for j in range(x_spr[0].shape[0]):
    #         plt.plot(x_spr[i][j]*info["c"],y_spr[i][j]*info["c"],"#008080")
    
    
    # plt.plot(np.array([info["tw"],info["tw"]])*info["c"],\
    #     np.array([0.1,-0.1])*info["c"])
    
    if vals[run]["show plot"]:
        if vals[run]["show legend"]:
            plt.legend()
        plt.xlabel("x/c")
        plt.ylabel("y/c")
        plt.axis("equal")
        if "close after" in vals[run]:
            plt.show(block=False)
            plt.pause(vals[run]["close after"])
            plt.close()
        else:
            plt.show()

    return xsvg, ysvg, ysvg*0., text
