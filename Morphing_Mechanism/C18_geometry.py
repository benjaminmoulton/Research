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
    dicti["ktt"] = vals[run]["kink thickness type"]
    
    dicti["Nres"]  = vals[run]["section resolution"]
    dicti["t"]  = vals[run]["shell thickness [mm]"] / dicti["c"] / in2mm
    dicti["f"] = vals[run]["flex start [x/c]"]
    dicti["to"] = dicti["f"] - np.abs(vals[run]["tongue start [in]"]) \
        / dicti["c"]
    dicti["mo"] = dicti["f"] - np.abs(vals[run]["mouth start [in]"]) \
        / dicti["c"]
    if dicti["to"] < dicti["mo"]:
        raise ValueError("Tongue cannot extend beyond mouth start")
    dicti["tw"] = vals[run]["TE wall start [x/c]"]
    dicti["mt"] = vals[run]["mouth clearance [mm]"] / dicti["c"] / in2mm
    dicti["fi"] = vals[run]["fillet side length [mm]"] / dicti["c"] / in2mm
    dicti["da"] = vals[run]["deflection angle [deg]"]
    dicti["dty"] = vals[run]["deflection type"]
    dicti["th"] = vals[run]["hinge at top"]
    dicti["sk"] = vals[run]["skip wing box"]

    # I beam values
    dicti["ith"] = vals[run]["I beam"]["thickness [mm]"] / dicti["c"] / in2mm
    dicti["ifw"] = vals[run]["I beam"]["flange width [x/c]"]
    dicti["hol"] = vals[run]["I beam"]["hollow [mm]"] / dicti["c"] / in2mm
    dicti["mhl"] = vals[run]["I beam"]["make hollow"]
    dicti["leth"] = vals[run]["leading edge thickness[mm]"] / dicti["c"] /in2mm

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
        # print(repeats); print(insct)
        removed = 0
        for index in repeats:
            if not index[0] == k:
                # print(index[0],k)
                insct = np.delete(insct,index[0]+removed,axis=0)
                removed -= 1

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
    if xy["dc"][it] == 0.0:
        xval = xp
    else:
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
def kincs(xy,info,shift=0.0):
    # initialize arcs arrays first index is 0 arc 1 skin
    x_kin = np.zeros((2,info["nk"]),dtype=np.ndarray)
    y_kin = np.zeros((2,info["nk"]),dtype=np.ndarray)
    # initialize kincs top and bottom coordinates arrays
    # [i][ ] 0 - arc on left (LE), 1 - arc on right (TE)
    # [ ][i] 0 -> 1 from LE to TE
    
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
        x_kin[0][i] = np.zeros((info["nkp"],))
        y_kin[0][i] = np.zeros((info["nkp"],))
        x_kin[1][i] = np.zeros((info["nkp"],))
        y_kin[1][i] = np.zeros((info["nkp"],))

        # set up x location
        if "k0" in info and "k1" in info:
            start = info["k0"]
        else:
            start = info["f"]
        x_loc = start + xy["arcdx"] * (i+1) + shift
        
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
        if info["ktt"] == "fixed" and len(info["kth"]) == 1:
            kink_thickness = info["kth"][0]
        elif info["ktt"] == "gradient" and len(info["kth"]) == 2:
            # caclulate fraction length
            frac_l = (l_loc-l_flex) / (l_end-l_flex)

            # calculate kink_thickness
            kink_thickness = info["kth"][0] + frac_l * \
                (info["kth"][1] - info["kth"][0])
        elif info["ktt"] == "specific" and len(info["kth"]) == info["nk"]:
            kink_thickness = info["kth"][i]
        else:
            raise ValueError("length of kink thickness array must be " + \
                "appropriate for type, {1 for 'fixed', 2 for 'gradient'," \
                    " num kinks for 'specific'")
            

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
            xb,yb = Bezier(xbpt,ybpt,k=3,n=info["nkp"])

            # add to array
            midx = np.append(midx,xb[:-1])
            midy = np.append(midy,yb[:-1])

        # add final point
        midx = np.append(midx,xb[-1])
        midy = np.append(midy,yb[-1])

        # create offset kink
        x_kin[0][i],y_kin[0][i] = offset(midx,midy,-kink_thickness/2.)
        x_kin[1][i],y_kin[1][i] = offset(midx,midy,kink_thickness/2.)

        # ##############################################333
        # if i==0:
        #     print()
        #     print("after offset")
        #     print("x kinc shape",0,i,"=",x_kin[0][i].shape)
        #     print("x kinc shape",1,i,"=",x_kin[1][i].shape)
        # ##########################################333333333333333333
        
        # find where the intersections occur
        it0 = intersections(x_kin[0][i],y_kin[0][i])
        it1 = intersections(x_kin[1][i],y_kin[1][i])

        # remove them
        x_kin[0][i],y_kin[0][i] = de_intersect_kink(x_kin[0][i],y_kin[0][i],\
            it0)
        x_kin[1][i],y_kin[1][i] = de_intersect_kink(x_kin[1][i],y_kin[1][i],\
            it1)

        # ##############################################333
        # if i==0:
        #     print("after intersect rem")
        #     print("x kinc shape",0,i,"=",x_kin[0][i].shape)
        #     print("x kinc shape",1,i,"=",x_kin[1][i].shape)
        # ##########################################333333333333333333

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
            xa = np.array([x_kin[a][i][ai[-1]]+shifty, x_kin[a][i][ai[2]],
                            x_kin[a][i][ai[1]], x_kin[a][i][ai[-1]]])
            ya = np.array([xy["t"](x_kin[a][i][ai[-1]]+shifty),
                            xy["t"](x_kin[a][i][ai[2]]),
                            y_kin[a][i][ai[1]], y_kin[a][i][ai[-1]]])
            xb = np.array([x_kin[b][i][bi[-1]]-shifty, x_kin[b][i][bi[2]],
                            x_kin[b][i][bi[1]], x_kin[b][i][bi[-1]]])
            yb = np.array([xy["t"](x_kin[b][i][bi[-1]]-shifty),
                            xy["t"](x_kin[b][i][bi[2]]),
                            y_kin[b][i][bi[1]], y_kin[b][i][bi[-1]]])
            xazt,yazt = Bezier(xa,ya,n=10)
            xbzt,ybzt = Bezier(xb,yb,n=10)

            # create final array
            x_kin[a][i] = np.insert(x_kin[a][i][ai[-1]+1:],0,xazt)
            y_kin[a][i] = np.insert(y_kin[a][i][ai[-1]+1:],0,yazt) 
            x_kin[b][i] = np.insert(x_kin[b][i][bi[-1]+1:],0,xbzt)
            y_kin[b][i] = np.insert(y_kin[b][i][bi[-1]+1:],0,ybzt)

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
            xa = np.array([x_kin[a][i][ai[-1]], x_kin[a][i][ai[1]],
                            x_kin[a][i][ai[2]], x_kin[a][i][ai[-1]]+shifty])
            ya = np.array([y_kin[a][i][ai[-1]], y_kin[a][i][ai[1]],
                            xy["b"](x_kin[a][i][ai[2]]),
                            xy["b"](x_kin[a][i][ai[-1]]+shifty)])
            xb = np.array([x_kin[b][i][bi[-1]], x_kin[b][i][bi[1]],
                            x_kin[b][i][bi[2]], x_kin[b][i][bi[-1]]-shifty])
            yb = np.array([y_kin[b][i][bi[-1]], y_kin[b][i][bi[1]],
                            xy["b"](x_kin[b][i][bi[2]]),
                            xy["b"](x_kin[b][i][bi[-1]]-shifty)])
                            
            xazt,yazt = Bezier(xa,ya,n=10)
            xbzt,ybzt = Bezier(xb,yb,n=10)

            # create final array
            x_kin[a][i] = np.append(x_kin[a][i][:ai[-1]],xazt)
            y_kin[a][i] = np.append(y_kin[a][i][:ai[-1]],yazt) 
            x_kin[b][i] = np.append(x_kin[b][i][:bi[-1]],xbzt)
            y_kin[b][i] = np.append(y_kin[b][i][:bi[-1]],ybzt)

        # ##############################################333
        # if i==0:
        #     print("after corner removal")
        #     print("x kinc shape",0,i,"=",x_kin[0][i].shape)
        #     print("x kinc shape",1,i,"=",x_kin[1][i].shape)
        #     print()
        # ##########################################333333333333333333
        
        # flip
        x_kin[0][i] = np.flip(x_kin[0][i]); y_kin[0][i] = np.flip(y_kin[0][i])
        x_kin[1][i] = np.flip(x_kin[1][i]); y_kin[1][i] = np.flip(y_kin[1][i])

        # # save top and bottom x y values to the points arrays in info
        # [i][ ][ ] 0 - top, 1 - bottom
        # [ ][i][ ] 0 -> 1 from LE to TE
        # [ ][ ][i] 0 - arc on left (LE), 1 - arc on right (TE)
        top = -1; bot = 0
        info["xpt"][0][i][0] = x_kin[0][i][top]
        info["xpt"][1][i][0] = x_kin[0][i][bot]
        info["xpt"][0][i][1] = x_kin[1][i][top]
        info["xpt"][1][i][1] = x_kin[1][i][bot]
        info["ypt"][0][i][0] = y_kin[0][i][top]
        info["ypt"][1][i][0] = y_kin[0][i][bot]
        info["ypt"][0][i][1] = y_kin[1][i][top]
        info["ypt"][1][i][1] = y_kin[1][i][bot]

    return x_kin,y_kin

# remove any arcs desired
def remove_kincs(x,y,info):
    # check to make sure all arcs to remove are available
    for i in range(len(info["ktr"])):
        if info["ktr"][i] < 1:
            raise ValueError("Kink to remove value {}".format(info["ktr"][i])+\
                " below allowable values (1)")
        elif info["ktr"][i] > info["nk"]:
            raise ValueError("Kink to remove value {}".format(info["ktr"][i])+\
                " above allowable value {}".format(info["nk"]))
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
        x["kin"][0][j] = np.array([x["rfl"][0][j+1][0],x["rfl"][0][j][-1]])
        y["kin"][0][j] = np.array([y["rfl"][0][j+1][0],y["rfl"][0][j][-1]])
        # create line at bottom to remove arcs
        if j == 0:
            x["kin"][1][j] = np.flip(np.array([x["ton"][0][-1],\
                x["rfl"][1][j+1][0]]))
            y["kin"][1][j] = np.flip(np.array([y["ton"][0][-1],\
                y["rfl"][1][j+1][0]]))
        else:
            x["kin"][1][j] = np.flip(np.array([x["rfl"][1][j][-1],\
                x["rfl"][1][j+1][0]]))
            y["kin"][1][j] = np.flip(np.array([y["rfl"][1][j][-1],\
                y["rfl"][1][j+1][0]]))

    return

# create outer mouth structure
def make_mouth(O,I,om,iui,ili,info):
    # initialize mouth array
    x_mouth = np.zeros((11,),dtype=np.ndarray)
    y_mouth = np.zeros((11,),dtype=np.ndarray)
    
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
    x_mouth[3] = np.linspace(info["mo"],info["f"]+info["ifw"]-0.5*info["ith"],\
        num=20)
    y_mouth[3] = om["b"](x_mouth[3])
    # find where flex starts as an interpolated point outer mouth
    x_fx_om = info["f"] + info["ifw"] - 0.5 * info["ith"]
    y_fx_om,p_fx_om,i_fx_om = get_y_loc(x_fx_om,om["bx"],om["by"])
    y_mouth[3][-1] = y_fx_om
    y_mouth[3][0] = y_ms_om

    # create circle around for ibeam bottom right corner
    r_corner = 0.5 * info["ith"]
    t_corner = np.linspace(-np.pi/2.,np.pi/2.,num=20)
    x_mouth[4] = r_corner * np.cos(t_corner) + x_fx_om
    y_mouth[4] = r_corner * np.sin(t_corner) + y_fx_om + 0.5 * info["ith"]

    # create inner lower i beam
    x_mouth[5] = np.linspace(info["f"] + info["ifw"] - 0.5 * info["ith"],\
        info["f"] + info["fi"],num=20)
    y_mouth[5] = ili["b"](x_mouth[5])
    # add in first point to ensure it stays connected
    x_mouth[5][0] = x_mouth[4][-1]; y_mouth[5][0] = y_mouth[4][-1]
    
    # create bezier curve to upper lip line
    # find where flex starts as an interpolated point outer mouth
    x_fx_ili = x_mouth[5][-1]
    y_fx_ili = y_mouth[5][-1]

    # find where the curve should end
    x_fx_cur = info["f"]
    y_fx_cur,_,_ = get_y_loc(x_fx_cur,ili["bx"],ili["by"])
    # add bezier curve shift
    y_fx_cur += info["fi"]

    # create curve input arrays
    x_points = np.array([x_fx_ili,(x_fx_ili+x_fx_cur)/2.,x_fx_cur,x_fx_cur])
    y_points = np.array([y_fx_ili,(y_fx_ili+y_fx_cur-info["fi"])/2.,\
        (y_fx_ili+y_fx_cur)/2.,y_fx_cur])

    # create fillet using a bezier curve maker
    x_mouth[6],y_mouth[6] = Bezier_Cubic(x_points,y_points)

    # create "upper lip" line
    # find where the flex starts as an interpolated point top inner
    x_fx_ti = info["f"]
    y_fx_ti,p_fx_ti,i_fx_ti = get_y_loc(x_fx_ti,iui["tx"],iui["ty"])
    # shift down due to bezier curve
    y_fx_ti -= info["fi"]

    # create upper lip
    x_mouth[7] = np.array([x_fx_ti,x_fx_cur])
    y_mouth[7] = np.array([y_fx_ti,y_fx_cur])

    # create bezier curve to nose bridge
    # find where the curve should end
    x_fx_ns = info["f"] + info["fi"]
    y_fx_ns,_,_ = get_y_loc(x_fx_cur,iui["tx"],iui["ty"])

    # create curve input arrays
    x_points = np.array([x_fx_ti,x_fx_ti,(x_fx_ti+x_fx_ns)/2.,x_fx_ns])
    y_points = np.array([y_fx_ti,(y_fx_ti+y_fx_ns)/2.,\
        (y_fx_ti+info["fi"]+y_fx_ns)/2.,y_fx_ns])

    # create fillet using a bezier curve maker
    x_mouth[8],y_mouth[8] = Bezier_Cubic(x_points,y_points)

    # create bottom of nose portion
    x_mouth[9] = np.linspace(x_fx_ns,\
        info["f"] + info["ifw"] - 0.5*info["ith"],num=20)
    y_mouth[9] = iui["t"](x_mouth[9])
    # fix first value
    y_mouth[9][0] = y_fx_ns

    # initialize values for end of bottom of nose portion
    x_ns_end = x_mouth[9][-1]; y_ns_end = y_mouth[9][-1]

    # create fillet from upper lip to inner airfoil
    # find where the fillet ends as an interpolated point top inner
    x_fi_ti = x_ns_end + info["ith"]
    y_fi_ti,_,_ = get_y_loc(x_fi_ti,I["tx"],I["ty"])

    # initialze point arrays to get the bezier curve
    x_points = np.array([x_ns_end,(x_ns_end+x_fi_ti)/2.,\
        (x_ns_end+x_fi_ti)/2.,x_fi_ti])
    y_points = np.array([y_ns_end,y_ns_end,y_fi_ti,y_fi_ti])

    # create fillet using a bezier curve maker
    x_mouth[10],y_mouth[10] = Bezier_Cubic(x_points,y_points)

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
    y_end,__,__ = get_y_loc(info["xpt"][1][0][0]-.03,I["bx"],I["by"])

    # create bezier curve
    xpts = np.array([x_fx_it,x_fx_it+.03,info["xpt"][1][0][0]-.03,\
        info["xpt"][1][0][0]])
    ypts = np.array([y_fx_it,y_start,y_end,\
        info["ypt"][1][0][0]])
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
    x_j_o = info["xpt"][1][0][0]
    y_j_o,p_j_o,i_j_o = get_y_loc(x_j_o,O["bx"],O["by"])

    # save to tongue tip
    x_tongue[5] = np.array([x_tongue[3][-1],x_j_o])
    y_tongue[5] = np.array([y_tongue[3][-1],y_j_o])
    # save value so info to be used in a later function
    info["outer"]["end"] = [x_j_o,y_j_o,i_j_o]

    return x_tongue,y_tongue

# create roofs and floors
def roofs_n_floors(I,info,x_fillet,y_tri_top,y_tri_bot):
    # initialize mouth array
    # [i][ ] 0 - roof, 1 - floor
    # [ ][i] 0 -> 1, LE -> TE
    # NOTE [1][0] will remain empty
    x_rofl = np.zeros((2,info["nk"]+1),dtype=np.ndarray)
    y_rofl = np.zeros((2,info["nk"]+1),dtype=np.ndarray)

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
                if j == info["nk"]:
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
def wing_box(Interps,info):
    # draw out Interps pieces
    I = Interps["I"]; im = Interps["im"]; iui = Interps["iui"]
    ili = Interps["ili"]; leth = Interps["leth"]

    # initialize wing box info
    n_box = 33 + 44 * info["mhl"]
    x_box = np.zeros((n_box,),dtype=np.ndarray)
    y_box = np.zeros((n_box,),dtype=np.ndarray)

    # ### determine the lines which form the x of the wing box
    x_line_up = [info["mo"],info["f"]-info["ith"]]
    y_line_up = [im["b"](info["mo"]),iui["t"](info["f"]-info["ith"])]
    x_line_down = [info["mo"],info["f"]-info["ith"]]
    y_line_down = [I["t"](info["mo"]),ili["b"](info["f"]-info["ith"])]

    # determine slopes and intercepts
    m_up = (y_line_up[1] - y_line_up[0]) / (x_line_up[1] - x_line_up[0])
    m_down = (y_line_down[1] - y_line_down[0]) / \
        (x_line_down[1] - x_line_down[0])
    b_up = y_line_up[0] - m_up * x_line_up[0]
    b_down = y_line_down[0] - m_down * x_line_down[0]

    # determine the y shift to the intercept for upper and lower shift of lines
    # the _# (u or l) indicates upper or lower shifted intercept
    y_shift_up = ((info["t"] / 2.)**2. / (1. + (m_up)**2.))**0.5
    y_shift_down = ((info["t"] / 2.)**2. / (1. + (m_down)**2.))**0.5
    b_up_u = b_up + y_shift_up
    b_up_l = b_up - y_shift_up
    b_down_u = b_down + y_shift_down
    b_down_l = b_down - y_shift_down

    # define triangle as north, east, south, and west.
    # determine non-center points of each triangle
    # n(orth) and s(outh) have l(eft) and r(ight)
    # e(ast) and w(est) have u(p) and d(own)
    # each as a c(enter) point
    x_n_l,y_n_l,i_n_l = crossing(I,m_down,b_down_u)
    x_n_r,y_n_r,_ = crossing(iui,m_up,b_up_u)
    x_e_u,y_e_u = info["f"]-info["ith"],m_up*(info["f"]-info["ith"])+b_up_l
    x_e_d,y_e_d = info["f"]-info["ith"],m_down*(info["f"]-info["ith"])+b_down_u
    x_s_l,y_s_l,_ = crossing(im,m_up,b_up_l,is_top=False)
    x_s_r,y_s_r,_ = crossing(ili,m_down,b_down_l,is_top=False)
    x_w_u,y_w_u = info["mo"], m_down * info["mo"] + b_down_l
    x_w_d,y_w_d = info["mo"],   m_up * info["mo"] + b_up_u
    # define triangle inward points
    x_n_c = (b_down_u - b_up_u) / (m_up - m_down)
    y_n_c = m_up * x_n_c + b_up_u
    x_e_c = (b_down_u - b_up_l) / (m_up - m_down)
    y_e_c = m_up * x_e_c + b_up_l
    x_s_c = (b_down_l - b_up_l) / (m_up - m_down)
    y_s_c = m_up * x_s_c + b_up_l
    x_w_c = (b_down_l - b_up_u) / (m_up - m_down)
    y_w_c = m_up * x_w_c + b_up_u

    # initialize k counter value
    k = 0

    # cycle through n e s w triangles, cw so 0th index can be gc

    # north triangle, start with "outer" line
    x_box[k] = np.linspace(x_n_l,info["f"]-1.5*info["ith"]-info["ifw"],num=10)
    y_box[k] = I["t"](x_box[k]); y_box[k][0] = y_n_l
    k += 1
    # north triangle, bezier curve to ibeam
    x_start = x_box[k-1][-1]; y_start = y_box[k-1][-1]
    x_end = info["f"]-0.5*info["ith"]-info["ifw"]; y_end = iui["t"](x_end)
    x_mid = (x_start + x_end) / 2.0
    ###### define bezier curve points and create curve
    x_bez = np.array([x_start,x_mid,x_mid,x_end])
    y_bez = np.array([y_start,y_start,y_end,y_end])
    x_box[k],y_box[k] = Bezier(x_bez,y_bez)
    k += 1
    # north triangle i beam left upper inner flange line
    x_box[k] = np.linspace(x_box[k-1][-1],x_n_r,num=10)
    y_box[k] = iui["t"](x_box[k])
    y_box[k][0] = y_box[k-1][-1]; y_box[k][-1] = y_n_r
    k += 1
    # north triangle right side
    x_box[k] = np.array([x_n_r,x_n_c])
    y_box[k] = np.array([y_n_r,y_n_c])
    k += 1
    # north triangle left side
    x_box[k] = np.array([x_n_c,x_n_l])
    y_box[k] = np.array([y_n_c,y_n_l])
    k += 1

    # east triangle, "outer" line (right side)
    x_box[k] = np.array([x_e_u,x_e_d])
    y_box[k] = np.array([y_e_u,y_e_d])
    k += 1
    # east triangle, lower line
    x_box[k] = np.array([x_e_d,x_e_c])
    y_box[k] = np.array([y_e_d,y_e_c])
    k += 1
    # east triangle, upper line
    x_box[k] = np.array([x_e_c,x_e_u])
    y_box[k] = np.array([y_e_c,y_e_u])
    k += 1

    # south triangle start on outer line, with forward running
    x_box[k] = np.linspace(x_s_r,info["f"]-0.5*info["ith"]-info["ifw"],num=10)
    y_box[k] = ili["b"](x_box[k])
    k += 1
    # south triangle bezier curve to ibeam
    x_start = x_box[k-1][-1]; y_start = y_box[k-1][-1]
    x_end = info["f"]-1.5*info["ith"]-info["ifw"]; y_end = im["b"](x_end)
    x_mid = (x_start + x_end) / 2.0
    ###### define bezier curve points and create curve
    x_bez = np.array([x_start,x_mid,x_mid,x_end])
    y_bez = np.array([y_start,y_start,y_end,y_end])
    x_box[k],y_box[k] = Bezier(x_bez,y_bez)
    k += 1
    # south triangle i beam left upper inner flange line
    x_box[k] = np.linspace(x_box[k-1][-1],x_s_l,num=10)
    y_box[k] = im["b"](x_box[k])
    y_box[k][0] = y_box[k-1][-1]; y_box[k][-1] = y_s_l
    k += 1
    # south triangle left upper side
    x_box[k] = np.array([x_s_l,x_s_c])
    y_box[k] = np.array([y_s_l,y_s_c])
    k += 1
    # south triangle right upper side
    x_box[k] = np.array([x_s_c,x_s_r])
    y_box[k] = np.array([y_s_c,y_s_r])
    k += 1

    # west triangle, "outer" line (left side)
    x_box[k] = np.array([x_w_d,x_w_u])
    y_box[k] = np.array([y_w_d,y_w_u])
    k += 1
    # west triangle, upper line
    x_box[k] = np.array([x_w_u,x_w_c])
    y_box[k] = np.array([y_w_u,y_w_c])
    k += 1
    # west triangle, lower line
    x_box[k] = np.array([x_w_c,x_w_d])
    y_box[k] = np.array([y_w_c,y_w_d])
    k += 1

    # create leading edge, with optional thickness
    # determine the upper and lower aft points of this piece
    x_upp = info["mo"]-info["t"]; x_low = x_upp
    y_upp = leth["t"](x_upp); y_low = leth["b"](x_low)
    # determine upper and lower aft first section points of this piece
    i_first = int(0.075 * info["Nres"])
    x_u_first = leth["tx"][i_first]
    # y_u_first = leth["ty"][i_first]
    x_l_first = leth["bx"][i_first]
    # y_l_first = leth["by"][i_first]

    # create aft line
    x_box[k] = np.array([x_upp,x_low])
    y_box[k] = np.array([y_upp,y_low])
    k += 1

    # initialize runaround arrays
    n_run = i_first*2+9
    x_run = np.zeros((n_run,))
    y_run = np.zeros((n_run,))

    # add first part
    x_run[0:5] = np.linspace(x_upp,x_u_first,num=5)
    y_run[0:5] = leth["t"](x_run[0:5])
    y_run[0] = y_upp

    # add next part (flipped top)
    x_run[5:5+i_first] = np.flip(leth["tx"][0:i_first])
    y_run[5:5+i_first] = np.flip(leth["ty"][0:i_first])

    # add next part (bottom)
    x_run[5+i_first:4+2*i_first] = leth["bx"][1:i_first]
    y_run[5+i_first:4+2*i_first] = leth["by"][1:i_first]

    # add final part
    x_run[4+2*i_first:] = np.linspace(x_l_first,x_low,num=5)
    y_run[4+2*i_first:] = leth["b"](x_run[4+2*i_first:])
    y_run[-1] = y_low

    # flip direction
    x_run = np.flip(x_run); y_run = np.flip(y_run)

    # initialize indices to pass to box array
    indices = [int(0.0625 * i * n_run) for i in range(16)]
    indices.append(n_run-1)

    # start at lower aft corner, circle clockwise
    for i in range(len(indices)-1):
        x_box[k] = x_run[indices[i]:indices[i+1]+1]
        y_box[k] = y_run[indices[i]:indices[i+1]+1]
        k += 1

    # create inner hollow I beam
    if info["mhl"]:

        # draw out Interps pieces
        huoi = Interps["huoi"]; huii = Interps["huii"]
        hloi = Interps["hloi"]; hlii = Interps["hlii"]

        # define various x values of import
        # x values at "end" of I beam (where semicircles start/end)
        x_left = info["f"] - info["ifw"] - 0.5*info["ith"]
        x_right = info["f"] + info["ifw"] - 0.5*info["ith"]
        # x values at mid of flange of I beam
        x_mid_left = info["f"] - 0.5*info["ifw"] - 0.5*info["ith"]
        x_mid_right = info["f"] + 0.5*info["ifw"]
        # x values at start left and right of mid of flange of I beam
        x_ml_left  = x_mid_left  - 0.5*info["hol"]
        x_ml_right = x_mid_left  + 0.5*info["hol"]
        x_mr_left  = x_mid_right - 0.5*info["hol"]
        x_mr_right = x_mid_right + 0.5*info["hol"]
        # x values near center of I beam
        x_cn_left  = info["f"] - 0.5*info["ith"] - 0.5*info["hol"]
        x_cn_right = info["f"] - 0.5*info["ith"] + 0.5*info["hol"]
        # x values at to the left and right of the center for the joint
        jnt_space = 0.5*info["hol"]
        x_jnt_left  = x_cn_left  - jnt_space
        x_jnt_right = x_cn_right + jnt_space

        # define various y values of import for i beam segments
        y_cn_left = huii["t"](x_cn_left); y_cn_lbottom = hlii["b"](x_cn_left)
        y_cn_right = huii["t"](x_cn_right); y_cn_rbottom =hlii["b"](x_cn_right)
        height_i_left  = y_cn_left - y_cn_lbottom
        height_i_right = y_cn_right - y_cn_rbottom
        y_up_left  = y_cn_left  - height_i_left  / 3. + 0.5*info["hol"]
        y_up_right = y_cn_right - height_i_right / 3. + 0.5*info["hol"]
        y_up_up = min(y_up_left,y_up_right)
        y_up_lo = y_up_up - info["hol"]
        y_lo_left  = y_cn_left  - 2.*height_i_left  / 3. + 0.5*info["hol"]
        y_lo_right = y_cn_right - 2.*height_i_right / 3. + 0.5*info["hol"]
        y_lo_up = max(y_lo_left,y_lo_right)
        y_lo_lo = y_lo_up - info["hol"]
        # y values at the top and bottom of the i beam for the joint
        y_up = min(y_cn_left - jnt_space,y_cn_right - jnt_space)
        y_lo = max(y_cn_lbottom + jnt_space,y_cn_rbottom + jnt_space)

        # run through points
        x_box[k] = np.linspace(x_left,x_ml_left,num=7)
        y_box[k] = huoi["t"](x_box[k]); y_left = y_box[k][0]
        k += 1

        # create wall
        x_box[k] = np.array([x_ml_left,x_ml_left])
        y_box[k] = np.array([y_box[k-1][-1],huii["t"](x_ml_left)])
        k += 1

        # run through points
        x_box[k] = np.linspace(x_ml_left,x_left,num=7)
        y_box[k] = huii["t"](x_box[k]); y_box[k][0] = y_box[k-1][-1]
        k += 1

        # define end of circle runaround
        x_end_circle = x_left
        y_mid_circle = (y_left + y_box[k-1][-1]) / 2.

        # define circle parameters
        r = 0.5 * info["hol"]
        t = np.linspace(3.*np.pi/2.,np.pi/2.,num=15)

        # create semi-circle
        x_box[k] = r * np.cos(t) + x_end_circle
        y_box[k] = r * np.sin(t) + y_mid_circle; y_box[k][0] = y_box[k-1][-1]
        y_box[k][-1] = y_left
        k += 1

        # create left of T

        # run through points
        x_box[k] = np.linspace(x_ml_right,x_jnt_left,num=7)
        y_box[k] = huoi["t"](x_box[k]); y_ml_right = y_box[k][0]
        k += 1

        # create wall
        x_box[k] = np.array([x_jnt_left,x_jnt_left])
        y_box[k] = np.array([y_box[k-1][-1],huii["t"](x_jnt_left)])
        k += 1

        # run through points
        x_box[k] = np.linspace(x_jnt_left,x_ml_right,num=7)
        y_box[k] = huii["t"](x_box[k]); y_box[k][0] = y_box[k-1][-1]
        k += 1

        # create wall
        x_box[k] = np.array([x_ml_right,x_ml_right])
        y_box[k] = np.array([y_box[k-1][-1],y_ml_right])
        k += 1

        # create right of T

        # run through points
        x_box[k] = np.linspace(x_jnt_right,x_mr_left,num=7)
        y_box[k] = huoi["t"](x_box[k]); y_jnt_right = y_box[k][0]
        k += 1

        # create wall
        x_box[k] = np.array([x_mr_left,x_mr_left])
        y_box[k]= np.array([y_box[k-1][-1],huii["t"](x_mr_left)])
        k += 1

        # run through points
        x_box[k] = np.linspace(x_mr_left,x_jnt_right,num=7)
        y_box[k] = huii["t"](x_box[k]); y_box[k][0] = y_box[k-1][-1]
        k += 1

        # create wall
        x_box[k] = np.array([x_jnt_right,x_jnt_right])
        y_box[k]= np.array([y_box[k-1][-1],y_jnt_right])
        k += 1

        # create bottom of T

        # add top
        x_box[k] = np.array([x_cn_left,x_cn_right])
        y_box[k] = np.array([y_up,y_up])
        k += 1

        # determine run down mid i beam 
        x_box[k] = np.array([x_cn_right,x_cn_right])
        y_box[k] = np.array([y_box[k-1][-1],y_up_up])
        k += 1

        # add bottom of T
        x_box[k] = np.array([x_cn_right,x_cn_left])
        y_box[k] = np.array([y_up_up,y_up_up])
        k += 1

        # determine run up mid i beam
        x_box[k] = np.array([x_cn_left,x_cn_left])
        y_box[k] = np.array([y_up_up,y_up])
        k += 1

        # create run over top of top right

        # run through points
        x_box[k] = np.linspace(x_mr_right,x_right,num=7)
        y_box[k] = huoi["t"](x_box[k]); y_mr_right = y_box[k][0]
        k += 1

        # define end of circle runaround
        x_end_circle = x_right
        y_mid_circle = (huii["t"](x_end_circle) + y_box[k-1][-1]) / 2.

        # define circle parameters
        t = np.linspace(np.pi/2.,-np.pi/2.,num=15)

        # create semi-circle
        x_box[k] = r * np.cos(t) + x_end_circle
        y_box[k] = r * np.sin(t) + y_mid_circle; y_box[k][0] = y_box[k-1][-1]
        k += 1

        # run through points
        x_box[k] = np.linspace(x_right,x_mr_right,num=7)
        y_box[k] = huii["t"](x_box[k]); y_box[k][0] = y_box[k-1][-1]
        k += 1

        # create wall
        x_box[k] = np.array([x_mr_right,x_mr_right])
        y_box[k] = np.array([y_box[k-1][-1],y_mr_right])
        k += 1

        # create center box

        # create top line
        x_box[k] = np.array([x_cn_left,x_cn_right])
        y_box[k] = np.array([y_up_lo,y_up_lo])
        k += 1

        # create right line
        x_box[k] = np.array([x_cn_right,x_cn_right])
        y_box[k] = np.array([y_up_lo,y_lo_up])
        k += 1

        # create bottom line
        x_box[k] = np.array([x_cn_right,x_cn_left])
        y_box[k] = np.array([y_lo_up,y_lo_up])
        k += 1

        # create left line
        x_box[k] = np.array([x_cn_left,x_cn_left])
        y_box[k] = np.array([y_lo_up,y_up_lo])
        k += 1

        # create bottom left piece

        # run through points
        x_box[k] = np.linspace(x_left,x_ml_left,num=7)
        y_box[k] = hlii["b"](x_box[k]); y_left = y_box[k][0]
        k += 1

        # create wall
        x_box[k] = np.array([x_ml_left,x_ml_left])
        y_box[k] = np.array([y_box[k-1][-1],hloi["b"](x_ml_left)])
        k += 1

        # run through points
        x_box[k] = np.linspace(x_ml_left,x_left,num=7)
        y_box[k] = hloi["b"](x_box[k]); y_box[k][0] = y_box[k-1][-1]
        k += 1

        # define end of circle runaround
        x_end_circle = x_left
        y_mid_circle = (y_left + y_box[k-1][-1]) / 2.

        # define circle parameters
        t = np.linspace(3.*np.pi/2.,np.pi/2.,num=15)

        # create semi-circle
        x_box[k] = r * np.cos(t) + x_end_circle
        y_box[k] = r * np.sin(t) + y_mid_circle; y_box[k][0] = y_box[k-1][-1]
        y_box[k][-1] = y_left
        k += 1

        # create bottom T part

        # create left of T

        # run through points
        x_box[k] = np.linspace(x_ml_right,x_jnt_left,num=7)
        y_box[k] = hlii["b"](x_box[k]); y_ml_right = y_box[k][0]
        k += 1

        # create wall
        x_box[k] = np.array([x_jnt_left,x_jnt_left])
        y_box[k] = np.array([y_box[k-1][-1],hloi["b"](x_jnt_left)])
        k += 1

        # run through points
        x_box[k] = np.linspace(x_jnt_left,x_ml_right,num=7)
        y_box[k] = hloi["b"](x_box[k]); y_box[k][0] = y_box[k-1][-1]
        k += 1

        # create wall
        x_box[k] = np.array([x_ml_right,x_ml_right])
        y_box[k] = np.array([y_box[k-1][-1],y_ml_right])
        k += 1

        # create right of T

        # run through points
        x_box[k] = np.linspace(x_jnt_right,x_mr_left,num=7)
        y_box[k] = hlii["b"](x_box[k]); y_jnt_right = y_box[k][0]
        k += 1

        # create wall
        x_box[k] = np.array([x_mr_left,x_mr_left])
        y_box[k]= np.array([y_box[k-1][-1],hloi["b"](x_mr_left)])
        k += 1

        # run through points
        x_box[k] = np.linspace(x_mr_left,x_jnt_right,num=7)
        y_box[k] = hloi["b"](x_box[k]); y_box[k][0] = y_box[k-1][-1]
        k += 1

        # create wall
        x_box[k] = np.array([x_jnt_right,x_jnt_right])
        y_box[k]= np.array([y_box[k-1][-1],y_jnt_right])
        k += 1

        # create bottom of T

        # add top
        x_box[k] = np.array([x_cn_left,x_cn_right])
        y_box[k] = np.array([y_lo_lo,y_lo_lo])
        k += 1

        # determine run down mid i beam 
        x_box[k] = np.array([x_cn_right,x_cn_right])
        y_box[k] = np.array([y_box[k-1][-1],y_lo])
        k += 1

        # add bottom of T
        x_box[k] = np.array([x_cn_right,x_cn_left])
        y_box[k] = np.array([y_lo,y_lo])
        k += 1

        # determine run up mid i beam
        x_box[k] = np.array([x_cn_left,x_cn_left])
        y_box[k] = np.array([y_lo,y_lo_lo])
        k += 1

        # create bottom right

        # run through points
        x_box[k] = np.linspace(x_mr_right,x_right,num=7)
        y_box[k] = hlii["b"](x_box[k]); y_mr_right = y_box[k][0]
        k += 1

        # define end of circle runaround
        x_end_circle = x_right
        y_mid_circle = (hloi["b"](x_end_circle) + y_box[k-1][-1]) / 2.

        # define circle parameters
        t = np.linspace(np.pi/2.,-np.pi/2.,num=15)

        # create semi-circle
        x_box[k] = r * np.cos(t) + x_end_circle
        y_box[k] = r * np.sin(t) + y_mid_circle; y_box[k][0] = y_box[k-1][-1]
        k += 1

        # run through points
        x_box[k] = np.linspace(x_right,x_mr_right,num=7)
        y_box[k] = hloi["b"](x_box[k]); y_box[k][0] = y_box[k-1][-1]
        k += 1

        # create wall
        x_box[k] = np.array([x_mr_right,x_mr_right])
        y_box[k] = np.array([y_box[k-1][-1],y_mr_right])
        k += 1

    # y = np.array([-0.5,1.0])
    # x = info["f"] - info["ifw"] - 0.5*info["ith"]
    # plt.plot(np.array([x,x])*info["c"],y,"b")
    # x = info["f"] + info["ifw"] - 0.5*info["ith"]
    # plt.plot(np.array([x,x])*info["c"],y,"b")

    # x = info["f"] - 0.5*info["ifw"] - 0.5*info["ith"]
    # plt.plot(np.array([x,x])*info["c"],y,"r")
    # x = info["f"] + 0.5*info["ifw"]
    # plt.plot(np.array([x,x])*info["c"],y,"r")

    # x = info["f"] - 0.5*info["ifw"] - 0.5*info["ith"] - 0.5*info["hol"]
    # plt.plot(np.array([x,x])*info["c"],y,"orange")
    # x = info["f"] - 0.5*info["ifw"] - 0.5*info["ith"] + 0.5*info["hol"]
    # plt.plot(np.array([x,x])*info["c"],y,"orange")
    # x = info["f"] + 0.5*info["ifw"] - 0.5*info["hol"]
    # plt.plot(np.array([x,x])*info["c"],y,"orange")
    # x = info["f"] + 0.5*info["ifw"] + 0.5*info["hol"]
    # plt.plot(np.array([x,x])*info["c"],y,"orange")

    # x = info["f"] - 0.5*info["ith"]
    # plt.plot(np.array([x,x])*info["c"],y,"g")

    # x = info["f"] - 0.5*info["ith"] - 0.5*info["hol"]
    # plt.plot(np.array([x,x])*info["c"],y,"y")
    # x = info["f"] - 0.5*info["ith"] + 0.5*info["hol"]
    # plt.plot(np.array([x,x])*info["c"],y,"y")

    # print("k:",k,", n:",n_box)

    return x_box,y_box

# create outer airfoil shape
def outer(O,info):
    # intialize how many sections to do
    number_extra = 11

    # initialize outer arrays
    x_outer = np.zeros((number_extra + 3,),dtype=np.ndarray)
    y_outer = np.zeros((number_extra + 3,),dtype=np.ndarray)

    # determine index to end at
    end_index = int(0.625 * info["Nres"])

    # create an indices array
    second_index = int(0.27 * info["Nres"])
    second_to_last_index = int(0.52 * info["Nres"])
    index_step = int((second_to_last_index - second_index) / (number_extra-3))
    
    indices = [0]
    for i in range(number_extra-1):
        indices.append(second_index + i * index_step)
    indices.append(end_index-1)

    # run through each index, create outer airfoil top runaround
    for i in range(len(indices)-1):
        x_outer[i] = O["x af"][indices[i]:indices[i+1]+1]
        y_outer[i] = O["y af"][indices[i]:indices[i+1]+1]

    # set min points to prevent running outside interpolation ranges
    b_min = O["x af"][end_index-1]; b_max = info["outer"]["start"][0]

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
        if group == "arc" or group == "rfl" or group == "kin":
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
        if group == "arc" or group == "rfl" or group == "kin":
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
        if group == "arc" or group == "rfl" or group == "kin":
            for i in range(x[group].shape[0]):
                x[group][i],y[group][i] = deflect.main(O,x[group][i],\
                    y[group][i],info)
        else:
            x[group],y[group] = deflect.main(O,x[group],y[group],info)
    
    return

# shift to quarter chord
def shift(x,y):
    # determine shift values
    x_shift = 0.25
    y_shift = 0.0#O["c"](x_shift)

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
        return x,y
    
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
        return x,y
    
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
    xsplit[j+1] = x["kin"][0][0]
    ysplit[j+1] = y["kin"][0][0]

    # set z array as zero
    zsplit = ysplit * 0.

    return xsplit,ysplit,zsplit

# determine the max and min of each "hole"
def guide_curves(x,y,I,vals,info,run):
    # initialize points matrix.
    # i 0 --> holes, with 0 = outer skin
    # j guide curve number
    # k 0 - x, 1 - y, 2 - z
    number_to_make = 6 + info["nk"] - (info["nk"] + 4) * info["sk"]
    guides = np.zeros((number_to_make,),dtype=np.ndarray)
    xpts = np.zeros((number_to_make,),dtype=np.ndarray)
    ypts = np.zeros((number_to_make,),dtype=np.ndarray)

    # initialize group and point tuples lists
    grp = [["mou","mou","mou","mou","mou","mou","mou","mou","mou","mou","mou",\
        "kin","kin","ton","ton","out","out"]]
    point = [
        [(0,-1),(1,0),(2,0),(3,-1),(5,0),(5,-1),(6,-1),(7,0),(8,-1),(9,-1),\
            (10,-1),((0,0),-1),((0,0),0),(2,0),(2,-1),(-2,-1),(-1,-1)]]

    # if not info["sk"]:
    #     grp = grp + [["box","box"],
    #         ["box","box","box"],
    #         ["box","box","box"],
    #         ["box","box","box"]
    #     ]
    #     point = point + [[(9,-1),(9,0)],
    #         [(7,0),(7,-1),(8,0)],
    #         [(5,-1),(4,-1),(4,0)],
    #         [(1,0),(1,-1),(2,0)]
    #     ]

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
        num_box = x["box"].shape[0]

        # run through each number and add that many to the to be added lists
        grp1 = []; point1 = []
        for i in range(num_box):
            grp1.append("box")
            point1.append((i,0))
        
        # add trailing edge triangle
        grp1.append("tri"); grp1.append("tri"); grp1.append("tri")
        point1.append((1,-1)); point1.append((0,0)); point1.append((0,-1))
        
        # add to 2d lists
        grp.append(grp1); point.append(point1)

        # # add a new index to the grp and point lists for each arc (minus 1)
        # j = 0
        # for i in range(info["nk"]-1):
        #     grp.append(["kin","kin","kin","kin"])
        #     point.append([((1,i),0),((0,i+1),0),((0,i+1),-1),((1,i),-1)])
        #     j += 1
        
        # # add for final arc and te triangle
        # grp.append(["kin","tri","tri","kin"])
        # point.append([((1,j),0),(3,-1),(3,0),((1,j),-1)])

        # # add trailing edge triangle
        # grp.append(["tri","tri","tri"])
        # point.append([(1,-1),(0,0),(0,-1)])

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
    
        # plt.plot(x_hole*info["c"],y_hole*info["c"],"bo")

    # if shift origin to quarter chord
    if "shift to c/4" in vals[run] and vals[run]["shift to c/4"]:
        # shift
        shift(xpts,ypts)

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
def split_part(x,y,I,vals,info,run):
    # create result dictionary
    results = {}
    
    # create dxf file of inner holes for loft cut
    # initialize inner dxf file
    num_inner = 4*info["nk"]
    x_inner = np.zeros((num_inner,),dtype=np.ndarray)
    y_inner = np.zeros((num_inner,),dtype=np.ndarray)

    # run through roofs and floors and add the required lines
    k = 0 # counter for inner arrays
    for i in range(info["nk"]):
        # add first wall
        x_inner[k] = x["kin"][1][i]
        y_inner[k] = y["kin"][1][i]
        k += 1

        # add roof
        x_inner[k] = x["rfl"][0][i+1]
        y_inner[k] = y["rfl"][0][i+1]
        k += 1

        # add second wall
        if i == info["nk"]-1:
            x_inner[k] = x["tri"][3]
            y_inner[k] = y["tri"][3]
        else:
            x_inner[k] = np.flip(x["kin"][0][i+1])
            y_inner[k] = np.flip(y["kin"][0][i+1])
        k += 1

        # add floor
        x_inner[k] = np.flip(x["rfl"][1][i+1])
        y_inner[k] = np.flip(y["rfl"][1][i+1])
        k += 1

    # save to results dictionary
    results["x inner"] = x_inner
    results["y inner"] = y_inner

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

    # save tongue tip midpoint
    xton0 = (x["ton"][2][0]+x["ton"][2][1])/2.
    yton0 = (y["ton"][2][0]+y["ton"][2][1])/2.
    xton1 = (info["f"] + info["to"])/2.
    yton1 = (I["it"]["b"](xton1) + I["ot"]["b"](xton1))/2.
    xton = np.array([xton0,xton1]); yton = np.array([yton0,yton1])

    # if shift origin to quarter chord
    if "shift to c/4" in vals[run] and vals[run]["shift to c/4"]:
        # shift
        shift(results["x inner"],results["y inner"])
        shift(results["x inn gc"],results["y inn gc"])
        shift(xton,yton)

    # resize chord of each segment
    resize(results["x inner"],results["y inner"],info["c"])
    resize(results["x inn gc"],results["y inn gc"],info["c"])
    resize(xton,yton,info["c"])

    # create inn guide curve dictionary element
    results["inn gc"] = np.zeros((num_inner,3))
    results["inn gc"][:,0] = results["x inn gc"]
    results["inn gc"][:,1] = results["y inn gc"]

    results["ton tip"] = np.zeros((2,3))
    results["ton tip"][:,0] = xton; results["ton tip"][:,1] = yton
    # plt.plot(xton,yton,"k")

    return results

# create actuation hole
def actuate_hole(x,y,I,vals,info,results,run):
    # create actuation hole # initialize array
    num = 4
    x_hole = np.zeros((num,),dtype=np.ndarray)
    y_hole = np.zeros((num,),dtype=np.ndarray)
    j = 0

    # create run along top of tongue
    x_hole[j] = np.linspace(info["mo"]+0.005,info["f"]+info["ifw"]+0.01,25)
    y_hole[j] = (I["it"]["b"](x_hole[j]) + I["om"]["b"](x_hole[j])) / 2.
    x_start = x_hole[j][0]; y_start = y_hole[j][0]
    j += 1

    # create aft wall
    x_hole[j] = np.array([x_hole[j-1][-1],x_hole[j-1][-1]])
    y_hole[j] = np.array([y_hole[j-1][-1],y_hole[j-1][-1] + 0.03])
    j += 1

    # create top
    x_hole[j] = np.array([x_hole[j-1][-1],x_start])
    y_hole[j] = np.array([y_hole[j-1][-1],y_hole[j-1][-1]])
    j += 1

    # create forward wall
    x_hole[j] = np.array([x_start,x_start])
    y_hole[j] = np.array([y_hole[j-1][-1],y_start])
    j += 1

    # for i in range(x_hole.shape[0]):
    #     plt.plot(x_hole[i]*info["c"],y_hole[i]*info["c"],"r")

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
    if "shift to c/4" in vals[run] and vals[run]["shift to c/4"]:
        # shift
        shift(results["x hole"],results["y hole"])
        shift(results["x hole gc"],results["y hole gc"])

    # resize chord of each segment
    resize(results["x hole"],results["y hole"],info["c"])
    resize(results["x hole gc"],results["y hole gc"],info["c"])

    # create inn guide curve dictionary element
    results["hole gc"] = np.zeros((num,3))
    results["hole gc"][:,0] = results["x hole gc"]
    results["hole gc"][:,1] = results["y hole gc"]

    return

# create servo box for visual check
def make_servo(O,it,info,vals,run):

    # initialize servo arrays
    x_servo = np.zeros((4,),dtype=np.ndarray)
    y_servo = np.zeros((4,),dtype=np.ndarray)

    # read in values
    servo = vals[run]["servo"]
    tilt_servo = vals[run][servo]["tilt"]
    width = vals[run][servo]["width [in]"] / info["c"]
    thickness = vals[run][servo]["thickness [in]"] / info["c"]
    arm_length = vals[run][servo]["arm length [in]"] / info["c"]
    x_cg = vals[run][servo]["x-cg [x/c]"]
    connector_height = vals[run][servo]["connector height [mm]"] / 25.4 \
        / info["c"]

    # start at forward upper corner of the servo, move clockwise
    # initialize index counter
    k = 0

    # find this point
    x_fp = x_cg - width / 2.
    y_fp = O["t"](x_fp) - info["t"]

    # find aft top point
    if tilt_servo:
        # create function to determine where the the servo hits upper surface
        def func(x):
            return O["t"](x) - info["t"] -((width**2. -(x-x_fp)**2.)**0.5+y_fp)
        x_ap = Bisect(func,x_fp+width*0.5,x_fp + width*1.5)
        y_ap = O["t"](x_ap) - info["t"]
        tilt = np.arctan2(y_ap-y_fp,x_ap-x_fp)
    else:
        tilt = 0.0
        x_ap = x_fp + width
        y_ap = y_fp + 0.0
    # print("tilt is ",np.rad2deg(tilt)," degrees")

    # add this line
    x_servo[k] = np.array([x_fp,x_ap])
    y_servo[k] = np.array([y_fp,y_ap])
    k += 1

    # add next line (aft vertical)
    x_servo[k] = np.array([x_ap,x_ap + np.sin(tilt) * thickness])
    y_servo[k] = np.array([y_ap,y_ap - np.cos(tilt) * thickness])
    k += 1

    # add next line (lower horizontal)
    x_servo[k] = np.array([x_ap + np.sin(tilt) * thickness,x_fp + \
        np.sin(tilt) * thickness])
    y_servo[k] = np.array([y_ap - np.cos(tilt) * thickness,y_fp - \
        np.cos(tilt) * thickness])
    k += 1

    # add final line (forward vertical)
    x_servo[k] = np.array([x_fp + np.sin(tilt) * thickness,x_fp])
    y_servo[k] = np.array([y_fp - np.cos(tilt) * thickness,y_fp])
    k += 1

    # report servo cg
    x_servo_pts = np.array([x[0] for x in x_servo])
    x_cg_avg = np.average(np.array(x_servo_pts))
    print("Servo CG at ",x_cg_avg," x/c")

    # initialize arm arrays
    x_arm = np.zeros((1,),dtype=np.ndarray)
    y_arm = np.zeros((1,),dtype=np.ndarray)

    # initialize servo arm x length along width
    x_servo_arm = vals[run][servo]["x servo arm [x/width]"] * width

    # determine x value and y value of the center of rotation
    x_rp = x_fp + np.cos(tilt) * x_servo_arm + np.sin(tilt) * thickness / 2.
    y_rp = y_fp + np.sin(tilt) * x_servo_arm - np.cos(tilt) * thickness / 2.

    # reinitialize k counter
    k = 0

    # add servo arm
    x_arm[k] = np.array([x_rp,x_rp])
    y_arm[k] = np.array([y_rp,y_rp - arm_length])
    k += 1

    # reinitialize counter
    k = 0

    # initialize arm arrays
    x_ahl = np.zeros((2,),dtype=np.ndarray)
    y_ahl = np.zeros((2,),dtype=np.ndarray)

    # find the point on the tongue where the actuation will occur
    x_hp = (info["f"] + info["to"]) / 2.
    y_hp = it["b"](x_hp)

    # create the first line
    x_ahl[k] = np.array([x_hp,x_hp])
    y_ahl[k] = np.array([y_hp,y_hp + connector_height])
    k += 1

    # create the second line
    x_ahl[k] = np.array([x_hp,0.1])
    y_ahl[k] = np.array([y_hp + connector_height,y_hp + connector_height])
    k += 1

    # report the overpass length
    print("overpass ",\
        ((y_hp + connector_height) - (y_rp - arm_length))*info["c"]," in")

    return x_servo,y_servo,x_arm,y_arm,x_ahl,y_ahl

# main function
def main(jsonfile,run="C18"):
    """
    This is an explanation of what this file does
    """
    # import json file if a dictionary or json file, OR raise error
    if type(jsonfile) == dict:
        vals = jsonfile
    elif os.path.isfile(jsonfile):
        vals = read_json(jsonfile)
    else:
        raise ValueError("C14 input must be dictionary or file path to a json file")

    # simplify dictionary
    info = simplify(vals,run)

    # create airfoil shape
    add = {}
    add["geometry"] = {}
    add["geometry"]["type"] = "outline_points"
    afname = vals[run]["airfoil dat file location"]
    if afname[:4] == "NACA":
        add["geometry"]["NACA"] = afname.split(" ")[-1]
        isNACA = True
    else:
        add["geometry"]["outline_points"] = afname
        isNACA = False

    # initialize airfoil database
    foil = adb.Airfoil("foil",add)
    coords = foil.get_outline_points(close_te=not(isNACA),top_first=True,\
        N=info["Nres"])
    x_af = coords[:,0]; y_af = coords[:,1]

    # enthicken the airfoil
    y_af *= vals[run]["enthicken"]

    # create inner airfoil
    polyx,polyy = inner(x_af,y_af,info["t"])

    # initialize interpolation dictionary
    I = {}

    # create interpolation of inner surface
    I["I"] = interpolate(polyx,polyy,make_camberline=True)

    # plt.axis('equal')
    # plt.plot(x_af*info["c"],y_af*info["c"])
    # plt.plot(polyx*info["c"],polyy*info["c"])
    # plt.show()

    # initialize x and y dictionaries
    x = {}; y = {}

    # determine where the arcs should start-ish
    if "kinks start [x/c]" in vals[run] and "kinks end [x/c]" in vals[run]:
        info["k0"] = vals[run]["kinks start [x/c]"]
        info["k1"] = vals[run]["kinks end [x/c]"]
        I["I"]["arcdx"] = (info["k1"] - info["t"] - info["k0"]) / (info["nk"]+1)
        if "straighten TE" in vals[run] and vals[run]["straighten TE"]:
            info["tw"] = vals[run]["kinks end [x/c]"]
    else:
        I["I"]["arcdx"] = (info["tw"] - info["t"] - info["f"]) / (info["nk"]+1)

    # find te triangle
    x["tri"],y["tri"] = te_triangle(I["I"],info)

    # determine kink shift value
    if "kink shift [x/c]" in vals[run]:
        kink_shift = vals[run]["kink shift [x/c]"]
    else:
        kink_shift = 0.0

    # find where each arc is
    x["kin"],y["kin"] = kincs(I["I"],info,kink_shift)

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

    # plt.axis('equal')
    # plt.plot(x_af*info["c"],y_af*info["c"])
    # plt.plot(polyx*info["c"],polyy*info["c"])
    # plt.plot(out_tx*info["c"],out_ty*info["c"])
    # plt.plot(inn_tx*info["c"],inn_ty*info["c"])
    # # out_mx,out_my = offset(x_af,y_af,2.*info["t"]+2.*info["mt"])
    # # plt.plot(out_mx*info["c"],out_my*info["c"])
    # # plt.plot(inn_mx*info["c"],inn_my*info["c"])
    # plt.show()

    # outer mouth
    out_mx,out_my = inner(x_af,y_af,2.*info["t"]+2.*info["mt"])
    I["om"] = interpolate(out_mx,out_my)
    # inner mouth
    inn_mx,inn_my = inner(x_af,y_af,3.*info["t"]+2.*info["mt"])
    I["im"] = interpolate(inn_mx,inn_my)
    # print(); print("here"); print()

    # plt.axis('equal')
    # plt.plot(x_af*info["c"],y_af*info["c"])
    # plt.plot(polyx*info["c"],polyy*info["c"])
    # plt.plot(out_tx*info["c"],out_ty*info["c"])
    # plt.plot(inn_tx*info["c"],inn_ty*info["c"])
    # plt.plot(out_mx*info["c"],out_my*info["c"])
    # plt.plot(inn_mx*info["c"],inn_my*info["c"])
    # plt.show()

    # inner upper i beam
    iui_mx,iui_my = inner(x_af,y_af,info["ith"])
    I["iui"] = interpolate(iui_mx,iui_my)
    # inner lower i beam
    ili_mx,ili_my = inner(x_af,y_af,2.*info["t"]+2.*info["mt"]+info["ith"])
    I["ili"] = interpolate(ili_mx,ili_my)
    # inner lower i beam
    leth_mx,leth_my = inner(x_af,y_af,info["leth"])
    I["leth"] = interpolate(leth_mx,leth_my)
    
    if info["mhl"]:
        # hollow upper outer i beam
        huoi_x,huoi_y = inner(x_af,y_af,0.5*info["ith"] - 0.5*info["hol"])
        I["huoi"] = interpolate(huoi_x,huoi_y)
        # hollow upper inner i beam
        huii_x,huii_y = inner(x_af,y_af,0.5*info["ith"] + 0.5*info["hol"])
        I["huii"] = interpolate(huii_x,huii_y)
        # hollow lower outer i beam
        hloi_x,hloi_y = inner(x_af,y_af,2.*info["t"]+2.*info["mt"] + \
            0.5*info["ith"] - 0.5*info["hol"])
        I["hloi"] = interpolate(hloi_x,hloi_y)
        # hollow lower inner i beam
        hlii_x,hlii_y = inner(x_af,y_af,2.*info["t"]+2.*info["mt"] + \
            0.5*info["ith"] + 0.5*info["hol"])
        I["hlii"] = interpolate(hlii_x,hlii_y)

        # plt.plot(huoi_x*info["c"],huoi_y*info["c"])
        # plt.plot(huii_x*info["c"],huii_y*info["c"])
        # plt.plot(hloi_x*info["c"],hloi_y*info["c"])
        # plt.plot(hlii_x*info["c"],hlii_y*info["c"])

    # create outer mouth structure
    x["mou"],y["mou"] =make_mouth(I["O"],I["I"],I["om"],I["iui"],I["ili"],info)

    # create tongue tip
    x["ton"],y["ton"] = tongue_tip(I["O"],I["I"],I["ot"],I["it"],info)

    # create roofs and floors
    x["rfl"],y["rfl"] = roofs_n_floors(I["I"],info,x["mou"][-1][-1],\
        y["tri"][3][0],y["tri"][3][1])

    # create wing box
    if not info["sk"]:
        x["box"],y["box"] = wing_box(I,info)

    # create outer run. start and end values are where break starts and ends
    x["out"],y["out"] = outer(I["O"],info)

    # remove kincs desired
    remove_kincs(x,y,info)

    # add servo box if desired
    if "show servo" in vals[run] and vals[run]["show servo"]:
        x["ser"],y["ser"],x["arm"],y["arm"],x["ahl"],y["ahl"] = \
            make_servo(I["O"],I["it"],info,vals,run)

    # report hinge point location
    if "report hinge point" in vals[run] and vals[run]["report hinge point"]:
        y_hp = (I["I"]["t"](info["f"]) + I["O"]["t"](info["f"])) / 2.
        print("Hinge Point at (x/c,y/c) : (",info["f"],",",y_hp,")")

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
            shift(u,v)

        # resize chord of each segment
        resize(u,v,info["c"])

        # create dxf file
        filename = vals[run]["dxf file path"]+vals[run]["dxf file name"]
        xvals,yvals,zvals = make_dxf(u,v,filename,vals[run]["write dxf"])

        # if split return desired
        if vals[run]["split return"]:
            xsplit,ysplit,zsplit = split(u,v)

        if vals[run]["guide curve return"]:
            guides = guide_curves( x,y,I,vals,info,run)
            splitdict = split_part(x,y,I,vals,info,run)
            
            # # test to ensure it is 
            # plt.axis("equal")
            # for i in range(splitdict["x inner"].shape[0]):
            #     plt.plot(splitdict["x inner"][i],splitdict["y inner"][i])
            #     plt.plot(splitdict["x inner"][i][1],splitdict["y inner"][i][1],"r*")
            #     plt.plot(splitdict["x inn gc"][i],splitdict["y inn gc"][i],"ko")
            # plt.show()
            # if actuation hole return
            if vals[run]["actuation hole return"]:
                actuate_hole(x,y,I,vals,info,splitdict,run)
            

        if vals[run]["show plot"]:
            # plt.plot(x_af*info["c"],y_af*info["c"])
            # plt.plot(polyx*info["c"],polyy*info["c"])
            # plt.plot(out_tx*info["c"],out_ty*info["c"])
            # plt.plot(inn_tx*info["c"],inn_ty*info["c"])
            # plt.plot(out_mx*info["c"],out_my*info["c"])
            # plt.plot(inn_mx*info["c"],inn_my*info["c"])
            # plt.plot(polyx*info["c"],polyy*info["c"])

            # initialize colors and labels
            g = {}
            g["tri c"] = "#6a0dad"; g["tri l"] = "Triangle"
            g["arc c"] = "#ff8c00"; g["arc l"] = "Arc"
            g["kin c"] = "#ff8c00"; g["kin l"] = "Kinc"
            g["mou c"] = "#0073ff"; g["mou l"] = "Mouth"
            g["ton c"] = "#ff0000"; g["ton l"] = "Tongue"
            g["rfl c"] = "#3cb371"; g["rfl l"] = "roofsNfloors"
            g["box c"] = "#aa6c39"; g["box l"] = "Wing Box"
            g["out c"] =       "k"; g["out l"] = "Outer"
            g["ser c"] =       "r"; g["ser l"] = "Servo"
            g["arm c"] =       "b"; g["arm l"] = "ServoArm"
            g["ahl c"] = "#6a0dad"; g["ahl l"] = "ActuationHole"

            # plot to be saved items
            for group in u:
                # color
                if vals[run]["all black"]:
                    c = "k"
                else:
                    c = g[group+" c"]
                
                # if roofs and floors or arcs, do i and j
                if group == "kin": 
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
                            elif not(i==1 and j==0):
                                plt.plot(u[group][i][j],v[group][i][j],c)
                elif group not in ["arm","ahl"]:
                    for i in range(u[group].shape[0]):
                        if i == 0:
                            l = g[group+" l"]
                        else:
                            l = ""
                        if group == "ser":
                            c = "b"
                        
                        if group not in ["box"] or (i < 33 or \
                            (i > 44 and i < 49) or (i > 52 and i < 57) or \
                                (i > 68 and i < 73)):
                            plt.plot(u[group][i],v[group][i],c,label=l)
    
    # # # plot the guides values
    # # for i in range(guides.shape[0]):
    # #     for j in range(guides[i].shape[0]):
    # #         plt.plot(guides[i][j,0],guides[i][j,1],"ok")

    # # plot holes
    # # for i in range(splitdict["x te"].shape[0]):
    # #     plt.plot(splitdict["x te"][i],splitdict["y te"][i],"r")
    # for i in range(splitdict["x inner"].shape[0]):
    #     plt.plot(splitdict["x inner"][i],splitdict["y inner"][i],"b")
    #     plt.plot(splitdict["x inn gc"][i],splitdict["y inn gc"][i],"ko")
    # for i in range(splitdict["x arcs"].shape[0]):
    #     plt.plot(splitdict["x arcs"][i],splitdict["y arcs"][i],"m")
    # # for i in range(info["nk"]):
    # #     for j in range(4):
    # #         plt.plot(splitdict["x arcs gc"][i,j],splitdict["y arcs gc"][i,j]\
    # #             ,"go")


    if vals[run]["show plot"]:
        if vals[run]["show legend"]:
            plt.legend()
        plt.xlabel("x/c")
        plt.ylabel("y/c")
        # plt.xlim((-0.5999992336796491,12.599999963508555))
        # plt.ylim((-3.389520272304698,3.831447030458203))
        plt.axis("equal")

        # import tikzplotlib
        # tikzplotlib.save("/home/ben/Desktop/kincs.tex")
        plt.show()
    
    # ######################################################
    # for i in range(splitdict["x inner"].shape[0]):
    #     plt.axis("equal")
    #     for j in range(splitdict["x inner"].shape[0]):
    #         plt.plot(splitdict["x inner"][j],splitdict["y inner"][j],"b")
    #     plt.plot(splitdict["x inn gc"][i],splitdict["y inn gc"][i],"ko")
    #     plt.show(block=False)
    #     plt.pause(0.6)
    #     plt.close()

    # ######################################################

    if vals[run]["split return"] and vals[run]["guide curve return"]:
        return xvals, yvals, zvals, xsplit, ysplit, zsplit, guides, splitdict
    elif vals[run]["split return"]:
        return xvals, yvals, zvals, xsplit, ysplit, zsplit
    elif vals[run]["guide curve return"]:
        return xvals, yvals, zvals, guides, splitdict
    else:
        return xvals, yvals, zvals

# # run file
# jsonfile = 'input_c18.json'
# main(jsonfile)