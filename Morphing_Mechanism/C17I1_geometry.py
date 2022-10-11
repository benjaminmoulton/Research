import airfoil_db as adb
import numpy as np
import json
import dxf
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
    dicti["t"]  = vals[run]["shell thickness [mm]"] / dicti["c"] / in2mm
    dicti["f"]  = vals[run]["flex start [x/c]"]
    dicti["to"] = dicti["f"] - np.abs(vals[run]["tongue start [in]"]) \
        / dicti["c"]
    dicti["te"] = dicti["f"] + np.abs(vals[run]["tongue end [in]"]) \
        / dicti["c"]
    dicti["mo"] = dicti["f"] - np.abs(vals[run]["mouth start [in]"]) \
        / dicti["c"]
    if dicti["to"] < dicti["mo"]:
        raise ValueError("Tongue cannot extend beyond mouth start")
    dicti["tw"] = vals[run]["TE wall start [x/c]"]
    dicti["ts"] = vals[run]["Shift TE wall for hole"]
    dicti["mt"] = vals[run]["mouth clearance [mm]"] / dicti["c"] / in2mm
    dicti["fi"] = vals[run]["fillet side length [mm]"] / dicti["c"] / in2mm
    dicti["sp"] = vals[run]["spar width [in]"] / dicti["c"]
    dicti["td"] = vals[run]["TE hole diameter [mm]"] / dicti["c"] / in2mm
    dicti["np"] = vals[run]["number of hole points"]

    dicti["da"] = vals[run]["deflection angle [deg]"]
    dicti["dty"] = vals[run]["deflection type"]
    dicti["th"] = vals[run]["hinge at top"]

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
    x_mouth[5],y_mouth[5] = Bezier_Cubic(x_points,y_points)

    return x_mouth,y_mouth

# create tongue tip geometry
def tongue_tip(O,I,ot,it,info):
    # initialize mouth array
    x_tongue = np.zeros((7,),dtype=np.ndarray)
    y_tongue = np.zeros((7,),dtype=np.ndarray)

    # create a bezier curve from the bottom of first arc to tongue tip
    # find interpolation point where flex starts on inner tongue
    x_fx_it = info["f"]
    y_fx_it,p_fx_it,i_fx_it = get_y_loc(x_fx_it,it["bx"],it["by"])

    # create points for the bezier curve
    y_start,_,_ = get_y_loc(x_fx_it+.03,it["bx"],it["by"])
    x_end = info["te"]
    y_end,__,__ = get_y_loc(  x_end-.03,I["bx"],I["by"])
    y_ten,p_ten,i_ten = get_y_loc(x_end,I["bx"],I["by"])

    # save x y of tongue end for future use
    info["tnd"] = [x_end,y_ten,p_ten,i_ten]

    # create bezier curve
    xpts = np.array([x_fx_it,x_fx_it+.03,x_end-.03,x_end])
    ypts = np.array([y_fx_it,y_start,y_end,y_ten])
    x_tongue[0],y_tongue[0] = Bezier_Cubic(xpts,ypts,n=60)

    # create top of tongue tip
    # find where the tongue starts inner tongue
    x_ts_it = info["to"]
    y_ts_it,p_ts_it,i_ts_it = get_y_loc(x_ts_it,it["bx"],it["by"])

    # initialize bottom of inner mouth
    x_tongue[1] = np.zeros((i_fx_it-i_ts_it+2,))
    y_tongue[1] = np.zeros((i_fx_it-i_ts_it+2,))

    # set first point
    x_tongue[1][0] = x_ts_it
    y_tongue[1][0] = y_ts_it

    # set in between points
    j = 1
    for i in range(i_ts_it,i_fx_it):
        x_tongue[1][j] = it["bx"][i]
        y_tongue[1][j] = it["by"][i]
        j += 1
    
    # set last point
    x_tongue[1][-1] = x_fx_it
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
    x_j_o = x_end
    y_j_o,p_j_o,i_j_o = get_y_loc(x_j_o,O["bx"],O["by"])

    # save to tongue tip
    x_tongue[6] = np.array([x_tongue[3][-1],x_j_o])
    y_tongue[6] = np.array([y_tongue[3][-1],y_j_o])
    # save value so info to be used in a later function
    info["outer"]["end"] = [x_j_o,y_j_o,i_j_o]

    return x_tongue,y_tongue

# create spar box
def spar(O,info):
    # initialize outer arrays
    x_spar = np.zeros((4,),dtype=np.ndarray)
    y_spar = np.zeros((4,),dtype=np.ndarray)

    # determine box center point and read in width
    x_center = 0.25 # define at quarter chord
    y_center = O["c"](x_center)
    half_width = info["sp"] / 2.

    # input bottom of box
    x_spar[0] = np.array([x_center-half_width,x_center+half_width])
    y_spar[0] = np.array([y_center-half_width,y_center-half_width])

    # input right of box
    x_spar[1] = np.array([x_center+half_width,x_center+half_width])
    y_spar[1] = np.array([y_center-half_width,y_center+half_width])

    # input top of box
    x_spar[2] = np.array([x_center+half_width,x_center-half_width])
    y_spar[2] = np.array([y_center+half_width,y_center+half_width])

    # input left of box
    x_spar[3] = np.array([x_center-half_width,x_center-half_width])
    y_spar[3] = np.array([y_center+half_width,y_center-half_width])


    return x_spar,y_spar

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

# determine where the hole can be placed inside the inner airfoil
def find_wall(O,foil,info):
    # define function to be used
    def fun(x):
        return foil.get_thickness(x) - info["td"]/2. - info["t"]
    
    # find where wall should start
    xl = 0.5; xu = np.max(O["tx"])
    info["tw"] = Bisect(fun,xl,xu)
    print("New TE wall value : {:<7.5} x/c".format(info["tw"]))
    
    return

# create te hole
def te_hole(O,I,info):
    # initialize outer arrays
    x_hole = np.zeros((1,),dtype=np.ndarray)
    y_hole = np.zeros((1,),dtype=np.ndarray)

    # determine hole center
    x_center = info["tw"] + 1.*info["td"]
    y_center = O["c"](x_center)

    # initialize theta, radius values
    theta = np.linspace(0,2.*np.pi,num=info["np"])
    radius = np.full(theta.shape,info["td"]/2.)

    # input array values
    x_hole[0] = radius * np.cos(theta) + x_center
    y_hole[0] = radius * np.sin(theta) + y_center

    return x_hole,y_hole

# create inner run shape
def inner_run(I,info):
    # initialize inner arrays
    x_inner = np.zeros((3,),dtype=np.ndarray)
    y_inner = np.zeros((3,),dtype=np.ndarray)

    # determine where to start the inner run
    # determine wall line location
    x_wall = info["tw"]
    y_wt,p_wt,i_wt = get_y_loc(x_wall,I["tx"],I["ty"])
    y_wb,p_wb,i_wb = get_y_loc(x_wall,I["bx"],I["by"])

    # initialize first arrays
    x_inner[0] = np.zeros((i_wb-info["tnd"][3]+2,))
    y_inner[0] = np.zeros((i_wb-info["tnd"][3]+2,))
    j = 0

    # add first value
    x_inner[0][j] = info["tnd"][0]
    y_inner[0][j] = info["tnd"][1]
    j += 1

    # add through values
    for i in range(info["tnd"][3],i_wb):
        x_inner[0][j] = I["bx"][i]
        y_inner[0][j] = I["by"][i]
        j += 1

    # add last value
    x_inner[0][j] = x_wall
    y_inner[0][j] = y_wb

    # wall line
    # input to arrays
    x_inner[1] = np.array([x_wall,x_wall])
    y_inner[1] = np.array([y_wb,y_wt])

    # initialize first arrays
    x_inner[2] = np.zeros((i_wt-info["fnd"][3]+2,))
    y_inner[2] = np.zeros((i_wt-info["fnd"][3]+2,))
    j = 0

    # add first value
    x_inner[2][j] = x_wall
    y_inner[2][j] = y_wt
    j += 1

    # add through values
    for i in range(i_wt-1,info["fnd"][3]-1,-1):
        x_inner[2][j] = I["tx"][i]
        y_inner[2][j] = I["ty"][i]
        j += 1

    # add last value
    x_inner[2][j] = info["fnd"][0]
    y_inner[2][j] = info["fnd"][1]

    return x_inner,y_inner

# create outer airfoil shape
def outer(O,info):
    # initialize outer arrays
    x_outer = np.zeros((2,),dtype=np.ndarray)
    y_outer = np.zeros((2,),dtype=np.ndarray)

    # initialize first run
    x_outer[0] = np.zeros((O["tx"].shape[0]+info["outer"]["start"][2]))
    y_outer[0] = np.zeros((O["ty"].shape[0]+info["outer"]["start"][2]))

    # run through top and add points
    j = 0
    for i in range(O["tx"].shape[0]-1,-1,-1):
        x_outer[0][j] = O["tx"][i]
        y_outer[0][j] = O["ty"][i]
        j += 1
    
    # add bottom points
    for i in range(1,info["outer"]["start"][2]):
        x_outer[0][j] = O["bx"][i]
        y_outer[0][j] = O["by"][i]
        j += 1
    
    # add final point
    x_outer[0][-1] = info["outer"]["start"][0]
    y_outer[0][-1] = info["outer"]["start"][1]

    # initialize second run
    x_outer[1] = np.zeros((O["bx"].shape[0]-info["outer"]["end"][2]+1,))
    y_outer[1] = np.zeros((O["bx"].shape[0]-info["outer"]["end"][2]+1,))

    # add first point
    x_outer[1][0] = info["outer"]["end"][0]
    y_outer[1][0] = info["outer"]["end"][1]

    # run through top and add points
    j = 1
    for i in range(info["outer"]["end"][2],O["bx"].shape[0]):
        x_outer[1][j] = O["bx"][i]
        y_outer[1][j] = O["by"][i]
        j += 1

    return x_outer,y_outer

# create dxf file
def make_dxf(x,y,filename,write_to_file=True):
    # determine number of arrays spots to make
    size = 0
    for group in x:
        shape = x[group].shape
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
    t["teh"] = [[0,True]]
    t["inn"] = [[2,False],[1,False],[0,False]]
    t["out"] = [[1,True],[0,True]]
    t["order"] = ["mou","inn","ton","out"]
    t["new order"] = ["spr","teh"]

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

def main(jsonfile,run="C17"):
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
    add = {}
    add["geometry"] = {}
    add["geometry"]["type"] = "outline_points"
    add["geometry"]["outline_points"] = vals[run]["airfoil dat file location"]

    # initialize airfoil database
    foil = adb.Airfoil("foil",add)
    coords = foil.get_outline_points()
    x_af = coords[:,0]; y_af = coords[:,1]

    # create inner airfoil
    polyx,polyy = inner(x_af,y_af,info["t"])

    # initialize interpolation dictionary
    I = {}

    # create further thicknesses inner airfoils, as well as interpolations
    # create interpolation of inner surface
    I["I"] = interpolate(polyx,polyy,make_camberline=True)
    # create interpolation of outer surface
    I["O"] = interpolate(x_af,y_af,make_camberline=True)
    # outer tongue
    out_tx,out_ty = inner(x_af,y_af,   info["t"]+   info["mt"])
    I["ot"] = interpolate(out_tx,out_ty)
    # inner tongue
    inn_tx,inn_ty = inner(x_af,y_af,2.*info["t"]+   info["mt"])
    I["it"] = interpolate(inn_tx,inn_ty)
    # outer mouth
    out_mx,out_my = inner(x_af,y_af,2.*info["t"]+2.*info["mt"])
    I["om"] = interpolate(out_mx,out_my)

    # initialize x and y dictionaries
    x = {}; y = {}

    # create outer mouth structure
    x["mou"],y["mou"] = make_mouth(I["O"],I["I"],I["om"],info)

    # create tongue tip
    x["ton"],y["ton"] = tongue_tip(I["O"],I["I"],I["ot"],I["it"],info)

    # create spar box
    x["spr"],y["spr"] = spar(I["O"],info)

    # if shift desired, shift hole
    if info["ts"]:
        find_wall(I["O"],foil,info)

    # create TE hole
    x["teh"],y["teh"] = te_hole(I["O"],I["I"],info)

    # create inner run
    x["inn"],y["inn"] = inner_run(I["I"],info)

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

        # create dxf file
        filename = vals[run]["dxf file path"]+vals[run]["dxf file name"]
        xvals,yvals,zvals = make_dxf(u,v,filename,vals[run]["write dxf"])
        
        # if points file desired, do so
        fname = vals[run]["points file path"]+vals[run]["points file name"]
        if vals[run]["write points"]:
            points(u,v,fname)

        if vals[run]["show plot"]:
            # initialize colors and labels
            g = {}
            g["mou c"] = "#0073ff"; g["mou l"] = "Mouth"
            g["ton c"] = "#ff0000"; g["ton l"] = "Tongue"
            g["spr c"] = "#aa6c39"; g["spr l"] = "Spar"
            g["teh c"] = "#ff8c00"; g["teh l"] = "TE Hole"
            g["inn c"] = "#3cb371"; g["inn l"] = "Inner"
            g["out c"] =       "k"; g["out l"] = "Outer"

            # plot to be saved items
            for group in u:
                # color
                if vals[run]["all black"]:
                    c = "k"
                else:
                    c = g[group+" c"]
                
                for i in range(u[group].shape[0]):
                    if i == 0:
                        l = g[group+" l"]
                    else:
                        l = ""
                    plt.plot(u[group][i],v[group][i],c,label=l)

    if vals[run]["show plot"]:
        if vals[run]["show legend"]:
            plt.legend()
        plt.xlabel("x/c")
        plt.ylabel("y/c")
        plt.axis("equal")
        plt.show()

    return xvals, yvals, zvals
