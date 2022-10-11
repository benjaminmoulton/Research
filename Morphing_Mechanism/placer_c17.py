import numpy as np
import json
import svg
import os
import C17_geometry as C17
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
    
    # put units for cpx,cpy,cpz,chord in terms of inches
    dist[:,:4] *= 12.
    
    # create dxf file of lines that can be used to define planes
    # determine where the "mid point" lies -- index of last negative y value
    index = np.max(np.argwhere(dist[:,1] < 0.0)[:,0])
    
    # input a distribution for 0 if not there
    if not dist[index+1,1] == 0.0:
        dist = np.insert(dist,index+1,np.array([dist[index,0],0.0,\
            dist[index,2],dist[index,3],0.0]),axis=0)

    # set start from index
    start = index + 1
    
    # create dimension specifier
    if vals["svg"]["2D"]:
        dim = "_2D"
    else:
        dim = "_3D"
    
    # create a folder in the desired directory for output files
    folder = vals["svg"]["file path"]+vals["svg"]["folder name"]+dim+"/"

    return dist,start,folder

# create rist dictionary
def make_rist(dist,vals):
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
    ismorphing = [True]
    name = ["bay"]
    spec = ["bay"]

    # determine lengths
    # determine where the fan portion starts
    b.append(vals["ducted fan y-cg and width [in]"][0]-\
        vals["ducted fan y-cg and width [in]"][1]/2.)
    ismorphing.append(False); name.append("fan"); spec.append("fan")
    # determine where the fan portion ends
    b.append(vals["ducted fan y-cg and width [in]"][0]+\
        vals["ducted fan y-cg and width [in]"][1]/2.)
    ismorphing.append(True); name.append("servo00"); spec.append("sv0")

    # determine the servo lengths
    n = vals["actuators per semispan"]
    for i in range(1,n):
        b.append((vals["wing tip y start [in]"] - b[2])/n*i + b[2])
        ismorphing.append(True)
        name.append("servo{}".format(str(i).zfill(2)))
        spec.append("sv{}".format(str(i).zfill(1)))
    
    # add the last two points, the start and end of the wing tips
    b.append(vals["wing tip y start [in]"])
    ismorphing.append(False); name.append("wingtip"); spec.append("wtp")
    b.append(dist[-1,1])
    # save to vals
    b = np.array(b)

    # initialize important things
    vals["b"] = b
    vals["is morphing"] = ismorphing
    vals["name"] = name
    vals["spec"] = spec


    # create rist dictionary as ribs lengths 
    y_max = np.max(dist[:,1])
    y_thk = vals["rib thickness [in]"]
    num_ribs = int(y_max/y_thk)

    # initialize rist array
    y = np.linspace(0,y_thk*num_ribs,num=num_ribs+1) + y_thk / 2.
    rist = np.zeros((y.shape[0],5))

    # # run through each rib and create a rib at the mid rib point
    rist[:,1] = y
    rist[:,0] = vals["x(y)"](y)
    rist[:,2] = vals["z(y)"](y)
    rist[:,3] = vals["c(y)"](y)
    rist[:,4] = vals["d(y)"](y)
    
    return rist

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

# create airframe based on airfoil 
def airfoil_interp(dist,vals):
    # create airfoil
    add = {}
    add["geometry"] = {}
    add["geometry"]["type"] = "outline_points"
    add["geometry"]["outline_points"] = \
        vals["C17"]["airfoil dat file location"]

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
    x_cam = 0.25; y_cam = foil.get_camber(x_cam)

    # save to vals array
    vals["x_af"] = x; vals["y_af"] = y
    vals["x_cam"] = x_cam; vals["y_cam"] = y_cam

    # save to airfoil shape
    x -= x_cam
    y -= y_cam

    # intialize airfoil array
    xs = np.zeros((dist.shape[0],x.shape[0]))
    ys = np.zeros((dist.shape[0],x.shape[0]))
    zs = np.zeros((dist.shape[0],x.shape[0]))

    # run through and add the points to each point based on 3d transform
    for i in range(dist.shape[0]):
        x_loc = x * dist[i,3]; y_loc = y * dist[i,3]
        xs[i],ys[i],zs[i] = transform_3D(x_loc,y_loc,y_loc*0.,dist,i,i)
        
    # run through each point and create an interpolation
    itpX_y = []; itpZ_y = []
    for i in range(x.shape[0]):
        itpX_y.append(interp.interp1d(ys[:,i],xs[:,i],kind="cubic"))
        itpZ_y.append(interp.interp1d(ys[:,i],zs[:,i],kind="cubic"))
    
    # create function to find the enthickening value
    def enthickness(yval):
        # determine local normal thickness
        y_loc = y*vals["c(y)"](yval)
        t_loc = np.max(y_loc) - np.min(y_loc)
        
        # determine local "airfoil" shape
        z_loc = np.zeros((x.shape[0],))
        for i in range(x.shape[0]):
            z_loc[i] = itpZ_y[i](yval)

        # calculate actual planar thickness
        t_enth = np.max(z_loc) - np.min(z_loc)
            
        # return thickness ratio
        return t_enth / t_loc
    
    return enthickness

# create shapes
def make_shapes(rist,vals,enth,folder):
    print()

    # save former values
    vals["c0"] = vals["C17"]["chord [in]"]
    vals["t0"] = vals["C17"]["tongue start [in]"]
    vals["m0"] = vals["C17"]["mouth start [in]"]

    # set "max" vals
    c_ymax = vals["c(y)"](vals["wing tip y start [in]"])
    vals["C17"]["chord [in]"] = c_ymax
    vals["C17"]["tongue start [in]"] = c_ymax * vals["t0"]/vals["c0"]
    vals["C17"]["mouth start [in]"]  = c_ymax * vals["m0"]/vals["c0"]

    # try:
    #     C17.main(vals)
    # except:
    #     raise ValueError("C17 geometry cannot be created at wing tip y start")
    
    # report
    print("Creating {} shape...".format(vals["distributions file"].split(".")\
        [0]))
    print()

    # create airfoil foil shape using airfoil database
    add = {}
    add["geometry"] = {}
    add["geometry"]["type"] = "outline_points"
    add["geometry"]["outline_points"] = \
        vals["C17"]["airfoil dat file location"]
    foil = adb.Airfoil("foil",add)
    coords = foil.get_outline_points()
    x_af = coords[:,0]; y_af = coords[:,1]
    O = C17.interpolate(x_af,y_af,make_camberline=True)

    # initialize shapes lists
    ax = np.zeros((rist.shape[0],),dtype=np.ndarray)
    ay = np.zeros((rist.shape[0],),dtype=np.ndarray)
    az = np.zeros((rist.shape[0],),dtype=np.ndarray)

    # create Failed bool
    failed = False

    # run through each rib and create it
    for i in range(rist.shape[0]):
        vals["C17"]["chord [in]"] = rist[i,3]
        vals["C17"]["tongue start [in]"] = rist[i,3] * vals["t0"]/vals["c0"]
        vals["C17"]["mouth start [in]"]  = rist[i,3] * vals["m0"]/vals["c0"]
        vals["C17"]["enthicken"] = enth(rist[i,1])
        vals["C17"]["svg text"] = str(i).zfill(3)

        # try to make C17 shape
        try:
            if not failed:
                # determine shift for TE hole 
                # find current hole
                info = {}
                info["t"]     = vals["C17"]["shell thickness [mm]"] / 25.4\
                    / rist[i,3]
                info["td"]    = vals["C17"]["TE hole diameter [mm]"] / 25.4\
                    / rist[i,3]
                info["tehxs"] = 0.0
                # determine new circle center
                curr_x = C17.find_wall(O,foil,info) + info["t"] + info["td"]/2.
                curr_y = O["c"](curr_x)
                print("curr ",curr_y,O["c"](0.5),end=" ")
                # shift to position
                curr_x -= vals["x_cam"]; curr_y -= vals["y_cam"]
                curr_x *= -rist[i,3]; curr_y *= rist[i,3]
                curr_x += rist[i,0]; curr_y += rist[i,2]
                # find next hole
                info["t"]     = vals["C17"]["shell thickness [mm]"] / 25.4\
                    / rist[i+1,3]
                info["td"]    = vals["C17"]["TE hole diameter [mm]"] / 25.4\
                    / rist[i+1,3]
                info["tehxs"] = 0.0
                # determine new circle center
                next_x = C17.find_wall(O,foil,info) + info["t"] + info["td"]/2.
                next_y = O["c"](next_x)
                print("next ",next_y,O["c"](0.5))
                next_x -= vals["x_cam"]; next_y -= vals["y_cam"]
                next_x *= -rist[i+1,3]; next_y *= rist[i+1,3]
                next_x += rist[i+1,0]; next_y += rist[i+1,2]

                # save to vals array
                vals["C17"]["TE hole x shift [in]"] = -(curr_x - next_x)
                vals["C17"]["TE hole y shift [in]"] = curr_y - next_y

                # run C17
                xsvg,ysvg,zsvg,text = C17.main(vals,"C17")
            else:
                raise ValueError               
        except ValueError:
            failed = True
            print("Airfoil c = {:>7.4f} [in] {:>3d}/{}".\
                format(rist[i,3],i,rist.shape[0]-1))
            # raise ValueError("FAIL")
            
            # save svg arrays
            xsvg = np.array([vals["x_af"] * rist[i,3]])
            ysvg = np.array([vals["y_af"] * rist[i,3]])
            zsvg = ysvg * 0.

            # create text dictionary
            text = {}
            text["text"] = vals["C17"]["svg text"]

            # pass in x y of start location
            text["start"] = [0.,0]
        else:
            print("C17 geo c = {:>7.4f} [in] {:>3d}/{}".\
                format(rist[i,3],i,rist.shape[0]-1))
        
        # create svg file
        filename = folder + vals["C17"]["svg text"]
        svg.svg(filename,xsvg,ysvg,geometry=vals["C17"]["svg type"],text=text)

        # create 3d arrays
        x3d = xsvg - vals["x_cam"] * rist[i,3]
        y3d = zsvg * 1.
        z3d = ysvg - vals["y_cam"] * rist[i,3]
        
        # shift back by control point value 
        x3d += rist[i,0]
        y3d += rist[i,1]
        z3d += rist[i,2]
        if i != 0:
            plt.axis("equal")
            for j in range(x3d.shape[0]):
                plt.plot(x3d[j],z3d[j],"k")
                plt.plot(px[j],pz[j],"b")
            plt.show()
        px = x3d; py = y3d; pz = z3d
        
        # add to shapes lists
        ax[i] = x3d; ay[i] = y3d; az[i] = z3d
    print()
    
    # save to list
    apts = [ax,ay,az]
    return apts

# plot half
def half_plotter(a,mins,maxs,vals,dist,axe,face,ymult=1.):
    # release values
    ax,ay,az = a; xmin,ymin,zmin = mins; xmax,ymax,zmax = maxs

    # run through each airfoil and plot it
    ind = 1; col = "k"
    for i in range(len(ax)):
        if vals["b"][ind] < dist[i,1]:
            if ind == (vals["b"].shape[0]-2):
                col = "0.3"
            elif ind == 1:
                col = "m"
            else:
                if col == "k":
                    col = "g"
                else:
                    col = "k"
            ind += 1
        # run through each line segment and plot it
        for j in range(ax[i].shape[0]):
            axe.plot(ax[i][j],ymult*ay[i][j],-az[i][j],c=col,zdir= face)

            # solve min and maxes
            axmin = np.min(ax[i][j]); axmax = np.max(ax[i][j])
            aymin = np.min(ymult*ay[i][j]); aymax = np.max(ymult*ay[i][j])
            azmin = np.min(-az[i][j]); azmax = np.max(-az[i][j])

            if axmin < xmin: xmin = axmin
            if axmax > xmax: xmax = axmax
            if aymin < ymin: ymin = aymin
            if aymax > ymax: ymax = aymax
            if azmin < zmin: zmin = azmin
            if azmax > zmax: zmax = azmax
    
    mins = [xmin,ymin,zmin]; maxs = [xmax,ymax,zmax]

    return mins, maxs

# plot full
def plotter(apts,vals,dist):
    print("Plotting...")

    # withdraw points
    ax,ay,az = apts

    # plot a 3d plot of the wing
    fig = plt.figure()
    axe = fig.add_subplot(111, projection='3d')
    face = 'z'
    
    # initialize min and max vars
    xmin = ax[0][0][0]; xmax = ax[0][0][0]
    ymin = ay[0][0][0]; ymax = ay[0][0][0]
    zmin = az[0][0][0]; zmax = az[0][0][0]

    # plot half
    a = [ax,ay,az]; mins = [xmin,ymin,zmin]; maxs = [xmax,ymax,zmax]
    mins,maxs = half_plotter(a,mins,maxs,vals,dist,axe,face)
    
    # if the whole design is to be plotted:
    if not vals["plot"]["halfspan"]:
        mins,maxs = half_plotter(a,mins,maxs,vals,dist,axe,face,ymult=-1)
    
    # release values
    xmin,ymin,zmin = mins; xmax,ymax,zmax = maxs
    
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
    print(x_lims,y_lims,z_lims)

    # set limits
    axe.set_xlim3d(x_lims[1], x_lims[0])
    axe.set_ylim3d(y_lims[0], y_lims[1])
    axe.set_zlim3d(z_lims[1], z_lims[0])

    # # plot where ducted fans go# plot a 3d plot of the wing
    # axe.plot([0,0],[ 4.5, 8.5],[0,0],c="m",zdir=face)
    # axe.plot([0,0],[-4.5,-8.5],[0,0],c="m",zdir=face)

    axe.view_init(20,-30)
    plt.show()

    return

def main(input_file):
    # read in json file
    all_vals = read_json(input_file)

    # find runner and run
    for run in all_vals["run"]:
        # set vals dict
        vals = all_vals[run]

        # get dist dictionary
        dist,start,folder = get_dist(vals)
        # 0 - cpx, 1 - cpy, 2 - cpz, 3 - chord, 4 - dihedral
        # body fixed coordinates

        # create airfoil interpolation
        enth = airfoil_interp(dist,vals)

        # create rist array for ribs
        rist = make_rist(dist,vals)

        # create airfoil interpolation
        enth = airfoil_interp(dist,vals)

        # create C17 shapes
        apts = make_shapes(rist,vals,enth,folder)
        
        # plot the data
        if vals["plot"][""]:
            plotter(apts,vals,rist)
            
    return

input_file = "placer_c17_input.json"
main(input_file)
