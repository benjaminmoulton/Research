import numpy as np
import json
import dxf
import os
import C14_geoPRECISE as C14
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

    return dist,start

# create plane lines
def make_planes(dist,start,vals):
    # initialize arrays for planar lines
    n = (dist.shape[0]-start) * 2 - 1
    xplane = np.zeros((n,),dtype=np.ndarray)
    yplane = np.zeros((n,),dtype=np.ndarray)
    zplane = np.zeros((n,),dtype=np.ndarray)

    # run through each point and create a line to be used to create the planes
    j = 0
    for i in range(start,dist.shape[0]):
        # initialize line
        xplane[j] = np.zeros((2,))
        yplane[j] = np.zeros((2,))
        zplane[j] = np.zeros((2,))

        # set last point as the planar control point
        yplane[j][1] = dist[i,1]
        xplane[j][1] = dist[i,2]

        # set first point as one unit away from the control point
        # at the dihedral angle
        line = vals["plane line length"]
        yplane[j][0] = yplane[j][1] - np.sin(dist[i,4]) * line
        xplane[j][0] = xplane[j][1] - np.cos(dist[i,4]) * line
        
        # initialize line
        xplane[j+1] = np.zeros((2,))
        yplane[j+1] = np.zeros((2,))
        zplane[j+1] = np.zeros((2,))

        if not j == 0:
            # set last point as the planar control point
            mn = -1/np.tan(dist[i,4])
            dy = ( (line/20.)**2./(1+mn**2.) )**0.5
            yplane[j+1][0] = yplane[j][0] + dy*mn
            xplane[j+1][0] = xplane[j][0] + dy

            # set first point as one unit away from the control point
            # at the dihedral angle
            yplane[j+1][1] = yplane[j][0] - dy*mn
            xplane[j+1][1] = xplane[j][0] - dy

            # increase index
            j += 2
        else:
            # increase index
            j += 1
    
    # reverse across the x axis
    xplane *= -1.
    
    # create dimension specifier
    if vals["dxf"]["2D"]:
        dim = "_2D"
    else:
        dim = "_3D"
    
    # create a folder in the desired directory for output files
    folder = vals["dxf"]["file path"]+vals["dxf"]["folder name"]+dim+"/"
    if not os.path.isdir(folder[:-1]):
        os.mkdir(folder[:-1])

    # save folder path
    vals["C14"]["dxf file path"] = folder
    
    # create a dxf file of the lines
    dxf.dxf(folder+"RightPlane2D",xplane,yplane,zplane,geometry="line")

    return

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

# create 3d, 2d, and guidecurves
def make_shapes(dist,start,vals):
    # save former values
    c0 = vals["C14"]["chord [in]"]
    t0 = vals["C14"]["tongue start [in]"]
    m0 = vals["C14"]["mouth start [in]"]

    # initialize counter index
    j = 0
    first_airfoil = True

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
    x_cam = 0.25; y_cam = foil.get_camber(x_cam)

    # initialize pointer list
    c14_guides = []; af_guides = []
    ax = []; ay = []; az = []
    # initialize first airfoil index
    index_saf = 1000

    # create a dxf file at each plane for the morphing shape
    for i in range(start,dist.shape[0]): # start+33,dist.shape[0]): # 
        # initilize info dictionary to calculate C14

        # set values
        vals["C14"]["chord [in]"] = dist[i,3]
        vals["C14"]["tongue start [in]"] = dist[i,3] * t0/c0
        vals["C14"]["mouth start [in]"]  = dist[i,3] * m0/c0
        vals["C14"]["dxf file name"] = "plane " + str(j).zfill(2)
        
        try:
            ix,iy,iz,isx,isy,isz,guides,splitdict = C14.main(vals)
        except:
            # initialize x y z arrays
            ix = (np.array([x]) -x_cam) * dist[i,3] # 
            iy = (np.array([y]) -y_cam) * dist[i,3] # 
            iz = iy * 0.

            # if this is the first failed airfoil, make an airfoil of the last
            if first_airfoil:
                isx = (np.array([x]) -x_cam) * dist[i-1,3] # 
                isy = (np.array([y]) -y_cam) * dist[i-1,3] # 
                isz = iy * 0.
                first_airfoil = False
                index_saf = i

                # initialize guide curve
                af_guide = np.zeros((2,3))
                imax = np.argwhere(np.max(isx[0]) == isx[0])[0][0]
                imin = np.argwhere(np.min(isx[0]) == isx[0])[0][0]
                af_guide[0,0] = isx[0][imax]; af_guide[1,0] = isx[0][imin]
                af_guide[0,1] = isy[0][imax]; af_guide[1,1] = isy[0][imin]

                # translate c14 guides
                # mirror across x axis
                af_guide[:,0] *= -1

                # flip iy and iz
                ia = af_guide[:,1] * -1.
                af_guide[:,1] = af_guide[:,2] * 1.
                af_guide[:,2] = ia * 1.

                # rotate by the dihedral
                ry = af_guide[:,1]*np.cos(dist[i-1,4]) + \
                    af_guide[:,2]*np.sin(dist[i-1,4])
                rz = af_guide[:,2]*np.cos(dist[i-1,4]) - \
                    af_guide[:,1]*np.sin(dist[i-1,4])
                af_guide[:,2] = rz; af_guide[:,1] = ry

                # shift to the c/4 point
                af_guide[:,0] += dist[i-1,0]
                af_guide[:,1] += dist[i-1,1]
                af_guide[:,2] += dist[i-1,2]
                
                # add to list
                af_guides.append(af_guide)
            else: # otherwise, make them zeros.
                isx = 0.; isy = 0; isz = 0
            
            # initialize guide curve
            af_guide = np.zeros((2,3))
            imax = np.argwhere(np.max(ix[0]) == ix[0])[0][0]
            imin = np.argwhere(np.min(ix[0]) == ix[0])[0][0]
            af_guide[0,0] = ix[0][imax]; af_guide[1,0] = ix[0][imin]
            af_guide[0,1] = iy[0][imax]; af_guide[1,1] = iy[0][imin]

            # translate c14 guides
            # mirror across x axis
            af_guide[:,0] *= -1

            # flip iy and iz
            ia = af_guide[:,1] * -1.
            af_guide[:,1] = af_guide[:,2] * 1.
            af_guide[:,2] = ia * 1.

            # rotate by the dihedral
            ry = af_guide[:,1]*np.cos(dist[i,4]) + \
                af_guide[:,2]*np.sin(dist[i,4])
            rz = af_guide[:,2]*np.cos(dist[i,4]) - \
                af_guide[:,1]*np.sin(dist[i,4])
            af_guide[:,2] = rz; af_guide[:,1] = ry

            # shift to the c/4 point
            af_guide[:,0] += dist[i,0]
            af_guide[:,1] += dist[i,1]
            af_guide[:,2] += dist[i,2]
            
            # add to list
            af_guides.append(af_guide)

            # run through each line and spit out how many points is in the file
            num_points = 0
            for n in range(ix.shape[0]):
                num_points += ix[n].shape[0]
            # report file study
            print("Airfoil c = {:>7.4f} [in] with {:>4d} points {:>2d}/{}".\
                format(dist[i,3],num_points,j,dist.shape[0]-start-1))
            extra = " airfoil"; out = " airfoil"
        else:
            # run through each line and spit out how many points is in the file
            num_points = 0
            for n in range(ix.shape[0]):
                num_points += ix[n].shape[0]
            # report file study
            print("C14 geo c = {:>7.4f} [in] with {:>4d} points {:>2d}/{}".\
                format(dist[i,3],num_points,j,dist.shape[0]-start-1))
            extra = ""; out = " out"

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
            
            # add to list
            c14_guides.append(guides)
        
        ### determine 2D values of the shape
        # mirror airfoil across y axis
        ix2d = -1. * ix; isx2d = -1. * isx

        # rotate by 90 degrees around origin
        rx2d = iy*-1.
        ry2d = ix2d*1.
        ix2d = rx2d; iy2d = ry2d
        rsx2d = isy*-1.
        rsy2d = isx2d*1.
        isx2d = rsx2d; isy2d = rsy2d

        # save index for future use
        index = i
        if not first_airfoil:
            index = i-1

        # determine how much to shift down due to dihedral
        # shift back x and y
        if not dist[i,4] == 0.0:
            # solve for shift value
            xshift  = get_x_shift(i,dist)
            xsshift = get_x_shift(index,dist)

            # add x shift and y shift to proper arrays
            ix2d  += xshift
            iy2d  += dist[i,0]
            isy2d += dist[index,0]
            isx2d += xsshift
                
        ### determine 3D values of the shape
        # mirror across x axis
        ix3d = -1. * ix; isx3d = -1. * isx

        # flip iy and iz
        ia3d = iy * -1.; iy3d = iz * 1.; iz3d = ia3d * 1.
        isa3d = isy * -1.; isy3d = isz * 1.; isz3d = isa3d * 1.

        # rotate by the dihedral
        ry3d = iy3d*np.cos(dist[i,4]) + iz3d*np.sin(dist[i,4])
        rz3d = iz3d*np.cos(dist[i,4]) - iy3d*np.sin(dist[i,4])
        iz3d = rz3d; iy3d = ry3d
        rsy3d = isy3d*np.cos(dist[index,4]) + isz3d*np.sin(dist[index,4])
        rsz3d = isz3d*np.cos(dist[index,4]) - isy3d*np.sin(dist[index,4])
        isz3d = rsz3d; isy3d = rsy3d

        # shift to the c/4 point
        ix3d += dist[i,0]; isx3d += dist[index,0]
        iy3d += dist[i,1]; isy3d += dist[index,1]
        iz3d += dist[i,2]; isz3d += dist[index,2]
        ax.append(ix3d); ay.append(iy3d); az.append(iz3d)

        # save the appropriate file         
        if vals["dxf"]["2D"]:
            ix = ix2d; iy = iy2d
            isx = isx2d; isy = isy2d
        else:
            ix = ix3d; iy = iy3d; iz = iz3d
            isx = isx3d; isy = isy3d; isz = isz3d

        # determine file location
        file_location = vals["C14"]["dxf file path"]+\
            vals["C14"]["dxf file name"]
        
        # specify dxf file type
        if first_airfoil:
            typ = "C14 "
        else:
            typ = "Airfoil "
        airfoil_type = vals["C14"]["airfoil dat file location"].split("/")[-1]\
            .split(".")[0]
        typ = airfoil_type + " " + typ

        # arrange file info for dxf file
        file_info = [
            typ + "geometry with chord length {:>7.4f} [in]".format(dist[i,3]),
            "control point at ({:>7.4f},{:>7.4f},{:>7.4f}) [in]".format(\
                dist[i,0],dist[i,1],dist[i,2])
        ]
        # write to dxf file if fisrt, last, or first airfoil
        if i == start or i == dist.shape[0]-1:
            dxf.dxf(file_location + extra,ix,iy,iz,file_info,geometry="spline")
        elif i == index_saf:
            dxf.dxf(prev_location+p_extra,px,py,pz,prev_info,geometry="spline")

        # if first airfoil, fix name
        if not first_airfoil:
            file_location = vals["C14"]["dxf file path"]+\
                "plane " + str(j-1).zfill(2)
            
        # write dxf file for outside shape
        if not (type(isx) == np.float64 and type(isy) == np.float64)\
            and not (type(isx) == float and type(isy) == float) and\
                (i == start or i == index_saf):
            dxf.dxf(file_location+out,isx,isy,isz,file_info,geometry="spline")
            if i == index_saf:
                dxf.dxf(prev_location+p_out,psx,psy,psz,prev_info,\
                    geometry="spline")

        # increase count
        j += 1

        # save previous values
        px = ix; py = iy; pz = iz; psx = isx; psy = isy; psz = isz
        prev_location = file_location; p_out = out; prev_info = file_info
        p_extra = extra
    
    # turn guides array into numpy array
    c14_guides = np.array(c14_guides)

    # turn airfoil guides array into numpy array
    af_guides = np.array(af_guides)

    return ax,ay,az,c14_guides,af_guides

# create c14 guidecurves dxf file
def c14_guides_dxf(c14_guides,vals):
    # NOTE: 
    # i -- plane #
    # j -- hole #
    # k -- guide curve
    # l -- x, y, or z value
    
    # run through each "hole" and make guide curves
    for j in range(c14_guides[0].shape[0]):
        # initialize extruder guide curves
        x_dxf = np.zeros((c14_guides[0][j].shape[0],),dtype=np.ndarray)
        y_dxf = np.zeros((c14_guides[0][j].shape[0],),dtype=np.ndarray)
        z_dxf = np.zeros((c14_guides[0][j].shape[0],),dtype=np.ndarray)

        # run through each guide curve
        for k in range(c14_guides[0][j].shape[0]): # 1): # 
            # initialize plane lists
            x_list = []; y_list = []; z_list = []

            # run through each plane
            for i in range(c14_guides.shape[0]):
                x_list.append(c14_guides[i][j][k][0])
                y_list.append(c14_guides[i][j][k][1])
                z_list.append(c14_guides[i][j][k][2])

            # add to dxf as numpy arrays
            x_dxf[k] = np.array(x_list)
            y_dxf[k] = np.array(y_list)
            z_dxf[k] = np.array(z_list)

        # save as a dxf file
        file_location = vals["C14"]["dxf file path"] + \
            "hole_" + str(j).zfill(2) + "_gc"
        dxf.dxf(file_location,x_dxf,y_dxf,z_dxf,geometry="spline")
    
    return

# create airfoil wing tip guidecurves dxf file
def af_guides_dxf(af_guides,vals):
    # NOTE: 
    # i -- plane #
    # k -- guide curve
    # l -- x, y, or z value
    
    # initialize extruder guide curves
    x_dxf = np.zeros((af_guides[0].shape[0],),dtype=np.ndarray)
    y_dxf = np.zeros((af_guides[0].shape[0],),dtype=np.ndarray)
    z_dxf = np.zeros((af_guides[0].shape[0],),dtype=np.ndarray)

    # run through each guide curve
    for k in range(af_guides[0].shape[0]): # 1): # 
        # initialize plane lists
        x_list = []; y_list = []; z_list = []

        # run through each plane
        for i in range(af_guides.shape[0]):
            x_list.append(af_guides[i][k][0])
            y_list.append(af_guides[i][k][1])
            z_list.append(af_guides[i][k][2])

        # add to dxf as numpy arrays
        x_dxf[k] = np.array(x_list)
        y_dxf[k] = np.array(y_list)
        z_dxf[k] = np.array(z_list)

    # save as a dxf file
    file_location = vals["C14"]["dxf file path"] + \
        "hole_af_gc"
    dxf.dxf(file_location,x_dxf,y_dxf,z_dxf,geometry="spline")

    return

# plot half
def half_plotter(a,mins,maxs,c14_guides,af_guides,vals,axe,face,ymult=1.0):
    # release values
    ax,ay,az = a; xmin,ymin,zmin = mins; xmax,ymax,zmax = maxs

    # run through each airfoil and plot it
    for i in range(len(ax)):
        # run through each line segment and plot it
        for j in range(ax[i].shape[0]):
            axe.plot(ax[i][j],ymult*ay[i][j],-az[i][j],c="k",zdir= face)

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
    
    if vals["plot"]["guide curves"]:
        # plot guide curves
        for j in range(c14_guides[0].shape[0]): # 1): # 
            for k in range(c14_guides[0][j].shape[0]):
                # initialize plane lists
                x_list = []; y_list = []; z_list = []

                # run through each plane
                for i in range(c14_guides.shape[0]):
                    x_list.append(c14_guides[i][j][k][0])
                    y_list.append(ymult*c14_guides[i][j][k][1])
                    z_list.append(-c14_guides[i][j][k][2])
                
                # set guidecurves blue if desired
                if vals["plot"]["guides solid color"]:
                    co = "b"
                else:
                    co = ""

                # plot
                axe.plot(x_list,y_list,z_list,c=co,zdir=face)
        
        for k in range(af_guides[0].shape[0]):
            # initialize plane lists
            x_list = []; y_list = []; z_list = []

            # run through each plane
            for i in range(af_guides.shape[0]):
                x_list.append(af_guides[i][k][0])
                y_list.append(ymult*af_guides[i][k][1])
                z_list.append(-af_guides[i][k][2])
            
            # set guidecurves blue if desired
            if vals["plot"]["guides solid color"]:
                co = "r"
            else:
                co = ""

            # plot
            axe.plot(x_list,y_list,z_list,c=co,zdir=face)
    
    mins = [xmin,ymin,zmin]; maxs = [xmax,ymax,zmax]

    return mins, maxs

# plot full
def plotter(ax,ay,az,c14_guides,af_guides,vals):
    print("Plotting...")

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
    mins,maxs = half_plotter(a,mins,maxs,c14_guides,af_guides,vals,axe,face)
    
    # if the whole design is to be plotted:
    if not vals["plot"]["halfspan"]:
        mins,maxs = half_plotter(a,mins,maxs,c14_guides,af_guides,\
            vals,axe,face,ymult=-1)
    
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

    # set limits
    axe.set_xlim3d(x_lims[1], x_lims[0])
    axe.set_ylim3d(y_lims[0], y_lims[1])
    axe.set_zlim3d(z_lims[1], z_lims[0])

    axe.view_init(20,-30)
    plt.show()

    return

def main(input_file):
    # read in json file
    vals = read_json(input_file)

    # get dist dictionary
    dist,start = get_dist(vals)
    # 0 - cpx, 1 - cpy, 2 - cpz, 3 - chord, 4 - dihedral
    # body fixed coordinates

    # create lines for planes
    make_planes(dist,start,vals)

    # create 3d, 2d, and guidecurves
    ax,ay,az,c14_guides,af_guides = make_shapes(dist,start,vals)

    # create c14 guidecurves dxf file
    c14_guides_dxf(c14_guides,vals)

    # create airfoil wing tip guidecurves dxf file
    af_guides_dxf(af_guides,vals)
    
    # plot the data
    if vals["plot"][""]:
        plotter(ax,ay,az,c14_guides,af_guides,vals)
        
    return


input_file = "placer_input.json"
main(input_file)