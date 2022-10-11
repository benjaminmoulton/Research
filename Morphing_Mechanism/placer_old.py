import numpy as np
import json
import dxf

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

def main(input_file):
    # read in json file
    vals = read_json(input_file)

    # read in distributions file
    dist = read_dist(vals["distributions file"])

    # remove all but control point, chord, and dihedral
    dist = dist[:,1:7]
    dist = np.delete(dist,4,axis=1)
    # 0 - cpx, 1 - cpy, 2 - cpz, 3 - chord, 4 - dihedral
    # body fixed coordinates

    # where the y value is negative (ie, left wing)
    # make dihedral negative
    for i in range(dist.shape[0]):
        if dist[i,1] < 0.0:
            dist[i,4] *= -1.

    # read in shape file
    name = vals["airfoil file"].split(".")
    if name[-1] == "c14":
        fx,fy,fz = read_c14(vals["airfoil file"])
    else:
        fx,fy,fz = read_text(vals["airfoil file"])
    # aero coordinates, ie.
    # x_b = - x_a
    # y_b =   z_a
    # z_b = - y_a
    
    # initialize dxf lists
    xd = []; yd = []; zd = []

    # get max x value to resize fx fy shape
    fx_max = 1.0

    # read through each point
    for i in range(dist.shape[0]):
        # read through each airfoil shape
        for j in range(fx.shape[0]):
            # initialize new points list
            x_new = np.zeros(fx[j].shape)
            y_new = np.zeros(fx[j].shape)
            z_new = np.zeros(fx[j].shape)

            # read through each line in the airfoil shape
            for k in range(fx[j].shape[0]):
                # x value
                # resize to unity chord
                x_pt = fx[j][k]/ fx_max
                # shift fx,fy shape values so 0 is at quarter chord
                x_pt -= 0.25
                # resize to new chord, negative for new coordinate frame
                x_pt *= -dist[i,3]
                # rotate with dihedral
                # not needed
                # shift to cp location
                x_pt += dist[i,0]
                
                # y value
                y_pt = fz[j][k]/ fx_max
                
                # z value
                # resize to unity chord
                z_pt = fy[j][k]/ fx_max
                # resize to new chord
                z_pt *= -dist[i,3]

                # rotate z and y with dihedral and save to points
                z_rot = z_pt*np.cos(dist[i,4]) - y_pt*np.sin(dist[i,4])
                y_rot = y_pt*np.cos(dist[i,4]) + z_pt*np.sin(dist[i,4])
                z_pt = z_rot; y_pt = y_rot

                # shift y and z to cp location
                y_pt += dist[i,1]
                z_pt += dist[i,2]

                # save to arrays
                x_new[k] = x_pt; y_new[k] = y_pt; z_new[k] = z_pt
            
            # save new point arrays to xd lists
            xd.append(x_new); yd.append(y_new); zd.append(z_new)
    
    # turn dxf lists into numpy arrays
    xd = np.array(xd); yd = np.array(yd); zd = np.array(zd)

    # turn into a dxf file
    dxf.dxf(vals["dxf file name"],xd,yd,zd,vals["write as"])
            
    return


input_file = "placer_input.json"
main(input_file)