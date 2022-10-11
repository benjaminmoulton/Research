import numpy as np
import dxf

# read in manta_values text file
filename = "/home/ben/Desktop/foilshapes/Manta_values.txt"
c = []; x = []; y = []; z = []
with open(filename,"r") as f:
    f.readline()
    lst = []
    for line in f:
        lst.append([ float(point) for point in line.split()])
    c = np.array([ point[0] for point in lst])
    x = np.array([ point[1] for point in lst])
    y = np.array([ point[2] for point in lst])
    z = np.array([ point[3] for point in lst])

# print(x)
# find number of E335 airfoils
num = x.shape[0] - 1
# initialize 2d array
x_dxf = []; y_dxf = []; z_dxf = []
epplerfilename = "/home/ben/Desktop/foilshapes/E335.txt"
for i in range(num):
    xe = []; ye = []; ze = []
    with open(epplerfilename,"r") as f:
        f.readline()
        lst = []
        for line in f:
            lst.append([ float(point) for point in line.split()])
        xe = np.array([ point[0]*c[i]+x[i] for point in lst])
        ye = np.array([ point[1]*c[i]+y[i] for point in lst])
        ze = np.array([z[i] for j in range(xe.shape[0])])
        f.close()
    x_dxf.append(xe)
    y_dxf.append(ye)
    z_dxf.append(ze)

# print(x_dxf,z_dxf)

# find nums of NACA 0012 airfoils
nacanum = 1
nacafilename = "/home/ben/Desktop/foilshapes/NACA 0012.txt"
for i in range(nacanum):
    xe = []; ye = []; ze = []
    with open(nacafilename,"r") as f:
        f.readline()
        lst = []
        for line in f:
            lst.append([ float(point) for point in line.split()])
        xe = np.array([ point[0]*c[num]+x[num] for point in lst])
        ye = np.array([ point[1]*c[num]+y[num] for point in lst])
        ze = np.array([z[num] for j in range(xe.shape[0])])
        f.close()
    x_dxf.append(xe)
    y_dxf.append(ye)
    z_dxf.append(ze)

# turn into numpy arrays
x_dxf = np.array(x_dxf)
y_dxf = np.array(y_dxf)
z_dxf = np.array(z_dxf)

# print(z_dxf)
# a = [0,1,2,3,4,5,6,7,8,9,10]
# print(a[-1])

dxfname = "MANTA"
dxf.dxf(dxfname,x_dxf,y_dxf,z_dxf,True)