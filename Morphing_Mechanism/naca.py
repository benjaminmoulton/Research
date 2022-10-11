import numpy as np
import airfoil_db as adb


# create airfoil shape
add = {}
add["geometry"] = {}
add["geometry"]["type"] = "NACA"
add["geometry"]["NACA"] = "2406"

# initialize airfoil database
foil = adb.Airfoil("foil",add)
coords = foil.get_outline_points(export="/home/ben/Desktop/NACA 2406.dat")
x = coords[:,0]; y = coords[:,1]
# print(x,y)

# with open("/home/ben/Desktop/NACA 0006.dat","w") as f:
#     for i in range(x.shape[0]):
#         f.write("{:<8.6f} {:<8.6f}\n".format(x[i],y[i]))