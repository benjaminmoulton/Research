import numpy as np
from datetime import datetime as dt

def eps_spline(filename,x,y,z,box = {},file_info = []):

    # create file
    f = open(filename + ".eps", "w+")

    # write eps file opening statements
    f.write("\%!PS-Adobe-3.0 EPSF-3.0\n")
    f.write("\%\%BoundingBox: ")
    f.write("{:<7.4f} {:<7.4f} ".format(box["llx"],box["lly"]))
    f.write("{:<7.4f} {:<7.4f}\n".format(box["urx"],box["ury"]))

    # run through file info and input it into the file
    for info in file_info:
        f.write("\%\%" + info + "\n")

    # write file creation time
    time = str(dt.now())
    f.write("\%\%CreationDate: " + time + "\n")

    # end comments
    f.write("\%\%EndComments\n")

    # # determine the number of splines to create
    # num_complines = x.shape[0]

    # # if the y and x spline counts are not equal in number, throw error
    # if y.shape[0] != num_complines:
    #     raise ValueError('y array must have same number of splines as x array')

    # # if the z and x spline counts are not equal in number, throw error
    # if z.shape[0] != num_complines:
    #     raise ValueError('z array must have same number of splines as x array')
    
    # # cycle through each compound line
    # for i in range(num_complines):

    #     # determine the number of points in the current spline
    #     num_lines = x[i].shape[0]

    #     # write dxf file spline opening statement
    #     f.write("0\nSPLINE\n8\nLayer\n100\nAcDbSpline\n71\n3\n74\n4\n")

    #     # if the y and x spline points are not equal in number, throw error
    #     if y[i].shape[0] != num_lines:
    #         raise ValueError('y array spline %d must have same number of points as x array' % i)

    #     # if the z and x spline points are not equal in number, throw error
    #     if z[i].shape[0] != num_lines:
    #         raise ValueError('z array spline %d must have same number of points as x array' % i)

    #     # determine start and end tangents (in components form)
    #     sx = x[i][1] - x[i][0]
    #     sy = y[i][1] - y[i][0]
    #     sz = z[i][1] - z[i][0]

    #     ex = x[i][-1] - x[i][-2]
    #     ey = y[i][-1] - y[i][-2]
    #     ez = z[i][-1] - z[i][-2]

    #     # cycle through line in each compound line
    #     for j in range(num_lines):

    #         # write next spline x coordinate
    #         f.write("11\n%.8f\n" % x[i][j])

    #         # write next spline y coordinate
    #         f.write("21\n%.8f\n" % y[i][j])

    #         # write next spline z coordinate
    #         f.write("31\n%.8f\n" % z[i][j])

    #         # # write line color (always black=0) # blue=5
    #         # f.write("62\n0\n")
        
    #     # add in start and end tangents
    #     # write start spline tangent x-component
    #     f.write("12\n%.8f\n" % sx)
    #     # write start spline tangent y-component
    #     f.write("22\n%.8f\n" % sy)
    #     # write start spline tangent z-component
    #     f.write("32\n%.8f\n" % sz)

    #     # write end spline tangent x-component
    #     f.write("13\n%.8f\n" % ex)
    #     # write end spline tangent y-component
    #     f.write("23\n%.8f\n" % ey)
    #     # write end spline tangent z-component
    #     f.write("33\n%.8f\n" % ez)

    # write eps file closing statement
    f.write("\%\%EOF\n\%\%EndDocument\n")
    
    # close file
    f.close()

    # # move file to desktop
    # pm.move_file_to_Desktop(filename + ".dxf")

    return

box = {}
box["llx"] = 0
box["lly"] = 0
box["urx"] = 12
box["ury"] = 12
eps_spline("/home/ben/Desktop/epser",0,0,0,box)