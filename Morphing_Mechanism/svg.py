import numpy as np
from datetime import datetime as dt
import os

def svg_polyline(filename,x,y,file_info,text,adobe):
    # create file
    f = open(filename + ".svg", "w+")

    # determine the number of lines to create
    num_complines = x.shape[0]

    # if the y and x spline counts are not equal in number, throw error
    if y.shape[0] != num_complines:
        raise ValueError('y array must have same number of lines as x array')

    # initialize max and min x and y values
    xmax = x[0][0]; xmin = x[0][0]; ymax = y[0][0]; ymin = y[0][0]

    # run through and determine the size of the file
    for i in range(num_complines):
        # determine the number of points in the current spline
        num_lines = x[i].shape[0]

        # if the y and x spline points are not equal in number, throw error
        if y[i].shape[0] != num_lines:
            raise ValueError('y array spline %d must have same' % i + \
                ' number of points as x array')
        
        # run through lines
        for j in range(num_lines):
            # determine local maxima and minima
            exmax = np.max(x[i][j]); exmin = np.min(x[i][j])
            eymax = np.max(y[i][j]); eymin = np.min(y[i][j])

            # compare and replace if necessary
            if exmax > xmax: xmax = exmax
            if exmin < xmin: xmin = exmin
            if eymax > ymax: ymax = eymax
            if eymin < ymin: ymin = eymin
    
    # move x and y such that they are a tenth of an inch within bounds
    ten = 0.1
    x += ten - xmin
    y -= ten + ymax

    # flip across y axis since coordinates in eps are opposite
    y *= -1.

    # translate to pixels
    if adobe:
        in_to_pixels = 72.0
    else:
        in_to_pixels = 96.0
    x *= in_to_pixels
    y *= in_to_pixels

    # save width and height of final shape
    wd = int((xmax - xmin + 2.*ten) * in_to_pixels)
    ht = int((ymax - ymin + 2.*ten) * in_to_pixels)

    # open svg file statement
    f.write("<svg width=\"{:d}px\" height=\"{:d}px\">\n".format(wd,ht))

    # create title, write to file
    title = filename.split("/")[-1]
    f.write("\t<title> {} </title>\n".format(title))
    
    # run through file info and input it into the file
    for info in file_info:
        f.write("\t<!--" + info + "-->\n")

    # write file creation time
    time = str(dt.now())
    f.write("\t<!--CreationDate: " + time + "-->\n")

    # create group for whole file
    f.write("\t<g>\n")
    
    # save values to svg file
    # run through complines
    for i in range(num_complines):
        # initialize new lines counter
        newlines = 0

        # determine the number of points in the current line
        num_lines = x[i].shape[0]

        # initialize path string
        path =  "\t\t<path fill=\"none\" stroke=\"#FF0000\" "
        path += "stroke-width=\"1\" d=\"M"
        
        # run through lines
        for j in range(num_lines):
            path += "{:>7.4f},{:<7.4f} ".format(x[i][j],y[i][j])
            if len(path)//100 > newlines:
                path += "\n\t\t"
                newlines += 1
        
        # write close of path to file
        path += "\"/>\n"

        f.write(path)
    
    if type(text) == dict:
        # determine where the text will start
        xloc =  (text["start"][0] + ten - xmin) * in_to_pixels
        yloc = -(text["start"][1] - ten - ymax) * in_to_pixels
        txt = text["text"]
        # yloc = int(ht/2)
        # xloc = int(0.15 * wd)

        # determine the font-size
        fontsize = int(ht/10)

        # turn into strings
        xloc = str(xloc)
        yloc = str(yloc)
        fontsize = str(fontsize)

        # set font family
        fontfamily = "'Arial'"

        # create line and write to file
        text_line = "\t\t<text x=\"" + xloc + "\" y=\"" + yloc + "\" "
        text_line += "fill=\"#000000\" font-size=\"" + fontsize + "\" "
        text_line += "font-family=\"" + fontfamily + "\">" + txt + "</text>\n"
        f.write(text_line)

    # create group for whole file
    f.write("\t</g>\n")

    # write eps file closing statement
    f.write("</svg>")
    
    # close file
    f.close()

    return

def svg_polygon(filename,x,y,file_info,text,adobe):
    # create file
    f = open(filename + ".svg", "w+")

    # determine the number of lines to create
    num_complines = x.shape[0]

    # if the y and x spline counts are not equal in number, throw error
    if y.shape[0] != num_complines:
        raise ValueError('y array must have same number of lines as x array')

    # initialize max and min x and y values
    xmax = x[0][0]; xmin = x[0][0]; ymax = y[0][0]; ymin = y[0][0]

    # run through and determine the size of the file
    for i in range(num_complines):
        # determine the number of points in the current spline
        num_lines = x[i].shape[0]

        # if the y and x spline points are not equal in number, throw error
        if y[i].shape[0] != num_lines:
            raise ValueError('y array spline %d must have same' % i + \
                ' number of points as x array')
        
        # run through lines
        for j in range(num_lines):
            # determine local maxima and minima
            exmax = np.max(x[i][j]); exmin = np.min(x[i][j])
            eymax = np.max(y[i][j]); eymin = np.min(y[i][j])

            # compare and replace if necessary
            if exmax > xmax: xmax = exmax
            if exmin < xmin: xmin = exmin
            if eymax > ymax: ymax = eymax
            if eymin < ymin: ymin = eymin
    
    # move x and y such that they are a tenth of an inch within bounds
    ten = 0.1
    x += ten - xmin
    y -= ten + ymax

    # flip across y axis since coordinates in eps are opposite
    y *= -1.

    # translate to pixels
    if adobe:
        in_to_pixels = 72.0
    else:
        in_to_pixels = 96.0
    x *= in_to_pixels
    y *= in_to_pixels

    # save width and height of final shape
    wd = int((xmax - xmin + 2.*ten) * in_to_pixels)
    ht = int((ymax - ymin + 2.*ten) * in_to_pixels)

    # open svg file statement
    f.write("<svg width=\"{:d}px\" height=\"{:d}px\">\n".format(wd,ht))

    # create title, write to file
    title = filename.split("/")[-1]
    f.write("\t<title> {} </title>\n".format(title))
    
    # run through file info and input it into the file
    for info in file_info:
        f.write("\t<!--" + info + "-->\n")

    # write file creation time
    time = str(dt.now())
    f.write("\t<!--CreationDate: " + time + "-->\n")

    # create group for whole file
    f.write("\t<g>\n")
    
    # save values to svg file
    # run through complines
    for i in range(num_complines):
        # initialize new lines counter
        newlines = 0

        # determine the number of points in the current line
        num_lines = x[i].shape[0]

        # initialize path string
        path =  "\t\t<polygon fill=\"none\" stroke=\"#FF0000\" "
        path += "stroke-width=\"1\" points=\""
        
        # run through lines
        for j in range(num_lines):
            path += "{:>7.4f},{:<7.4f} ".format(x[i][j],y[i][j])
            if len(path)//100 > newlines:
                path += "\n\t\t"
                newlines += 1
        
        # write close of path to file
        path += "\"/>\n"

        f.write(path)
    
    if type(text) == dict:
        # determine where the text will start
        xloc =  (text["start"][0] + ten - xmin) * in_to_pixels
        yloc = -(text["start"][1] - ten - ymax) * in_to_pixels
        txt = text["text"]
        # yloc = int(ht/2)
        # xloc = int(0.15 * wd)

        # determine the font-size
        fontsize = int(ht/10)

        # turn into strings
        xloc = str(xloc)
        yloc = str(yloc)
        fontsize = str(fontsize)

        # set font family
        fontfamily = "'Arial'"

        # create line and write to file
        text_line = "\t\t<text x=\"" + xloc + "\" y=\"" + yloc + "\" "
        text_line += "fill=\"#000000\" font-size=\"" + fontsize + "\" "
        text_line += "font-family=\"" + fontfamily + "\">" + txt + "</text>\n"
        f.write(text_line)

    # finish group for whole file
    f.write("\t</g>\n")

    # write eps file closing statement
    f.write("</svg>")
    
    # close file
    f.close()

    return

def svg(filename,x,y,file_info = [],geometry = "polygon",text="",adobe = True):

    # create directory if not existant yet
    directory = filename.split("/")[1:]
    path = "/"
    for i in range(len(directory)-1):
        path += directory[i] + "/"
        if not os.path.isdir(path[:-1]):
            os.mkdir(path[:-1])

    if geometry == "polyline":
        svg_polyline(filename,x,y,file_info,text,adobe)
    elif geometry == "polygon":
        svg_polygon(filename,x,y,file_info,text,adobe)
    else:
        raise ValueError("Geometry incorrect")
    return



# # house example svg
# x = np.array([
#     np.array([0,1]),
#     np.array([1,1]),
#     np.array([1,0.5,0]),
#     np.array([0,0]),
# ])
# y = np.array([
#     np.array([0,0]),
#     np.array([0,1]),
#     np.array([1,1.5,1]),
#     np.array([1,0.5]),
# ])
# texter = "000"

# svg("/home/ben/Desktop/house_polygon",x,y,geometry="polyline",text=texter)