#!/usr/bin/env python3   ###### rename file. constructor?
# -*- coding: utf-8 -*-
"""
Created on Fri Jun  7 09:37:00 2019

@author: benjamin
"""

import numpy as np

class point:
    
#    x = 0; y = 0; z = 0;

    def __init__(self, x = 0, y = 0, z = 0):
        self.x = x
        self.y = y
        self.z = z
    
    
    def fix(self, x, y, z = 0):
        import copy
        self.x = copy.deepcopy(x)
        self.y = copy.deepcopy(y)
        self.z = copy.deepcopy(z)
#        
#        self.x = x
#        self.y = y
#        self.z = z
    
    
    def set_point(self, point):
        import copy
        self.x = copy.deepcopy(point.x)
        self.y = copy.deepcopy(point.y)
        self.z = copy.deepcopy(point.z)
#        import copy
#        
#        self = copy.deepcopy(point)
##        self = selfcopy
##        self.coord()
    
    def move(self, x_shift = 0, y_shift = 0, z_shift = 0):
        self.x += x_shift
        self.y += y_shift
        self.z += z_shift
    
    def coord(self):
        x_string = '%.3f' % self.x
        y_string = '%.3f' % self.y
        z_string = '%.3f' % self.z
        if(self.x < 0 and self.x > -1):
            x_string = x_string[:1] + x_string[2:]
        if(self.y < 0 and self.y > -1):
            y_string = y_string[:1] + y_string[2:]
        if(self.z < 0 and self.z > -1):
            z_string = z_string[:1] + z_string[2:]
        coordinate = '( %s , %s , %s )' % (x_string,y_string,z_string)
        print(coordinate)
    
    def provemove(self,x_shift,y_shift,z_shift = 0):
        self.coord()
        self.move(x_shift,y_shift,z_shift)
        self.coord()
        print();print()

class points: # rename
    
#    points = [];
#    
#    x = [];
#    
#    x_min_index = 0
#    
#    y = [];
#    
#    z = [];
    
    
    def __init__(self, x = [], y = [], z = [],ccw = True):
        self.points = np.array([])
        """
        
        figure out how to deal with given a only a set of x's or so on
        
        """
        #####################################################stopped numpy transition here#######
        if(len(z) == len(x)):
            
            for i in range(len(x)):
                
                self.points.append(point(x[i],y[i],z[i]))
            
        elif(len(z) == 0 and not len(z) == len(x) and not len(y) == 0):
            
            for i in range(len(x)):
                
                self.points.append(point(x[i],y[i],0))
            
        elif(len(z) == 0 and not len(z) == len(x) and len(y) == 0):
            
            for i in range(len(x)):
                
                self.points.append(point(x[i],0,0))
            
        else:
            
            for i in range(len(x)):
                
                self.points.append(point(x[i],y[i],0))
        
        if(len(self.points) != 0):
            
            self.x_min_index = self.x().index(min(self.x()))
            self.x_max_index = self.x().index(max(self.x()))
        
        self.ccw = ccw
        
        return
    
    def x(self):
        
        x_list = []
        
        for i in range(len(self.points)):
            x_list.append(self.points[i].x)
        
        return x_list
    
    def y(self):
        
        y_list = []
        
        for i in range(len(self.points)):
            y_list.append(self.points[i].y)
        
        return y_list
    
    def z(self):
        
        z_list = []
        
        for i in range(len(self.points)):
            z_list.append(self.points[i].z)
        
        return z_list
    
    def x_max(self):
        
        return self.points[self.x_max_index].x
    
    def x_min(self):
        
        return self.points[self.x_min_index].x
    
    def chord(self):
        
        max_x = self.points[self.x_max_index].x
        min_x = self.points[self.x_min_index].x
        chord = round(max_x - min_x,2)
        
        return chord
    
    def unityChord(self):
        
        c = self.chord()
        
        for i in range(len(self.points)):
            self.points[i].x /= c
            self.points[i].y /= c
        
        self.old_c = c
    
    def change_chord(self,newChord,oldChord = 0):
        if(oldChord == 0):
            oldChord = self.chord()
        
        for u in range(len(self.points)):
            
            self.points[u].x *= ( newChord / oldChord )
            
            self.points[u].y *= ( newChord / oldChord )
    

    def thickness(self):
        
        y = self.y()
        max_y = max(y)
        max_y_index = y.index(max_y)
        min_y_at_max_y_index = self.interp(self.points[max_y_index].x,False)
        thickness = round(max_y - min_y_at_max_y_index,2)
        
        return thickness
    
    def percentThickness(self):
        
        t = self.thickness() / self.chord()
        
        return t
    
    #a function that is still to be developed, but is in use.
    def quarter_chord(self): # finds the cg of a given NACA 4 digit airfoil
        qc = point()
        chord = self.chord()
        qc.x = chord / 4
        qc.x += self.points[self.x_min_index].x
        #some code that finds the exact position of y of cg on wing using NACAfour number
        avg_y = 0;
        
        for i in range(len(self.points)):
            avg_y += self.points[i].y
        
        qc.y = avg_y / len(self.points)
        return qc
    
    def shift(self, x_shift = 0, y_shift = 0, z_shift = 0):
        
        for i in range(len(self.points)):
            self.points[i].x += x_shift
            self.points[i].y += y_shift
            self.points[i].z += z_shift
            
    def add(self,data,insert_at_beginning = False,location = 9000):
        
        if(insert_at_beginning):
            location = 0
        
        if(location == 9000):
            
            if(type(data) == point):
                
                self.points.append(data)
                
                self.x_min_index = self.x().index( min( self.x() ) )
                self.x_max_index = self.x().index( max( self.x() ) )
                
                return
            
            elif(type(data) == points):
                
                self.points.extend(data.points)
                
                self.x_min_index = self.x().index(min(self.x()))
                self.x_max_index = self.x().index(max(self.x()))
                
                return
            
        else:
            
            if(type(data) == point):
                
                self.points.insert(location,data)
                
                self.x_min_index = self.x().index( min( self.x() ) )
                self.x_max_index = self.x().index( max( self.x() ) )
                
                return
            
            elif(type(data) == points):
                
                for i in range(len(data.points)-1,-1,-1):
                    
                    self.points.insert(location,data.points[i])
                
                self.x_min_index = self.x().index(min(self.x()))
                self.x_max_index = self.x().index(max(self.x()))
                
                return
    
    def xmax_in_range(self,start_percentage = 0.1,end_percentage = 0.8, isTop = True):
        c = self.chord()
        start = start_percentage * c + self.points[self.x_min_index].x
        end = end_percentage * c + self.points[self.x_min_index].x
        
        max_ = -1000
        ind = -100
        add = 0
        
        if(not isTop):
            add += self.x_min_index
        
        for i in range(len(self.points) // 2 ):
            if( self.points[i + add].x >= max_ and self.points[i + add].x >= start and self.points[i + add].x <= end):
                max_ = self.points[i + add].x
                ind = i + add - 1
        
        return max_, ind
 
    def ymax_in_range(self,start_percentage = 0.1,end_percentage = 0.8, isTop = True):
        c = self.chord()
        start = start_percentage * c + self.points[self.x_min_index].x
        end = end_percentage * c + self.points[self.x_min_index].x
        
        max_ = -1000
        ind = -100
        add = 0
        
        if(not isTop):
            add += self.x_min_index - 1
        
        for i in range((len(self.points) // 2) - 1 ):
            ydist = abs(self.points[i+add].y - self.points[i+add+1].y)
            if( ydist >= max_ and self.points[i + add].x >= start and self.points[i + add].x <= end):
                max_ = ydist
                ind = i + add; 
                if(not isTop): 
                    ind += 1
        
        return max_, ind  
    
    # if a list is not whole, ( a closed loop ) make it so
    def makeWhole(self, start_at_end = True): # verses ending at the start ie sectionfoil, start on bottom or top
        
#        for i in range(len(self.points)):
#            if(type ( self.points[i].x ) == np.ndarray): 
#                print( type ( self.points[i].x )); print(self.points[i].x)
#            if(type ( self.points[i].y ) == np.ndarray): 
#                print( type ( self.points[i].y )); print(self.points[i].y)
#            if(type ( self.points[i].z ) == np.ndarray): 
#                print( type ( self.points[i].z )); print(self.points[i].z)
        if(type(self.points[0].x ) == np.ndarray or type(self.points[len(self.points)-1].x ) == np.ndarray or
           type(self.points[0].y ) == np.ndarray or type(self.points[len(self.points)-1].y ) == np.ndarray or
           type(self.points[0].z ) == np.ndarray or type(self.points[len(self.points)-1].z ) == np.ndarray):
            if (np.round(self.points[len(self.points)-1].x,5) != np.round(self.points[0].x,5) or 
                np.round(self.points[len(self.points)-1].y,5) != np.round(self.points[0].y,5) or 
                np.round(self.points[len(self.points)-1].z,5) != np.round(self.points[0].z,5)):
                
                pt = point()
                
                if(start_at_end):
                    pt.set_point(self.points[len(self.points)-1])
                    self.add(pt,True)
                    
                else:# end at start
                    pt.set_point(self.points[0])
                    self.add(pt)    
            
        else:
            if (round(self.points[len(self.points)-1].x,5) != round(self.points[0].x,5) or 
                round(self.points[len(self.points)-1].y,5) != round(self.points[0].y,5) or 
                round(self.points[len(self.points)-1].z,5) != round(self.points[0].z,5)):
                
                pt = point()
                
                if(start_at_end):
                    pt.set_point(self.points[len(self.points)-1])
                    self.add(pt,True)
                    
                else:# end at start
                    pt.set_point(self.points[0])
                    self.add(pt)
        
    
    def remove(self,index_to_remove = 0):
        
        self.points.pop(index_to_remove)
        
        self.x_min_index = self.x().index(min(self.x()))
        self.x_max_index = self.x().index(max(self.x()))
        
        return
        
    
    def split(self,start_index = 0,end_index = 0):
        
        if(end_index == 0):
            end_index = len(self.points)
        
        import copy
        
        # create copy of self
        
        new_points = copy.deepcopy(self)
        
        # split the lists within the object based on given indices.
        
        new_points.points = new_points.points[start_index:end_index]
        
        if(len(new_points.points) > 0):
            new_points.x_min_index = new_points.x().index(min(new_points.x()))
            new_points.x_max_index = new_points.x().index(max(new_points.x()))
        
        # return the new points object
        
        return new_points
    
    
    def interp(self,x_value,isTop = True):
        
        # bring in interpolation solver
        from scipy.interpolate import interp1d as itp
        
        # inialize top and bottom point lists
        # and return the interpolated value
        
        if(isTop):
            top = self.split(end_index = self.x_min_index)
            
            top_interp = itp(top.x(),top.y())#,kind = 'cubic')
            
            return top_interp(x_value)
        
        else:
            bottom = self.split(start_index = self.x_min_index)
            
            bottom_interp = itp(bottom.x(),bottom.y())#,kind = 'cubic')
            
            return bottom_interp(x_value)
        
        
    
    def size(self):
        
        return len(self.points)
    
    def show(self):
        
        print()
        
        for i in range(len(self.points)):
            
            self.points[i].coord()
        
        print()
        
        return
    
    def plot(self,other_plotters = [],axiss = 'equal',start_split = 0, end_split = 0,display = [''],Subplot = 0):
        
        import matplotlib.pyplot as plt
        
        if(len(other_plotters) > 0):
            disp = len(display)
            for i in range(len(other_plotters) - disp + 1):
                
                ind = disp + i - 1
                
                if(type(other_plotters[ind]) == point):
                    display.append( 'o' )
                    
                else:
                    display.append( '' )
        
        if(Subplot == 1):
            plt.subplot(211)
            plt.title('root foil')
        elif(Subplot == 2):
            plt.subplot(212)
            plt.title('tip foil')
        
        
        plt.axis(axiss)
        splitted = self.split(start_split,end_split)
        plt.plot(splitted.x(),splitted.y(),display[0])
        
        for i in range(len(other_plotters)):
            if(type(other_plotters[i]) == points):
                plt.plot(other_plotters[i].x(),other_plotters[i].y(),display[i+1])
            elif(type(other_plotters[i]) == point):
                plt.plot(other_plotters[i].x,other_plotters[i].y,display[i+1])
                
        
        plt.show()
        
        return
        
    def print_lengths(self):
        
#        print('  x length %d' % len(self.x))
#        print('  y length %d' % len(self.y))
#        print('  z length %d' % len(self.z))
        print('pts length %d' % len(self.points))
        
        return

class airfoil(points):
    
    def camber_line(self,camberline):
        self.camberline = camberline;
    
    def get_chord(self):
        
        maxx = max(self.x())
        
        minx = min(self.x())
        
        c = maxx - minx
        
        return c
    

class wing:
    
#    xyz = points()
#    
#    uvw = points()
    
    
    def __init__(self, xyz_points = airfoil(), uvw_points = airfoil() ):
        self.xyz = points(xyz_points.x,xyz_points.y,xyz_points.z)
        
        self.uvw = points(uvw_points.x,uvw_points.y,uvw_points.z)
        
        return
    
    
    def swap_sides(self):
        
        temp = self.xyz
        
        self.xyz = self.uvw
        
        self.uvw = temp
        
        return
    

#pt1 = point(5,1)
#pt2 = point(6,2)
#pt3 = point(0.0001,0.110100)
#lis = points()
#
#lis.add(pt1)
#lis.add(pt2)
#lis.add(pt3)
#lis.add(point(3,5))
#lis.add(point(5,7))
#lis.add(point(6,7))
#
#
#
#lis.remove(1)
##print()
#
##print(len(lis.x))
##print(lis.x)
#half = lis.split(0,3)
#
##print(half.x)
##print(lis.x)
#
#new_x = lis.x
#
#new_y = lis.y
#
#new_points = points(new_x,new_y)
#
#
#new_points.show()
#print()
#print()
#print()
#new_points.points[0].x = 100
#new_points.show()




#w = wing()
#
#w.xyz = lis
#
#w.uvw = lis
#
#w.xyz.show()
#w.uvw.show()




#
#p = point(5,1)
#p.fix(4,4)
#p.x = 5
#print(p.x)
#
#print(p.y)
#
#p.coord()
#
#
#
#print(type(p))





