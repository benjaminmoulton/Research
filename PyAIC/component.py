import numpy as np
import json
import matplotlib.pyplot as plt

class Component:
    """A default class for calculating and containing the mass properties of an object.

    Parameters
    ----------
    input_vars : dict , optional
        Must be a python dictionary
    """
    def __init__(self,input_dict={}):

        # initialize cg location, inertia tensor, mass, density
        self.density = 0.0
        self.mass = 0.0
        self.volume = 0.0
        self.inertia_tensor = np.zeros((3,3))
        self.cg_location = np.zeros((3,1))
        self.locations = {
            "root" : np.zeros((3,1)),
            "tip"  : np.zeros((3,1))
        }
        self._components = []

        # retrieve info
        self._retrieve_info(input_dict)


    def _retrieve_info(self,input_dict):
        """A function which retrieves the information and stores it globally.
        """
        
        # store variables from file input dictionary

        # store input values
        self.type = input_dict.get("type","box")


    def update_densities(self):
        return 0


    def get_density(self):
        return self.density

  
    def shift_properties_to_location(self,input_location):
        """Method which determines the mass properties of the given component
        about a given location.
        
        Parameters
        ----------
        input_location : array
            The location about which to determine the mass properties. Must be
            formatted as input_location = [[x],[y],[z]].
        """

        # if input_location is a list, turn into numpy array
        if isinstance(input_location, list):
            input_location = np.array(input_location)
        
        # determine new location
        s = input_location - self.cg_location
        # print(self.cg_location)

        # calculate mass shift from parallel axis theorem
        inner_product = np.matmul(s.T,s)[0,0]
        outer_product = np.matmul(s,s.T)
        I_shift = self.mass*( inner_product * np.eye(3) - outer_product )

        # calculate inertia tensor about new location
        Inew = self.inertia_tensor + I_shift

        # return dictionary of values
        output_dict = {
            "mass" : self.mass,
            "cg_location" : self.cg_location,
            "inertia_tensor" : Inew
        }
        return output_dict


    def get_mass_properties(self):
        """Method which returns mass, cg, I about cg rotated to total cframe"""

        output_dict = {
            "mass" : self.mass,
            "cg_location" : self.cg_location,
            "inertia_tensor" : self.inertia_tensor
        }
        return output_dict


class PseudoPrismoid(Component):
    """A default class for calculating and containing the mass properties of a
    Psuedo Prismoid.

    Parameters
    ----------
    input_vars : dict , optional
        Must be a python dictionary
    """
    def __init__(self,input_dict={}):

        # invoke init of parent
        Component.__init__(self,input_dict)

        # retrieve additional info
        self._retrieve_info(input_dict)

        # initialize the volume
        cr,ct = self._cr,self._ct
        tr,tt = self._tr,self._tt
        cr2, ct2 = cr**2., ct**2.
        crct = cr * ct
        ka = tr * (3.*cr2 + 2.*crct + ct2) + tt *(cr2 + 2.*crct + 3.*ct2)
        u0 = 1.0
        self.volume = self._b / 12. * ka * u0


    def _retrieve_info(self,input_dict):
        """A function which retrieves the information and stores it globally.
        """
        
        # store variables from file input dictionary

        # store input values
        self._side = input_dict.get("side","right")
        if self._side == "right": self._delta = 1.0
        else: self._delta = -1.0
        # store mass / density variables
        if "density" in input_dict:
            self.density = input_dict.get("density")
            self._given_density_not_mass = True
        elif "mass" in input_dict:
            self.mass = input_dict.get("mass")
            self._given_density_not_mass = False
        else: # raise error if neither given
            raise ValueError("No mass / density property given")
        
        # store other variables
        geometry = input_dict.get("geometry",{})
        self._b = geometry.get("span",1.0)
        self._cr,self._ct = geometry.get("chord",[1.0,1.0])
        self._tr,self._tt = geometry.get("thickness",[0.1,0.1])
        self._ur,self._ut = geometry.get("camber",[0.01,0.01])
        self._xmt = geometry.get("max_thickness_location",0.1)
        self._Lambda = np.deg2rad(input_dict.get("sweep",0.0))
        self._Gamma = np.deg2rad(input_dict.get("dihedral",0.0))
        root_location = input_dict.get("root_location",[0.0,0.0,0.0])
        tip_location  = input_dict.get( "tip_location",[0.0,0.0,0.0])
        self.locations = {
            "root" : np.array(root_location)[:,np.newaxis],
            "tip"  : np.array( tip_location)[:,np.newaxis]
        }


    def _get_kappa_values(self):

        # initialize certian values for ease of calculation
        cr,ct = self._cr,self._ct
        tr,tt = self._tr,self._tt
        # times values
        cr2, ct2 = cr**2., ct**2.
        crct = cr * ct
        cr3, ct3 = cr**3., ct**3.
        cr2ct, crct2 = cr2 * ct, cr * ct2
        cr4, ct4 = cr**4., ct**4.
        cr3ct, cr2ct2, crct3 = cr3 * ct, cr2 * ct2, cr * ct3

        # calculate kappa values
        self._ka = tr * (3.*cr2 + 2.*crct + ct2) + tt *(cr2 + 2.*crct + 3.*ct2)
        
        one = 4.*cr3 + 3.*cr2ct + 2.*crct2 +    ct3
        two =    cr3 + 2.*cr2ct + 3.*crct2 + 4.*ct3
        self._kb = tr * one + tt * two

        one = 3.*cr2 + 4.*crct + 3.*ct2
        two =    cr2 + 3.*crct + 6.*ct2
        self._kc = tr * one + 2. * tt * two

        one = (cr + ct) * (2.*cr2 + crct + 2.*ct2)
        two = cr3 + 3.*cr2ct + 6.*crct2 + 10.*ct3
        self._kd = tr * one + tt * two

        one = 5.*cr4 + 4.*cr3ct + 3.*cr2ct2 + 2.*crct3 +    ct4
        two =    cr4 + 2.*cr3ct + 3.*cr2ct2 + 4.*crct3 + 5.*ct4
        self._ke = tr * one + tt * two

        one = cr2 + 2.*crct +  2.*ct2
        two = cr2 + 4.*crct + 10.*ct2
        self._kf = tr * one + tt * two

        one = 35.*cr4 + 20.*cr3ct + 10.*cr2ct2 +  4.*crct3 +     ct4
        two = 15.*cr4 + 20.*cr3ct + 18.*cr2ct2 + 12.*crct3 +  5.*ct4
        thr =  5.*cr4 + 12.*cr3ct + 18.*cr2ct2 + 20.*crct3 + 15.*ct4
        fou =     cr4 +  4.*cr3ct + 10.*cr2ct2 + 20.*crct3 + 35.*ct4
        self._kg = tr**3. *one + tr**2.*tt *two + tr*tt**2. *thr + tt**3. *fou


    def _get_upsilon_values(self):

        # calculate upsilon values
        self._u0 = 1.0
        self._u1 = 1.0
        self._u2 = 1.0
        self._u3 = 1.0


    def get_mass_properties(self):
        """Method which returns mass, cg, I about cg rotated to total cframe"""

        # calculate kappa values
        self._get_kappa_values()

        # calculate upsilon values
        self._get_upsilon_values()

        # calculate mass
        # self.volume = self._b / 12. * self._ka
        if self._given_density_not_mass:
            self.mass = self.density * self.volume
        
        # calculate center of gravity values
        num1 = 3. * self._kb * self._u1
        num = num1 + 4. * self._b * self._kc * self._u0 * np.tan(self._Lambda)
        xbar = - num / 20. / self._ka / self._u0
        ybar = self._delta * self._b * self._kc * self._u0/5./self._ka/self._u0
        self.cg_location = np.array([xbar,ybar,0.0])[:,np.newaxis]

        # calculate moments and products of inertia about the wing root c/4
        num = 56. * self._b**2. * self._kf * self._u0 + self._kg * self._u3
        Ixxo = self.mass * num / 280. / self._ka / self._u0
        
        one = 2. * self._b * self._kf * self._u0 * np.tan(self._Lambda)**2.
        two = self._kd * self._u1 * np.tan(self._Lambda)
        num1 = 84. * self._b * (one + two) + 49. * self._ke * self._u2
        num = num1 + 3. * self._kg * self._u3
        Iyyo = self.mass * num / 840. / self._ka / self._u0

        one = self._b * ( np.tan(self._Lambda)**2. + 1. ) * self._kf * self._u0
        two = self._kd * self._u1 * np.tan(self._Lambda)
        num = 12. * self._b * (2. * one + two) + 7. * self._ke * self._u2
        Izzo = self.mass * num / 120. / self._ka / self._u0

        num1 = 4. * self._b * self._kf * self._u0 * np.tan(self._Lambda)
        num = num1 + self._kd * self._u1
        Ixyo = -self._delta * self._b * self.mass * num / 20./self._ka/self._u0

        Ixzo = 0.0
        Iyzo = 0.0

        # create inertia tensor
        Io = np.array([
            [ Ixxo,-Ixyo,-Ixzo],
            [-Ixyo, Iyyo,-Ixzo],
            [-Ixzo,-Ixzo, Izzo]
        ])

        # calculate mass shift from parallel axis theorem
        s = self.cg_location
        inner_product = np.matmul(s.T,s)[0,0]
        outer_product = np.matmul(s,s.T)
        I_shift = self.mass*( inner_product * np.eye(3) - outer_product )

        # calculate inertia tensor about the cg
        Icg = Io - I_shift
        # rows = [0,1,2,0,0,1]
        # cols = [0,1,2,1,2,2]
        # I_cg = Icg[rows,cols]
        # lbm_slug = 32.17405082
        # m_wng = 2.66974726 / lbm_slug
        # cg_wng = np.array([ 
        #     -0.13712589,
        #     0.72910392,
        #     0.00000000000001
        # ])
        # I_wng = np.array([
        #     2.16737992,
        #     0.15240118,
        #     2.31726753,
        #     0.23551551,
        #     -0.00000002,
        #     -0.00000100
        # ]) / lbm_slug
        # I_wng_abt_cg = np.array([
        #     0.74816221,
        #     0.10220056,
        #     0.84784920,
        #     -0.03140323,
        #     -0.00000015,
        #     -0.00000032
        # ]) / lbm_slug
        # for i in range(6):
        #     print(I_wng[i])
        #     print(I_cg[i])
        #     print(I_wng_abt_cg[i])
        #     print()
        # # print( Icg[0,0])
        # # print( Icg[1,1])
        # # print( Icg[2,2])
        # # print(-Icg[0,1])
        # # print(-Icg[0,2])
        # # print(-Icg[1,2])
        # # print()

        # define a rotation matrix
        CG = np.cos(-self._delta*self._Gamma)
        SG = np.sin(-self._delta*self._Gamma)
        Rx = np.array([
            [1.0, 0.0, 0.0],
            [0.0,  CG, -SG],
            [0.0,  SG,  CG]
        ])

        # rotate cg location and inertia for dihedral
        # print(self.cg_location)
        self.cg_location = np.matmul(Rx,self.cg_location)
        # print(self.cg_location)
        self.inertia_tensor = np.matmul(Rx,np.matmul(Icg,Rx.T))

        # shift cg by root location given
        self.cg_location += self.locations["root"]

        output_dict = {
            "mass" : self.mass,
            "cg_location" : self.cg_location,
            "inertia_tensor" : self.inertia_tensor
        }
        return output_dict
        


class SymmetricAirfoil(PseudoPrismoid):
    """A default class for calculating and containing the mass properties of a
    Symmetric Airfoil.

    Parameters
    ----------
    input_vars : dict , optional
        Must be a python dictionary
    """
    def __init__(self,input_dict={}):

        # invoke init of parent
        PseudoPrismoid.__init__(self,input_dict)

        # save a coefficients for NACA 4-digit thickness distribution
        self._a0 = 2.969
        self._a1 = -1.260
        self._a2 = -3.516
        self._a3 = 2.843
        self._a4 = -1.036

        # initialize the volume
        cr,ct = self._cr,self._ct
        tr,tt = self._tr,self._tt
        cr2, ct2 = cr**2., ct**2.
        crct = cr * ct
        ka = tr * (3.*cr2 + 2.*crct + ct2) + tt *(cr2 + 2.*crct + 3.*ct2)
        a0 = self._a0
        a1 = self._a1
        a2 = self._a2
        a3 = self._a3
        a4 = self._a4
        u0 = 1./ 60.*( 40.*a0+ 30.*a1+ 20.*a2+ 15.*a3+ 12.*a4)
        self.volume = self._b / 12. * ka * u0


    def _get_upsilon_values(self):

        # initialize certain values for ease of calculation
        a0 = self._a0
        a1 = self._a1
        a2 = self._a2
        a3 = self._a3
        a4 = self._a4
        a03, a13, a23, a33, a43 = a0**3., a1**3., a2**3., a3**3., a4**3.
        a02, a12, a22, a32, a42 = a0**2., a1**2., a2**2., a3**2., a4**2.

        # calculate upsilon values
        self._u0 = 1./ 60.*( 40.*a0+ 30.*a1+ 20.*a2+ 15.*a3+ 12.*a4)
        self._u1 = 1./ 60.*( 56.*a0+ 50.*a1+ 40.*a2+ 33.*a3+ 28.*a4)
        self._u2 = 1./980.*(856.*a0+770.*a1+644.*a2+553.*a3+484.*a4)
        one1 = 2./5.*a03 + 1*a02*a1 + 3./4.*a02*a2 + 3./5.*a02*a3
        one2 = 1./2.*a02*a4 + 6./7.*a0*a12 + 4./3.*a0*a1*a2
        one3 = 12./11.*a0*a1*a3 + 12./13.*a0*a1*a4
        two1 = 6./11*a0*a22 + 12./13.*a0*a2*a3 + 4./5.*a0*a2*a4 +2./5.*a0*a32
        two2 = 12./17.*a0*a3*a4 + 6./19.*a0*a42 + 1./4.*a13 + 3./5.*a12*a2
        thr1 = 1./2.*a12*a3 + 3./7.*a12*a4 + 1./2.*a1*a22 + 6./7.*a1*a2*a3
        thr2 = 3./4.*a1*a2*a4 + 3./8.*a1*a32 + 2./3.*a1*a3*a4 + 3./10.*a1*a42
        fou1 = 1./7.*a23 + 3./8.*a22*a3 + 1./3.*a22*a4 + 1./3.*a2*a32
        fou2 = 3./5.*a2*a3*a4 + 3./11.*a2*a42 + 1./10.*a33 + 3./11.*a32*a4
        fou3 = 1./4.*a3*a42 + 1./13.*a43
        one = one1 + one2 + one3
        self._u3 = one + two1 + two2 + thr1 + thr2 + fou1 + fou2 + fou3



class DiamondAirfoil(PseudoPrismoid):
    """A default class for calculating and containing the mass properties of a
    Diamond Airfoil.

    Parameters
    ----------
    input_vars : dict , optional
        Must be a python dictionary
    """
    def __init__(self,input_dict={}):

        # invoke init of parent
        PseudoPrismoid.__init__(self,input_dict)

        # initialize the volume
        cr,ct = self._cr,self._ct
        tr,tt = self._tr,self._tt
        cr2, ct2 = cr**2., ct**2.
        crct = cr * ct
        ka = tr * (3.*cr2 + 2.*crct + ct2) + tt *(cr2 + 2.*crct + 3.*ct2)
        u0 = 1. / 2.
        self.volume = self._b / 12. * ka * u0


    def _get_upsilon_values(self):

        # calculate upsilon values
        self._u0 = 1. /  2.
        self._u1 = ( 4. * self._xmt     + 1. ) / 6.
        self._u2 = ( 8. * self._xmt**2. + 3. ) / 14.
        self._u3 = 1. / 4.

