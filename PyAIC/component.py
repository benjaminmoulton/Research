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
        self.cg_location = np.zeros((3,))

        # retrieve info
        self._retrieve_info(input_dict)


    def _retrieve_info(self,input_dict):
        """A function which retrieves the information and stores it globally.
        """
        
        # store variables from file input dictionary

        # store input values
        self.type = input_dict.get("type","box")


    def get_density(self):
        return self.density


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
        self._Lambda = np.deg2rad(input_dict.get("sweep",0.0))
        self._Gamma = np.deg2rad(input_dict.get("dihedral",0.0))


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
        self._kb = tr * one + tr * two

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


    def get_mass_properties(self):

        # calculate kappa values
        self._get_kappa_values()

        # calculate mass
        self._V = self._b / 12. * self._ka
        if self._given_density_not_mass:
            self.mass = self.density * self._V
        
        # calculate center of gravity values
        num = 3. * self._kb + 4. * self._b * self._kc * np.tan(self._Lambda)
        xbar = - num / 20. / self._ka
        ybar = self._delta * self._b * self._kc / 5. / self._ka
        self.cg_location = np.array([xbar,ybar,0.0])

        # calculate moments and products of inertia about the wing root c/4
        num = 56. * self._b**2. * self._kf + self._kg
        Ixxo = self.mass * num / 280. / self._ka
        
        one = 2. * self._b * self._kf * np.tan(self._Lambda)**2.
        two = self._kd * np.tan(self._Lambda)
        num = 84. * self._b * (one + two) + 49. * self._ke + 3. * self._kg
        Iyyo = self.mass * num / 840. / self._ka

        one = self._b * self._kf * ( np.tan(self._Lambda)**2. + 1. )
        two = self._kd * np.tan(self._Lambda)
        num = 12. * self._b * (one + two) + 7. * self._ke
        Izzo = self.mass * num / 120. / self._ka

        num = 4. * self._b * self._kf * np.tan(self._Lambda) + self._kd
        Ixyo = self._delta * self._b * self.mass * num / 20. / self._ka

        Ixzo = 0.0
        Iyzo = 0.0

        # create inertia tensor
        Io = np.array([
            [ Ixxo,-Ixyo,-Ixzo],
            [-Ixyo, Iyyo,-Ixzo],
            [-Ixzo,-Ixzo, Izzo]
        ])

        # calculate mass shift from parallel axis theorem
        s = self.cg_location[:,np.newaxis]
        inn_prod = np.matmul(s.T,s)[0,0]
        out_prod = np.matmul(s,s.T)
        I_shift = self.mass*( inn_prod * np.eye(3) - out_prod )

        # calculate inertia tensor about the cg
        self.inertia_tensor = Io - I_shift

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