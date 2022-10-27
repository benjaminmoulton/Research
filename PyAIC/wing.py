import numpy as np
import json
import matplotlib.pyplot as plt

from component import Component, PseudoPrismoid, SymmetricAirfoil

class Wing:
    """A class calculating and containing the mass properties of an aircraft.

    Parameters
    ----------
    input_vars : string or dict , optional
        Must be a .json file or python dictionary.

    Raises
    ------
    TypeError
        If the input_vars type is not a dictionary or the file path to 
        a .json file
    """
    def __init__(self,input_dict={}):

        # retrieve info
        self._retrieve_info(input_dict)


    def _retrieve_info(self,input_dict):
        """A function which retrieves the information and stores it globally.
        """
        
        # store variables from file input dictionary
        wing_type = input_dict.get("type","symmetric_airfoil")

        side = input_dict.get("side","both")

        # store mass / density variables
        if "density" in input_dict:
            self.density = input_dict.get("density")
            self._given_density_not_mass = True
        elif "mass" in input_dict:
            self.mass = input_dict.get("mass")
            self._given_density_not_mass = False
        else: # raise error if neither given
            raise ValueError("No mass / density property given")
        
        # import connect_to info
        connect = input_dict.get("connect_to",{})
        connect_ID = connect.get("ID",0)
        location = connect.get("location", "tip")
        connect_dx = connect.get("dx", 0.0)
        connect_dy = connect.get("dy", 0.0)
        connect_dz = connect.get("dz", 0.0)
        connect_y_offset = connect.get("y_offset", 0.0)

        # import geometry info
        b = input_dict.get("semispan",1.0)
        dihedral = input_dict.get("dihedral",1.0)
        sweep = input_dict.get("sweep",0.0)
        chord = input_dict.get("chord",1.0)
        thickness = input_dict.get("thickness",0.10)
        camber = input_dict.get("camber",0.02)


    def get_mass_properties(self):

        print("yeah")