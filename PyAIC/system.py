import numpy as np
import json
import matplotlib.pyplot as plt

from wing import Wing
from component import Component, PseudoPrismoid, SymmetricAirfoil

class AircraftSystem:
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
    def __init__(self,input_vars={},verbose=True):

        # get info or raise error
        self._get_input_vars(input_vars)

        # retrieve info
        self._retrieve_info()


    def _get_input_vars(self,input_vars):
        # get info or raise error

        # determine if the input_vars is a file or a dictionary
        input_vars_type = type(input_vars)

        # dictionary
        if input_vars_type == dict:
            self.input_dict = input_vars
        
        # json file
        elif input_vars_type == str and input_vars.split(".")[-1] == "json":
            self.input_dict = self._get_json(input_vars)

        # raise error
        else:
            raise TypeError("input_vars must be json file path, or " + \
                "dictionary, not {0}".format(input_vars_type))
    

    def _get_json(self,file_path):
        # import json file from file path
        json_string = open(file_path).read()

        # save to vals dictionary
        input_dict = json.loads(json_string)
        
        return input_dict


    def _retrieve_info(self):
        """A function which retrieves the information and stores it globally.
        """
        
        # store variables from file input dictionary

        # store component input values
        components = self.input_dict.get("components",{})

        wing_types = ["pseudo_prismoid","symmetric_airfoil","diamond_airfoil"]

        # initialize components, and each
        self.components = {}
        for component in components:
            input_dict = components[component]
            id_number = input_dict.get("ID")

            # initialize PseudoPrismoid
            if input_dict["type"] in wing_types:
                self.components[id_number] = Wing(input_dict)



    def get_mass_properties(self):

        # determine properties of each component
        for component in self.components:
            self.components[component].get_mass_properties()


if __name__ == "__main__":
    AS = AircraftSystem("test_input.json")
    AS.get_mass_properties()