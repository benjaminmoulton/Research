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

        # initialize get each id number
        ids = []; attach_ids = []; name = []
        for component in components:
            input_dict = components[component]
            ids.append(input_dict.get("ID"))
            attach_ids.append(input_dict["connect_to"].get("ID",0))
            name.append(component)
        
        # throw error if an id given is repeated
        unique_ids,counts = np.unique(ids,return_counts=True)
        repeated_ids = unique_ids[counts > 1]
        if len(unique_ids) != len(ids):
            raise TypeError("Repeated ID values: {}".format(repeated_ids))
        
        # throw error if attach order is circular
               
        # reorganize attach_ids in order of run

        quit()

        self.components = {}
        for component in components:

            # initialize wings
            input_dict = components[component]
            id_number = input_dict.get("ID")
            if input_dict["type"] in wing_types:
                self.components[id_number] = Wing(input_dict)


    def report_as_SolidWorks_report(self,info,positive_tensor=True):
        """Method which reports the mass and inertia properties as given in 
        SolidWorks.
        
        Parameters
        ----------
        info : dictionary
            The dictionary with 'mass', 'cg_location' and 'inertia_tensor'
            information.
        
        positive_tensor : boolean, optional
            Whether to report in positive tensor formulation. Defaults to true.
        """


        print("Mass = {:> 10.8f} slugs".format(info["mass"]))
        print()

        print("Center of mass: (feet)")
        print("\tX = {:> 10.8f}".format(info["cg_location"][0,0]))
        print("\tY = {:> 10.8f}".format(info["cg_location"][1,0]))
        print("\tZ = {:> 10.8f}".format(info["cg_location"][2,0]))
        print()
        I = info["inertia_tensor"] * 1.0
        if positive_tensor:
            I[[0,0,1,1,2,2],[1,2,0,2,0,1]] *= -1.0
        [[Ixx,Ixy,Ixz],[Iyx,Iyy,Iyz],[Izx,Izy,Izz]] = I
        print("Moment of inertia:( slugs * square feet)")
        if positive_tensor:
            print("\t\tPositive Tensor Formulation")
        else:
            print("\t\tNegative Tensor Formulation")
        print("\tIxx = {:> 10.8f}\tIxy = {:> 10.8f}\tIxz = {:> 10.8f}".format(\
            Ixx,Ixy,Ixz))
        print("\tIyx = {:> 10.8f}\tIyy = {:> 10.8f}\tIyz = {:> 10.8f}".format(\
            Iyx,Iyy,Iyz))
        print("\tIzx = {:> 10.8f}\tIzy = {:> 10.8f}\tIzz = {:> 10.8f}".format(\
            Izx,Izy,Izz))
        print()


    def get_mass_properties(self,report=False,positive_tensor=True):

        # determine properties of each component
        for i in self.components:
            self.components[i].get_mass_properties()
        
        # determine total mass
        self.mass = 0.0
        for i in self.components:
            self.mass += self.components[i].mass
        
        # determine total cg location
        self.cg_location = np.zeros((3,1))
        for i in self.components:
            self.cg_location += self.components[i].mass * \
                self.components[i].cg_location
        self.cg_location /= self.mass
        
        # determine total inertia
        self.inertia_tensor = np.zeros((3,3))
        location = self.cg_location
        for i in self.components:
            self.inertia_tensor += \
                self.components[i].shift_properties_to_location(\
                location)["inertia_tensor"]

        # return dictionary of values
        output_dict = {
            "mass" : self.mass,
            "cg_location" : self.cg_location,
            "inertia_tensor" : self.inertia_tensor
        }

        # report
        if report:
            self.report_as_SolidWorks_report(output_dict,positive_tensor)
        
        return output_dict


if __name__ == "__main__":
    # AS = AircraftSystem("test_input.json")
    # AS.get_mass_properties()
    AS = AircraftSystem("horizon.json")
    AS.get_mass_properties(report=True)