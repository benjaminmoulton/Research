import numpy as np
import json
import matplotlib.pyplot as plt

from wing import Wing
from component import Component, Cuboid, Cylinder, Sphere, Rotor

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
        self._initialize_objects()


    def _get_input_vars(self,input_vars):
        # get info or raise error

        # determine if the input_vars is a file or a dictionary
        input_vars_type = type(input_vars)

        # dictionary
        if input_vars_type == dict:
            self.input_dict = input_vars
            self.file_name = "Aircraft"
        
        # json file
        elif input_vars_type == str and input_vars.split(".")[-1] == "json":
            self.file_name = input_vars
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


    def _initialize_objects(self):

        # store component input values
        components = self.input_dict.get("components",{})

        # check if mass in dictionary
        given_total_mass    = "mass"    in self.input_dict
        given_total_density = "density" in self.input_dict

        # check if we are given english units
        self.given_english_units = self.input_dict.get("units","English") == \
            "English"

        if given_total_mass:
            self.mass = self.input_dict.get("mass")
        elif given_total_density:
            self.density = self.input_dict.get("density")

        wing_types = ["prismoid","pseudo_prismoid","symmetric_airfoil",
        "diamond_airfoil"]
        comp_types = ["cuboid","cylinder","sphere","rotor"]

        # initialize get each id number
        ids = []; attach_ids = []; name = []
        for component in components:
            input_dict = components[component]
            ids.append(input_dict.get("ID"))
            input_dict["connect_to"] = input_dict.get("connect_to",{})
            attach_ids.append(input_dict["connect_to"].get("ID",0))
            name.append(component)
        
        # throw error if an id given is repeated
        unique_ids,counts = np.unique(ids,return_counts=True)
        repeated_ids = unique_ids[counts > 1]
        if len(unique_ids) != len(ids):
            raise TypeError("Repeated ID values: {}".format(repeated_ids))
               
        # reorganize attach_ids in order of run
        run_order = np.argsort(attach_ids).tolist()
        name_order = np.array(name)[run_order].tolist()
        attach_id_order = np.array(attach_ids)[run_order].tolist()

        # initialize components, save zero component
        self.components = {}
        self.components[0] = Component({"mass":0.0})


        for i in range(len(name_order)):

            component = name_order[i]

            # define component input dictionary
            input_dict = components[component]
            id_number = input_dict.get("ID")

            # overwrite mass and or densities if given
            if given_total_mass:
                input_dict["density"] = 1.0
                if "mass" in input_dict: input_dict.pop("mass")
            elif given_total_density:
                input_dict["density"] = self.density
                if "mass" in input_dict: input_dict.pop("mass")
            
            # shift to root or tip values if attached to a wing
            input_dict["connect_to"] = input_dict.get("connect_to",{})
            attach_loc = input_dict["connect_to"].get("location","tip")
            loc = self.components[attach_id_order[i]].locations[attach_loc]
            input_dict["connect_to"]["dx"] = \
                input_dict["connect_to"].get("dx",0.0) + loc[0,0]
            input_dict["connect_to"]["y_offset"] = \
                input_dict["connect_to"].get("y_offset",0.0) + loc[1,0]
            input_dict["connect_to"]["dz"] = \
                input_dict["connect_to"].get("dz",0.0) + loc[2,0]


            # initialize wings
            if input_dict["type"] in wing_types:
                self.components[id_number] = Wing(input_dict)
            elif input_dict["type"] == "cuboid":
                self.components[id_number] = Cuboid(input_dict)
            elif input_dict["type"] == "cylinder":
                self.components[id_number] = Cylinder(input_dict)
            elif input_dict["type"] == "sphere":
                self.components[id_number] = Sphere(input_dict)
            elif input_dict["type"] == "rotor":
                self.components[id_number] = Rotor(input_dict)
            
            # save "name"
            if input_dict["type"] in wing_types + comp_types:
                self.components[id_number].name = component
        
        # calculate total volume
        self.volume = 0.0
        for comp_id in self.components:
            self.volume += self.components[comp_id].volume
        
        # rewrite densities if given total mass
        if given_total_mass:
            # calculate total density, apply to each
            self.density = self.mass / self.volume
            for comp_id in self.components:
                self.components[comp_id].density = self.density * 1.
                self.components[comp_id].update_densities()


    def report_as_SolidWorks_report(self,info,positive_tensor=True,use_Lanham=False,name=""):
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

        if name == "":
            partin = ""
        else:
            partin = name + " in "
        if use_Lanham:
            symbol = "-*"
            intro = "Lanham "
        else:
            symbol = "=="
            intro = ""
        print((symbol * 50)[:-3])
        print(intro+"Mass properties of",partin + self.file_name)
        print()

        # fix units if not given english units
        if not self.given_english_units:
            info["mass"] /= 14.59390
            info["cg_location"] /= 0.3048
            info["volume"] /= (0.3048**3.)
            info["inertia_tensor"] /= 14.59390 * (0.3048)**2.

        print(intro+"Mass = {:> 10.8f} slugs".format(info["mass"]))
        print()

        print(intro+"Volume = {:> 10.8f} cubic feet".format(info["volume"]))
        print()

        print(intro+"Center of mass: (feet)")
        print("\tX = {:> 14.8f}".format(info["cg_location"][0,0]))
        print("\tY = {:> 14.8f}".format(info["cg_location"][1,0]))
        print("\tZ = {:> 14.8f}".format(info["cg_location"][2,0]))
        print()
        
        print(intro+"Angular momentum: (slugs * square feet / seconds)")
        print("\tX = {:> 14.8f}".format(info["angular_momentum"][0,0]))
        print("\tY = {:> 14.8f}".format(info["angular_momentum"][1,0]))
        print("\tZ = {:> 14.8f}".format(info["angular_momentum"][2,0]))
        print()
        
        I = info["inertia_tensor"] * 1.0
        if positive_tensor:
            I[[0,0,1,1,2,2],[1,2,0,2,0,1]] *= -1.0
        [[Ixx,Ixy,Ixz],[Iyx,Iyy,Iyz],[Izx,Izy,Izz]] = I
        print(intro+"Moment of inertia about the CG:( slugs * square feet)")
        if positive_tensor:
            print("\t\tPositive Tensor Formulation")
        else:
            print("\t\tNegative Tensor Formulation")
        print("\tIxx = {:> 19.8f}\tIxy = {:> 19.8f}\tIxz = {:> 19.8f}".format(\
            Ixx,Ixy,Ixz))
        print("\tIyx = {:> 19.8f}\tIyy = {:> 19.8f}\tIyz = {:> 19.8f}".format(\
            Iyx,Iyy,Iyz))
        print("\tIzx = {:> 19.8f}\tIzy = {:> 19.8f}\tIzz = {:> 19.8f}".format(\
            Izx,Izy,Izz))
        print()
        
        I = info["origin_inertia_tensor"] * 1.0
        if positive_tensor:
            I[[0,0,1,1,2,2],[1,2,0,2,0,1]] *= -1.0
        [[Ixx,Ixy,Ixz],[Iyx,Iyy,Iyz],[Izx,Izy,Izz]] = I
        print(intro+"Moment of inertia about the origin:( slugs * square feet)")
        if positive_tensor:
            print("\t\tPositive Tensor Formulation")
        else:
            print("\t\tNegative Tensor Formulation")
        print("\tIxx = {:> 19.8f}\tIxy = {:> 19.8f}\tIxz = {:> 19.8f}".format(\
            Ixx,Ixy,Ixz))
        print("\tIyx = {:> 19.8f}\tIyy = {:> 19.8f}\tIyz = {:> 19.8f}".format(\
            Iyx,Iyy,Iyz))
        print("\tIzx = {:> 19.8f}\tIzy = {:> 19.8f}\tIzz = {:> 19.8f}".format(\
            Izx,Izy,Izz))
        print()
        print((symbol * 50)[:-3])


    def get_mass_properties(self,report=False,individual=False,use_Lanham=False,positive_tensor=True):

        # determine properties of each component
        for i in self.components:
            self.components[i].get_mass_properties()
        
        # determine total mass
        self.mass = 0.0
        for i in self.components:
            self.mass += self.components[i].get_mass()
        
        # determine lanham mass
        self.mass_lanham = 0.0
        for i in self.components:
            self.mass_lanham += self.components[i].get_mass(True)
        
        # determine lanham volume
        self.volume_lanham = 0.0
        for i in self.components:
            self.volume_lanham += self.components[i].get_volume(True)
        
        # determine total cg location
        self.cg_location = np.zeros((3,1))
        for i in self.components:
            self.cg_location += self.components[i].mass * \
                self.components[i].cg_location
        self.cg_location /= self.mass
        
        # determine lanham cg location
        self.cg_location_lanham = np.zeros((3,1))
        for i in self.components:
            self.cg_location_lanham += self.components[i].get_mass(True) * \
                self.components[i].get_cg_location(True)
        self.cg_location_lanham /= self.mass_lanham
        
        # determine total angular momentum
        self.angular_momentum = 0.0
        for i in self.components:
            self.angular_momentum += self.components[i].angular_momentum
        
        # determine total inertia
        self.inertia_tensor = np.zeros((3,3))
        location = self.cg_location
        for i in self.components:
            self.inertia_tensor += \
                self.components[i].shift_properties_to_location(\
                location)["inertia_tensor"]
        
        # determine lanham inertia
        self.inertia_tensor_lanham = np.zeros((3,3))
        location = self.cg_location_lanham
        for i in self.components:
            self.inertia_tensor_lanham += \
                self.components[i].shift_properties_to_location(\
                location,True)["inertia_tensor"]
        
        # determine origin inertia
        self.origin_inertia_tensor = np.zeros((3,3))
        location = np.zeros((3,1))
        for i in self.components:
            self.origin_inertia_tensor += \
                self.components[i].shift_properties_to_location(\
                location)["inertia_tensor"]
        
        # determine lanham origin inertia
        self.origin_inertia_tensor_lanham = np.zeros((3,3))
        location = np.zeros((3,1))
        for i in self.components:
            self.origin_inertia_tensor_lanham += \
                self.components[i].shift_properties_to_location(\
                location,True)["inertia_tensor"]

        # return dictionary of values
        self.properties_dict = {
            "mass" : self.mass,
            "volume" : self.volume,
            "cg_location" : self.cg_location,
            "angular_momentum" : self.angular_momentum,
            "origin_inertia_tensor" : self.origin_inertia_tensor,
            "inertia_tensor" : self.inertia_tensor
        }

        # Lanham dictionary of values
        lanham = {
            "mass" : self.mass_lanham,
            "volume" : self.volume_lanham,
            "cg_location" : self.cg_location_lanham,
            "angular_momentum" : self.angular_momentum,
            "origin_inertia_tensor" : self.origin_inertia_tensor_lanham,
            "inertia_tensor" : self.inertia_tensor_lanham
        }

        # report
        if report:
            if not individual:
                self.report_as_SolidWorks_report(self.properties_dict,\
                    positive_tensor)
                if use_Lanham:
                    self.report_as_SolidWorks_report(lanham,positive_tensor,\
                        use_Lanham)
            else:
                for i in self.components:
                    if i != 0:
                        info = self.components[i].properties_dict
                        info["origin_inertia_tensor"] = \
                            self.components[i].shift_properties_to_location(\
                            np.zeros((3,1)))["inertia_tensor"]
                        name = self.components[i].name
                        self.report_as_SolidWorks_report(info,positive_tensor,False,name)
                        if use_Lanham:
                            info = self.components[i].shift_properties_to_location(\
                                self.cg_location_lanham,True)
                            info["origin_inertia_tensor"] = \
                                self.components[i].shift_properties_to_location(\
                                np.zeros((3,1)),True)["inertia_tensor"]
                            self.report_as_SolidWorks_report(info,positive_tensor,use_Lanham,name)
        
        if use_Lanham:
            return self.properties_dict, lanham
        else:
            return self.properties_dict


    def get_mass_properties_about_point(self,point,report=False,individual=False,positive_tensor=True):

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
        
        # determine total angular momentum
        self.angular_momentum = 0.0
        for i in self.components:
            self.angular_momentum += self.components[i].angular_momentum
        
        # determine total inertia
        self.inertia_tensor = np.zeros((3,3))
        location = point
        for i in self.components:
            self.inertia_tensor += \
                self.components[i].shift_properties_to_location(\
                location)["inertia_tensor"]

        # return dictionary of values
        self.properties_dict = {
            "mass" : self.mass,
            "volume" : self.volume,
            "cg_location" : self.cg_location,
            "angular_momentum" : self.angular_momentum,
            "inertia_tensor" : self.inertia_tensor
        }

        # report
        if report:
            if not individual:
                self.report_as_SolidWorks_report(self.properties_dict,\
                    positive_tensor)
            else:
                for i in self.components:
                    if i != 0:
                        info = self.components[i].properties_dict
                        info["origin_inertia_tensor"] = \
                            self.components[i].shift_properties_to_location(\
                            np.zeros((3,1)))["inertia_tensor"]
                        name = self.components[i].name
                        self.report_as_SolidWorks_report(info,positive_tensor,name)
        
        return self.properties_dict


if __name__ == "__main__":
    # AS = AircraftSystem("test_input.json")
    # AS.get_mass_properties()
    # AS = AircraftSystem("hunsaker_test.json")
    # AS.get_mass_properties(report=True)
    # AS = AircraftSystem("simple_foam_wings.json")
    # AS.get_mass_properties(report=True,individual=True)
    # AS = AircraftSystem("CRM.json")
    # AS.get_mass_properties(report=True)
    AS = AircraftSystem("horizon.json")
    AS.get_mass_properties(report=True,use_Lanham=True)#,individual=True)
    # AS = AircraftSystem("propeller.json")
    # AS.get_mass_properties(report=True)#,individual=True)
    # AS = AircraftSystem("test_untwisted_propeller.json")
    # AS.get_mass_properties(report=True)#,individual=True)
    # AS = AircraftSystem("test_alternate_untwisted_propeller.json")
    # AS.get_mass_properties(report=True)#,individual=True)