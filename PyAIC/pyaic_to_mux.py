import machupX as mx
import numpy as np
import json

from wing import Wing

def pyaic_to_mux(filename,force_symmetric=False,untwist=False,straighten_c_4=False,tag=""):
    # read in json file
    json_string = open(filename).read()
    airplane_dict = json.loads(json_string)

    # rename wings dictionary
    airplane_dict["wings"] = airplane_dict.pop("components")
    key0 = list(airplane_dict["wings"].keys())[0]
    airplane_dict["wings"][key0]["is_main"] = True
    
    # force airfoils to be symmetric airfoils of specified thickness
    if force_symmetric:

        # initialize x values
        n = 200 # airfoil points
        theta = np.linspace(np.pi,0.0,n)
        x = - 0.5 * np.cos(theta) + 0.5
        x_lo = 1.0 * np.flip(x)[1:]

        for wing in airplane_dict["wings"]:
        
            # read in a coefficients
            avals = airplane_dict["wings"][wing].get(\
                "thickness_distribution_coefficients","open_trailing_edge")
            if avals == "open_trailing_edge":
                avals = [2.969, -1.260, -3.516, 2.843, -1.015]
            elif avals == "closed_trailing_edge":
                avals = [2.969, -1.260, -3.516, 2.843, -1.036]
            elif avals == "closed_hunsaker":
                avals = [2.980, -1.320, -3.286, 2.441, -0.815]
            a0 = avals[0] * 1.0
            a1 = avals[1] * 1.0
            a2 = avals[2] * 1.0
            a3 = avals[3] * 1.0
            a4 = avals[4] * 1.0

            # determine airfoil shape
            y_up = (a0*x**0.5 + a1*x + a2*x**2. + a3*x**3. + a4*x**4.)/2.
            y_lo = -1. * np.flip(y_up)[1:]

            # create airfoil shape numpy array
            xs =    x.tolist() + x_lo.tolist()
            ys = y_up.tolist() + y_lo.tolist()
            xys = np.vstack((xs,ys))
            
            thickness = airplane_dict["wings"][wing].get("thickness",0.12)
            if isinstance(thickness,(float,int)):
                thickness = [[0.0,thickness*1.0],[1.0,thickness*1.0]]
            thickness = np.array(thickness)

            # run through airfoils, create symmetric airfoil, replace in airfoils
            af0 = list(airplane_dict["airfoils"].keys())[0]
            airfoil = airplane_dict["wings"][wing].get("airfoil",af0)
            if isinstance(thickness,(float,int)):
                airfoil = [[0.0,airfoil],[1.0,airfoil]]
            for i in range(len(airfoil)):
                # determine local thickness
                tk = np.interp(airfoil[i][0],thickness[:,0],thickness[:,1])
                
                # determine airfoil
                afi = xys * 1.0
                afi[:,1] *= tk

                # replace airfoil in airfoils
                airplane_dict["airfoils"][airfoil[i][1]].pop("NACA","0008")
                airplane_dict["airfoils"][airfoil[i][1]]["outline_points"] = \
                    afi.tolist()

    # remove twist from each wing
    if untwist:
        for wing in airplane_dict["wings"]:
            airplane_dict["wings"][wing].pop("twist",0)
    
    # straighten quarter chord
    if straighten_c_4:
        wing_keys = list(airplane_dict["wings"].keys())
        for wing in wing_keys:
            old_wing = airplane_dict["wings"].pop(wing)

            # add mass if not in dictionary
            old_wing["mass"] = 1.0
            class_wing = Wing(old_wing)
    
    quit()

    # add weight to dictionary
    airplane_dict["weight"] = 1.0

    # make scene dict
    if tag == "":
        tag = filename.split(".")[0]
    scene_dict = {
        "tag" : tag,
        "scene" : {
            "aircraft" : {
                tag : {
                    "file" : airplane_dict,
                    "state" : {
                            "position": [0.0,0.0,0.0],
                        "velocity" : 100.0,
                        "alpha" : 0.0,
                        "beta" : 0.0
                    }
                }
            }
        }
    }

    my_scene = mx.Scene(scene_dict)
    my_scene.display_wireframe(show_vortices=False)
    # my_scene.export_dxf(number_guide_curves=6)


if __name__ == "__main__":
    # pyaic_to_mux("simple_foam_wings.json")
    # pyaic_to_mux("CRM.json",tag="CRM_OML") # OML
    # pyaic_to_mux("CRM.json",force_symmetric=True,tag="CRM_symm")
    # pyaic_to_mux("CRM.json",force_symmetric=True,untwist=True,tag="CRM_notwist")
    pyaic_to_mux("CRM.json",force_symmetric=True,untwist=True,straighten_c_4=True,tag="CRM_strait")