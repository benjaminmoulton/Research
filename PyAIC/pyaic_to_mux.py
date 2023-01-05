import machupX as mx
import numpy as np
import json
import copy
from matplotlib import pyplot as plt

from system import AircraftSystem

def pyaic_to_mux(filename,force_symmetric=False,untwist=False,straight_c_4=False,English_units=True,tag=""):
    # report
    print("initializing",tag,"from",filename)

    # read in json file
    json_string = open(filename).read()
    airplane_dict = json.loads(json_string)

    # rename wings dictionary
    system_dict = copy.deepcopy(airplane_dict)
    airplane_dict["wings"] = airplane_dict.pop("components")
    
    # force airfoils to be symmetric airfoils of specified thickness
    if force_symmetric:
        print("\tmaking airfoils symmetric...")

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
            y_up = -(a0*x**0.5 + a1*x + a2*x**2. + a3*x**3. + a4*x**4.)/2.
            y_lo = -1. * np.flip(y_up)[1:]

            # create airfoil shape numpy array
            xs =    x.tolist() + x_lo.tolist()
            ys = y_up.tolist() + y_lo.tolist()
            xys = np.vstack((xs,ys)).T
            
            thickness = airplane_dict["wings"][wing].get("thickness",0.12)
            if isinstance(thickness,(float,int)):
                thickness = [[0.0,thickness*1.0],[1.0,thickness*1.0]]
            thickness = np.array(thickness)

            # run through airfoils, create symmetric airfoil, replace in airfoils
            af0 = list(airplane_dict["airfoils"].keys())[0]
            airfoil = airplane_dict["wings"][wing].get("airfoil",af0)
            if isinstance(airfoil,str):
                airfoil = [[0.0,airfoil],[1.0,airfoil]]
            for i in range(len(airfoil)):
                # determine local thickness
                tk = np.interp(airfoil[i][0],thickness[:,0],thickness[:,1])
                
                # determine airfoil
                afi = xys * 1.0
                afi[:,1] *= tk

                # plt.plot(afi[:,0],afi[:,1],"k")
                # plt.axis("equal")
                # plt.show()

                # replace airfoil in airfoils
                airplane_dict["airfoils"][airfoil[i][1]]["geometry"].pop("NACA","0008")
                # airplane_dict["airfoils"][airfoil[i][1]]["geometry"].pop("outline_points",0)
                # airplane_dict["airfoils"][airfoil[i][1]]["geometry"]["NACA"] = "0012"
                airplane_dict["airfoils"][airfoil[i][1]]["geometry"]["outline_points"] = \
                    afi

    # remove twist from each wing
    if untwist:
        print("\tremoving twist...")
        for wing in airplane_dict["wings"]:
            airplane_dict["wings"][wing].pop("twist",0)
    
    # straighten quarter chord
    if straight_c_4:
        print("\tstraightening quarter-chord line...")
        System = AircraftSystem(system_dict)
        wing_keys = list(airplane_dict["wings"].keys())
        id_counter = \
            max([airplane_dict["wings"][wing]["ID"] for wing in wing_keys]) + 1
        for wing in wing_keys:
            old_wing = airplane_dict["wings"].pop(wing)
            # add mass if not in dictionary
            old_wing["mass"] = 1.0
            curr_id = old_wing["ID"]
            class_wing = System.components[curr_id]

            # run through sub-wings
            counter = 0
            for subwing in class_wing._component_inputs:
                subwing_dict = class_wing._component_inputs[subwing].copy()
                subwing_name = str(wing) + "_" + str(counter)

                # connect_to
                root_loc = np.copy(subwing_dict.pop("root_location"))
                dict1 = { 
                    "connect_to" : {
                        "dx" : root_loc[0] * 1.0,
                        "dy" : root_loc[1] * 1.0,
                        "dz" : root_loc[2] * 1.0
                    }
                }

                # geometry
                geom_dict = subwing_dict.pop("geometry")
                dict2 = {
                    "ID" : id_counter,
                    "semispan" : geom_dict["span"] * 1.0,
                    "chord" : [ [0.0,geom_dict["chord"][0] * 1.0],
                                [1.0,geom_dict["chord"][1] * 1.0]],
                    "thickness" : [ [0.0,geom_dict["thickness"][0] * 1.0],
                                    [1.0,geom_dict["thickness"][1] * 1.0]],
                    "airfoil" : [   [0.0,geom_dict["airfoil"][0]],
                                    [1.0,geom_dict["airfoil"][1]]]
                }

                # remove other pieces
                tip_loc = subwing_dict.pop("tip_location")
                subwing_dict.pop("thickness_distribution_coefficients")
                subwing_dict.pop("density")

                # new dictionary
                new_wing_dict = {**dict2,**subwing_dict,**dict1}.copy()

                airplane_dict["wings"][subwing_name] = new_wing_dict
                
                # if filename.split(".")[0] == "CRM":
                #     pre = "& "
                #     ratio = 1.0 / 0.3048
                # else:
                #     pre = ""
                #     ratio = 1.0
                # if counter <= class_wing._num_spans-2:
                #     q = 0
                #     b = geom_dict["span"]
                #     cr = geom_dict["chord"][q] * ratio
                #     tr = geom_dict["thickness"][q] * 100.0
                #     mr = geom_dict["camber"][q] * 100.0
                #     xr = root_loc[0] * ratio
                #     yr = root_loc[1] * ratio
                #     zr = root_loc[2] * ratio
                #     Gr = subwing_dict["dihedral"]
                #     Lr = subwing_dict["sweep"]
                #     print(pre+"{:>4} & {:>9.5f} & {:>9.5f} & {:>9.5f} & ".format(counter+1,b,cr,tr),end="")
                #     print("{:>10.5f} & {:>10.5f} & {:>10.5f} & ".format(xr,yr,zr),end="")
                #     print("{:>10.5f} & {:>10.5f} \\\\".format(Gr,Lr))#,end="")
                # if counter == class_wing._num_spans-2:
                #     q = 1
                #     cr = geom_dict["chord"][q] * ratio
                #     tr = geom_dict["thickness"][q] * 100.0
                #     mr = geom_dict["camber"][q] * 100.0
                #     xr = tip_loc[0] * ratio
                #     yr = tip_loc[1] * ratio
                #     zr = tip_loc[2] * ratio
                #     Gr = subwing_dict["dihedral"]
                #     Lr = subwing_dict["sweep"]
                #     print(pre+"{:>4} & {:>9} & {:>9.5f} & {:>9.5f} & ".format("tip","-",cr,tr),end="")
                #     print("{:>10} & {:>10} & {:>10} & ".format("-","-","-"),end="")
                #     print("{:>10} & {:>10} \\\\".format("-","-"))#,end="")
                #     print()

                id_counter += 1
                counter += 1
            
    #         print("----------------------")
    # quit()
    # make one wing main
    key0 = list(airplane_dict["wings"].keys())[0]
    airplane_dict["wings"][key0]["is_main"] = True

    # add weight to dictionary
    airplane_dict["weight"] = 1.0

    # new_dict = json.loads(open(filename).read())
    # new_dict["wings"] = new_dict.pop("components")
    # key0 = list(new_dict["wings"].keys())[0]
    # new_dict["wings"][key0]["is_main"] = True
    # for wing in new_dict["wings"]:
    #     new_dict["wings"][wing].pop("twist",0)
    # new_dict["weight"] = 1.0

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
                }#,
                # "perfect" : {
                #     "file" : new_dict,
                #     "state" : {
                #             "position": [0.0,0.0,0.0],
                #         "velocity" : 100.0,
                #         "alpha" : 0.0,
                #         "beta" : 0.0
                #     }
                # }
            }
        }
    }

    if not English_units:
        scene_dict["units"] = "SI"
    

    # print(json.dumps(airplane_dict,sort_keys=True, indent=4))

    print("\tinitializing MachUpX model...")
    my_scene = mx.Scene(scene_dict)
    # print("\tplotting wireframe...")
    # my_scene.display_wireframe(show_vortices=False,show_legend=True)
    print("\texporting DXF files...")
    my_scene.export_dxf(number_guide_curves=12)


if __name__ == "__main__":
    # pyaic_to_mux("curvy.json",force_symmetric=True,straight_c_4=True,tag="curvy_straight")
    # pyaic_to_mux("simple_foam_wings.json")
    # pyaic_to_mux("CRM.json",English_units=False,tag="CRM_OML") # OML # redo this one...
    # pyaic_to_mux("CRM.json",force_symmetric=True,English_units=False,tag="CRM_symm")
    # pyaic_to_mux("CRM.json",force_symmetric=True,untwist=True,English_units=False,tag="CRM_notwist")
    # pyaic_to_mux("CRM.json",force_symmetric=True,untwist=True,straight_c_4=True,English_units=False,tag="CRM_straight")
    # pyaic_to_mux("horizon.json",tag="horizon_OML")
    # pyaic_to_mux("horizon.json",force_symmetric=True,tag="horizon_symm")
    # pyaic_to_mux("horizon.json",force_symmetric=True,straight_c_4=True,tag="horizon_straight")
    # pyaic_to_mux("propeller_for_mux.json",tag="propeller_OML")
    pyaic_to_mux("propeller_for_mux.json",force_symmetric=True,tag="propeller_symm")
