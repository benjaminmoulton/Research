import machupX as mx
import numpy as np
import json

def pyaic_to_mux(filename):
    # read in json file
    json_string = open(filename).read()
    airplane_dict = json.loads(json_string)

    # rename wings dictionary
    airplane_dict["wings"] = airplane_dict.pop("components")
    key0 = list(airplane_dict["wings"].keys())[0]
    airplane_dict["wings"][key0]["is_main"] = True

    # add weight to dictionary
    airplane_dict["weight"] = 1.0

    # make scene dict
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
    # my_scene.display_wireframe(show_vortices=False)
    my_scene.export_dxf(number_guide_curves=6)


if __name__ == "__main__":
    pyaic_to_mux("simple_foam_wings.json")