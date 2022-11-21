import machupX as mx
import numpy as np
import json

def pyaic_to_mux(filename):
    # read in json file
    json_string = open(filename).read()
    input_dict = json.loads(json_string)

    # save components as wings
    wings = input_dict["components"]

    # add airfoils, this is airplane

    # make scene dict
    tag = filename.split(".")[0]
    scene_dict = {
        "tag" : tag,
        "scene" : {
            "aircraft" : {
                tag : {
                    "file" : wings,
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


if __name__ == "__main__":
    pyaic_to_mux("simple_foam_wings.json")