import machupX as mx
import os
import json

path0 = "/home/ben/Desktop/VO_input.json"
path1 = "V1_input.json"

# # import json file
# json_string=open(path0).read()
# vals = json.loads(json_string)
THIS_FOLDER = os.path.dirname(os.path.abspath(__file__))
my_file = os.path.join(THIS_FOLDER, 'V2_input.json')

v0 = mx.Scene(my_file)

# see wireframe
v0.display_wireframe(show_vortices=False)

v1 = mx.Scene(path1)

# see wireframe
v1.display_wireframe(show_vortices=False)