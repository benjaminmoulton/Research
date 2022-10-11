import numpy as np
import json
import C18_geometry as c14

json_string = open("input_c18.json").read()
vals = json.loads(json_string)

for i in range(31):
    d = -15 + i
    print(d,"\n")
    vals["C18"]["deflection angle [deg]"] = [d]
    if d == 0:
        vals["C18"]["deflection angle [deg]"] = d

    c14.main(vals)