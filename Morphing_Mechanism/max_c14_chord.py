import C14
import numpy as np
import json

# function that reads in a json file as a dictionary
def read_json(filename):
    # import json file
    json_string=open(filename).read()
    vals = json.loads(json_string)
    return vals

def main(json_file):
    # read in vals
    vals = read_json(json_file)

    chord = np.linspace(4.545,4.55,1000)

    c0 = vals["C14"]["chord [in]"]
    t0 = vals["C14"]["tongue start [in]"]
    m0 = vals["C14"]["mouth start [in]"]

    for i in range(chord.shape[0]):

        vals["C14"]["chord [in]"] = chord[i]
        vals["C14"]["tongue start [in]"] = chord[i] * t0/c0
        vals["C14"]["mouth start [in]"]  = chord[i] * m0/c0

        try:
            C14.main(vals)
        except:
            print("chord length {:<10.8f} failed".format(chord[i]))
        else:
            print("SUCCESS!!! chord length {:<5.4f}".format(chord[i]))
            break
    
    return

jsonfile = 'input_c14.json'
main(jsonfile)

"""

Max chord length is on the order of 4.5 in

"""