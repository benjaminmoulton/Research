import numpy as np
from matplotlib import pyplot as plt

def get_info(filename):
    with open(filename,"r") as f:

        # initialize points array
        data = []
        # run through each line
        for line in f:
            # split the line
            data.append( float(line.split(",")[0]) )
        # close file
        f.close()
    
    # turn points list into numpy array
    data = np.array(data)

    return data

if __name__ == "__main__":

    t = np.arange(1800.0,2020.0,step = 1.0)

    burble = get_info("burble.txt")
    
    plt.plot(t,burble)
    plt.show()