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

    ngram_files = [ # check all these...
        "burble.txt", "burbled.txt", "burbles.txt", "burbling.txt",
        "compressibility_burble.txt", # only _INF version
        "flow_separated.txt", "flow_separation.txt", "flow_separations.txt",
        "separation_of_flow.txt", "separation_of_flowing.txt", 
        "separation_of_flows.txt",
    ]

    ngrams = {}
    folder = "./ngrams/"

    for file in ngram_files:
        ngram = file.split(".")[0]
        ngrams[ngram] = get_info(folder + file)
        plt.plot(t,ngrams[ngram],label = ngram.replace("_"," "))
    
    plt.xlabel("Year")
    plt.ylabel("Word Use as Percentage of Annual Words Published")
    plt.legend()
    plt.show()
    
    ngrams_simp = {}
    ngrams_simp["burble"] = (ngrams["burble"] + ngrams["burbled"] 
        + ngrams["burbles"] + ngrams["burbling"] 
        - ngrams["compressibility_burble"])
    ngrams_simp["compressibility_burble"] = ngrams["compressibility_burble"]
    ngrams_simp["flow_separation"] = (ngrams["flow_separated"] 
        + ngrams["flow_separation"] + ngrams["flow_separations"] 
        + ngrams["separation_of_flow"] + ngrams["separation_of_flowing"]
        + ngrams["separation_of_flows"])

    for ngram in ngrams_simp:
        plt.plot(t,ngrams_simp[ngram],label = ngram.replace("_"," "))
    
    plt.xlabel("Year")
    plt.ylabel("Word Use as Percentage of Annual Words Published")
    plt.legend()
    plt.show()