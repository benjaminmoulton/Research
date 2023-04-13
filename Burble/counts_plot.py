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

    # find from "timeseries" in file after ctrl + shift + I
    ngram_files = [ # check all these...
        "burble.txt", "burbled.txt", "burbles.txt", "burbling.txt",
        "Burble_.txt", "BURBLE__.txt", 
        "compressibility_burble.txt", # only _INF version
        "Compressibility_burble_.txt", "Compressibility_Burble__.txt", 
        "COMPRESSIBILITY_BURBLE___.txt", 
        "flow_separate.txt", "flow_separated.txt", 
        "flow_separates.txt", "flow_separating.txt", 
        "flow_separation.txt", "flow_separations.txt",
        "separation_of_flow.txt", "separation_of_flowing.txt", 
        "separation_of_flows.txt",
    ]

    ngrams = {}
    folder = "./ngrams/"
    plots_folder = "./plots/"
    show_plots = True

    for file in ngram_files:
        ngram = file.split(".")[0]
        ngrams[ngram] = get_info(folder + file)
        plt.plot(t,ngrams[ngram],label = ngram.replace("_"," "))
    
    plt.xlabel("Year")
    plt.ylabel("Word Use as Percentage of Annual Words Published")
    plt.legend()
    plt.show()
    
    ngrams_simp = {}
    ngrams_simp["compressibility_burble"] = (
        + ngrams["compressibility_burble"] + ngrams["Compressibility_burble_"]
        + ngrams["Compressibility_Burble__"] + ngrams["COMPRESSIBILITY_BURBLE___"])
    ngrams_simp["burble"] = (ngrams["burble"] + ngrams["burbled"] 
        + ngrams["burbles"] + ngrams["burbling"]
        + ngrams["Burble_"] + ngrams["BURBLE__"]
        - ngrams_simp["compressibility_burble"])
    ngrams_simp["flow_separation"] = (
        ngrams["flow_separate"] + ngrams["flow_separated"] 
        + ngrams["flow_separates"] + ngrams["flow_separating"]
        + ngrams["flow_separation"] + ngrams["flow_separations"]
        + ngrams["separation_of_flow"] + ngrams["separation_of_flowing"]
        + ngrams["separation_of_flows"])

    for ngram in ngrams_simp:
        plt.plot(t,ngrams_simp[ngram],label = ngram.replace("_"," "))
    
    plt.xlabel("Year")
    plt.ylabel("Word Use as Percentage of Annual Words Published")
    plt.legend()
    plt.savefig(fname=plots_folder+"condensed.png",dpi=300.0,transparent=True)
    if show_plots:
        plt.show()
    else:
        plt.close()