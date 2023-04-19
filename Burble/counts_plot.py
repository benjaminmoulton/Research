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
    it = np.argwhere(t>=1900.0)[0,0]

    # find from "timeseries" in file after ctrl + shift + I
    ngram_files = [ # check all these...
        "burble.txt", "burbled.txt", "burbles.txt", "burbling.txt",
        "Burble_.txt", "BURBLE__.txt", "Burbles_.txt", "BURBLES__.txt", 
        "Burbling_.txt", 

        "burble_fence.txt", 

        "burble_point.txt", 

        "compressibility_burble.txt", # only _INF version
        "Compressibility_burble_.txt", "Compressibility_Burble__.txt", 
        "COMPRESSIBILITY_BURBLE___.txt", 

        "flow_separate.txt", "flow_separated.txt", 
        "flow_separates.txt", "flow_separating.txt", 
        "flow_separation.txt", "flow_separations.txt",
        "separation_of_flow.txt", "separation_of_flowing.txt", 
        "separation_of_flows.txt",

        "babble.txt", 
        "bubble.txt", 
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
        + ngrams["compressibility_burble"] 
        # + ngrams["Compressibility_burble_"]
        # + ngrams["Compressibility_Burble__"] + ngrams["COMPRESSIBILITY_BURBLE___"]
    )
    ngrams_simp["burble_point"] = (
        + ngrams["burble_point"]
    )
    ngrams_simp["burble_fence"] = (
        + ngrams["burble_fence"]
    )
    ngrams_simp["burble"] = (ngrams["burble"] 
        # + ngrams["burbled"] 
        # + ngrams["burbles"] + ngrams["burbling"]
        # + ngrams["Burble_"] + ngrams["BURBLE__"]
        # + ngrams["Burbles_"] + ngrams["BURBLES__"]
        # + ngrams["Burbling_"]
        - ngrams_simp["compressibility_burble"]
        - ngrams_simp["burble_point"]
        - ngrams_simp["burble_fence"]
    )
    ngrams_simp["flow_separation"] = (
        + ngrams["flow_separation"] 
        # + ngrams["flow_separate"] + ngrams["flow_separated"] 
        # + ngrams["flow_separates"] + ngrams["flow_separating"]
        #  + ngrams["flow_separations"]
        # + ngrams["separation_of_flow"] + ngrams["separation_of_flowing"]
        # + ngrams["separation_of_flows"]
    )
    ngrams_simp["babble"] = (
        + ngrams["babble"]
    )
    ngrams_simp["bubble"] = (
        + ngrams["bubble"]
    )

    # color-blind colors
    cs = ["#F5793A","#A95AA1","#85C0F9","#0F2080","k"]

    # plot condensed plot
    order = ["burble","flow_separation",
    "compressibility_burble","burble_point","burble_fence"]

    for i in range(len(order)):
        ngram = order[i]
        lbl = ngram.replace("_"," ")
        plt.plot(t[it:],ngrams_simp[ngram][it:],c=cs[i],label=lbl)
    
    plt.xlabel("Year")
    plt.ylabel("Word Use / Annual Words Published")
    plt.legend()
    plt.savefig(fname=plots_folder+"condensed.png",dpi=300.0,transparent=True)
    if show_plots:
        plt.show()
    else:
        plt.close()

    # plot comparison plot
    order = ["burble","babble","bubble"]

    for i in range(len(order)):
        ngram = order[i]
        lbl = ngram.replace("_"," ")
        plt.plot(t[it:],ngrams_simp[ngram][it:],c=cs[i],label=lbl)
    
    plt.xlabel("Year")
    plt.ylabel("Word Use / Annual Words Published")
    plt.legend()
    plt.savefig(fname=plots_folder+"comparison.png",dpi=300.0,transparent=True)
    if show_plots:
        plt.show()
    else:
        plt.close()

    # plot normalized comparison plot
    order = ["burble","babble","bubble"]

    for i in range(len(order)):
        ngram = order[i]
        lbl = ngram.replace("_"," ")
        div_2019 = ngrams_simp[ngram][it:] / ngrams_simp[ngram][-1]
        plt.plot(t[it:],div_2019,c=cs[i],label=lbl)
    
    plt.xlabel("Year")
    plt.ylabel("Word Use / Latest Word Use")
    plt.legend()
    plt.savefig(fname=plots_folder+"norm_comparison.png",dpi=300.0,\
        transparent=True)
    if show_plots:
        plt.show()
    else:
        plt.close()