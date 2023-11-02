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

    # look for "timeseries" in graph?content.... under Sources after ctrl + shift + I

    t = np.arange(1800.0,2020.0,step = 1.0)
    it = np.argwhere(t>=1900.0)[0,0]

    # find from "timeseries" in file after ctrl + shift + I
    ngram_files = [ # check all these...
        "burble.txt", "burbled.txt", "burbles.txt", "burbling.txt",
        "Burble_.txt", "BURBLE__.txt", "Burbles_.txt", "BURBLES__.txt", 
        "Burbling_.txt", 

        "burble_fence.txt", 

        "burble_point.txt", 

        "separation_region.txt",

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
        "gurgle.txt",
    ]

    ngrams = {}
    folder = "./ngrams/"
    plots_folder = "./plots/"
    plots_folder_inv = "./plots_inv/"
    show_plots = False

    # color blind colors
    cs = ["#F5793A","#A95AA1","#85C0F9","#0F2080",
    '#377eb8', '#4daf4a',
    '#f781bf', '#a65628', '#984ea3',
    '#999999', '#e41a1c', '#dede00']
    cs_inv = ["#0a86c5","#56a55e","#7a3f06","#f0df7f",
    '#c88147', '#b250b5',
    '#087e40', '#59a9d7', '#67b15c',
    '#666666', '#1be5e3', '#2121ff']

    for i,file in enumerate(ngram_files):
        ngram = file.split(".")[0]
        ngrams[ngram] = get_info(folder + file)
        c_i = i % len(cs)
        plt.plot(t,ngrams[ngram],c=cs[c_i],label = ngram.replace("_"," "))
    
    plt.xlabel("Year")
    plt.ylabel("Word Use as Percentage of Annual Words Published")
    plt.legend()
    if show_plots:
        plt.show()
    else:
        plt.close()
    
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
    ngrams_simp["separation_region"] = (
        + ngrams["separation_region"]
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
    ngrams_simp["gurgle"] = (
        + ngrams["gurgle"]
    )

    # plot condensed plot
    order = ["burble","flow_separation",
    "compressibility_burble","separation_region","burble_point","burble_fence"]

    for run_dark in [True,False]:
        if run_dark:
            save_folder = plots_folder_inv
            plt.style.use(['dark_background'])
            csr = cs_inv
        else:
            save_folder = plots_folder
            plt.style.use(['default'])
            csr = cs
        
        # run
        for i in range(len(order)):
            ngram = order[i]
            lbl = ngram.replace("_"," ")
            plt.plot(t[it:],ngrams_simp[ngram][it:],c=csr[i],label=lbl)
        
        plt.xlabel("Year")
        plt.ylabel("Word Use / Annual Words Published")
        plt.legend()
        plt.savefig(fname=save_folder+"condensed.pdf",dpi=300.0,\
            transparent=True)
        if show_plots:
            plt.show()
        else:
            plt.close()

    # plot comparison plot
    order = ["burble","babble","bubble","gurgle"]

    for run_dark in [True,False]:
        if run_dark:
            save_folder = plots_folder_inv
            plt.style.use(['dark_background'])
            csr = cs_inv
        else:
            save_folder = plots_folder
            plt.style.use(['default'])
            csr = cs
        
        # run
        for i in range(len(order)):
            ngram = order[i]
            lbl = ngram.replace("_"," ")
            plt.plot(t[it:],ngrams_simp[ngram][it:],c=csr[i],label=lbl)
        
        plt.xlabel("Year")
        plt.ylabel("Word Use / Annual Words Published")
        plt.legend()
        plt.savefig(fname=save_folder+"comparison.pdf",dpi=300.0,\
            transparent=True)
        if show_plots:
            plt.show()
        else:
            plt.close()

    # plot normalized comparison plot
    order = ["burble","babble","bubble","gurgle"]

    for run_dark in [True,False]:
        if run_dark:
            save_folder = plots_folder_inv
            plt.style.use(['dark_background'])
            csr = cs_inv
        else:
            save_folder = plots_folder
            plt.style.use(['default'])
            csr = cs
        
        # run
        for i in range(len(order)):
            ngram = order[i]
            lbl = ngram.replace("_"," ")
            div_2019 = ngrams_simp[ngram][it:] / ngrams_simp[ngram][-1]
            plt.plot(t[it:],div_2019,c=csr[i],label=lbl)
        
        plt.xlabel("Year")
        plt.ylabel("Word Use / Latest Word Use")
        plt.legend()
        plt.savefig(fname=save_folder+"norm_comparison.pdf",dpi=300.0,\
            transparent=True)
        if show_plots:
            plt.show()
        else:
            plt.close()