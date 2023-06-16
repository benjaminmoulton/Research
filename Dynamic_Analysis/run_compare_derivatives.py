import compare_derivatives as cd

if __name__ == "__main__":

    # compare general aviation, RC glider, and fighter aircraft
    run_files = ["F_16.json",
    "NT_33A.json","F_104A.json","F_4C.json","X_15.json","HL_10.json",
    "Lockheed_Jetstar.json","Convair_880M.json","boeing_747.json","C_5A.json",
    "XB_70A.json",
    "A_7A.json", "A_4D.json", "F_105B.json", "navion.json", "DC_8.json",
    "Cessna_172.json"]
    run_files = ["aircraft_database/" + i for i in run_files]
    cd.Comparison(run_files)
