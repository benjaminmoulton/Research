import numpy as np

def organize_equation(equation,main_coeff = "B",main_subcoeff = "a"):
    first,equation = equation.split("*(")
    equation,last  = equation.split(")/")
    test = "+ " + equation

    print("splitting...")
    test_list = test.replace("+ ","+").replace("- ","-").split(" ")

    # initialize strings to assign to
    val_dict = {}
    orig_vals = ["Xle","Xte","ts","twb","tan(Lr)","b_2","xcg","ycg","zcg"]
    key_vals = [k + "**2" for k in orig_vals] + orig_vals
    key_vals = [k + "**3" for k in orig_vals] + key_vals
    key_vals = [k + "**4" for k in orig_vals] + key_vals

    for value in test_list:
        
        dict_folder = ""
        keys = []

        for key in key_vals:
            if key in value and key + "**" not in value:

                dict_folder += key + "*"
                keys.append(key)
        
                # print(dict_folder,keys, value)
        
        # replace components of string
        piece = value + ""
        for key in keys:
            piece = piece.replace("*" + key,"")
        
        # rename folder if for additional pieces
        if len(dict_folder) == 0:
            dict_folder = "1"
        else:
            dict_folder = dict_folder[:-1]
        
        if dict_folder not in val_dict:
            val_dict[dict_folder] = []
        
        # append to array
        val_dict[dict_folder].append(piece)

    # run through and "flip" signs
    val_dict_keys = list(val_dict.keys())
    for key in val_dict_keys:
        # check if all negative
        if all([piece[0]=="-" for piece in val_dict[key]]):
            # remove the "-", add to key name
            pieces = [piece[1:] for piece in val_dict[key]]
            new_key = "-" + key

            # remove old, replace
            val_dict[new_key] = pieces
            val_dict.pop(key)
            
        elif all([piece[0]=="+" for piece in val_dict[key]]):
            # remove the "+", add to key name
            pieces = [piece[1:] for piece in val_dict[key]]
            new_key = "+" + key

            # remove old, replace
            val_dict[new_key] = pieces
            val_dict.pop(key)
        else:
            raise Exception("Mixed list of + and - signs. Figure it out.")

    # find greatest common divisor
    val_dict_keys = list(val_dict.keys())
    for key in val_dict_keys:

        # create list of  divisors
        divisors = []
        safe = True
        for i in range(len(val_dict[key])):

            number = val_dict[key][i].split("*")[0]

            try:
                divisors.append(int(number))
            except:
                safe = False
                break
        
        if safe:
            gcd = np.gcd.reduce(divisors)
            
            reduced = [int(divisor / gcd) for divisor in divisors]

            for i in range(len(val_dict[key])):
                # make strings
                strdiv = str(divisors[i])
                strred = str(reduced[i])

                # remove old multiplier
                val_dict[key][i] = val_dict[key][i].replace(strdiv,strred)

                if val_dict[key][i][:2] == "1*": 
                    val_dict[key][i] = val_dict[key][i][2:]
            
            # modify key name
            new_key = key[:1] + str(gcd) + "*" + key[1:]
            val_dict[new_key] = val_dict.pop(key)

    # # find ones that should be split up
    # val_dict_keys = list(val_dict.keys())
    # te = ["tr","tt"]
    # for k in val_dict_keys:
    #     # check if all have a tr or tt term or none have it
    #     if all([te[0] in val_dict[k][i] or te[1] in val_dict[k][i] \
    #         for i in range(len(val_dict[k]))]):
    #         pass
    #     elif all([te[0] not in val_dict[k][i] and te[1] not in val_dict[k][i] \
    #         for i in range(len(val_dict[k]))]):
    #         pass
    #     else:
            
    #         # split up the family
    #         have_t = []
    #         havent = []
    #         havent_key = "1*"

    #         for i in range(len(val_dict[k])):
    #             if te[0] in val_dict[k][i] or te[1] in val_dict[k][i]:
    #                 have_t.append(val_dict[k][i])
    #             else:
    #                 havent.append(val_dict[k][i])
            
    #         # create a second setting in the dictionary, add havent's
    #         new_key = k[0] + "1*" + k[1:]
    #         val_dict[new_key] = havent
    #         val_dict[k] = have_t

    # find similar multiplied values
    val_dict_keys = list(val_dict.keys())
    for key in val_dict_keys:
        pop_keys = []
        # determine keys currently in dictionary
        val_dict_subkeys = list(val_dict.keys())
        for subkey in val_dict_subkeys:
            if key != subkey and key in val_dict:
                if val_dict[key] == val_dict[subkey]:
                    pop_keys.append(subkey)
        
        # modify key name
        if pop_keys:
            new_key = key
            for popkey in pop_keys: new_key += " " + popkey
            val_dict[new_key] = val_dict.pop(key)
            for popkey in pop_keys: val_dict.pop(popkey)

    # split off coefficients
    val_dict_keys = list(val_dict.keys())
    coeff = 0
    coeffs = {}
    for key in val_dict_keys:
        # create string from list
        longun = val_dict[key][0]
        for i in range(len(val_dict[key])-1):
            longun += " + " + val_dict[key][i+1]
        
        # initialize new coeffs name
        if coeff > 25:
            main_subcoeff = "A"; coeff = 0
        B_n = main_coeff + chr(ord(main_subcoeff) + coeff)
        coeffs[B_n] = longun
        coeff += 1

        # change in old array
        val_dict[key] = B_n


    # find greatest common divisor of non-coefficients
    val_dict_keys = list(val_dict.keys())
    for key in val_dict_keys:

        # create the list
        components = key.split()

        # create list of  divisors
        divisors = []
        safe = True
        for i in range(len(components)):

            number = components[i].split("*")[0]

            try:
                divisors.append(int(number))
            except:
                safe = False
                break
        
        if safe:
            gcd = np.gcd.reduce(divisors)
            if max(divisors) < 0:
                gcd = -gcd
            
            reduced = [int(divisor / gcd) for divisor in divisors]

            for i in range(len(components)):
                # make strings
                strdiv = str(divisors[i])
                strred = str(reduced[i])

                if gcd < 0:
                    strred = "+" + strred

                # remove old multiplier
                components[i] = components[i].replace(strdiv,strred)

                if components[i][1:3] == "1*": 
                    components[i] = components[i][:1] + components[i][3:]
            
            # create new key
            new_key = ""
            for i in range(len(components)):
                new_key += components[i] + " "
            new_key = new_key[:-1]
            
            # modify key name
            val_dict[new_key] = val_dict.pop(key)

            # modify coefficient
            gcd = str(gcd)
            if gcd[0].isnumeric():
                gcd = "+" + gcd
            val_dict[new_key] = gcd + "*" + val_dict[new_key]

    # create a new string equation
    equation = ""
    for key in val_dict:
        equation += val_dict[key] + " * ( " + key + " ) "
    equation = equation[:-1]

    # add spacing
    equation = equation.replace(" -"," - ").replace(" +"," + ")


    # report
    eq1 = "equation = " + first + "*( "
    print(eq1)
    print(" " * len(eq1) + equation)
    print(" )/" + last)
    print("\nwhere\n")
    for key in coeffs:
        print(key,"=",coeffs[key])
    print()