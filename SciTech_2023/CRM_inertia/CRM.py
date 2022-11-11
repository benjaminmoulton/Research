import machupX as mx
import json
import numpy as np

# read in json files
f = open('CRM_MUX_aircraft.json', 'r')
craftDict = json.load(f)
f.close()
f = open('CRM_MUX_scene.json', 'r')
sceneDict = json.load(f)
f.close()

one = "CRM_OML"
sceneDict['scene']['aircraft'][one] = sceneDict['scene']['aircraft'].pop("CRM")
sceneDict['scene']['aircraft'][one]['file'] = dict(craftDict)
CRM_OML = mx.Scene(sceneDict)
# CRM_OML.display_wireframe(show_vortices=False)
# CRM_OML.export_dxf(aircraft=one,number_guide_curves=8)

for i in CRM_OML._airplanes:
    plane = CRM_OML._airplanes[i]
    for j in plane.wing_segments:
        segment = plane.wing_segments[j]
        max_cambers = np.zeros((segment._num_airfoils,))
        max_thicknesses = np.zeros((segment._num_airfoils,))
        for k in range(segment._num_airfoils):
            max_cambers[k] = segment._airfoils[k].get_max_camber()
            max_thicknesses[k] = segment._airfoils[k].get_max_thickness()
        print(max_thicknesses)
        print(max_cambers)
        print()
quit()

for aircraft in CRM_OML._airplanes:
    print(CRM_OML._airplanes[aircraft].wing_segments)
    for _,segment in CRM_OML._airplanes[aircraft].wing_segments.items():
        print(segment)
        print("thickness")
        print("[\t",end="")
        for i in range(len(segment._airfoils)):
            if i != 0:
                print("\t",end="")
            print("[span,\t{}\t]".format(segment._airfoils[i].get_max_thickness()),end="")
            if i != len(segment._airfoils) - 1:
                print(",")
            else:
                print("],")
        print("camber")
        print("[\t",end="")
        for i in range(len(segment._airfoils)):
            if i != 0:
                print("\t",end="")
            print("[span,\t{}\t]".format(segment._airfoils[i].get_max_camber()),end="")
            if i != len(segment._airfoils) - 1:
                print(",")
            else:
                print("],")

two = "CRM_rectangular"
sceneDict['scene']['aircraft'][two] = sceneDict['scene']['aircraft'].pop(one)
sceneDict['scene']['aircraft'][two]['file'] = dict(craftDict)
CRM_rectangular = mx.Scene(sceneDict)
# CRM_rectangular.display_wireframe(show_vortices=False)
CRM_rectangular.export_dxf(aircraft=two,export_as_prismoid=True)

three = "CRM_no_twist"
sceneDict['scene']['aircraft'][three] = sceneDict['scene']['aircraft'].pop(two)
sceneDict['scene']['aircraft'][three]['file'] = dict(craftDict)
sceneDict['scene']['aircraft'][three]['file']["wings"]["main_wing"].pop("twist")
CRM_no_twist = mx.Scene(sceneDict)
# CRM_rectangular.display_wireframe(show_vortices=False)
CRM_no_twist.export_dxf(aircraft=three,export_as_prismoid=True)

# make prismoids object
IDnum = 2

# create constant sweep sections which meet up with horizon other sections
def split_wing(craftDict,name,IDnumber):
    # report
    print("splitting", name,"...")

    # save certain values
    bw_2 = craftDict['wings'][name]["semispan"]
    sweep_old = craftDict['wings'][name]["sweep"]
    chord_old = craftDict['wings'][name]["chord"]
    airfoil_old = craftDict['wings'][name]["airfoil"]
    dihedral_old = craftDict['wings'][name]["dihedral"]
    spans = [i[0] for i in chord_old]

    sweep_is_list = type(sweep_old) == list
    dihedral_is_list = type(dihedral_old) == list
    
    sections = len(spans) - 1
    chord = np.zeros((2*sections,2))

    # cycle through and make sections
    bw_2_new = 0.0
    for i in range(sections):
        j = int(i*2)
        chord[j  ,0] = 0.0 # spans[i  ]
        chord[j+1,0] = 1.0 # spans[i+1]

        # determine local span
        b_loc = bw_2 * (spans[i+1] - spans[i])
        
        # determine sweep angles
        if sweep_is_list:
            sweep_diff  = sweep_old[i+1][1]-sweep_old[i][1]
            sweep_start_deg = spans[i  ] * sweep_diff + sweep_old[i][1]
            sweep_end_deg   = spans[i+1] * sweep_diff + sweep_old[i][1]
        else:
            sweep_diff = 0.0
            sweep_start_deg = sweep_old * 1.0
            sweep_end_deg = sweep_old * 1.0
        
        # determine dihedral angles
        if dihedral_is_list:
            dihedral_diff  = dihedral_old[i+1][1]-dihedral_old[i][1]
            dihedral_start_deg = spans[i  ] * dihedral_diff + dihedral_old[i][1]
            dihedral_end_deg   = spans[i+1] * dihedral_diff + dihedral_old[i][1]
        
            # determine "equivalent" sweep angles for calculation
            c_sw_deg = sweep_start_deg
            a_sw_deg = sweep_end_deg - c_sw_deg
            c_di_deg = dihedral_start_deg
            a_di_deg = dihedral_end_deg - c_di_deg

            # switch to radians
            c_sw_rad = np.deg2rad(c_sw_deg)
            a_sw_rad = np.deg2rad(a_sw_deg)
            c_di_rad = np.deg2rad(c_di_deg)
            a_di_rad = np.deg2rad(a_di_deg)

            # calculate shift distance
            x_shift = b_loc/a_sw_rad * \
                np.log(np.abs(np.cos(c_sw_rad)/np.cos(a_sw_rad+c_sw_rad)))
            z_shift = -b_loc/a_di_rad * \
                ( np.cos(c_di_rad) - np.cos(a_di_rad + c_di_rad) )
            y_shift = b_loc/a_di_rad * \
                ( np.sin(a_di_rad + c_di_rad) - np.sin(c_di_rad) )
            b_loc = np.sqrt(z_shift**2. + y_shift**2.)
            
            # calculate corresponding constant sweep angle
            L_rad = np.arctan(x_shift / b_loc)
            L_deg = np.rad2deg(L_rad)
            sweep_new = L_deg
            G_rad = np.arctan(-z_shift / y_shift)
            G_deg = np.rad2deg(G_rad)
            dihedral_new = G_deg

        else:
            sweep_new = sweep_old
            dihedral_new = dihedral_old
        
        chord[j  ,1] = chord_old[i  ][1]
        chord[j+1,1] = chord_old[i+1][1]
        
        # create new wing
        section = dict(craftDict["wings"][name])
        section.pop("control_surface")
        section.pop("grid")
        if "twist" in section:
            section.pop("twist")

        # add new semispan sweep, dihedral, chord
        section["semispan"] = b_loc
        section["sweep"] = sweep_new
        section["dihedral"] = dihedral_new
        section["chord"] = np.array(chord_old[i:i+2]).tolist()
        section["chord"][0][0] = 0.0
        section["chord"][1][0] = 1.0
        section["airfoil"] = np.array(airfoil_old[i:i+2]).tolist()
        section["airfoil"][0][0] = 0.0
        section["airfoil"][1][0] = 1.0
        IDnumber += 1
        section["ID"] = IDnumber
        if i == 0:
            if "connect_to" in craftDict['wings'][name]:
                section["connect_to"] = \
                    dict(craftDict['wings'][name]["connect_to"])
        else:
            section["connect_to"] = dict(ID = IDnumber - 1)
        
        # save to dictionary
        craftDict["wings"][name + "_" + str(i)] = section

        # return IDnumber


# make constant
split_wing(craftDict,"main_wing",IDnum)
split_wing(craftDict,"horizontal_tail",len(craftDict["wings"]))

with open("CRM_MUX_discretized.json", "w") as outfile:
    json.dump(craftDict, outfile, indent = 4)

## prismoids file

four = "CRM_prismoids"
sceneDict['scene']['aircraft'][four] = sceneDict['scene']['aircraft'].pop(three)
sceneDict['scene']['aircraft'][four]['file'] = dict(craftDict)

CRM_prismoids = mx.Scene(sceneDict)
# CRM_prismoids.display_wireframe(show_vortices=False)
CRM_prismoids.export_dxf(aircraft=four,export_as_prismoid=True)

