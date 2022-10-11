import json
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
import machupX as mx
import numpy as np

## constants
#####################################################################
#####################################################################
#####################################################################

## span frac locations
s_bay = 0.07560964056743188
s_taper = 0.1613005665438547
s_0 = 0.42047294 #0.32904045
s_1 = 0.65427572 #0.49678034
s_2 = 0.83982266 #0.66452023
s_3 = 0.95895113 #0.83226011

## optional quick reference keyword argument inputs
aircraft = {'aircraft': 'Horizon'}
sacs = {'control_state': {'sym': None, 'asym':None}}
sas  = {'state': {  "velocity": None,
                    "alpha": None,
                    # "beta": None,
                    "angular_rates": [None]*3,
                    "angular_rate_frame": "stab"}}
forcesOptions = {'body_frame':False, 'stab_frame':True, 'wind_frame':False,
                    'dimensional':False}#, 'verbose':True}

## readin in json files
f = open('HorizonAircraft_BiggerWinglet.json', 'r')
horizonDict = json.load(f)
f.close()
f = open('HorizonScene.json', 'r')
sceneDict = json.load(f)
f.close()

## extract needed data
Sw = horizonDict['reference']['area']
bw = horizonDict['reference']['lateral_length']
cbar = horizonDict['reference']['longitudinal_length']

## useful functions
#####################################################################
#####################################################################
#####################################################################

def interpolate(x1,x2,x,y1,y2):
    return (y2-y1)/(x2-x1)*(x-x1)+y1

def isIterable(val):
    try:
        iter(val)
        return True
    except:
        return False

def chord(s):
    if s < 0. or s > 1.: raise ValueError('Unexpected value for span fraction')
    if s <= s_bay: return 2.75      ## bay chord is 2.75 ft
    if s >= s_taper: return interpolate(s_taper, 1., s,
                                1.1458333333333333, 0.9166666666666666)
    return -3.16641665e+00 + s * (1.85435652e+02 + s * (-1.80033656e+03 +
                                                        5.06167082e+03 * s))
def Chord(S):
    if isIterable(S):
        return np.array([chord(s) for s in S])
    else:
        return chord(S)

def generateControlFunctions(sd, ad):
    def sym(s, symDefl=sd):
        if s < 0. or s > 1.: raise ValueError('Unexpected value for span fraction')
        if s < s_bay: return symDefl[0]
        if s < s_taper: return 0.
        if s < s_0: return symDefl[1]
        if s < s_1: return symDefl[2]
        if s < s_2: return symDefl[3]
        if s < s_3: return symDefl[4]
        return symDefl[5]
    def asym(s, asymDefl=ad):
        if s < 0. or s > 1.: raise ValueError('Unexpected value for span fraction')
        if s < s_taper: return 0.
        if s < s_0: return asymDefl[0]
        if s < s_1: return asymDefl[1]
        if s < s_2: return asymDefl[2]
        if s < s_3: return asymDefl[3]
        return asymDefl[4]
    def Sym(S):
        if isIterable(S):
            return np.array([sym(s) for s in S])
        else:
            return sym(S)
    def Asym(S):
        if isIterable(S):
            return np.array([asym(s) for s in S])
        else:
            return asym(S)
    return Sym, Asym

def updateControls(symDefl, asymDefl, sc):
    ## generate control functions and update scene
    Sym, Asym = generateControlFunctions(symDefl, asymDefl)
    sacs['control_state'][ 'sym'] =  Sym
    sacs['control_state']['asym'] = Asym
    sc.set_aircraft_control_state(**sacs,**aircraft)

def updateState(V, aoa, beta, omega, sc):
    ## update state
    sas['state']['velocity'] = V
    sas['state']['alpha'] = aoa
    sas['state']['beta'] = beta
    sas['state']['angular_rates'] = omega
    sc.set_aircraft_state(**sas,**aircraft)

## update the input json files with necessary items
#####################################################################
#####################################################################
#####################################################################

## chord
horizonDict['wings']['wing']['chord'] = Chord

## horizon file
sceneDict['scene']['aircraft']['Horizon']['file'] = horizonDict

## control surfaces
Sym, Asym = generateControlFunctions([0.]*6, [0.]*5)
sceneDict['scene']['aircraft']['Horizon']['control_state'][ 'sym'] =  Sym
sceneDict['scene']['aircraft']['Horizon']['control_state']['asym'] = Asym

# create constant sweep sections which meet up with horizon other sections
def make_constant(horizonDict,name):
    # horizonDict['wings'][name + "2"] = dict(horizonDict['wings'][name])
    bw_2 = horizonDict['wings'][name]["semispan"]
    sweep_old = horizonDict['wings'][name]["sweep"]
    dihedral_old = horizonDict['wings'][name]["dihedral"]
    if name == "wing":
        spans = [0.0,s_bay,s_taper,s_0,s_1,s_2,s_3,1.0]
    else:
        spans = [0.0,1.0]
    sections = len(spans) - 1
    sweep = np.zeros((2*sections,2))
    dihedral = np.zeros((2*sections,2))
    bw_2_new = 0.0
    for i in range(sections):
        j = int(i*2)
        dihedral[j  ,0]  = sweep[j  ,0] = spans[i  ]
        dihedral[j+1,0]  = sweep[j+1,0] = spans[i+1]

        # determine local span
        b_loc = bw_2 * (spans[i+1] - spans[i])
        
        # determine sweep angles
        sweep_diff  = sweep_old[1][1]-sweep_old[0][1]
        sweep_start_deg = spans[i  ] * sweep_diff + sweep_old[0][1]
        sweep_end_deg   = spans[i+1] * sweep_diff + sweep_old[0][1]
        dihedral_diff  = dihedral_old[1][1]-dihedral_old[0][1]
        dihedral_start_deg = spans[i  ] * dihedral_diff + dihedral_old[0][1]
        dihedral_end_deg   = spans[i+1] * dihedral_diff + dihedral_old[0][1]
        
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
        if name != "BiggerWinglet":
            z_shift = -b_loc/a_di_rad * \
                ( np.cos(c_di_rad) - np.cos(a_di_rad + c_di_rad) )
            y_shift = b_loc/a_di_rad * \
                ( np.sin(a_di_rad + c_di_rad) - np.sin(c_di_rad) )
            bw_2_new += np.sqrt(z_shift**2. + y_shift**2.)
        
        # calculate corresponding constant sweep angle
        L_rad = np.arctan(x_shift / b_loc)
        L_deg = np.rad2deg(L_rad)
        sweep[j  ,1] = sweep[j+1,1] = L_deg
        if name != "BiggerWinglet":
            G_rad = np.arctan(-z_shift / y_shift)
            G_deg = np.rad2deg(G_rad)
            dihedral[j  ,1] = dihedral[j+1,1] = G_deg

    horizonDict['wings'][name]["sweep"] = sweep.tolist()
    if name != "BiggerWinglet":
        horizonDict['wings'][name]["dihedral"] = dihedral.tolist()
        horizonDict['wings'][name]["semispan"] = bw_2_new

# make constant
make_constant(horizonDict,"wing")
make_constant(horizonDict,"BiggerWingletBase")
make_constant(horizonDict,"BiggerWinglet")

# chord array
spans = [0.0,s_bay,s_taper,s_0,s_1,s_2,s_3,1.0]
names = ["bay","fan","s0","s1","s2","s3","s4"]
chords = np.zeros((len(spans),2))
for i in range(len(spans)-1):
    chords[i,0] = spans[i]
    chords[i+1,0] = spans[i+1]
    chords[i  ,1] = chord(spans[i  ])
    chords[i+1,1] = chord(spans[i+1])

horizonDict['wings']["wing"]["chord"] = chords.tolist()

# split main wing into sections
horizonDict["wings"]["wing"].pop("grid")
horizonDict["wings"]["wing"].pop("control_surface")
bwb = horizonDict["wings"].pop("BiggerWingletBase")
bwb["ID"] = 8; bwb["connect_to"] = dict(ID = bwb["ID"] - 1)
bw  = horizonDict["wings"].pop("BiggerWinglet")
bw["ID"] = 9; bw["connect_to"] = dict(ID = bw["ID"] - 1)

mainwing = horizonDict["wings"].pop("wing")
for i in range(len(spans)-1):
    section = dict(mainwing)
    section["ID"] += i
    section["semispan"] = mainwing["semispan"] * (spans[i+1] - spans[i])
    j = 2*i
    section["sweep"] = mainwing["sweep"][j][1]
    section["dihedral"] = mainwing["dihedral"][j][1]
    section["chord"] = np.array(mainwing["chord"][i:i+2]).tolist()
    section["chord"][0][0] = 0.0
    section["chord"][1][0] = 1.0
    section["airfoil"] = "NACA_2412_fcf56_linearAlpha"
    section["connect_to"] = dict(ID = section["ID"] - 1)


    # save section
    horizonDict["wings"][names[i]] = section

horizonDict["wings"]["BiggerWingletBase"] = bwb
horizonDict["wings"]["BiggerWinglet"] = bw


with open("disc_Horizon.json", "w") as outfile:
    json.dump(horizonDict, outfile, indent = 4)

## MUX
#####################################################################
#####################################################################
#####################################################################
## create scene
scene = mx.Scene(sceneDict)

scene.export_dxf(export_as_prismoid=True)

# scene.display_wireframe(show_vortices=False, show_legend=False)


# updateState(54., 20., 0.,[0]*3, scene)
# updateControls([-5]*6, [0]*5, scene)
# fm = scene.solve_forces(**forcesOptions)['Horizon']['total']
# print(fm)

'''
                Setting the Control Surfaces
#####################################################################
#####################################################################
#####################################################################

Currently the control surfaces are all set to 0 degrees deflection. To set
the control surfaces, you can use the function 'updateControls'

updateControls( [symBay,  sym0,  sym1,  sym2,  sym3,  sym4],
                [        asym0, asym1, asym2, asym3, asym4],
                scene )

where 'sym' represents the elevator deflection of that servo and 'asym'
represents the aileron deflection of that servo. Note, that

abs(symI) + abs(asymI) <= 20.     for I = 0 to 4

This represents the deflection limits for flap deflection in the airfoil
database


                Setting the Aircraft State
#####################################################################
#####################################################################
#####################################################################

Currently the state is set to:
    velocity = 54.     ft/sec
    alpha = 0.
    beta = 0.
    omega = [0., 0., 0.]      [p, q, r]

To set the aircraft state you can use the function 'updateState'

updateState(V, aoa, beta, omega, scene)

'''

