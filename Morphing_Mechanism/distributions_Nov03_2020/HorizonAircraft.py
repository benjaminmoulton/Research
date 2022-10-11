import json
import numpy as np
# import warnings
# warnings.simplefilter(action='ignore', category=FutureWarning)
import machupX as mx

def interpolate(x1,x2,x,y1,y2):
    return (y2-y1)/(x2-x1)*(x-x1)+y1

def isIterable(val):
    try:
        iter(val)
        return True
    except:
        return False

f = open('HorizonAircraft.json', 'r')
horizonDict = json.load(f)
f.close()
f = open('HorizonScene.json', 'r')
sceneDict = json.load(f)
f.close()

## update the chord function in the main wing of the aircraft
#####################################################################
def chord(s):
    if s < 0. or s > 1.: raise ValueError('Unexpected value for span fraction')
    
    s_bay = 0.07560964056743188
    c_bay = 2.75                ## ft
    if s <= s_bay: return c_bay
    
    s_taper = 0.1613005665438547
    
    c_tip = 0.9166666666666666      ## 11. / 12.
    c_taper = 1.1458333333333333     ## c_tip / 0.8
    
    if s >= s_taper: return interpolate(s_taper, 1., s, c_taper, c_tip)
    
    A, B, C, D = 5.06167082e+03, -1.80033656e+03,  1.85435652e+02, -3.16641665e+00
    
    return A * s**3. + B * s**2. + C * s + D

def Chord(S):
    if isIterable(S):
        return np.array([chord(s) for s in S])
    else:
        return chord(S)

horizonDict['wings']['wing']['chord'] = Chord

## update scene with horizon dict
sceneDict['scene']['aircraft']['Horizon']['file'] = horizonDict

## update control functions in the control state of the scene
######################################################################
s_bay = 0.07560964056743188
s_taper = 0.1613005665438547
s_0 = 0.42047294 #0.32904045
s_1 = 0.65427572 #0.49678034
s_2 = 0.83982266 #0.66452023
s_3 = 0.95895113 #0.83226011
def sym(s, symCenter=0., symDefl=[0,0,0,0,0]):
    if s < 0. or s > 1.: raise ValueError('Unexpected value for span fraction')
    sym0, sym1, sym2, sym3, sym4 = symDefl
    if s < s_bay: return symCenter
    if s < s_taper: return 0.
    if s < s_0: return sym0
    if s < s_1: return sym1
    if s < s_2: return sym2
    if s < s_3: return sym3
    return sym4
def asym(s, asymDefl=[0,0,0,0,0]):
    if s < 0. or s > 1.: raise ValueError('Unexpected value for span fraction')
    asym0, asym1, asym2, asym3, asym4 = asymDefl
    if s < s_taper: return 0.
    if s < s_0: return asym0
    if s < s_1: return asym1
    if s < s_2: return asym2
    if s < s_3: return asym3
    return asym4
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
sceneDict['scene']['aircraft']['Horizon']['control_state'][ 'sym'] =  Sym
sceneDict['scene']['aircraft']['Horizon']['control_state']['asym'] = Asym

## MUX
######################################################################


scene = mx.Scene(sceneDict)

