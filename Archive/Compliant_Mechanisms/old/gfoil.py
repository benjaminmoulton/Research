#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 20 16:12:54 2018

@author: benjamin
This creates a gfoil

The user may edit inputs in the getParams() function such that the .txt file can be located 
to input point values (as from a seligdat file).
Else, parameters can be inserted for NACA 4 digit airfoils.
"""

def getParams():#               
    # stm # what vars will modify which foilshape output: self(rootairfoil etc), .txt, .mu
    
      #     # root foil shape # SHOULD be a naca four digit
    rootAirfoil = "NACA 0016"
    
      #     # tip foil shape  # SHOULD be a naca four digit
    tipAirfoil = "NACA 0016"  
    
      ###   # offset value in inches of the foam from the XY axis (z value) 
            # as well as how tall the piece of foam is (y height)( for use in the movements of the HotwireCNC) thickness
            # and the length of the foam peice (x height) ( from Trailing Edge to Leading Edge) chord 
            # foam must be positioned 1 inch out from the wire.
    offset = 12.
    depth = 2
    length = 10
    foamPositionOffset = [offset,depth,length]
    
      ###   # length of the wire from the attachment at the spring to the opposite eyescrew
    wireLength = 31.5 
         
      ##    # This is the chord length of the root foil
    rootChord = 1. 
           
      ##    # This is the chord length of the tip foil
    tipChord = 1. 
            
      ##    # This is the degree value of the rotation of the tipfoil downward, or CW on a right wing.
    washout = 0.  
    
      ##    # This is the degree value, where + is CW(right wing) toward the tail.             
    sweep = 1. 
                
      ##    # semispan of wing (ie. left or right only)
    span = 7  
               
      ###   # whether the values are to be in inches. the other option is milimeters
    isInInches = True  
       
      ##    # Whether the wing to cut is the left wing or not(right)
    isLeftWing = False  
      
      ###     # only to be modified if user wants the wire to cut differently than normal
    rootWireKerf = 0.01
      
      ###     # only to be modified if user wants the wire to cut differently than normal
    tipWireKerf = 0.01         
    
    # offset for where the foilshape is? u axis is along the baseshape. determine how to find where the foam begins/ends
    
      ###    # whether the airfoil desired will be a fishBac design !!!!!! only for symmetrical airfoils !!!!!!
    isFishbone = False
    
      ###   # whether the fishbone design will be a default setting, ie. (startperc = location of max thickness,endperc = 0.85,
              # boneWidth = 0.0225,gapWidth = 0.0225,spineWidth = 0.25,isTlocal = False) or not (params to be used are those following)
    isDefaultFishbone = False 
    
      ###    # if Fishbone, at what % chord should the fishbone start from the Leading Edge
    startperc = .3 
    
      ###    # if Fishbone, at what % chord should the fishbone end from the Leading Edge
    endperc = .85  
    
      ###   # if Fishbone, how wide should the bones be in % chord       
    boneWidth = 0.035 
    
      ###   # if Fishbone, how wide should the gaps between bones be in % chord     
    gapWidth = 0.035  
    
      ###   # if Fishbone, how wide should the spine be in % thickness     
    spineWidth = .5 
    
      ###   # if Fishbone, should the Thickness be based on the Max thickness, or the local thickness at each point    
    isTlocalFish = True 
          
    FBVars = [isDefaultFishbone,startperc,endperc,boneWidth,gapWidth,spineWidth,isTlocalFish]
    
    isOscispine = False
    
    isDefaultOscispine = False
    
    startPerc = .35
    
    endPerc = .85
    
    period = 0.1
    
    spineWidth = 0.2
    
    firstUp = True
    
    isTlocalOsc = False
    
    OSVars = [isOscispine,isDefaultOscispine,startPerc,endPerc,period,spineWidth,firstUp,isTlocalOsc]
    
    # selig datfile location if data points are to be used instead. for tip
    fileLocationTip = 'none'  
    
    # selig datfile location if data points are to be used instead. for root
    fileLocationRoot = 'none'   
    # /home/benjamin/Desktop/foilshapes/NACA 64A210.txt
    # /home/benjamin/Desktop/foilshapes/wing_one.mu
    # /home/aerolab/Desktop/naca.txt
    # /home/benjamin/Desktop/oscispine.txt
    
    return rootAirfoil,tipAirfoil,foamPositionOffset,wireLength,rootChord,tipChord,washout,sweep,span,isInInches,isFishbone,FBVars,OSVars,isLeftWing,rootWireKerf,tipWireKerf,fileLocationTip,fileLocationRoot










def gFoil(): 
    
    # files are imported for use in the program ::::
    
    
    # gcode file consists of functions used in the manipulation of files
    #  and file formats.
    
    import gcode as gc                   # seligdatfile.txt
    
    # NACApoints file consists of functions used in the creation, modification,
    #  and use of data points for naca four digit airfoils
    
    import NACApoints as nc
    
    # geometry file consists of functions used in creating the geometry specified
    # by the user in the getParams() function
    
    import geometry as geo
    
    # the pyplot class is imported for use in plotting airfoil data to show the user
    #  how foilshapes are constructed after input.
    
    from matplotlib import pyplot as plt
    
    # 
    rootAirfoil,tipAirfoil,offset,wireLength,rootChord,tipChord,washout,sweep,span,isInInches,isFishbone,fishvars,osvars,isLeftWing,rootWireKerf,tipWireKerf,fileLocationTip,fileLocationRoot = getParams()
    
    fileType = ''
    
    if(fileLocationRoot[-4:] == '.txt' and fileLocationTip[-4:] == '.txt'): # != 'none'
        fileType = '.txt file of data points'
        x,y = gc.NACAFileReader(fileLocationRoot)
        u,v = gc.NACAFileReader(fileLocationTip)
        x,y = nc.enlarge(x,y,rootChord)
        u,v = nc.enlarge(u,v,tipChord)
        isDefault = fishvars[0]
        if(osvars[0] == False):
            if(isDefault):
                fx,fy = nc.FishBonePoints(x,y,geo.getThickness(x,y),geo.getChord(x))
                fu,fv = nc.FishBonePoints(u,v,geo.getThickness(u,v),geo.getChord(u))
            else:
                start = fishvars[1]
                end = fishvars[2]
                bone = fishvars[3]
                gap= fishvars[4]
                spine = fishvars[5]
                islocal = fishvars[6]
                fx,fy = nc.FishBonePoints(x,y,geo.getThickness(x,y),geo.getChord(x),start,end,bone,gap,spine,islocal)
                fu,fv = nc.FishBonePoints(u,v,geo.getThickness(u,v),geo.getChord(u),start,end,bone,gap,spine,islocal)
        else:
            if(osvars[1] == True):
                ox,oy = nc.OsciSpine(x,y,geo.getChord(x),geo.getThickness(x,y))
                ou,ov = nc.OsciSpine(u,v,geo.getChord(u),geo.getThickness(u,v))
            else:
                start = osvars[2]
                end = osvars[3]
                period = osvars[4]
                spineWidth = osvars[5]
                firstUp = osvars[6]
                isTlocalOsc = osvars[7]
                ox,oy = nc.OsciSpine(x,y,geo.getChord(x),geo.getThickness(x,y),start,end,period,spineWidth,firstUp,isTlocalOsc)
                ou,ov = nc.OsciSpine(u,v,geo.getChord(u),geo.getThickness(u,v),start,end,period,spineWidth,firstUp,isTlocalOsc)
    
    elif(fileLocationRoot[-3:] == '.mu' and fileLocationTip[-3:] == '.mu'): # != 'none'
        fileType = '.mu file'
        isLeft_mu,span_mu,sweep_mu,dihedral,mountingAngle,washout_mu,rchord_mu,tchord_mu,rfoil_mu,tfoil_mu = gc.MachUpFileReader(fileLocationRoot)
        # we won't use dihedral or mounting angle. Those parameters can be physically modified by the user after the wing is cut
        x,y,fx,fy,ox,oy = nc.infoCompile(rfoil_mu,rchord_mu,isFishbone,fishvars,osvars)
        u,v,fu,fv,ou,ov = nc.infoCompile(tfoil_mu,tchord_mu,isFishbone,fishvars,osvars)
        isLeftWing = isLeft_mu
        washout = washout_mu
        sweep = sweep_mu
        span = span_mu
        
    else: # reading in input values  # if(fileLocationRoot == 'none' and fileLocationTip == 'none'):
        fileType = '%s root airfoil and %s tip airfoil' % (rootAirfoil,tipAirfoil)
        x,y,fx,fy,ox,oy = nc.infoCompile(rootAirfoil,rootChord,isFishbone,fishvars,osvars)
        u,v,fu,fv,ou,ov = nc.infoCompile(tipAirfoil,tipChord,isFishbone,fishvars,osvars)

    
    



    #aligns outer quarter chord with root quarter chord
    if(isFishbone):
        x = fx
        y = fy
        u = fu
        v = fv
    if(osvars[0]):
        x = ox
        y = oy
        u = ou
        v = ov
    u = geo.alignQuarterChord(x,u)
    
    if(washout != 0):
        u,v = geo.giveWashout(u,v,washout)
    if(sweep != 0):
        u = geo.giveSweep(u,sweep,span)
    if(isLeftWing == True): # this might need to be swapped becuase of the way the hotwireCutter positions 4 axes
        tempx = x
        tempy = y
        x = u; y = v; u = tempx; v = tempy;
    
    if(rootWireKerf != 0 and tipWireKerf != 0):
        kx,ky = nc.kerfMaker(x,y,rootWireKerf)
        ku,kv = nc.kerfMaker(u,v,tipWireKerf)
    if(rootWireKerf == 0 or tipWireKerf == 0):
        kx,ky = x,y
        ku,kv = u,v
    
    wx,wy,wu,wv = geo.fixAxisPoints(kx,ky,ku,kv,span,offset[0],wireLength)
    
    
    plt.subplot(211)
    plt.axis('equal')
    plt.title("root foil")
    plt.plot(x,y)
    if(rootWireKerf != 0 and tipWireKerf != 0):
        plt.plot(kx,ky)
    plt.show()
    plt.subplot(212)
    plt.title("tip foil")
    plt.axis('equal')
    plt.plot(u,v)
    if(rootWireKerf != 0 and tipWireKerf != 0):
        plt.plot(ku,kv)
    plt.show()
    
    from mpl_toolkits.mplot3d import Axes3D
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    face = 'y'
    ax.plot(wx,wy,0,zdir= face)
    ax.plot(x,y,offset[0],zdir= face)
    ax.plot(u,v,offset[0]+span,zdir= face)
    ax.plot(wu,wv,wireLength,zdir= face)
    ax.view_init(20,-30)
    for i in range(0,len(wx),4):
        ax.plot([wx[i],wu[i]],[wy[i],wv[i]],[0,wireLength],zdir= face)
    
    wy,wv = geo.shiftUp(wy,wv,offset[1])
    wx,wu = geo.shiftRight(wx,wu)
    
    
    wx = geo.roundValuesG(wx)
    wy = geo.roundValuesG(wy)
    wu = geo.roundValuesG(wu)
    wv = geo.roundValuesG(wv)
    
    
    
    
    #insert into the first values repositioning so that the machine can cut the desired shape
    
    
    fishTEXT = ''
    if(isFishbone):
        fishTEXT = 'Fishbone VCCW,'
    
    osciTEXT = ''
    if(osvars[0] == True):
        osciTEXT = 'Oscillatory Spine VCCW,'
    
    leftTEXT = 'right wing,'
    if(isLeftWing):
        leftTEXT = 'left wing,'
    
    inch = 'in'
    if(not isInInches):
        inch = 'mm'
    
    fileExplanation = """(This file is a G-Code generated wing shape for use in the Aerolab CNCHotwire)
(Wing shape is based off a %s)
(Wing Characteristics::)\n(%s %s %s root chord %.1f %s, tip chord %.1f %s,)
(%.1f degrees washout, %.1f degrees sweep,span %.1f %s)\n""" % (fileType,fishTEXT, osciTEXT,
    leftTEXT,rootChord,inch,tipChord,inch,washout,sweep,span,inch)                                                     
#    print(fileExplanation)
    
    
    gc.taperedFileWriter(wx,wy,wu,wv,isInInches,fileExplanation,offset[1],offset[2])
    return 0


gFoil()




















