



#the function below takes a set of coordinates for a non-cambered wing and determines the 
# coordinates for a fishBAC design on the wing through cubic splines, and returns the 
# new coordinates to the user.
def FishboneMaker(x,y,xbone,ybone,x1,x2,y1,y2):
    fishx = []
    fishy = []
    nextBoneIndex = 0

    i = 0
    while(i < len(x1)):      
        if((xbone[nextBoneIndex] >= x1[i]) and nextBoneIndex < int(len(xbone)/2)):
            for h in range(4):
                fishx.append(xbone[nextBoneIndex])
                fishy.append(ybone[nextBoneIndex])
                nextBoneIndex += 1
        elif((nextBoneIndex != 0 and xbone[nextBoneIndex-1] > x1[i]) or nextBoneIndex == 0): # (xbone[nextBoneIndex] < x1[i]):
            fishx.append(x1[i])
            fishy.append(y1[i]) 
            i += 1
        else:
            i += 1
    turnpt = nextBoneIndex
    j = 0    
    while(j < len(x2)):
        if(nextBoneIndex < len(ybone) and (xbone[nextBoneIndex] <= x2[j])):
             for h in range(4):
                 fishx.append(xbone[nextBoneIndex])
                 fishy.append(ybone[nextBoneIndex])
                 nextBoneIndex += 1
        elif((xbone[nextBoneIndex-1] < x2[j]) or nextBoneIndex == turnpt): # (xbottombones[nextBoneIndex] > x2[j]):
            fishx.append(x2[j])
            fishy.append(y2[j])
            j += 1
        else:
            j += 1
    
    return fishx,fishy


#the function below gets the points to be used by the FishboneMaker
# as the locations for the fishbone points
# uses the FishboneMaker function to create the fish bone, and return the x and y values
# to the user
def FishBonePoints(x,y,thickness,chord,startperc = 0,endperc = 0.85,boneWidth = 0.0225,gapWidth = 0.0225,spineWidth = 0.25,isTlocal = False):
    from scipy.interpolate import interp1d as interp
    if(startperc == 0):
        startperc = x[y.index(max(y))] / chord
    xpost = endperc*chord
    topbonex = []
    while(xpost >= ((startperc + gapWidth + boneWidth)*chord)):
        topbonex.append(xpost)
        topbonex.append(xpost)
        xpost -= gapWidth*chord
        topbonex.append(xpost)
        topbonex.append(xpost)
        xpost -= boneWidth*chord
    
    bottombonex = list(reversed(topbonex))
    xone = x[:len(x)//2]
    xtwo = x[len(x)//2:]
    yone = y[:len(x)//2]
    ytwo = y[len(x)//2:]
    xonef = xone[::-1]
    yonef = yone[::-1]
    fone = interp(xonef,yonef)
    ftwo = interp(xtwo,ytwo)
    topboney = []
    bottomboney = []
    if(isTlocal == True):
        for g in range(0,len(topbonex),4):
            topboney.append(float(fone(topbonex[g])))
            topboney.append(spineWidth*float(fone(topbonex[g+1])) )#/ 2)
            topboney.append(spineWidth*float(fone(topbonex[g+2])) )#/ 2)
            topboney.append(float(fone(topbonex[g+3])))
            bottomboney.append(float(ftwo(bottombonex[g])))
            bottomboney.append(spineWidth*float(ftwo(bottombonex[g+1])) )#/ 2)
            bottomboney.append(spineWidth*float(ftwo(bottombonex[g+2])) )#/ 2)
            bottomboney.append(float(ftwo(bottombonex[g+3])))
    else:
        for g in range(0,len(topbonex),4):
            topboney.append(float(fone(topbonex[g])))
            topboney.append(spineWidth*thickness / 2.0)
            topboney.append(spineWidth*thickness / 2.0)
            topboney.append(float(fone(topbonex[g+3])))
            bottomboney.append(float(ftwo(bottombonex[g])))
            bottomboney.append(-spineWidth*thickness / 2.0)
            bottomboney.append(-spineWidth*thickness / 2.0)
            bottomboney.append(float(ftwo(bottombonex[g+3])))
    bonex = []
    bonex.extend(topbonex)
    bonex.extend(bottombonex)
    boney = []
    boney.extend(topboney)
    boney.extend(bottomboney)
    fx,fy = FishboneMaker(x,y,bonex,boney,xone,xtwo,yone,ytwo)
    
    return fx,fy






