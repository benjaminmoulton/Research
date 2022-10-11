# -*- coding: utf-8 -*-
"""
Created on Fri Sep  4 15:51:16 2020

@author: Dallin
"""

# prac prog for Guide Curve handpicking for Horizon motor smoothing


def prelim():
    alone = []
    behind = []
    over = []
    front = []
    
    for i in range(0,23,1):
        alone.append(i)
    for i in range(23,37,1):
        behind.append(i)
    for i in range(37,69,1):
        over.append(i)
    for i in range(69,89,1):
        front.append(i)
    for i in range(89,200,1):
        alone.append(i)
    
    return(alone,behind,over,front)


def listcom(pop,kid):
    new = []
    while len(kid)> 0:
        if pop[0] < kid[0]:
            new.append(pop[0])
            pop.remove(pop[0])
        else:
            new.append(kid[0])
            kid.remove(kid[0])
    return new


# Input the desired number of hand picked guide curves
# Returns the hand picked guide curves
def handpick(num,alone,behind,over,front):
    # alone,behind,over,front = prelim()
    """ Initial hpicks entries ---------------------------"""
    hpicks = []
    hpicks.append(min(alone))
    hpicks.append(min(behind))
    hpicks.append(max(behind))
    hpicks.append(min(over))
    hpicks.append(max(over))
    hpicks.append(min(front))
    hpicks.append(max(front))
    mid = int( (min(front)+max(front))/2)
    diff = (round(max(alone)/2)) - mid
    hpicks.append( (round(max(alone)/2)) + diff )
    hpicks.append(max(alone))
    #print('\nhpicks:',hpicks)
    """ Leftover determination and verification ----------"""
    percs = [.02,.07,.25,.15,.45,.06]
    #percs = [.18,.18,.18,.18,.18,.10]
    leftover = num - (len(hpicks)-1)   
    #print(leftover)
    lefts = []
    for elem in percs:
        lefts.append(round(leftover*elem))
    #print('\nlefts is:',lefts)
    #print('sum is: ',sum(lefts))
    if (leftover != (sum(lefts))):
        diff = leftover - (sum(lefts))
        lefts[4] += diff
    #print('diff is: ',diff)
    #print('\nlefts is now: ',lefts)
    #print('sum is now: ',sum(lefts))
    """ In between indeces determination -----------------"""
    bump = 0
    adder = []
    for i in range(len(hpicks)-1):
        if (hpicks[i+1]-hpicks[i]) > 1:
            step = (hpicks[i+1]-hpicks[i]) /(lefts[bump]+1)
            for j in range(lefts[bump]):
                adder.append(round(hpicks[i]+(step*(j+1))))
            bump += 1
    #print('\nadder is:',adder)
    #print(len(hpicks)+len(adder))
    hpicks = listcom(hpicks,adder)
    #print('\nhpicks: ',hpicks,'\n',len(hpicks))
    """---------------------------------------------------"""
    return hpicks


# handpick(48)
    
    
    
