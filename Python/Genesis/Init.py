#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 22 10:23:54 2017

@author: justinpringle

This is the initialization script

Used to numerically model shoreline change due to longshore transport:
    init Grid
    init Shoreline
    init Structures
    init Boundaries
    Waves
    Morphology
    
"""
from grid_and_shoreline import gridGen as grid
from morphology import morphology as morph
from wave import wave as wave
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as anim

def init():
    line.set_data([],[])
    return line,

def animate(i):
    

    ax.set_title('num %d, dir = %.2f'%(i,dStore[i]))
    line.set_data(test.x,s0[i+1])
    
    u = h[i]*np.cos(an[i])
    v = h[i]*np.sin(an[i])
    q.set_UVC(u,v,h[i])
    
    
    
    return line, q, 

def interpParam(time,parm,dtParm):
    '''
    time = global time
    interpolates the parameter at said time.
    Param can have own time step.
    Parm must be an array of len N
    '''
    
#    time = i*dt
    rem = np.remainder(time,dtParm)/float(dtParm)
    ind = int((time-rem)/dtParm)
#    print(rem)
    P0 = parm[ind]
    if ind<len(parm):
        P1 = parm[ind+1]
        value = P0+(P1-P0)*rem
        
    else:
        value = P0
    
    return value

def defaultWave(H0,T0,D0,dt,timeTot):
    '''
    returns wave climate for default model.
    '''
    N = int(timeTot/dt)
#    check
    if N*dt<timeTot:
        N+=1
        
    H = np.asarray([H0 for i in range(N+1)])
    T = np.asarray([T0 for i in range(N+1)])
    D = np.asarray([D0 for i in range(N+1)])
#    D[-3::] = 130
    
    return H,T,D

def defaultFill(Qmax,dt,timeStart,dur,timeTot):
    '''
    adds fill to the shoreline.
    Volume is distributed as Gauss func. along shoreline.
    Qmax is the max flow rate -> assume linear distribution.
    timeStart is the start time of the nourishment
    dt is the time interval of nourishment
    dur is the duration of the nourishment, Qmax at dur/2
    '''
    N = int(timeTot/dt)
    
    if N*dt<timeTot:
        N+=1
        
#    indS = int(timeStart/dt)
    qNourish  = np.zeros(N)
    for i in range(N):
        tm = dt*i
        if tm>= timeStart and tm<=timeStart + float(dur)/2:
            qNourish[i] = (tm-timeStart)/(float(dur)/2)*Qmax
        elif tm>timeStart+float(dur)/2 and tm <= timeStart+dur:
#            print(Qmax - Qmax*(timeStart+float(dur)/2 - tm)/(float(dur)/2))
            qNourish[i] = Qmax - Qmax*(tm-timeStart-float(dur)/2)/(float(dur)/2)
            
    return qNourish
            
    

if __name__ == '__main__':
    
    
    # initialize grid
    # calcs shoreline and shore normals
    smooth = False
    dt = 5 #seconds
    T = 1000*15
    rot = 200
    test = grid('',rot)
    test.defaultGrid(smooth,dX=10,N=200,perturbation=True)
    # read in any seawalls
    test.readSWL('',default=True)
    # read in any structures
    test.readStr('',default=True,on=True)
    ## some storage lists
    s0 = [[i for i in test.y]]
    an = []
    h=[]
    
    H0 = 5
    T0 = 16
    D0 = 100
    dtWave = 0.1*3600 #time step in seconds
    Hs,Tp,D = defaultWave(H0,T0,D0,dtWave,T)
    dStore = []
    
    #nourishment
    dtNourish = 180
    qNour = defaultFill(1,dtNourish,500,10*dtNourish,T+dtNourish)
    
    # initialize time steps
    # loop through and update wave accordingly
    t = 0
    while t<T:
        
        #interpolate waves
        hs = interpParam(t,Hs,dtWave)
        tp = interpParam(t,Tp,dtWave) 
        d = interpParam(t,D,dtWave)
        qN = interpParam(t,qNour,dtNourish)
        dStore.append(d)
#        print(d)
        #use the offshore contour to calculate wave refraction
        wv = wave(hs,tp,d,test.normalsOff)
        #wave angle relative to offshore contour
        wv.orient()
        #refract and shoal waves -> still need to add difraction
        hb = wv.iteration()
        #init morphology 
        bed = morph(test.x,test.y,test.yInit,test.xPoints,test.qPoints,test.swl,
                    test.yStr,test.xStr,hs,tp,1)
        #calculate potential LST
        Q = bed.CERC(wv.thetaB,hb,wv.Cgb,wv.sgn,rot,test.normals)
        #quantify erosion and accretion
        Q,qAn = bed.orient(wv.thetaB,test.normals,rot,Q,wv.sgn)
        #add source and sinks
        bed.source_sinks(qi=qN)
        Dg = bed.A*test.yAdd**(2/float(3))        
        dyPot,dy = bed.update(6,2,dt,Q,Q,hb,rot,d,Dg,0.8)
        
        test.update(bed.y,bed.qPoints,smooth)
        s0.append([i for i in test.y])
        an.append([np.deg2rad(270-test.rot)-np.deg2rad(test.normalsOff[i])+
                   wv.sgn[i]*wv.thetaB[i] for i in range(len(wv.thetaB))])#([np.deg2rad(270)-i for i in wv.sgn*wv.thetaB])
        h.append([i for i in hb])
        
        t+=dt
    
    
    fig,ax = plt.subplots()
    
    line, = ax.plot([],[])
    ax.plot(test.x,s0[0])
    ax.plot(test.x,s0[-1],alpha=0.5)
    ax.plot(test.x,test.yOff+np.max(s0[0]),alpha=0.5,linestyle='--')
    ax.vlines(test.xPoints[1:-1],0,50,alpha=0.25,linewidth=0.5)
    q = ax.quiver(test.xPoints[1:-1],[50 for i in range(len(h[0]))],
                  np.asarray(h[0])*np.sin(np.asarray(an[0])),
                  np.asarray(h[0])*np.cos(np.asarray(an[0])),width=0.005)
#    ax.plot(test.)
    aniM = anim.FuncAnimation(fig,animate,frames=len(s0)-1,blit=False,init_func=init)
    plt.show()    

    
    
