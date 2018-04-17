#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 22 21:11:04 2017

@author: justinpringle

morphology calcs
"""

import numpy as np

class morphology():
    
    def __init__(self,x,y,yInit,xPoints,qPoints,swl,yStr,xStr,h0,t0,d50):
        
        self.x = x
        self.y=y
        self.xPoints = xPoints
        self.qPoints= qPoints
        
        ## sed props
        self.rhos = 2.65e3 #(kg/m3)
        self.rho  = 1.03e3 #(kg/m3)
        self.p = 0.4 #porosity
        self.d50 = d50 #units in mm
        self.h0 = h0
        self.t0 = t0
        self.yInit = yInit
        
        self.swl = swl
        self.yStr = yStr
        self.xStr = xStr
        
         # average bottom slope from shoreline to depth of active 
        #LST - default empirical calc
        
    def kamphuis(self,theta,hs,):
        '''
        calcs the kaphuis longshore transport
        '''
        
    def CERC(self,theta,hs,cg,sgn,rot,normals):
        '''
        calcs the longshore transport per CERC eqn.
        these are Q potentials, then checks for structures and seawalls etc.
        theta: angle of breaking rel to local shoreline
        hs: sig wave height
        cg: group velocity
        dhdx: longshore gradient in brekaing wave height
        '''
        ##profile characteristics (taken from empirical relations in Manual)
        L0 = 9.81*self.t0**2/(2*np.pi)
        Dtlo = (2.3-10.9*self.h0)*(self.h0/L0)
        ## subject to interpretation
        ## NOTE units of A are in 1/m**3
        if self.d50<0.4:
            self.A = 0.41*(self.d50)**0.94
        elif self.d50>=0.4 and self.d50<10:
            self.A = 0.23*(self.d50)**0.32
        elif self.d50>=10 and self.d50 <40:
            self.A = 0.23*(self.d50)**0.28
        elif self.d50 >=40:
            self.A = 0.46*(self.d50)**0.11            
            
        self.m = ((self.A**3)/Dtlo)**(1/2)
        
        thetaR = np.abs(np.deg2rad(360-rot-normals)-sgn*theta)#np.deg2rad(theta)
#        for i in thetaR:
#            print(np.rad2deg(i))
        self.K1 = 0.77
        self.K2 = 0.5*self.K1
        
        a1 = self.K1/(16*(self.rhos/self.rho - 1)*(1-self.p)*(1.416**(5/2)))
        a2 = self.K2/(8*(self.rhos/self.rho-1)*(1-self.rho)*self.m*1.416**(7/2))
        
        dhdx = np.zeros(len(hs))
        dhdx[0:-1] = (hs[1::]-hs[0:-1])/(self.x[2::]-self.x[1:-1])
        
        Q = hs**2*cg*(a1*np.sin(2*thetaR)-a2*np.cos(thetaR)*dhdx)
        
        #check swl
        if len(self.swl) >0:
            #do something here
            pass
        
        #check str
        #set Q to zero at structures
        if len(self.xStr)>0:
            for i in range(len(self.yStr)):
#                if self.y[i]<self.yStr[i]:
                Q[self.xStr] = 0
#                else:
#                    continue
        
        return Q
    
    def relativeAngles(self,theta,rot,normals):
        '''
        calculates the relative angles
        '''
        
        return None
    
    def orient(self,theta,normals,rot,Q,sgn):
        '''
        erosion or accretion depending on shoreline.
        '''
        angDiff = 360-rot-normals - sgn*np.rad2deg(theta)
        sgn=[]
        for i in range(len(angDiff)):
            if angDiff[i]>0 or angDiff[i]<0:
                sgn.append(angDiff[i]/np.abs(angDiff[i])) #this keeps sin(theta) positive whn calculating Q
            else:
                sgn.append(1)
        sgn = np.asarray(sgn)
        
        return Q*sgn,angDiff
        
        
    def byPass(self,QPot,hs,rot,D0,perm):
        '''
        calculates bypassing
        rot tells me what 'updrift' side is on
        Dg depth of groyne -> can be array of groynes but must be same length as
        hs
        '''
        Dlt = 1.6*hs
#        for i in Dlt:
#            print i
        yDiff = np.round(self.y,2)
        Dg = []
        normal = 360-rot
        sgn = normal-D0
        if len(self.yStr) >0:
            for i in range(len(self.yStr)):
            
                Dg.append(self.A*self.yStr[i]**(2/float(3)))
                Dlt = 1.6*hs[self.xStr[i]]
                if Dg[i]<Dlt:
                    byp = 1 - Dg[i]/Dlt
                else:
                    byp = 0
                F = perm*(1-byp)+byp
#                print(F,sgn)
                if sgn>0:
#                    print(dyPot[self.xStr[i]-1])
#                    dyPot[self.xStr[i]-1]*byp
#                    print(dyPot[self.xStr[i]])
                    if self.y[self.xStr[i]]<self.yStr[i]:
                        QPot[self.xStr[i]]=QPot[self.xStr[i]-1]*F
                    else:
                        QPot[self.xStr[i]]=QPot[self.xStr[i]-1]
#                        QPot[self.xStr[i]+1]=QPot[self.xStr[i]-1]*F
#                    print(QPot[self.xStr[i]-1],QPot[self.xStr[i]],QPot[self.xStr[i]+1])
#                    print(dyPot[self.xStr[i]-1]*(1-byp))
                else:
#                    dyPot[self.xStr[i]+1]*=byp
                    QPot[self.xStr[i]]=QPot[self.xStr[i]+1]*F
            
                
#        Dg = self.A*self.y**(2/float(3))
        
#        BYP = 1-Dg/Dlt
#        for i in BYP:
#            print i
        #check incoming wave angle if + or -
        
#        print len(BYP)
#        if sgn>0:
#            #updrift is left side
                
#            dy = [dyPot[0]]
#            dy.extend([dyPot[i]*(1-BYP[i])+dyPot[i-1]*BYP[i] for i in range(1,len(dyPot)-1)])
#            dy.append(dyPot[-1])
#            
#        else:
#            #updrift is right side            
#            dy=[dyPot[i]*(1-BYP[i])+dyPot[i+1]*BYP[i] for i in range(len(dyPot)-1)]
#            dy.append(dyPot[-1])
#        print(len(dy))    
        return np.asarray(QPot)

    def checkStr(self,dy,rot,D0):
        '''
        checks if sand level past end of structures.
        '''
        normal = 360-rot
        sgn = normal-D0
        
        if len(self.yStr)>0:
            for i in range(len(self.swl)):
                if sgn>0:
                    if self.y[i]+dy[i] < self.swl[i]:
#                        dy[i+1] = dy[i]
                        dy[i] = 0
#                        print(i)
                elif sgn<0:
                    if self.y[i]+dy[i] < self.swl[i]:
#                        dy[i-1] = dy[i]
                        dy[i] = 0
        return dy
    
    def source_sinks(self,qso = None,qi = None,default = 'gauss'):
        '''
        adds local sources and sinks
        qso and qsi must be arrays qso -> source, qsi->sink
        '''
        if qso == None:
            qso = np.zeros(len(self.x))
        if qi == None or qi == 0:
            qsi = np.zeros(len(self.x))
        
        if qi !=None and qi>0 and default == 'gauss':
            #apply gaussian distributed fill at coords
            X0 = 140*10
            qsi = qi/10*np.exp(-(self.x-X0)**2/(2*(qi)**2))
        self.q = qso-qsi
        
    def boundaries(self,Q,m=0,n=0):
        '''
        user defined boundaries
        Neumann -> gradient boundary, user to define gradient
        m = gradient at start, n = gradient at end -> default zero 
        '''
#        q = [i for i in Q]
#        q.extend(Q)
        q = [Q[0]+m*Q[0]]
#        q=[]
        q.extend(Q)
        q.append(Q[-1]+n*Q[-1])
        
        return np.asarray(q)
        
    def update(self,db,dc,dt,Q,Qt,hs,rot,D0,Dg,perm):
        '''
        updates the shoreline positions based on Q
        does this for each point
        dY is an array to speed up computation
        db,dc are also arrays for each grid point
        Q = longshore transport rates at current time step
        Qt = longshore transport rates and next time step (estimated)
        '''
        dx = self.xPoints[1::] - self.xPoints[0:-1]
        Q = self.byPass(Q,hs,rot,D0,perm)
        Qt = self.byPass(Qt,hs,rot,D0,perm)
        
        Q = self.boundaries(Q,m=0)
        Qt = self.boundaries(Qt,m=0)
        
        dQ0 = Q[1::]-Q[0:-1]
        dQ1 = Qt[1::]-Qt[0:-1]
#        for i in dQ0:
#            print i
#        dQ0 = self.boundaries(dQ0)
#        dQ1 = self.boundaries(dQ1)
        
        dQdX = 0.5*(dQ0/dx+dQ1/dx) # take the time mean
#        for i in dQdX:
#            print i
#        print(dQdX[-1])
        
        #calc the potential change in Y
        #run bypass and alter change in Y accordingly
        dyPot = float(dt)/(db+dc)*(dQdX-self.q)
#        for i in dyPot:
#            print i
        dY = self.checkStr(dyPot,rot,D0)#self.byPass(dyPot,hs,rot,D0)
#        for i in dY:
#            print i
        self.y += dY
#        for i in dY:
#            print i
        ### could change this
#        self.qPoints[0] += dY[0] #effectively Neumann boundary
#        self.qPoints[1::] += dY
        return Q,dY
        
        
        
        
        
        
        