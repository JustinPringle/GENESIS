#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 25 11:32:02 2017

@author: justinpringle
"""

import numpy as np
#from grid_and_shoreline import gridGen


class wave():
    '''
    wave module of GENESIS.
    includes refraction, diffraction, shoaling
    Child of Grid
    '''
    
    def __init__(self,H0,T0,D0,shoreNormals):
        '''
        bind deep water conditions
        '''
        self.H0 = H0
        self.T0 = T0
        self.D0 = D0
        self.normals = shoreNormals
        self.angles = None
        
    def orient(self):
        '''
        calcs the angle between incoming wave angle and shore normal.
        Must take note of what is pos and what is neg direction
        
        '''
        
        self.angDiff = self.normals-self.D0 #- self.normals
        self.sgn=[]
        for i in range(len(self.angDiff)):
            if self.angDiff[i]>0 or self.angDiff[i]<0:
                self.sgn.append(self.angDiff[i]/np.abs(self.angDiff[i])) #this keeps sin(theta) positive whn calculating Q
            else:
                self.sgn.append(1)
        self.sgn = np.asarray(self.sgn)
    
    def estL(self,db):
        '''
        see Reeve et al 2004
        '''
        
        y = 4*np.pi**2/(9.81*self.T0**2)*db
        L = self.T0*np.sqrt(9.81*db*(y+(1+0.6522*y+0.4622*y**2+0.0864*y**4+0.0675*y**5)**-1)**-1)
        
        return L
    
    def snell(self,L):
        '''
        calculates angle at breaking.
        '''
        
        self.thetaB = np.arcsin(L/self.L0*np.sin(self.angles))
        
    def iteration(self):
        '''
        iterates to find wave angle, height and length at breaking
        '''
        #Init wave conditions
        angles = np.deg2rad(np.abs(self.angDiff))
        #mask to prevent >90 angle diff
        self.angles = angles#np.ma.masked_where(angles==np.isnan,angles)
        
        
        self.L0 = 9.81*self.T0**2/(2*np.pi)
        self.C0 = self.L0/self.T0
        self.Cg0 = 0.5*self.C0
        
        self.gamma = 0.78
        db = self.H0/self.gamma
        L = self.estL(db)
#        print(L)
        
        self.updateCg(db,L)
        self.snell(L)
        
        self.refraction()
        self.shoaling()
        
        HbOld = self.H0
        HbNew = self.Kr*self.Ks*self.H0
#        print(HbOld)
#        print(self.Kr)

        err = np.abs(HbNew-HbOld)/HbOld*100
#        for i in range(len(angles)):
#            print(np.rad2deg(angles[i]),self.normals[i])
#        print('start')

        
        while np.max(err)>1:
#        for i in range(2):
            db = HbNew/self.gamma
            L = self.estL(db)
#            print(db)
            self.updateCg(db,L)
            self.snell(L)
        
            self.refraction()
            self.shoaling()
            
            HbOld = HbNew
            HbNew = self.Kr*self.Ks*self.H0
            err = np.abs(HbNew-HbOld)/HbOld*100
#            for i in range(len(db)):
#                print(self.Ks[i],self.Cg0,self.Cgb[i])
#            print('new')
#            print(err)
#            print(HbNew)
            
#            print(HbNew)
#        HbNew = self.boundaries(HbNew)
#        self.thetaB = self.boundaries(self.thetaB)
#        self.Cgb = self.boundaries(self.Cgb)
#        self.sgn = self.boundaries(self.sgn)
        return HbNew
        
    def updateCg(self,d,L):
        
        n = 0.5*(1+(2*np.pi*d/L)/np.sinh(2*np.pi*d/L))
#        print(n)
        C = L/self.T0
        self.Cgb = n*C
        
    
    def refraction(self):
        '''
        calculates the refraction coefficient along wave rays at each grid point
        '''
#        print('ref')
#        for i in self.angles:
#            print(i)
        self.Kr = np.sqrt(np.abs(np.cos(self.angles))/np.abs(np.cos(self.thetaB)))       
        
    def diffraction(self):
        '''
        calculates the diffraction around structures
        '''
        
        self.Kd = 1
        
    def shoaling(self):
        '''
        calculates wave shoaling.
        '''        
        self.Ks = np.sqrt(self.Cg0/self.Cgb)
        
    def boundaries(self,H,m=0,n=0):
        '''
        user defined boundaries
        Neumann -> gradient boundary, user to define gradient
        m = gradient at start, n = gradient at end -> default zero 
        '''
        h =[]
        h.append(H[0]+m*H[0])
        h.extend(H)
        
        return np.asarray(h)
        
        
if __name__ == '__main__':
    
    
    test = wave('',200)
    test.bind(2,16,160)
    test.defaultGrid()
    test.orient()
    
    ho,hb=test.iteration()
        
        
        
        
        
        
        
        