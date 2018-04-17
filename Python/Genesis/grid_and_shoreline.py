#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 22 10:26:03 2017

@author: justinpringle
"""

import numpy as np

class gridGen():
    '''
    this class initializes a grid based user settings to be stored in a file.
    If no file is given then default shoreline with default spacing is 
    initialized.
    
    if N shoreline (X,Y) points then N+1 Q points
    '''
    def __init__(self,directory,rot):
        self.x = None
        self.y = None
        self.xPoints = None
        self.qPoints = None
        self.Path = directory
        self.rot = rot
        self.smoothY = None
        self.yFill = None
    
    def defaultGrid(self,smooth,dX = 5, N = 100,perturbation=True):
        '''
        sets the default grid with triangular/gaussian perdubation due to a fill.
        
        NB:
            Q points are calculated at cell walls. To avoid numerical error 
            Q points must be placed at strategic positions i.e. at any changes 
            in shoreline angle, groynes etc.
            Therefore default grid frst calcs the y positions of Q points, and 
            then interpolates morphological Y points
            Note offshore contour // to initial shoreline and not perturbation,
            therefore wave refraction only due to angle of offshore contour.
        '''
        ####################################
#        N = 200 #
        if np.remainder(N,2) == 0:
            half = int(N/2)
        else:
            half = int(N+1)/2 #check if odd, this is halfway @ Y Points
#        dX = 2.5 #metres
        amp = 50 #amplitude of perturbation
        grad = 45 #gradient on perturbation
        rad = np.deg2rad(grad)
        ###################################
        
        self.x = np.linspace(0,dX*(N-1),num=N) #because 0 is included
        self.y = np.zeros(N)
        #offshore contour for wave refraction
        self.yOff = np.asarray([i for i in self.y])
        nGrid = int(amp/np.tan(rad)/dX)
        m = amp/nGrid
        start = half-nGrid
        end = half+nGrid
        
        #gauss
        x0 = self.x[half]
        if perturbation:
            self.yFill = amp/2*np.exp(-(self.x-x0)**2/(2*(amp/2)**2))
            yoffGauss = amp/10*np.exp(-(self.x-x0)**2/(2*(amp)**2))
            self.yOff+=yoffGauss
            self.y+=self.yFill
        else:
            self.y+=0
        #additional amount of sand... "headland"
        self.yAdd = np.max(self.y)
#        for i in range(int(start),int(half)):
#            self.y[i] = m*(i-start)
##            print(m*(i-start))
#        for i in range(int(half),int(end)):
#            self.y[i] = amp - m*(i-half)
            
        self.calcQpoints() 
        self.normals = self.calcNormals(self.x,self.y,smooth=smooth)
        self.normalsOff = self.calcNormals(self.x,self.yOff,smooth=smooth)
        
        self.yInit = np.asarray([i for i in self.y])
        
            
    def calcQpoints(self):
        '''
        calcs Q points from X,Y points. Y points located at cell centres.
        '''
        N = len(self.y)+1
        self.xPoints = np.zeros(N)
        self.qPoints = np.zeros(N)
        
        for i in range(N-2):       
            self.qPoints[i+1] = (self.y[i]+self.y[i+1])/2
            self.xPoints[i+1] = (self.x[i]+self.x[i+1])/2
#        print(type(self.qPoints))
    
    def calcSmooth(self):
        
        self.smoothY = [self.y[0],self.y[1]]
        
        for i in range(2,len(self.y)-2):
            self.smoothY.append((self.y[i-2]+self.y[i+1]+self.y[i]+self.y[i-1]+self.y[i+2])/5)
            
        self.smoothY.append(self.y[-2])
        self.smoothY.append(self.y[-1])
        self.smoothY = np.asarray(self.smoothY)
        
        
    def calcNormals(self,x,y,smooth=False):
        '''
        calculates the shorenormals
        works from left to right.
        first calcs in local coords and then applies world rotation, i.e:
            ^
            |
            |_____> (x,y)
        
        angles are calculated between Q points and applied at Y points.
        Works in radians
        returns back in degrees
        
        NOTE:
            need to say something about shoreline orientation in real world
            coords vs cartesian etc
        '''
#        N = len(self.y)
#        angles = np.zeros(N)
        
#        YDiff = self.qPoints[1::]-self.qPoints[0:-1]
#        XDiff = self.xPoints[1::]-self.xPoints[0:-1]
        if not smooth:
            YDiff = y[1::]-y[0:-1]
            XDiff = x[1::]-x[0:-1]
        else:
            self.calcSmooth()
            YDiff = self.smoothY[1::]-self.smoothY[0:-1]
            XDiff = x[1::]-x[0:-1]
        R = np.sqrt(YDiff**2+XDiff**2)
#        print()
        
        #cartesian coords
        angles = np.rad2deg(np.arcsin(YDiff/R)) #+ np.deg2rad(90) #avoids divide by zero with tan
        
        #now return angles to real world coords
        normals = 360-self.rot-angles
        
        return normals
        
            
    def readShoreline(self,xFile,yFile,smooth):
        '''
        reads the shoreline in from user defined.
        then calcs the Y points
        '''
        xPoints = []
        yPoints = []
        
        with open('%s%s'%(self.Path,xFile),'r') as f:
            for line in f.readlines():
                xPoints.append(float(line))
                
        with open('%s%s'%(self.Path,yFile),'r') as f:
            for line in f.readlines():
                yPoints.append(float(line))
                
        self.x = np.asarray(xPoints)
        self.y = np.asarray(yPoints)
        
        self.calcQpoints()
        self.calcNormals(smooth=smooth)
        
        self.yInit = np.asarray([i for i in self.y])
    
    def readSWL(self,file,default=True):
        '''
        reads the seawall.
        default is none
        '''
        
        if default:
            self.swl = np.asarray([-2 for i in range(len(self.y))])
        else:
            self.swl =[]
            with open('%s%s'%(self.Path,file),'r') as f:
                for ln in f.readlines():
                    self.swl.append(float(ln))
                    
            self.swl = np.asarray(self.swl)
            
    def readStr(self,file,default=True,on=True):
        '''
        reads in structures -> groynes.
        default places groyne at centre.
        Note bypassing NB now.
        must supply depths at groyne tips
        '''
        
        if not on:
            self.xStr = []
            self.yStr = []
        else:
            self.yStr = []
            self.xStr = []
            if default:
                self.yStr = [20,20,2]
                self.xStr = [150,180,165]
            
            else:
                with open('%s%s'%(self.Path,file),'r') as f:
                    for ln in f.readlines():
                        temp = ln.split(',')
                        self.xStr.append(temp[0])
                        self.yStr.append(temp[1])
            
    def update(self,y,qPoints,smooth):
        
        self.y = y
        self.calcQpoints()
#        print(type(self.qPoints))
#        self.calcNormals(smooth=smooth)
        self.normals = self.calcNormals(self.x,self.y,smooth=smooth)
#        self.normalsOff = self.calcNormals(self.x,self.yOff,smooth=smooth)
            
if __name__ == '__main__':
    
    import matplotlib.pyplot as plt
    
    test = gridGen('',200)
    test.defaultGrid()
    
    plt.plot(test.x,test.y)
    plt.vlines(test.xPoints,0,50,alpha=0.5,linewidth=0.5)
    plt.scatter(test.x,test.y)
    
    plt.show()
    
    
    
    
    
    
    
    
    
    
    
    
            
        
        