#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 14 03:12:19 2024

@author: star
"""
import random

class Lattice(object):
    """A single atomic sheet with sizes: nx and ny; onsite potential energy: eps; hopping: t; magnetic flux 1/Q;"""

    def __init__(self, Nx, Ny, eps, t, Q):
        """
        Args:
            Nx (int  ): number of sites on x-direction
            Ny (int  ): number of sites on y-direction
            eps(float): onsite energy
            t  (float): nearest neighbor hopping potential
            1/Q(int  ): magnetic flux
            w  (float): impurity width
        """
        self.Nx  = Nx                             
        self.Ny  = Ny        
        self.eps = eps                           
        self.t   = t                              
        self.Q   = Q
         
    def getEps(self): 
        return self.eps
    def setEps(self,eps):
        self.eps = eps
    
    def getT(self):
        return self.t
    def setT(self, t):
        self.t = t
    
    def getNx(self):
        return self.Nx
    def setNx(self, Nx):
        self.Nx = Nx
    
    def getNy(self):
        return self.Ny
    def setNy(self, Ny):
        self.Ny = Ny
    
    def getQ(self):
        return self.Q
    def setQ(self, Q):
        self.Q = Q

    def sites(self):
        """
        Calculates the total number of sites on the Lattice sheet. Called by the 
        setHamiltonian() in Hamiltonian class. 

        Returns:
            An integer: total number of sites, Nx * Ny
        """
        return self.Nx * self.Ny
    
    def __str__(self):
        return "This lattice model has a size: "            \
            + self.description()
    
    def description(self):
        return str(self.getNx()) + " x " + str(self.getNy()) \
            + ", eps = "                                    \
            + str(self.getEps())                            \
            + ", t = " + str(self.getT()) + ", Q = " + str(self.getQ())


class DoppedLattice(Lattice):
    tag = 1
    def __init__(self, Nx, Ny, eps, t, Q, w):    
        Lattice.__init__(self, Nx, Ny, eps, t, Q)
        self.w = w
        self.dID = DoppedLattice.tag
        DoppedLattice.tag += 1
        
    def getW(self):
        return self.w
    def setW(self, w):
        self.w = w    
        
    def impurityPotentials(self):
        if self.w != 0:
            p = []
            start = self.eps - self.w/2
            for i in range(self.sites()):
                p.append(random.random() * self.w + start)
        return p    
    
    def __str__(self):
        return "This dopped model has a size: "      \
            + self.description()                     \
            + ", w = " + str(self.getW())

        
# Test instance 
smallSheet = DoppedLattice(10, 10, 0, 1, 10, 1)
print (smallSheet)
#print (smallSheet.impurityPotentials())
  