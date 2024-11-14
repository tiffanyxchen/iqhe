#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 14 03:12:19 2024

@author: star
"""
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
        """
        self.Nx  = Nx                             
        self.Ny  = Ny        
        self.eps = eps                           
        self.t   = t                              
        self.Q   = Q
        
    def getEps(self):
        return self.eps
    
    def getT(self):
        return self.t
    
    def getNx(self):
        return self.Nx
    
    def getNy(self):
        return self.Ny
    
    def getQ(self):
        return self.Q
    
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
            + str(self.getNx()) + " x " + str(self.getNy()) \
            + ", eps = "                                    \
            + str(self.getEps())                            \
            + ", t = " + str(self.getT()) + ", Q = " + str(self.getQ())  

# Test instance 
smallSheet = Lattice(10, 10, 0, 1, 10)
print (smallSheet)  
 