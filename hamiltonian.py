
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 14 03:15:29 2024

@author: star
"""
import numpy as np
from  math import pi
from cmath import e
import matplotlib.pyplot as plt
from lattice import Lattice

class Hamiltonian(object):
    """ 
    A Hamiltonian matrix representing a square lattice under a perpendicular magnetic field.
    The Hamiltonian includes onsite energy, hopping energy, finite lattice size, and magnetic flux.
    """
    def __init__(self, lat):
        self.matrix = self.setMatrix(lat)

    def getMatrix(self):
        return self.matrix
    
    def onsite(self, eps, N, mat):
        # On site energy eps on the diagonal
        for i in range(N):
            mat[i][i] = eps  
        return mat
    
    def onsiteImpurities(self, mat, lat):
        # On site energy eps on the diagonal
        p = lat.impurityPotentials()
        for i in range(mat.shape[0]):
            mat[i][i] = p[i]
        return mat

    def xHopping(self, t, Ny, N, mat):        
        #kinetic energy in x-direction
        for i in range(Ny):
            n = i + Ny
            while (n < N):
                mat[i][n] = t
                mat[n][i] = t
                i =  n 
                n += Ny
        return mat

    def yHopping(self, t, Ny, N, theda, mat):
        #kinetic energy in y-direction
        p = 0
        m = 1
        while (p < N):
            for r in range (p, p + Ny - 1):  
                theda_m = theda * m
                phase   = e**( 1j * theda_m)
                phaseC  = e**(-1j * theda_m)
                mat[r][r + 1] = t * phase
                mat[r + 1][r] = t * phaseC
            m = m + 1    
            p = p + Ny
        return mat

    def xBoundary(self, t, Nx, Ny, mat):    
        #kinetic energy in x-direction, Boundary    
        for c in range (Nx):
            mat[c][c + Ny * (Nx - 1)] = t
            mat[c + Ny * (Nx - 1)][c] = t    
        return mat
 
    def yBoundary(self, t, Nx, Ny, theda, mat):
        #kinetic energy in y-direction, Boundary
        for c in range (Nx):
            m = c + 1
            theda_m = theda * m
            phase   = e**( 1j * theda_m)
            phaseC  = e**(-1j * theda_m)
            x = c * Ny
            mat[x] [(c + 1) * Ny - 1] = t * phaseC   
            mat[(c + 1) * Ny - 1] [x] = t * phase
        return mat
    
    def setMatrix(self, lat):
        
        N   = lat.sites ()
        mat = np.full((N,N),0j)
        
        eps = lat.getEps()
        w   = lat.getW  ()
        t   = lat.getT  ()
 
        if w != 0:
            mat = self.onsiteImpurities(mat, lat)
        else:
            if eps != 0:
                mat = self.onsite(eps, N, mat)

        
        if t != 0:
            Nx    = lat.getNx() 
            Ny    = lat.getNy()
            phi   = 1/lat.getQ()
            theda = 2 * pi * phi
            
            mat   = self.xHopping (t, Ny, N,         mat)
            mat   = self.yHopping (t, Ny, N,  theda, mat)
            mat   = self.xBoundary(t, Nx, Ny,        mat)
            mat   = self.yBoundary(t, Nx, Ny, theda, mat)
        
        return mat
        
 
# Example Usage        
bigger   = Lattice(10, 10, 0, 1, 10, 1)
squareH  = Hamiltonian(bigger)
#biggerSquareH.graphEnergyLevels()
print (squareH)
'''
small = Lattice(2,2,0,1,2,1)
smallSquareH = Hamiltonian(small)
print(smallSquareH.getMatrix())
'''