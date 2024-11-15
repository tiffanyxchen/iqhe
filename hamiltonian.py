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
        self.lat           = lat
        self.matrix        = self.setMatrix      (               )   
        self.eigenvals     = self.eigenvalues    (self.matrix    )
        self.matrix_imp    = self.setMatrix_imp  (               )   
        self.eigenvals_imp = self.eigenvalues    (self.matrix_imp)
        self.eDegeneracies = self.degeneracies()[0]
        self.nDegeneracies = self.degeneracies()[1]
        self.levels        = len(self.eDegeneracies)

    def getMatrix(self):
        return self.matrix
    
    def getEigenvals(self):
        return self.eigenvals
    
    def getLevels(self):
        return self.levels
    
    def onsite(self, eps, N, mat):
        # On site energy eps on the diagonal
        for i in range(N):
            mat[i][i] = eps  
        return mat
    
    def onsiteImpurities(self, mat):
        # On site energy eps on the diagonal
        for i in range(mat.shape[0]):
            mat[i][i] = self.lat.impurityPotentials()[i]
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
    
    def setMatrix(self):
        
        N   = self.lat.sites ()
        mat = np.full((N,N),0j)
        
        eps = self.lat.getEps()
        #w   = self.lat.getW  ()
        t   = self.lat.getT  ()
        
        if eps != 0:
            mat = self.onsite(eps, N, mat)

        #if w != 0:
        #    mat = self.onsiteImpurites(eps, N, mat, w)

        
        if t != 0:
            Nx    = self.lat.getNx() 
            Ny    = self.lat.getNy()
            phi   = 1/self.lat.getQ()
            theda = 2 * pi * phi
            
            mat   = self.xHopping (t, Ny, N,         mat)
            mat   = self.yHopping (t, Ny, N,  theda, mat)
            mat   = self.xBoundary(t, Nx, Ny,        mat)
            mat   = self.yBoundary(t, Nx, Ny, theda, mat)
        
        return mat
    
    
    def setMatrix_imp(self):
        w   = self.lat.getW()
        mat = self.matrix    
        if w != 0:
            mat = self.onsiteImpurities(mat)
        return mat
       
        
    def eigenvalues(self, matrix):
        """Calculates the eigenvalues of the H"""
        
        eigenvalues = np.linalg.eigvals(matrix)
        realE = (np.round(eigenvalues.real,2))
        realE_sorted = sorted(realE)
        return realE_sorted
                                                                                                                                                                         
    def degeneracies(self):
        #count energy levels and degeneracies
        eigenvals = self.getEigenvals()
        N         = len(self.eigenvals)
        
        nLevels = 0
        energyLevels  = []
        degeneracies  = []

        c = 0
        while (c < N):
            print( str(nLevels) + "st energy level is: " + str(eigenvals[c]))
            energyLevels.insert(nLevels, eigenvals[c])
            
            i = c + 1
            numDeg = 1
            while (i < N and eigenvals[c] == eigenvals[i]):
                numDeg += 1
                i += 1
                
            print("  Its degeneracy is: " + str(numDeg))
            degeneracies.insert(nLevels, numDeg)
            print("    ")
            
            if (i < N and eigenvals[c] != eigenvals[i]):
                    nLevels = nLevels + 1
            c = i
        return energyLevels, degeneracies

    
    def traceEigenvals(self):
        return sum(self.eigenvals)

    def printTrace(self):
        print("Adding eigenvalues: " + str(np.round(self.traceEigenvals(), 1)))
    
    def totEigenstates(self):
        return sum(self.nDegeneracies)
    
    def printTotEigenstates(self):
        print("\nAdding density of states: " + str(self.totEigenstates()))
        
    
    def __str__(self):
        description = (
            f"\nThe {self.lat.getNx()} x {self.lat.getNy()} square lattice has \n"
            f"On site energy = {self.lat.getEps()}; t = {self.lat.getT()}; theta = 2*pi/{self.lat.getQ()}\n\n"
            f"Its Hamiltonian is represented by: the {self.matrix.shape[0]} x {self.matrix.shape[0]} matrix with B-field:\n\n"
            f"The eigenvalues of the matrix with B-field:\n{self.getEigenvals()}\n\n"
            f"# of unique energy levels: {len(self.eDegeneracies)}.\n\nThese are:\n{self.eDegeneracies}\n\n"
            f"The corresponding degeneracies are:\n{self.nDegeneracies}\n"
                        )
        return description


    def graphEnergyLevels(self):
        
        num_bins = self.levels * 12
        plt.hist(self.eigenvals    , num_bins, facecolor='black', alpha=0.8, range=[-8, 8])
        plt.hist(self.eigenvals_imp, num_bins, facecolor='blue' , alpha=0.4, range=[-8, 8])
        plt.ylim(0, 1)  # Set the y-axis range here

        plt.xlabel('Energy Levels')
        plt.ylabel('Degeneracy')
        title = f"{self.lat.getNx()}x{self.lat.getNy()} Square Lattice: eps={self.lat.getEps()}, t={self.lat.getT()}, Q={self.lat.getQ()}, w={self.lat.getW()}"
        plt.title(title)
        plt.show()
        plt.savefig(title)

    
# Example Usage        
bigger        = Lattice(10, 10, 0, 2, 10, 2)
biggerSquareH = Hamiltonian(bigger)
biggerSquareH.graphEnergyLevels()
print (biggerSquareH)
