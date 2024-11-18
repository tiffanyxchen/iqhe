#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 15 17:29:28 2024

@author: star
"""
import numpy as np
import matplotlib.pyplot as plt
from lattice import Lattice
from lattice import DoppedLattice
from hamiltonian import Hamiltonian

class Analysis(object):
    """ 
    A Hamiltonian matrix representing a square lattice under a perpendicular magnetic field.
    The Hamiltonian includes onsite energy, hopping energy, finite lattice size, and magnetic flux.
    """
    def __init__(self, h_matrix, h_matrix_imp):
        self.eigenvals     = self.eigenvalues    (h_matrix    )
        self.eigenvals_imp = self.eigenvalues    (h_matrix_imp)
        self.eDegeneracies = self.degeneracies()[0]
        self.nDegeneracies = self.degeneracies()[1]
        self.levels        = len(self.eDegeneracies)
    
    def getEigenvals(self):
        return self.eigenvals
    
    def getLevels(self):
        return self.levels       
        
    def eigenvalues(self, h_matrix):
        """Calculates the eigenvalues of the H"""
        
        eigenvalues = np.linalg.eigvals(h_matrix)
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


    def graphEnergyLevels(self, lat):
        
        num_bins = self.levels * 12
        plt.hist(self.eigenvals    , num_bins, facecolor='black', alpha=0.8, range=[-8, 8])
        plt.hist(self.eigenvals_imp, num_bins, facecolor='blue' , alpha=0.4, range=[-8, 8])
        plt.ylim(0, 1)  # Set the y-axis range here

        plt.xlabel('Energy Levels')
        plt.ylabel('Degeneracy')
        title = f"{lat.getNx()}x{lat.getNy()} Square Lattice: eps={lat.getEps()}, t={lat.getT()}, Q={lat.getQ()}, w={lat.getW()}"
        plt.title(title)
        plt.show()
        plt.savefig(title)

    
# Example Usage        
bigger        = Lattice(10, 10, 0, 2, 10)
biggerSquareH = Hamiltonian(bigger)

biggerD        = DoppedLattice(10, 10, 0, 2, 10, 2)
biggerSquareHD = Hamiltonian(biggerD)

result = Analysis(biggerSquareH.getMatrix(), biggerSquareHD.getMatrix())
result.graphEnergyLevels(biggerD)
 