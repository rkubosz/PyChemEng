#!/usr/bin/python

from numpy import *
from collections import defaultdict

from Phase import *

R = 8.314

################################################################################
################################################################################
class IdealGas(Phase):

    def __init__(self, name):
        self.name = name
        self.T = 298.15
        self.x = {}
        self.conc = {}

    def methodName(self):
        print "Ideal Gas"

    def chemicalPotential(self):
        Ntot = 0.0
        for i in self.conc:
            Ntot += self.conc[i]
        for i in self.conc:
            self.x[i] = self.conc[i]/Ntot

        self.setReferenceState()
        mu = {}
        for i in self.conc:
            mu[i] = self.mu_ref[i] + R*self.T*log(self.x[i])
        return mu

    def setReferenceState(self):
        self.mu_ref = {}
        for i in self.conc:
            self.mu_ref[i] = self.moleculeDict[i].gRef()
