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
        self.set_x()
        self.setReferenceState()
        mu = {}
        for i in self.conc:
            mu[i] = -1.0e30
            if (self.x[i] > 0.0):
                mu[i] = self.mu_ref[i] + R*self.T*log(self.x[i])
        return mu

    def setReferenceState(self):
        self.mu_ref = {}
        for i in self.conc:
            self.mu_ref[i] = self.moleculeDict[i].gRef()
            pvap = Phase.moleculeDict[i].vaporPressure(self.T)
            self.mu_ref[i] -= R*self.T*log(pvap/self.p)
