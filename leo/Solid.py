#!/usr/bin/python

from numpy import *
from collections import defaultdict

from Phase import *


R = 8.314

################################################################################
################################################################################
class Solid(Phase):

    def __init__(self, name):
        self.name = name
#        self.molecule_dict = {}
        self.x = {}
        self.conc = {}
        self.T = 298.15

    def methodName(self):
        print "Solid"

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
