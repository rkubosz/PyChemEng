#!/usr/bin/python

from numpy import *
from collections import defaultdict

import Molecule
from Phase import Phase

R = 8.314

################################################################################
################################################################################
class IdealSolution(Phase):

    def __init__(self, name, group_list, A, B):
        self.name = name
        self.group_list = group_list
#        self.molecule_dict = {}
        self.conc = {}
        self.x = {}
        self.T = 298.15
        self.p = 1.0e5

    def methodName(self):
        print 'Ideal Solution'

    def chemicalPotential(self):
        mu = {}
        mu_ref = {}
        for i in self.conc:
            pvap = Phase.moleculeDict[i].vaporPressure(self.T)
            mu_ref[i] = R*self.T*log(pvap/self.p)
            mu[i] = -1.0e30
            if (self.x[i] > 0.0):
                mu[i] = mu_ref[i] + R*self.T*log(self.x[i])
        return mu

