#!/usr/bin/python


class Phase:

    moleculeDict = {}
    
    def __init__(self, name):
        self.name = name
        
    def get_name(self):
        return self.name

    def molecules(self):
        return self.moleculeDict
    
    def moleculeList(self):
        return self.conc.keys()
    
    def moleNumbers(self):
        return self.conc

    def getMoleNumber(self, a):
        N = 0.0
        if (a in self.conc):
            N = self.conc[a]
        return N
    
    def get_T(self):
        return self.T

    def setTemperature(self, T):
        self.T = T

    def setPressure(self, p):
        self.p = p
    
    def set_conc(self, name, conc):
        self.conc[name] = conc
        
    def set_x(self):
        sum = 0.0
        for name in self.conc.keys():
            sum += self.conc[name]
        self.x = {}
        for name in self.conc.keys():
            self.x[name] = self.conc[name] / sum

    def get_x(self):
        return self.x

    def add_molecule(self, mol):
        name = mol.get_name()
        self.moleculeDict[name] = mol
        self.conc[name] = 0.0

    def print_molecule_list(self):
        for name in self.moleculeDict.keys():
            print name, self.conc[name]

    def print_conc(self):
        for name in self.moleculeDict.keys():
            print name, self.conc[name]

    def methodName(self):
        pass

    def chemicalPotential(self):
        pass

    def gibbsFreeEnergy(self):
        mu = self.chemicalPotential()
        G = 0.0
        for i in self.conc:
            G += self.conc[i]*mu[i]
        return G

