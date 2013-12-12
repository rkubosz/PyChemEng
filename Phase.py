#!/usr/bin/python


class Phase:
    def __init__(self, name):
        self.name = name
        self.conc = {}
        
    def get_name(self):
        return self.name

    def molecules(self):
        return self.molecule_dict
    
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

    def set_T(self, T):
        self.T = T
    
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
        self.molecule_dict[name] = mol
        self.conc[name] = 0.0

    def print_molecule_list(self):
        for name in self.molecule_dict.keys():
            print name, self.conc[name]

    def print_conc(self):
        for name in self.molecule_dict.keys():
            print name, self.conc[name]

    def methodName(self):
        pass

    def chemicalPotential(self):
        pass

    def gibbsFreeEnergy(self):
        mu = self.chemicalPotential()
        G = 0.0
        for i in self.molecule_dict:
            G += self.conc[i]*mu[u]
        return G

