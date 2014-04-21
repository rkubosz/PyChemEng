#!/usr/bin/python

import UNIFACGroup

class Molecule:

    def __init__(self, name, id):
        self.name = name
        self.id = id
        self.group_list = []
        self.nu = {}
        self.r = 0.0
        self.q = 0.0
        self.G = {}
        self.antoineA = 0.0
        self.antoineB = 0.0
        self.antoineC = 0.0

    def get_name(self):
        return self.name

    def add_group(self, group, nu):
        self.group_list.append( group )
        name = group.get_name()
        self.nu[name] = nu

    def print_group_list(self):
        for group in self.group_list:
            name = group.get_name()
            print name, self.nu[name]

    def set_parameters(self):
        r = 0.0
        q = 0.0
        self.G = {}
        for group in self.group_list:
            name = group.get_name()
            r += self.nu[name] * group.get_R()
            q += self.nu[name] * group.get_Q()
            self.G[name] = self.nu[name] * group.get_Q()
        self.r = r
        self.q = q

    def setAntoineCoefficients(self, coeff):
        self.antoineA = coeff[0]
        self.antoineB = coeff[1]
        self.antoineC = coeff[2]

    def vaporPressure(self, T):
        p = 10.0 ** (self.antoineA - self.antoineB/(T+self.antoineC))
        return p * 1.0e5
        
    def get_r(self):
        return self.r

    def get_q(self):
        return self.q

    def get_G(self):
        return self.G

    def gRef(self):
        return 0.0

