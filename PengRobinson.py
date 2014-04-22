#!/usr/bin/python



class PengRobinson(Phase):

    def __init__(self):
        self.conc = {}

    def get_pressure(self):
        amix = 0.0
        bmix = 0.0
        for i in self.conc:
            for j in self.conc:
                amix += a[i][j]*self.x[i]*self.x[j]
            bmix += b[i]*self.x[j]

        p = R*self.T/(self.V-bmix) \
            - amix/(self.V*self.V+2.0*bmix*self.V-bmix*bmix)

        return p

