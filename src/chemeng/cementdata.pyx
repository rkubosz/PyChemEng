#!/usr/bin/env python
from chemeng import *
import chemeng.NASAdata
import csv
import math

################ NIST data source

registerSpecies("CaAl2SiO6", Components({'Ca':1, 'Al':2, 'Si':1, 'O':6}))
registerSpecies("CaAl2Si2O8", Components({'Ca':1, 'Al':2, 'Si':2, 'O':8}))
registerSpecies("Ca2Al2Si3O10(OH)2", Components({'Ca':2, 'Al':2, 'Si':3, 'O':12, 'H':2}))
registerSpecies("Ca2Al2SiO7", Components({'Ca':2, 'Al':2, 'Si':1, 'O':7}))
registerSpecies("Ca3Al2Si3O12", Components({'Ca':3, 'Al':2, 'Si':3, 'O':12}))
registerSpecies("Al2Si4O10(OH)2", Components({'Al':2, 'Si':4, 'O':12, 'H':2}))
registerSpecies("Al2Si2O5(OH)4", Components({'Al':2, 'Si':2, 'O':9, 'H':4}))
registerSpecies("Ca2Al3Si3O12(OH)", Components({'Ca':2, 'Al':3, 'Si':3, 'O':13, 'H':1}))
registerSpecies("CaAl4Si2O10(OH)2", Components({'Ca':1, 'Al':4, 'Si':2, 'O':12, 'H':2}))
registerSpecies("CaSiO3", Components({'Ca':1, 'Si':1, 'O':3}))
registerSpecies("Ca2SiO4", Components({'Ca':2, 'Si':1, 'O':4}))
registerSpecies("Ca3SiO5", Components({'Ca':3, 'Si':1, 'O':5}))
registerSpecies("Ca3Si2O7", Components({'Ca':3, 'Si':2, 'O':7}))


class CementThermoData(ThermoConstantsType):
    def __init__(self, Tmin, Tmax, a, notes=""):
        ThermoConstantsType.__init__(self, Tmin, Tmax, notes) #Required 
        self.a = a
    
    def Cp0(self, T):
        return self.a[0] * T**(-2) + self.a[2] * T**(-0.5) + self.a[4] +2 * self.a[5] * T + self.a[6] * T ** 2

    def Hf0(self, T):
        return -self.a[0] / T + self.a[1] + 2 * self.a[2] * T**(0.5) + self.a[4] * T + self.a[5] * T**2 + self.a[6] * T ** 3 / 3

    def S0(self, T):
        return -self.a[0] / (2 * T**2) - 2 * self.a[2] * T**(-0.5) + self.a[3] + self.a[4] * math.log(T)+ 2 *  self.a[5] * T + self.a[6] * T**(2) / 2

    def __str__(self):
        retval = "CementPolynomial{Tmin="+str(self.Tmin)+", Tmax="+str(self.Tmax)+", notes='"+self.comments+"', a=["
        for i in range(7):
            retval+=str(self.a[i])+", "
        retval = retval[:-2] + "]}"
        return retval

    def __repr__(self):
        return self.__str__()

import chemeng.config
import os
with open(os.path.join(chemeng.config.datadir, 'Cement.csv'), 'rb') as datafile:
    reader = csv.reader(filter(lambda row: row[0]!='!', datafile), delimiter=',', quotechar='"')
    reader.next() #Skip the header
    for row in reader:
        #Load data
        species = row[0]
        phase = row[1]
        Tmin = float(row[6])
        Tmax = float(row[7])
        a = map(float, row[8:15])
        if len(row[3].strip()) != 0:
            Hf0 = float(row[3]) * 1000.0
            a[1] += Hf0
        notes = row[15]
        speciesData[species].registerPhase(phase)
        speciesData[species].registerPhaseCoeffs(CementThermoData(Tmin, Tmax, a, notes), phase)
        if len(row[4].strip()) != 0:
            V0 = float(row[4])



################ Data gathered from multiple sources

registerSpecies("Ca3Al2O6", Components({'Ca':3, 'Al':2, 'O':6}))
registerSpecies("Ca12Al14O33", Components({'Ca':12, 'Al':14, 'O':33}))
registerSpecies("CaAl2O4", Components({'Ca':1, 'Al':2, 'O':4}))
registerSpecies("CaAl4O7", Components({'Ca':1, 'Al':4, 'O':7}))
registerSpecies("CaFe2O4", Components({'Ca':1, 'Fe':2, 'O':4}))
registerSpecies("Ca2Fe2O5", Components({'Ca':2, 'Fe':2, 'O':5}))
registerSpecies("MgFe2O4", Components({'Mg':1, 'Fe':2, 'O':4}))
registerSpecies("Al2TiO5", Components({'Al':2, 'Ti':1, 'O':5}))
registerSpecies("Fe2TiO4", Components({'Fe':2, 'Ti':1, 'O':4}))
registerSpecies("Zn2TiO4", Components({'Zn':2, 'Ti':1, 'O':4}))


class Cement2ThermoData(ThermoConstantsType):
    def __init__(self, Tmin, Tmax, a, notes=""):
        ThermoConstantsType.__init__(self, Tmin, Tmax, notes) #Required 
        self.a = a
    
    def Cp0(self, T):
        return self.a[0] + 2* self.a[1] * T - self.a[2]/(T**2)

    def Hf0(self, T):
        return self.a[0]*T + self.a[1] * (T**2) + self.a[2]/T + self.a[3]

    def S0(self, T):
        return self.a[0]*math.log(T) + 2 * self.a[1] * T + self.a[2]/(2*T**2) + self.a[4] 

    def __str__(self):
        retval = "CementPolynomial{Tmin="+str(self.Tmin)+", Tmax="+str(self.Tmax)+", notes='"+self.comments+"', a=["
        for i in range(7):
            retval+=str(self.a[i])+", "
        retval = retval[:-2] + "]}"
        return retval

    def __repr__(self):
        return self.__str__()

import chemeng.config
import os
with open(os.path.join(chemeng.config.datadir, 'Cement2.csv'), 'rb') as datafile:
    reader = csv.reader(filter(lambda row: row[0]!='!', datafile), delimiter=',', quotechar='"')
    reader.next() #Skip the header
    for row in reader:
        #Load data
        species = row[0]
        phase = row[1]
        Tmin = float(row[4])
        Tmax = float(row[5])
        a = map(float, row[6:11])
        a[:] = [x*4.184 for x in a] 
        if len(row[3].strip()) != 0:
            Hf0 = float(row[3])*1000.0 
            a[3] += Hf0
        if len(row[2].strip()) != 0:
            S0 = float(row[2])*4.184
            a[4] += S0 
        notes = row[11]
        speciesData[species].registerPhase(phase)
        speciesData[species].registerPhaseCoeffs(Cement2ThermoData(Tmin, Tmax, a, notes), phase)
