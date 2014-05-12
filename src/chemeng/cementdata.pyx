from chemeng import *
import csv
import math

print "\n\n\n\n\n"

class CementThermoData(ThermoConstantsType):
    def __init__(self, Tmin, Tmax, a, notes=""):
        ThermoConstantsType.__init__(self, Tmin, Tmax, notes) #Required 
        self.a = a
    
    def Cp0(self, T):
        return self.a[0] * T**(-2) + self.a[2] * T**(-0.5) + self.a[4] +2 * self.a[5] * T + self.a6 * T ** 2

    def Hf0(self, T):
        return -self.a[0] / T + self.a[1] + 2 * self.a[2] * T**(0.5) + self.a[4] * T + self.a[5] * T**2 + self.a[6] * T ** 3 / 3

    def S0(self, T):
        return -self.a[0] / (2 * T**2) - 2 * self.a[2] * T**(-0.5) + self.a[3] + self.a[4] * math.log(T)+ 2 *  self.a[5] * T + a[6] * T**(2) / 2

with open('/usr/local/PyChemEng/data/Cement.csv', 'rb') as datafile:
    reader = csv.reader(datafile, delimiter=',', quotechar='"')
    reader.next() #Skip the header
    for row in reader:
        species = row[0]
        phase = row[1]
        #S0 = float(row[2])
        #Hf0 = float(row[3])
        #V0 = float(row[4])
        #Gf0 = float(row[5])
        Tmin = float(row[6])
        Tmax = float(row[7])
        coeffs = map(float, row[8:15])
        notes = row[15]
        speciesData[species].registerPhase(phase)
        speciesData[species].registerPhaseCoeffs(CementThermoData(Tmin, Tmax, coeffs, notes), phase)
