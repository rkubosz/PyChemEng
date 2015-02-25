#!/usr/bin/env python
import math
import os

from chemeng.speciesdata cimport SpeciesDataType,AntioneConstants
from chemeng.speciesdata import speciesData
from chemeng.components import Components
import chemeng.NASAdata

####################################################################
# Antoine polynomials
# Functions should all take T as Kelvin and return P as Pascals, but
# vary due to coefficients.
####################################################################

#Original expression takes kelvin and produces bar, using the form 10^(A-B / (T+C))
cdef class AntionePolynomial(AntioneConstants):
    cdef double C[3]
    def __init__(AntionePolynomial self, double Tmin, double Tmax, C, comments):
        AntioneConstants.__init__(self, Tmin, Tmax, comments)
        cdef int i
        for i in range(3):
            self.C[i] = C[i]

    cpdef double Pvap(self, double T):
        return 1e5 * (10.0 ** (self.C[0] - self.C[1] / (T + self.C[2])))


#Not sure about this one
cdef class AntioneCnRPolynomial(AntioneConstants):
    cdef double C[3]
    def __init__(AntionePolynomial self, double Tmin, double Tmax, C, comments):
        AntioneConstants.__init__(self, Tmin, Tmax, comments)
        cdef int i
        for i in range(3):
            self.C[i] = C[i]

    cpdef double Pvap(self, double T):
        return (1.01325e5 / 760.0) * math.exp(self.C[0] - self.C[1] / (T + self.C[2]))

#Takes celcius and produces mmHg, using the form 10^(A-B / (T+C))
cdef class Antoine10CmmHgPolynomial(AntioneConstants):
    cdef double C[3]
    def __init__(AntionePolynomial self, double Tmin, double Tmax, C, comments):
        AntioneConstants.__init__(self, Tmin, Tmax, comments)
        cdef int i
        for i in range(3):
            self.C[i] = C[i]

    cpdef double Pvap(self, double T):
        return (1.01325e5 / 760.0) * (10 ** (self.C[0] - self.C[1] / (T - 273.15 + self.C[2])))

####################################################################
# Datafile loading
####################################################################
def parseTemperature(string):
    if string[-1] == 'C':
        return float(string[:-1])+273.15
    return float(string)

def initDataDir(directory):
    datafile = open(os.path.join(directory, 'antoine.inp'), 'r')
    for line in datafile:
        datafields = line.split()
    
        #Skip comments or other problems
        if line[0] == "#" or len(datafields) == 0:
            continue

        Tmin = parseTemperature(datafields[1])
        Tmax = parseTemperature(datafields[2])
        fitFunction = datafields[3]
        constants = map(float, datafields[4:])
        speciesData[datafields[0]].registerPhase("Liquid")

        if datafields[3] == "Antoine":
            speciesData[datafields[0]].registerAntoineCoeffs(AntionePolynomial(Tmin, Tmax, constants, "antione.inp"), "Liquid")
        elif datafields[3] == "AntoineCnR":
            speciesData[datafields[0]].registerAntoineCoeffs(AntioneCnRPolynomial(Tmin, Tmax, constants, "antione.inp"), "Liquid")
        elif datafields[3] == "Antoine10CmmHg":
            speciesData[datafields[0]].registerAntoineCoeffs(Antoine10CmmHgPolynomial(Tmin, Tmax, constants, "antione.inp"), "Liquid")
        else:
            raise Exception("Unrecognised Antione polynomial type"+datafields[3])

import chemeng.config
initDataDir(chemeng.config.datadir)
