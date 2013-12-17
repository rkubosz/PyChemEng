#!/usr/bin/env python
####################################################################
# Physical constants
####################################################################
#Constant used in the NASA data set! (If changed, will need to rescale the data)
R = 8.31451 

#STP is 25 celcius at 1 bar
T0 = 273.15 + 25
P0 = 1.0e5

from Elements import elements

####################################################################
# Fit functions
####################################################################
#fitFunctions is a dictionary of functions used for data fitting. Each
#function must take two arguments, T and Coeffs.
fitFunctions={}

def registerCpFitFunction(name, function, integratedfunction, integratedfunctionOverT):
    if name in fitFunctions:
        raise Exception("This function name is already in use!")
    fitFunctions[name] = function
    fitFunctions[name+"Integrated"] = integratedfunction
    fitFunctions[name+"IntegratedOverT"] = integratedfunctionOverT

####################################################################
# Helper functions
####################################################################
#Functions used throughout the code
def relativeError(val, ref):
    import math
    return math.fabs((val - ref) / (ref + (ref==0)))

####################################################################
# Species data
####################################################################
#Species data is a dictionary of thermodynamic data on different
#species/components and their phases. This is addressed using the
#chemical formula of the species to be studied.  
speciesData={}

class SpeciesDataType:
    """
    This class represents the isobaric (P=P0) data for a species, and may include multiple phases
    """
    from collections import namedtuple
    ThermoDataType = namedtuple('ThermoDataType', ['Tmin', 'Tmax', 'fitFunction', 'constants', 'HConst', 'SConst'])

    def __init__(self, name, mass, elementalComposition):
        self.name = name
        self.phaseCoefficients = {}

    def inDataRange(self, T, phase = 0):
        for Tmin, Tmax, func, C, Hconst, Sconst in self.phaseCoefficients[phase]:
            if T >= Tmin and T <= Tmax:
                return True
        return False
        
    def registerPhaseCoeffs(self, Coeffs, phase):
        if phase not in self.phaseCoefficients:
            self.phaseCoefficients[phase] = []
        self.phaseCoefficients[phase].append(Coeffs)
        #Ensure that the data is sorted from lowest to highest temperature range
        self.phaseCoefficients[phase].sort(key = lambda x : (x.Tmin, x.Tmax))
        
    def Cp0(self, T, phase=0):
        for Tmin, Tmax, func, C, Hconst, Sconst in self.phaseCoefficients[phase]:
            if T >= Tmin and T <= Tmax:
                return R * fitFunctions[func](T, C)
        raise Exception("Cannot find valid Cp0 expression for "+self.name+" at "+str(T)+"K")

    def Hf0(self, T, phase=0):
        for Tmin, Tmax, func, C, Hconst, Sconst in self.phaseCoefficients[phase]:
            if T >= Tmin and T <= Tmax:
                return R * (fitFunctions[func+"Integrated"](T, C) + Hconst)
        raise Exception("Cannot find valid Hf0 expression for "+self.name+" at "+str(T)+"K")

    def S0(self, T, phase=0):
        for Tmin, Tmax, func, C, Hconst, Sconst in self.phaseCoefficients[phase]:
            if T >= Tmin and T <= Tmax:
                return R * (fitFunctions[func+"IntegratedOverT"](T, C) + Sconst)
        raise Exception("Cannot find valid S0 expression for "+self.name+" at "+str(T)+"K")

    def Gibbs0(self, T, phase=0):
        return Hf0(T, phase) - T * S0(T, phase)

def registerSpecies(name, elementalComposition, mass=None):
    calcMass = 0
    for element, amount in elementalComposition.iteritems():
        if element not in elements:
            raise Exception("Species "+name+" with elemental composition "+str(elementalComposition)+" has an unknown element, "+element)
        else:
            calcMass += elements[element].mass * amount
    if mass is None:
        mass = calcMass
    else:
        if relativeError(calcMass, mass) > 0.01:
            raise Exception("Calculated species mass is significantly different when compared to passed value. Is the elemental composition or molMass correct?\n" + name + ", " + str(elementalComposition) + ", " + str(mass) + ", " + str(calcMass))
    if name not in speciesData:
        speciesData[name] = SpeciesDataType(name, mass, elementalComposition)
