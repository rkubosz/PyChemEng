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

from collections import namedtuple
ConstantsType = namedtuple('ConstantsType', ['Tmin', 'Tmax', 'fitFunction', 'constants', 'HConst', 'SConst'])

class SpeciesDataType:
    """
    This class represents the isobaric (P=P0) data for a species, and may include multiple phases
    """
    class PhaseData:        
        def __init__(self, name, comments):
            self.name = name
            self.comments = comments
            self.constants = []
            
        def __str__(self):
            output = "Phase{"+self.name+", \""+str(self.comments)+"\", "+str(len(self.constants))+" constants, "
            for data in self.constants:
                output += "["+str(data.Tmin)+", "+str(data.Tmax)+"] "
            return output + "}"

    def __init__(self, name, mass, elementalComposition):
        self.name = name
        self.mass = mass
        self.elementalComposition = elementalComposition
        self.phases = {}

    def __str__(self):
        output = "Species{"+self.name+", ["
        for num, phase in self.phases.iteritems():
            output += str(num)+":"+phase.name+", "
        return output[:-2] +"]}"

    def inDataRange(self, T, phase):
        for Tmin, Tmax, func, C, Hconst, Sconst in self.phases[phase].constants:
            if T >= Tmin and T <= Tmax:
                return True
        return False
        
    def registerPhase(self, phasenumber, phasename, comments):
        #Check if this phase has been registered before, if not,
        #create it
        if phasenumber not in self.phases:
            self.phases[phasenumber] = self.PhaseData(phasename, comments)

        #Check that the data is consistent (tests for consistency if
        #the phase existed before)
        if phasename != self.phases[phasenumber].name:
            raise Exception("Trying to register phase "+phasename+":"+str(phasenumber)+" but we have "+str(self.phaseNames))
        
    def registerPhaseCoeffs(self, Coeffs, phase):
        self.phases[phase].constants.append(ConstantsType(*Coeffs))
        #Ensure that the data is sorted from lowest to highest temperature range
        self.phases[phase].constants.sort(key = lambda x : (x.Tmin, x.Tmax))
        
    def Cp0(self, T, phase):
        for Tmin, Tmax, func, C, Hconst, Sconst in self.phases[phase].constants:
            if T >= Tmin and T <= Tmax:
                return R * fitFunctions[func](T, C)
        raise Exception("Cannot find valid Cp0 expression for "+self.name+" at "+str(T)+"K")

    def Hf0(self, T, phase):
        for Tmin, Tmax, func, C, Hconst, Sconst in self.phases[phase].constants:
            if T >= Tmin and T <= Tmax:
                return R * (fitFunctions[func+"Integrated"](T, C) + Hconst)
        raise Exception("Cannot find valid Hf0 expression for "+self.name+" at "+str(T)+"K")

    def S0(self, T, phase):
        for Tmin, Tmax, func, C, Hconst, Sconst in self.phases[phase].constants:
            if T >= Tmin and T <= Tmax:
                return R * (fitFunctions[func+"IntegratedOverT"](T, C) + Sconst)
        raise Exception("Cannot find valid S0 expression for "+self.name+" at "+str(T)+"K")

    def Gibbs0(self, T, phase):
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
