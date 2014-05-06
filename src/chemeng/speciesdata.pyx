#!/usr/bin/env python
# distutils: language = c++
# cython: profile=True

from chemeng.components cimport Components
from chemeng.elementdata import elements

cdef class ThermoConstantsType:
    def __init__(ThermoConstantsType self, double Tmin, double Tmax):
        self.Tmin = Tmin
        self.Tmax = Tmax

    cpdef double Cp0(self, double T):
        return 0.0

    cpdef double S0(self, double T):
        return 0.0

    cpdef double Hf0(self, double T):
        return 0.0
    
cdef class PhaseData:
    cdef public string name
    cdef public string comments
    cdef public list constants
    def __init__(self, name, comments):
        self.name = name
        self.comments = comments
        self.constants = []
        
    def __str__(self):
        output = "Phase{'"+self.name+"', T=["
        if len(self.constants) > 0:
            for data in self.constants:
                output += str(data.Tmin)+'->'+str(data.Tmax)+"K, "
            output = output[:-2]
        return output+"], comments='"+self.comments+"'}"
    
    def __repr__(self):
        return self.__str__()

cdef class SpeciesDataType:
    """
    This class represents the isobaric (P=P0) data for a species, and may include multiple phases
    """

    def __init__(self, name, mass, elementalComposition):
        self.name = name
        self.mass = mass
        self.elementalComposition = Components(elementalComposition)
        self.phases = {}
        self.antoineData = []

    def __str__(self):
        output = "Species{'"+self.name+"', phases=["
        for name in self.phases:
            output += name+", "
        return output[:-2] +"], elementalComposition="+str(self.elementalComposition)+"}"

    def __repr__(self):
        return self.__str__()

    def registerPhase(self, phasename, comments):
        if phasename not in self.phases:
            self.phases[phasename] = PhaseData(phasename, comments)
        if self.phases[phasename].comments != comments:
            raise Exception("Comments on phase ("+phasename+") are changing current:an"+self.phases[phasename].comments+"\nNew:\n"+comments)
        
    def registerPhaseCoeffs(self, coeffs, string phasename):
        self.phases[phasename].constants.append(coeffs)
        
    cpdef bint inDataRange(SpeciesDataType self, double T, string phase):
        cdef PhaseData data = self.phases[phase]
        cdef ThermoConstantsType datum
        for datum in data.constants:
            if T >= datum.Tmin and T <= datum.Tmax:
                return True
        return False
        
    cpdef double Cp0(SpeciesDataType self, double T, string phase):
        cdef PhaseData data = self.phases[phase]
        cdef ThermoConstantsType datum
        for datum in data.constants:
            if T >= datum.Tmin and T <= datum.Tmax:
                return datum.Cp0(T)
        raise Exception("Cannot find valid Cp0 expression for "+self.name+" at "+str(T)+"K\n"+self.validRanges(phase))

    cpdef double Hf0(SpeciesDataType self, double T, string phase):
        cdef PhaseData data = self.phases[phase]
        cdef ThermoConstantsType datum
        for datum in data.constants:
            if T >= datum.Tmin and T <= datum.Tmax:
                return datum.Hf0(T)
        raise Exception("Cannot find valid Hf0 expression for "+self.name+" at "+str(T)+"K\n"+self.validRanges(phase))

    cpdef double S0(SpeciesDataType self, double T, string phase):
        cdef PhaseData data = self.phases[phase]
        cdef ThermoConstantsType datum
        for datum in data.constants:
            if T >= datum.Tmin and T <= datum.Tmax:
                return datum.S0(T)
        raise Exception("Cannot find valid S0 expression for "+self.name+" at "+str(T)+"K\n"+self.validRanges(phase))

    cpdef double Gibbs0(self, double T, string phase):
        cdef PhaseData data = self.phases[phase]
        cdef ThermoConstantsType datum
        for datum in data.constants:
            if T >= datum.Tmin and T <= datum.Tmax:
                return datum.Hf0(T) - T * datum.S0(T) 
        raise Exception("Cannot find valid Gibbs0 expression for "+self.name+" at "+str(T)+"K\n"+self.validRanges(phase))

    cpdef str validRanges(self, phase):
        cdef str msg = "Valid ranges:\n"
        cdef ThermoConstantsType datum
        for datum in self.phases[phase].constants:
            msg += msg+"["+str(datum.Tmin)+", "+str(datum.Tmax)+"]\n"
        return msg
    
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

from libcpp.string cimport string
from libcpp.map cimport map

speciesData = {}

cpdef list findSpeciesData(str search):
    cdef list retval = []
    for key,species in speciesData.iteritems():
        if search in key:
            retval.append(species)
            continue
        for phasekey, phase in species.phases.iteritems():
            if search in phasekey or search in phase.name or search in phase.comments:
                retval.append(species)
                continue
    return retval

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
