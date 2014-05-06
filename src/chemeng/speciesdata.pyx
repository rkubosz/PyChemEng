#!/usr/bin/env python
# distutils: language = c++
# cython: profile=True

from chemeng.components cimport Components
from chemeng.elementdata import elements

cdef class ThermoConstantsType:
    def __init__(ThermoConstantsType self, double Tmin, double Tmax, str comments):
        self.Tmin = Tmin
        self.Tmax = Tmax
        self.comments = comments

    cpdef double Cp0(self, double T):
        return 0.0

    cpdef double S0(self, double T):
        return 0.0

    cpdef double Hf0(self, double T):
        return 0.0
    
cdef class PhaseData:
    cdef public string name
    cdef public list constants
    def __init__(self, name):
        self.name = name
        self.constants = []
        
    def __str__(self):
        output = "Phase{'"+self.name+"', T=["
        if len(self.constants) > 0:
            for data in self.constants:
                output += str(data.Tmin)+'->'+str(data.Tmax)+"K, "
            output = output[:-2]
        return output+"]}"
    
    def __repr__(self):
        return self.__str__()

cdef class SpeciesDataType:
    """
    This class represents the isobaric (P=P0) data for a species, and may include multiple phases
    """

    def __init__(self, name, Components elementalComposition):
        self.name = name
        self.mass = 0
        for key,amount in elementalComposition.iteritems():
            self.mass += elements[key].mass * amount
        self.elementalComposition = elementalComposition
        self.phases = {}
        self.antoineData = []

    def __str__(self):
        output = "Species{'"+self.name+"', phases=["
        for name in self.phases:
            output += name+", "
        return output[:-2] +"], elementalComposition="+str(self.elementalComposition)+"}"

    def __repr__(self):
        return self.__str__()

    def registerPhase(self, phasename):
        if phasename not in self.phases:
            self.phases[phasename] = PhaseData(phasename)
        
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

cpdef registerSpecies(name, Components elementalComposition):
    for element in elementalComposition.keys():
        if element not in elements:
            raise Exception("Species "+name+" with elemental composition "+str(elementalComposition)+" has an unknown element, "+element)

    if name not in speciesData:
        speciesData[name] = SpeciesDataType(name, elementalComposition)
