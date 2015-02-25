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

cdef class AntioneConstants:
    def __init__(AntioneConstants self, double Tmin, double Tmax, str comments):
        self.Tmin = Tmin
        self.Tmax = Tmax
        self.comments = comments
        
    cpdef double Pvap(self, double T):
        return 0.0

    
cdef class PhaseData:
    cdef public string name
    cdef public list constants
    cdef public list antioneconstants

    def __init__(self, name):
        self.name = name
        self.constants = []
        self.antioneconstants = []
        
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

    def registerAntoineCoeffs(self, coeffs, string phasename):
        if not str(phasename) in self.phases.keys():
            raise Exception("No phase "+phasename+" in species "+self.name)
        self.phases[phasename].antioneconstants.append(coeffs)
        
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

    cpdef double Pvap(SpeciesDataType self, double T, string phase):
        cdef PhaseData data = self.phases[phase]
        cdef AntioneConstants datum
        for datum in data.antioneconstants:
            if T >= datum.Tmin and T <= datum.Tmax:
                return datum.Pvap(T)
        raise Exception("Cannot find valid Pvap expression for "+self.name+" at "+str(T)+"K\n"+self.validRanges(phase))

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

cpdef list findSpeciesData(str keyword = "", Components composition = Components({}), list elements = []):
    cdef list retval = []

    for key,species in speciesData.iteritems():
        #Check that the requested elements are there
        fail = False
        for element in elements:
            if element not in species.elementalComposition:
                fail = True
                break
        for element, amount in composition.iteritems():
            if element not in species.elementalComposition:
                fail = True
                break
            if amount != species.elementalComposition[element]:
                fail = True
                break
        if fail:
            continue

        fail = True
        #Now check for keywords
        if keyword in key:
            fail = False
        for phasekey, phase in species.phases.iteritems():
            if keyword in phasekey or keyword in phase.name:
                fail = False
                break
            for constants in phase.constants:
                if keyword in constants.comments:
                    fail = False
                    break
                
        if not fail:
            retval.append(species)

    return retval

cpdef registerSpecies(name, Components elementalComposition):
    for element in elementalComposition.keys():
        if element not in elements:
            raise Exception("Species "+name+" with elemental composition "+str(elementalComposition)+" has an unknown element, "+element)

    if name not in speciesData:
        speciesData[name] = SpeciesDataType(name, elementalComposition)
