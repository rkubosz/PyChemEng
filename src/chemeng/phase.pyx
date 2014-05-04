#!/usr/bin/env python
# distutils: language = c++
# cython: profile=True

from chemeng.components cimport Components
from chemeng.speciesdata import speciesData
from chemeng.speciesdata cimport SpeciesDataType
#Ensure that the thermodata is loaded, so the data is available
from chemeng.thermodata import T0, P0
import math

cdef public double R = 8.31451

####################################################################
# Phase base class
####################################################################
#This class implements all of the helper functions and data types
#common to a single phase.
cdef class Phase:
    """A base class which holds fundamental methods and members of a single phase which may contain multiple components"""

    def __init__(self, components, T, P, phase):
        """The constructor for a stream"""
        #The temperature of the phase
        self.T = T
        #The component species of the phase
        self.components = Components(components)
        #The pressure of the phase
        self.P = P
        #An identifier used to fetch thermodynamic properties for the phase
        self.phase = phase

    def __str__(self):
        return "<"+self.__class__.__name__+", %g mol, %g K, %g bar, " % (self.components.total(), self.T, self.P / 1e5) + str(self.components)+">"

#The thermodynamic functions below do not include any effect of
#pressure (they assume the phase is at the reference pressure P0), and
#assumes that the phase is an ideal mixture. Classes which derive from
#Phase can override these methods to add corrections to these
#assumptions.
    cpdef double Cp(self) except + : # Units are J
        """A calculation of the isobaric heat capacity (Cp) of the phase"""
        cdef double sum = 0.0
        cdef SpeciesDataType sp
        for entry in self.components._list:
            sp = speciesData[entry.first]
            sum += entry.second * sp.Cp0(self.T, phase=self.phase)
        return sum

    cpdef double enthalpy(self) except + : # J
        """A calculation of the enthalpy (H) of the Phase"""
        cdef double sum = 0.0
        cdef SpeciesDataType sp
        for entry in self.components._list:
            sp = speciesData[entry.first]
            sum += entry.second * sp.Hf0(self.T, phase=self.phase)
        return sum

    cpdef double entropy(Phase self) except + : # J / K
        """A calculation of the entropy S of the phase."""
        cdef double total = self.components.total()
        #Individual component entropy
        cdef double sumEntropy = 0.0
        cdef SpeciesDataType sp
        for entry in self.components._list:
            sp = speciesData[entry.first]
            sumEntropy += entry.second * sp.S0(self.T, phase=self.phase)

        #Mixing entropy
        for entry in self.components._list:
            if entry.second > 0.0:
                sumEntropy -= R * entry.second * math.log(entry.second / total)
        return sumEntropy

    cpdef double gibbsFreeEnergy(Phase self): #J
        """A calculation of the Gibbs free energy (G) of the phase"""
        return self.enthalpy() - self.T * self.entropy()

    cpdef double internalEnergy(Phase self):# Units are J
        return self.enthalpy() - self.P * self.volume()

    cpdef double helmholtzFreeEnergy(Phase self): #J
        """A calculation of the Gibbs free energy (G) of the Stream at a temperature T. If T is not given, then it uses the stream temperature"""
        return self.gibbsFreeEnergy() - self.P * self.volume()

#Helper functions to manipulate the thermodynamic properties of a phase
    def setInternalEnergy(self, double U):
        import scipy
        import scipy.optimize
        def worker(double T):
            self.T = T
            return self.internalEnergy() - U
        self.T = scipy.optimize.fsolve(worker, self.T+0.0)[0]

    def setEnthalpy(self, double enthalpy):
        import scipy
        import scipy.optimize
        def worker(double T):
            self.T = T
            return self.enthalpy() - enthalpy
        self.T = scipy.optimize.fsolve(worker, self.T)[0]

    cpdef list createStateVector(Phase self):
        return [moles for moles in self.components.values()] + [self.T, self.P]

    cpdef restoreStateVector(Phase self, list vec):
        for key, val in zip(self.components.keys(), vec[:-2]):
            self.components[key] = val
        self.T = vec[-2]
        self.P = vec[-1]

#Functions which must be overridden by derived classes
    cpdef double volume(self):
        raise Exception("Function missing from "+self.__class__.__name__+"!")

####################################################################
# Ideal gas class
####################################################################
cdef class IdealGasPhase(Phase):
    def __init__(self, components, T, P):
        Phase.__init__(self, components, T=T, P=P, phase="Gas")

    cpdef IdealGasPhase copy(IdealGasPhase self):
        cdef IdealGasPhase retval = IdealGasPhase.__new__(IdealGasPhase)
        retval.T = self.T
        retval.P = self.P
        retval.phase = self.phase
        retval.components = self.components.copy()        
        return retval

    #As the enthalpy of an ideal gas is constant with pressure, we do
    #not need to override the base class definition

    cpdef double entropy(self): # J / K
        cdef double pressureEntropy =  - R * self.components.total() * math.log(self.P / P0)
        return Phase.entropy(self) + pressureEntropy

    cpdef double volume(self):
        return self.components.total() * R * self.T / self.P

    cpdef double Cv(self) except + :
        return self.Cp() - R * self.total()

    def __add__(self, IdealGasPhase other):
        cdef IdealGasPhase output = self.copy()
        output.components += other.components
        output.setEnthalpy(self.enthalpy() + other.enthalpy())
        return output

####################################################################
# Incompressible Phase
####################################################################
cdef class IncompressiblePhase(Phase):
    def __init__(self, components, phaseID, T, P, molarvolume=0.0):
        Phase.__init__(self, components, T=T, P=P, phase=phaseID)
        self.molarvolume = molarvolume

    cpdef IncompressiblePhase copy(IncompressiblePhase self):
        cdef IncompressiblePhase retval = IncompressiblePhase.__new__(IncompressiblePhase)
        retval.T = self.T
        retval.P = self.P
        retval.phase = self.phase
        retval.components = self.components.copy()
        retval.molarvolume = self.molarvolume
        return retval

    #Entropy is constant with pressure for incompressible phases, so keep Phase definition

    cpdef double enthalpy(IncompressiblePhase self) except +:
        cdef double base = Phase.enthalpy(self)
        return Phase.enthalpy(self) + self.volume() * (self.P - P0)
    
    cpdef double volume(IncompressiblePhase self):
        return self.components.total()  * self.molarvolume

    cpdef double Cv(self) except + :
        return self.Cp()

    def __add__(self, IncompressiblePhase other):
        cdef IncompressiblePhase output = self.copy()
        output.components = output.components + other.components
        output.setEnthalpy(self.enthalpy() + other.enthalpy())
        return output
