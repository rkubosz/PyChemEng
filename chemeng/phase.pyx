#!/usr/bin/env python
from chemeng.components cimport Components
from chemeng.speciesdata import R, T0, P0, speciesData
#Ensure that the thermodata is loaded, so the data is available
import chemeng.thermodata
import math

####################################################################
# Phase base class
####################################################################
#This class implements all of the helper functions and data types
#common to a single phase.
cdef class Phase:
    """A base class which holds fundamental methods and members of a single phase which may contain multiple components"""

    def __init__(self, T, components, P, phase):
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
        for entry in self.components._list:
            sum += entry.second * speciesData[entry.first].Cp0(self.T, phase=self.phase)
        return sum

    cpdef double enthalpy(self) except + : # J
        """A calculation of the enthalpy (H) of the Phase"""
        cdef double sum = 0.0
        for entry in self.components._list:
            sum += entry.second * speciesData[entry.first].Hf0(self.T, phase=self.phase)
        return sum

    cpdef double entropy(self) except + : # J / K
        """A calculation of the entropy S of the phase."""
        cdef double total = self.components.total()
        #Individual component entropy
        cdef double sumEntropy = 0.0
        for entry in self.components._list:
            sumEntropy += entry.second * speciesData[entry.first].S0(self.T, phase=self.phase)

        #Mixing entropy
        for entry in self.components._list:
            if entry.second > 0.0:
                sumEntropy -= R * entry.second * math.log(entry.second / total)
        return sumEntropy

    cpdef double gibbsFreeEnergy(self): #J
        """A calculation of the Gibbs free energy (G) of the phase"""
        return self.enthalpy() - self.T * self.entropy()

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
    
#Functions which must be overridden by derived classes
    cpdef double internalEnergy(self):
        raise Exception("Function missing from "+self.__class__.__name__+"!")

    cpdef double volume(self):
        raise Exception("Function missing from "+self.__class__.__name__+"!")

####################################################################
# Ideal gas class
####################################################################
cdef class IdealGasStream(Phase):
    def __init__(self, T, components, P):
        Phase.__init__(self, T, components, P, 0)

    cpdef IdealGasStream copy(self):
        cdef IdealGasStream retval = IdealGasStream.__new__(IdealGasStream)
        retval.T = self.T
        retval.P = self.P
        retval.phase = 0
        retval.components = self.components.copy()        
        return retval

    #As the enthalpy of an ideal gas is constant with pressure, we do
    #not need to override the base class definition
    #def enthalpy(self):

    cpdef double internalEnergy(self):# Units are J
        cdef double PV = self.components.total() * R * self.T
        return self.enthalpy() - PV

    cpdef double entropy(self): # J / K
        """A calculation of the entropy (S) of the phase"""
        #We need to include the effect of pressure on the entropy
        cdef double pressureEntropy =  - R * self.components.total() * math.log(self.P / P0)
        return Phase.entropy(self) + pressureEntropy

    cpdef double helmholtzFreeEnergy(self): #J
        """A calculation of the Gibbs free energy (G) of the Stream at a temperature T. If T is not given, then it uses the stream temperature"""
        PV = self.components.total() * R * self.T
        return self.gibbsFreeEnergy() - PV

    cpdef double volume(self):
        return self.components.total() * R * self.T / self.P

    def __add__(self, IdealGasStream other):
        """An operator to allow mixing of phases (with calculation of exit temperature)"""
        cdef IdealGasStream output = IdealGasStream((self.T + other.T) / 2.0, {}, P=min(self.P, other.P))
        output.components = self.components + other.components
        output.setEnthalpy(self.enthalpy() + other.enthalpy())
        return output

####################################################################
# Incompressible Solid
####################################################################
class IncompressibleSolid(Phase):
    def __init__(self, T, components, P, molardensity, phaseID):
        Phase.__init__(self, T, components, P, phaseID)
        self.molardensity = molardensity

    #As the enthalpy of an ideal gas is constant with pressure, we do
    #not need to override the base class definition
    #def enthalpy(self):

    def internalEnergy(self):# Units are J
        return self.enthalpy() - self.P * self.volume()

    def helmholtzFreeEnergy(self): #J
        """A calculation of the Gibbs free energy (G) of the Stream at a temperature T. If T is not given, then it uses the stream temperature"""
        return self.gibbsFreeEnergy() - self.P * self.volume()

    def volume(self):
        return self.components.total() / self.molardensity

    def __add__(self, other):
        """An operator to allow mixing of phases (with calculation of exit temperature)"""
        if type(self) != type(other):
            raise Exception("Cannot mix two phases of different types, "+str(self)+"+"+str(other))
        #Make an output stream with the correct concentration
        output = self.__class__((self.T + other.T) / 2.0, self.components + other.components, P=min(self.P, other.P))
        output.molardensity = output.components.total() / (self.volume() + other.volume())
        output.setEnthalpy(self.enthalpy() + other.enthalpy())
        return output
