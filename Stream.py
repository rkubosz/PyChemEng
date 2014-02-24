#!/usr/bin/env python
from Components import Components
from Data import R, T0, P0, speciesData
#Ensure that the thermodata is loaded
import ThermoData
import math

####################################################################
# Phase base class
####################################################################
#This class implements all of the helper functions and data types
#common to a single phase.
class Phase(object):
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

    def __add__(self, other):
        """An operator to allow mixing of phases (with calculation of exit temperature)"""
        if type(self) != type(other):
            raise Exception("Cannot mix two phases of different types, "+str(self)+"+"+str(other))
        #Make an output stream with the correct concentration
        output = self.__class__((self.T + other.T) / 2.0, self.components + other.components, P=min(self.P, other.P))
        output.setEnthalpy(self.enthalpy() + other.enthalpy())
        return output

    def __str__(self):
        return "<"+self.__class__.__name__+", %g mol, %g K, %g bar, " % (self.components.total(), self.T, self.P / 1e5) + str(self.components)+">"

#The thermodynamic functions below do not include any effect of
#pressure (they assume the phase is at the reference pressure P0), and
#assumes that the phase is an ideal mixture. Classes which derive from
#Phase can override these methods to add corrections to these
#assumptions.
    def Cp(self): # Units are J
        """A calculation of the isobaric heat capacity (Cp) of the phase"""
        return sum([flow * speciesData[x].Cp0(self.T, phase=self.phase) for x, flow in self.components.iteritems()])

    def enthalpy(self): # J
        """A calculation of the enthalpy (H) of the Phase"""
        return sum([flow * speciesData[x].Hf0(self.T, phase=self.phase) for x, flow in self.components.iteritems()])

    def entropy(self): # J / K
        """A calculation of the entropy S of the phase."""
        total = self.components.total()
        #Individual component entropy
        componentEntropy = sum(flow * speciesData[x].S0(self.T, phase=0) for x, flow in self.components.iteritems())
        #Assuming ideal mixing entropy!
        mixingEntropy = - R * sum([flow * math.log((flow + 0.0) / total) for x, flow in self.components.iteritems() if flow != 0])
        return componentEntropy + mixingEntropy

    def gibbsFreeEnergy(self): #J
        """A calculation of the Gibbs free energy (G) of the phase"""
        return self.enthalpy() - self.T * self.entropy()
#Helper functions to manipulate the thermodynamic properties of a phase
    def setInternalEnergy(self, U):
        import scipy
        import scipy.optimize
        def worker(T):
            self.T = T
            return self.internalEnergy() - U
        self.T = scipy.optimize.fsolve(worker, self.T+0.0)[0]
    def setEnthalpy(self, enthalpy):
        import scipy
        import scipy.optimize
        def worker(T):
            self.T = T
            return self.enthalpy() - enthalpy
        self.T = scipy.optimize.fsolve(worker, self.T)[0]
    
#Functions which must be overridden by derived classes
    def internalEnergy(self):
        raise Exception("Function missing from "+self.__class__.__name__+"!")

    def internalEnergy(self):
        raise Exception("Function missing from "+self.__class__.__name__+"!")

    def volume(self):
        raise Exception("Function missing from "+self.__class__.__name__+"!")

####################################################################
# Ideal gas class
####################################################################
class IdealGasStream(Phase):
    def __init__(self, T, components, P):
        super(IdealGasStream, self).__init__(T, components, P, 0)

    #As the enthalpy of an ideal gas is constant with pressure, we do
    #not need to override the base class definition
    #def enthalpy(self):

    def internalEnergy(self):# Units are J
        PV = self.components.total() * R * self.T
        return self.enthalpy() - PV

    def entropy(self): # J / K
        """A calculation of the entropy (S) of the phase"""
        #We need to include the effect of pressure on the entropy
        pressureEntropy =  - R * self.components.total() * math.log(self.P / P0)
        return super(IdealGasStream, self).entropy() + pressureEntropy

    def helmholtzFreeEnergy(self): #J
        """A calculation of the Gibbs free energy (G) of the Stream at a temperature T. If T is not given, then it uses the stream temperature"""
        PV = self.components.total() * R * self.T
        return self.gibbsFreeEnergy() - PV

    def volume(self):
        return self.components.total() * R * self.T / self.P

####################################################################
# Incompressible Solid
####################################################################
class IncompressibleSolid(Phase):
    def __init__(self, T, components, P, molardensity, phaseID):
        super(IncompressibleSolid, self).__init__(T, components, P, phaseID)
        self.molardensity = molardensity

    #As the enthalpy of an ideal gas is constant with pressure, we do
    #not need to override the base class definition
    #def enthalpy(self):

    def internalEnergy(self):# Units are J
        return self.enthalpy() - P * self.volume()

    def entropy(self): # J / K
        """A calculation of the entropy (S) of the phase"""
        return super(IdealGasStream, self).entropy()

    def helmholtzFreeEnergy(self): #J
        """A calculation of the Gibbs free energy (G) of the Stream at a temperature T. If T is not given, then it uses the stream temperature"""
        return self.gibbsFreeEnergy() - P * self.volume()

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
