#!/usr/bin/env python
from Components import Components
from ThermoData import Cp, Hf, S
from Data import R, T0, P0
import math
	  
####################################################################
# Stream class
####################################################################
class IdealGasStream():
    """A class which represents an ideal gas stream of Components and
    their temperature"""
    components = Components()
    T = 0
    P = 0
    def __init__(self, T, components, P=1.01325e5):
        """The constructor for a stream"""
        self.T = T
        self.components = Components(components)
        self.P = P

    def Cp(self, T=None): # J
        """A calculation of the heat capacity (Cp) of the Stream at a temperature T. If T is not given, then it uses the stream temperature"""
        if T == None:
            T = self.T
        return sum([flow * Cp(x, T) for x, flow in self.components.iteritems()])

    def internalEnergy(self, T=None):# J
        if T == None:
            T = self.T
        PV = self.components.total() * R * T
        return self.enthalpy(T) - PV

    def setInternalEnergy(self, U): # J
        import scipy
        import scipy.optimize
        self.T = scipy.optimize.fsolve(lambda T : self.internalEnergy(T) - U, self.T)[0]

    def enthalpy(self, T=None): # J
        """A calculation of the enthalpy (H) of the Stream at a temperature T. If T is not given, then it uses the stream temperature"""
        if T == None:
            T = self.T
        return sum([flow * Hf(x, T) for x, flow in self.components.iteritems()])

    def setEnthalpy(self, enthalpy): # J
        import scipy
        import scipy.optimize
        self.T = scipy.optimize.fsolve(lambda T : self.enthalpy(T) - enthalpy, self.T)[0]

    def entropy(self, T=None): # J / K
        """A calculation of the entropy (S) of the Stream at a temperature T. If T is not given, then it uses the stream temperature"""
        if T == None:
            T = self.T
        total = self.components.total()
        #Individual component entropy
        componentEntropy = sum(flow * S(x, T) for x, flow in self.components.iteritems())
        #Assuming ideal mixing entropy!
        mixingEntropy = - R * sum([flow * math.log((flow + 0.0) / total) for x, flow in self.components.iteritems() if flow != 0])
        pressureEntropy =  - R * total * math.log(self.P / P0)
        return componentEntropy + mixingEntropy + pressureEntropy

    def gibbsFreeEnergy(self, T=None): #J
        """A calculation of the Gibbs free energy (G) of the Stream at a temperature T. If T is not given, then it uses the stream temperature"""
        if T == None:
            T = self.T
        return self.enthalpy(T) - T * self.entropy(T)

    def helmholtzFreeEnergy(self, T=None): #J
        """A calculation of the Gibbs free energy (G) of the Stream at a temperature T. If T is not given, then it uses the stream temperature"""
        if T == None:
            T = self.T
        PV = self.components.total() * R * T
        return self.gibbsFreeEnergy(T) - PV

    def volume(self):
        return self.components.total() * R * self.T / self.P
                
    def __add__(self, other):
        """An operator to allow mixing streams (with calculation of exit temperature)"""
        #Make an output stream with the correct concentration
        output = IdealGasStream((self.T + other.T) / 2.0, self.components + other.components, P=min(self.P, other.P))
        output.setEnthalpy(self.enthalpy() + other.enthalpy())
        return output

    def __str__(self):
        return "<%g mol/s, %g K, %g bar, " % (self.components.total(), self.T, self.P / 1e5) + str(self.components)+">"
