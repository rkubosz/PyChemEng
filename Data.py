#!/usr/bin/env python
####################################################################
# Physical constants
####################################################################
#Constant used in the NASA data set! (If changed, will need to rescale the data)
R = 8.31451 

#STP is 25 celcius at 1 bar
T0 = 273.15 + 25
P0 = 1.0e5

####################################################################
# Species data
####################################################################
#Species data is a dictionary of thermodynamic data on different
#species/components. This is addressed using the chemical formula of
#the species to be studied.

speciesData={}

class PhaseData:
    def Cp0(T):
        pass

    def Hf0(T):
        pass

    def S0(T):
        pass

    def Gibbs0(T):#AKA chemical potential
        return Hf0(T) - T * S0(T)
    
class IdeaPhaseData(PhaseData):
