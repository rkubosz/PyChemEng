#!/usr/bin/env python

#We need to import ThermoData to initialise the thermodynamic database
import ThermoData

##################################################################
# Components
##################################################################
from Components import Components

#Components are used in mass balances. They are dictionaries of
#species and their molar amounts/flows.
Air = Components({"N2":78.084, "O2":20.9476,  "Ar":0.934, "CO2":0.0314, "Ne":0.001818, "He": 0.000524, "CH4":0.0002})

#Each species can be accessed by its formula
molesO2 = Air["O2"]

#We can mix and multiply Components as well
Gas=Components({"CH4":1})
StoichiometricFuel = Air + 0.5 * molesO2 * Gas

#We can also normalise the compositions using the total molar count
StoichiometricFuel /= StoichiometricFuel.total()

#We can get the elemental composition of the components, which is also
#a Components class.
print StoichiometricFuel.elementalComposition()
#C{'C':0.0950944, 'Ar':0.00845453, 'H':0.379241, 'He':4.74323e-06, 'Ne':1.64565e-05, 'O':0.379802, 'N':1.41363}

##################################################################
# Streams
##################################################################
from Stream import IdealGasStream

#Streams are components with an associated temperature and pressure
fuelstream = IdealGasStream(298.15, StoichiometricFuel, P=1.01325e5)

#Once these are fixed, we can calculate a range of parameters:
print fuelstream.Cp(), fuelstream.enthalpy(), fuelstream.entropy(), fuelstream.gibbsFreeEnergy(), fuelstream.helmholtzFreeEnergy()
#29.7267355968 -7184.68907772 200.139970908 -66856.4214039 -69335.3925604 668.375700264

#If we add together streams, the output has the lowest pressure of the
#input streams, and enthalpy conservation is used to calculate the
#exit temperature
outletstream = fuelstream + IdealGasStream(800, {"CO2":1, "H2O":1}, P=1.01325e5)
print outletstream.T
#668.375700264

##################################################################
# Reaction
##################################################################
import Reaction

#We can calculate reaction equilibria just by specifying an input
#stream, any possible product species and if it is either at constant
#Pressure or Volume and if it is either at constant Temperature or
#Adiabatic.
OutStream = Reaction.react(fuelstream, {"N2", "H2O", "CO2", "CO", "O2", "OH", "H2", "O", "H"}, constP=True, constT=False)

print OutStream
#<1.00736 mol/s, 2230.04 K, 1.01325 bar, C{'CO2':0.0863113, 'CO':0.00878318, 'OH':0.00334971, 'H2':0.00345241, 'H':0.000384071, 'Ne':1.64565e-05, 'O':0.000234096, 'H2O':0.184301, 'Ar':0.00845453, 'CH4':0, 'N2':0.706813, 'O2':0.00525572, 'He':4.74323e-06}>
