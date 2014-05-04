#!/usr/bin/python

from chemeng import *

print elements['Ca']
#Element{Ca, Z=20, AW=40.078, 24 isotopes}
print elements['Li'].mass
#6.941 (g/mol)

print elements['C'].isotopes[6]
#Isotope{C, Z=6, N=6, M=12.0(0.0), P=0.9893}

print elements[(1,2)] #Tritium
#Element{T, Z=1, AW=3.0160492777}

air = Components({'O2':0.21, 'N2':0.79})
print air
#C{'N2':0.79, 'O2':0.21}

#Specify fuel and air streams
fuel = Components({'CH4':1.0})
air = Components({'N2':0.79, 'O2':0.21})

#Determine the required O2 for combustion
elementalfuel = fuel.elementalComposition()
requiredO2 = elementalfuel['C'] + elementalfuel['H'] / 4

#Create a stochiometric mixture
stochiometricMix = fuel + air * (requiredO2 / air['O2'])

print stochiometricMix 
#C{'CH4':1, 'N2':7.52381, 'O2':2}

print stochiometricMix.normalised() 
#C{'CH4':0.0950226, 'N2':0.714932, 'O2':0.190045}

print stochiometricMix.total() 
#10.5238095238 (mol)

print stochiometricMix.totalMass() 
#290.807545714 (g)

print stochiometricMix.elementalComposition()
#C{'C':1, 'H':4, 'N':15.0476, 'O':4}

print speciesData['CO2']
#Species{CO2, phases=[Gas], elementalComposition=C{'C':1, 'O':2}}

print speciesData['CO2'].phases['Gas']
#Phase{Gas, T=[200.0-&gt;1000.0K, 1000.0-&gt;6000.0K, 6000.0-&gt;20000.0K], comments='Gurvich,1991 pt1 p27 pt2 p24.'}

print findSpeciesData("C8")
#[Species{'C8H18,isooctane', phases=[Gas], elementalComposition=C{'C':8, 'H':18}},
# Species{'C8H18,n-octane', phases=[Gas], elementalComposition=C{'C':8, 'H':18}},
# Species{'C8H8,styrene', phases=[Gas], elementalComposition=C{'C':8, 'H':8}},
# Species{'C8H16,1-octene', phases=[Gas], elementalComposition=C{'C':8, 'H':16}},
# Species{'C8H17,n-octyl', phases=[Gas], elementalComposition=C{'C':8, 'H':17}},
# Species{'C8H10,ethylbenz', phases=[Gas], elementalComposition=C{'C':8, 'H':10}}]


print speciesData['H2O'].Cp0(298.15, 'Gas')#Heat capacity
#33.5877103224 (J/(mol K))

print speciesData['H2O'].Hf0(298.15, 'Gas')#Enthalpy of formation
#-241826.00034 (J/mol)

print speciesData['SO'].Hf0(298.15, 'Gas') #Entropy
#188.829115517 (J/(mol K))


vapour=IdealGasPhase({'H2O':1.0}, T=179.9+273.15, P=10.e5)
print vapour.Cp()
#34.7503667276 (J/ K)
print vapour.enthalpy()
#-236543.780767 (J)

print speciesData['H2O']
#Species{'H2O', phases=[1cr, Gas, Liquid], elementalComposition=C{'H':2, 'O':1}}
liquid=IncompressiblePhase({'H2O':1.0}, T=179.9+273.15, P=10.e5, phaseID="Liquid", molarvolume=0.018 / 998.0)
print vapour.enthalpy() - liquid.enthalpy()
#37427.074787 (J/mol)
print (vapour.enthalpy() - liquid.enthalpy()) / 18.0
#2079.28193261 (J/g)

liquid=IncompressiblePhase({'H2O':1.0}, T=179.9+273.15, P=10.e5, phaseID="Liquid")
print (vapour.enthalpy() - liquid.enthalpy()) / 18.0
#2080.18373622 (J/g)
