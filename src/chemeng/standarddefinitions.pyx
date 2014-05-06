#!/usr/bin/env python
from chemeng.components import Components

#Composition of dry air (to 1ppmv)
DryAir = Components({"N2":78.084, "O2":20.946,  "Ar":0.934, "CO2":0.0397, "Ne":0.001818, "He": 0.000524, "CH4":0.000179, "Kr":0.000114})

#Definitions taken from the GasEq program
StandardHydrocarbonCombustionComponents = Components({"N2":0, "H2O":0, "CO2":0, "CO":0, "O2":0, "OH":0, "H":0, "O":0, "H2":0, "NO":0})
ExtendedHydrocarbonCombustionComponents = StandardHydrocarbonCombustionComponents + Components({"HCO":0, "CH4":0, "CH3":0, "HO2":0, "NO2":0, "NH3":0, "NH2":0, "N":0, "HCN":0, "CN":0, "N2O":0, "C2":0, "CH":0}) #Removed until the Burcat DB is in "CH2O":0
RichSootingProducts = ExtendedHydrocarbonCombustionComponents + Components({"C2H4":0, "C2H2":0, "NH":0, "CH":0, "CH2":0, "C(S)":0})
