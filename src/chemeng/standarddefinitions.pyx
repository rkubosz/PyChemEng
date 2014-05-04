#!/usr/bin/env python
from chemeng.components import Components

#Composition of dry air (to 1ppmv)
DryAir = Components({"N2":78.084, "O2":20.946,  "Ar":0.934, "CO2":0.0397, "Ne":0.001818, "He": 0.000524, "CH4":0.000179, "Kr":0.000114})

#Definitions taken from the GasEq program
StandardHydrocarbonCombustionComponents = {"N2", "H2O", "CO2", "CO", "O2", "OH", "H", "O", "H2", "NO"}
ExtendedHydrocarbonCombustionComponents = StandardHydrocarbonCombustionComponents.union({"HCO", "CH2O", "CH4", "CH3", "HO2", "NO2", "NH3", "NH2", "N", "HCN", "CN", "N2O", "C2", "CH"})
RichSootingProducts = ExtendedHydrocarbonCombustionComponents.union({"C2H4", "C2H2", "NH", "CH", "CH2", "C(S)"})
