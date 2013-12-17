#!/usr/bin/env python
import math
from Components import Components
import ThermoData 
from Data import T0, speciesData, registerFitFunction

####################################################################
# Antoine polynomials
####################################################################
#Functions should take T as Kelvin and return P as Pascals
registerFitFunction("Antoine", lambda T, C : 1e5 * 10.0 ** (C[0] - C[1] / (T + C[2])))


####################################################################
# Datafile loading
####################################################################
datafile = open(os.path.join(os.path.dirname(__file__), 'datafiles/antoine.inp'), 'r')
for line in datafile:
    datafields = line.split()

    #Skip comments or other problems
    if line[0] == "#" or len(datafields) == 0:
        continue
    speciesData[datafields[0]].registerAntoineData(Tmin = float(datafields[1]),
                                                   Tmax = float(datafields[2]),
                                                   fitFunction = datafields[2], 
                                                   constants = map(float, datafields[3:]))

#Examples of manual loading
#Taken from the NIST webbook
#speciesData["CO2"].registerAntoineData(292.77, 366.63, "Antoine", [5.24677, 1598.673, -46.424])
#speciesData["CO2"].registerAntoineData(273.0, 351.7, "Antoine", [5.37229, 1670.409, -40.191])
#speciesData["CO2"].registerAntoineData(364.8, 513.91, "Antoine", [4.92531, 1432.526, -61.819])
