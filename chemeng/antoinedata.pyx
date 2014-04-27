#!/usr/bin/env python
import math
import os
from chemeng.components import Components
from chemeng.speciesdata import T0, speciesData, registerFitFunction

####################################################################
# Antoine polynomials
####################################################################
#Functions should take T as Kelvin and return P as Pascals
registerFitFunction("Antoine", lambda T, C : 1e5 * (10.0 ** (C[0] - C[1] / (T + C[2]))))
registerFitFunction("AntoineCnR", lambda T, C : 133.2895 * math.exp(C[0] - C[1] / (T + C[2])))

####################################################################
# Datafile loading
####################################################################
def parseTemperature(string):
    if string[-1] == 'C':
        return float(string[:-1])+273.15
    return float(string)

def initDataDir(directory):
    datafile = open(os.path.join(directory, 'antoine.inp'), 'r')
    for line in datafile:
        datafields = line.split()
    
        #Skip comments or other problems
        if line[0] == "#" or len(datafields) == 0:
            continue
    
        speciesData[datafields[0]].registerAntoineData(Tmin = parseTemperature(datafields[1]),
                                                       Tmax = parseTemperature(datafields[2]),
                                                       fitFunction = datafields[3], 
                                                       constants = map(float, datafields[4:]))

#Examples of manual loading
#Taken from the NIST webbook
#speciesData["CO2"].registerAntoineData(292.77, 366.63, "Antoine", [5.24677, 1598.673, -46.424])
#speciesData["CO2"].registerAntoineData(273.0, 351.7, "Antoine", [5.37229, 1670.409, -40.191])
#speciesData["CO2"].registerAntoineData(364.8, 513.91, "Antoine", [4.92531, 1432.526, -61.819])

import sys
import os.path
initDataDir('/usr/local/PyChemEng/data')
