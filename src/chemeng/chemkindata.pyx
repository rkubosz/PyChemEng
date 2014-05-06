#!/usr/bin/env python
import math
#Ensure that the elemental database has been loaded first
import chemeng.elementdata
from chemeng.components import Components
from chemeng.speciesdata import speciesData, registerSpecies, registerCpFitFunction, relativeError
from chemeng.speciesdata cimport SpeciesDataType

####################################################################
# Physical constants
####################################################################
cdef public double R = 8.31451
T0 = 273.15 + 25.0
P0 = 1.0e5

registerCpFitFunction("ChemKin",
                      #The function (f)
                      lambda T, C : C[0] + C[1] * T + C[2] * T **2 + C[3] + T,
                      #The integrated function (f) (without the integration constant)
                      lambda T, C : R * sum([constant * math.log(T) if order == -1 else constant * (T ** (order+1)) / (order+1) for constant, order in C]),
                      #The integrated function over temperature (f/T)
                      lambda T, C : R * sum([constant * math.log(T) if order == 0 else constant * (T ** (order)) / (order) for constant, order in C]))

cpdef parseCHEMKINDataFile(filename, quiet=True):
    cdef SpeciesDataType sp
    file = open(filename, "r")
    lineit = iter(file)

    #First line is THERMO or THERMO ALL
    line = lineit.next()
    if line.strip()[:6] != "THERMO":
        raise Exception("Does not look like a CHEMKIN database")

    #Next line is temperature range and common temperature T
    line = lineit.next()
    commonT = float(line[9:19])
    #Begin parsing records
    linecount=0
    while True:
        #Keep reading from the file until we run out of lines
        try:
            linecount += 1
            line = lineit.next()
        except StopIteration:
            break

        if line.strip() == "END":
            break

        species = line[0:18].split()[0].strip()
        comments = line[18:24].strip()
        if not quiet:
            print ">>>Parsing:'"+species+"' '"+comments+"'"

        ########################################
        try:
            linecount += 1
            line = lineit.next()
        except StopIteration:
            break

        if not quiet:
            print "  Chemical Formula:",line    
        MolecularFormula = {}
        for entry in range(4):
            shift = 5 * entry
            atom = line[24+shift:26+shift]
            natom= int(line[26+shift:28+shift])
            MolecularFormula[atom] = natom

        phase=line[44]
        if phase=="G":
            phase = "Gas"
        elif phase=="L":
            phase = "Liquid"
        elif phase=="S":
            phase = "Solid"
        
        lowT = float(line[45:55])
        highT = float(line[55:65])
        midT = commonT
        if line[65:73].strip() != "":
            midT = float(line[65:73])

        try:
            linecount += 1
            line = lineit.next()
        except StopIteration:
            break
        
        C=[]
        for i in range(5):
            C.append(float(line[i*15:(i+1)*15]))

        try:
            linecount += 1
            line = lineit.next()
        except StopIteration:
            break

        for i in range(5):
            C.append(float(line[i*15:(i+1)*15]))

        try:
            linecount += 1
            line = lineit.next()
        except StopIteration:
            break

        for i in range(4):
            C.append(float(line[i*15:(i+1)*15]))

        coeffs=[]
        coeffs.append([lowT, midT, "ChemKin", C[], HConst, Sconst])
            
            phasename = "Gas"
            if species[-1] == ")":
                #This has a phase qualifier at the end in parentheses, grab it
                import re
                m = re.match('(.*?)\(([^()]*?)\)$', species)
                species = m.group(1)
                phasename = m.group(2)
                if phasename == "L":
                    phasename = "Liquid"
                else:
                    phasename = str(phase)+phasename
            registerSpecies(species, MolecularFormula, MW)

            sp = speciesData[species]

            sp.registerPhase(phasename, comments=comments)
            for C in coeffs:
                sp.registerPhaseCoeffs(C, phasename)
                    
            if quiet:
                continue
            
            if sp.inDataRange(T0, phasename):
                HfCalc = sp.Hf0(T0, phasename)
                error = relativeError(HfCalc, HfRef)
                if error > HfMaxError:
                    print "Warning: Species \""+species+"\" and phase "+phasename+" in file:"+filename+" at line "+str(linecount)+" has a Hf of "+str(HfRef)+" but a calculated value of "+str(HfCalc)+" a relative error of "+str(error)+" for Hf at "+str(T0)

        except Exception as e:
            if not quiet:
                import traceback
                print traceback.print_exc()
                print "Error parsing record for",species,"in file:",filename,"at line",linecount
                print "   ",e.message
                raise

def initDataDir(directory):
    import os
    parseNASADataFile(os.path.join(directory, 'NASA_CEA.inp'), quiet=True)
    print "Loaded NASA dataset (database at",len(speciesData),"species)"

import sys
import os.path
initDataDir('/usr/local/PyChemEng/data')
