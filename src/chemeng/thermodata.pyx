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
#Constants used in the NASA data set! (If changed, will need to
#rescale the data)
cdef public double R = 8.31451
T0 = 273.15 + 25.0
P0 = 1.0e5

##Add the NASA polynomial
registerCpFitFunction("Poly",
                      #The function (f)
                      lambda T, C : R * sum([constant * (T ** order) for constant, order in C]),
                      #The integrated function (f) (without the integration constant)
                      lambda T, C : R * sum([constant * math.log(T) if order == -1 else constant * (T ** (order+1)) / (order+1) for constant, order in C]),
                      #The integrated function over temperature (f/T)
                      lambda T, C : R * sum([constant * math.log(T) if order == 0 else constant * (T ** (order)) / (order) for constant, order in C]))

####################################################################
# Test functions
####################################################################
HfMaxError=0.01 #1% error
MWMaxError=0.0003 #0.03% error

####################################################################
# NASA Glenn Thermodynamic Database
####################################################################
cpdef parseFortanFloat(string):
    #First, try to parse the exponentiated formats
    data = string.split('D')
    if len(data) != 2:
        data = string.split('E')
    if len(data) != 2:
        data = string.split('d')
    if len(data) == 2:
        return float(data[0].strip() + "e" + data[1].strip())
    #Try to parse it as a standard number
    return float(string.strip())

cpdef parseNASADataFile(filename, quiet=True):
    cdef SpeciesDataType sp     

    file = open(filename, "r")
    lineit = iter(file)
    commentlines=2
    linecount=0
    while True:
        linecount += 1
        #Keep reading from the file until we run out of lines
        try:
            line = lineit.next()
        except StopIteration:
            break

        #Skip comments
        if line[0] == "!": continue
        
        #Skip first two record lines
        if commentlines > 0:
            commentlines -= 1
            continue
        
        #Skip blank lines
        if line.strip() == "":
            continue
        
        #Find the exit conditions
        if line[0:12] == "END PRODUCTS":
            break
        
        #component=line.split(" ",1)[0]
        startline=linecount

        ######NEED A SMART COMPONENT PARSER#####
        firstline = line.strip()
        ########################################
        linecount += 1
        line = lineit.next()
        intervals=int(line[1])
        try:
            refdatecode=line[3:9].strip()
            ####Parse Chemical formula?
            phase=int(line[50:52])
            MW=parseFortanFloat(line[52:65])#g/mol
            HfRef=parseFortanFloat(line[65:80])# @298.15K , in J/mol
            coeffs = []
            MolecularFormula = {}
            for offset in range(5):
                startatom=10+8*offset
                endatom=startatom+3
                startnumber=13+8*offset
                endnumber=startnumber+5
                atom = line[startatom:endatom].strip()
                if len(atom) == 1 and atom[0] == "E":
                    atom = atom[0].lower()+"-"
                if len(atom)==2:
                    atom = atom[0]+atom[1].lower()
                if len(atom)==3:
                    atom = atom[0]+atom[1].lower()+atom[2].lower()
                natom = float(line[startnumber:endnumber].strip())
                if natom != 0:
                    MolecularFormula[atom] = natom
            for entry in range(intervals):
                linecount += 1
                line = lineit.next()
                Tmin = float(line[1:12])
                Tmax = float(line[12:22])

                Ncoeffs=int(line[22])
                orders = []
                for offset in range(Ncoeffs):
                    start = 24 + 5 * offset
                    orders.append(float(line[start:start + 5].strip()))
                if Ncoeffs > 8:
                    raise Exception("Cannot handle more than 8 Cp coefficients")
                C=[]
                C2 = []

                linecount += 1
                line = lineit.next()
                for offset in range(5):
                    C.append(parseFortanFloat(line[16 * offset: 16 * (offset + 1)]))

                linecount += 1
                line = lineit.next()
                for offset in range(Ncoeffs - 5):
                    C.append(parseFortanFloat(line[16 * offset: 16 * (offset + 1)]))

                HConst = parseFortanFloat(line[48:48 + 16])
                Sconst = parseFortanFloat(line[64:64 + 16])

                fitFunction = "Poly"
                constants = [[C[i], orders[i]] for i in range(Ncoeffs)]
                coeffs.append([Tmin, Tmax, fitFunction, constants, R * HConst, R * Sconst])

            species, comments = firstline.split(" ", 1)
            comments = comments.strip()
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
            #print firstline
            #print "species =",species," phase =",phasetype
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
            while linecount < startline + 1 + intervals * 3:
                linecount += 1
                lineit.next()
    print ""

def initDataDir(directory):
    import os
    print "Loading thermodynamic data"
    parseNASADataFile(os.path.join(directory, 'thermo.inp'), quiet=False)
    print "Loaded",len(speciesData),"thermodynamic species"
    #The Burcat database is very inconsistent which makes it hard to
    #parse in Hf just due to its non-standard formatting
    #parseNASADataFile(os.path.join(directory, '/NEWNASA.inp'), quiet=False)

import sys
import os.path
initDataDir('/usr/local/PyChemEng/data')
