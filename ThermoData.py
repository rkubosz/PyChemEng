#!/usr/bin/env python
import math
from Elements import elements
from Components import Components
from Data import speciesData, R, T0, P0

####################################################################
# Fit function registration
####################################################################
#fitFunctions is a dictionary of functions used for data fitting. Each
#function must take two arguments, T and Coeffs. For each fit
#function, you will need to provide two integrated variants.
fitFunctions={}

def registerFitFunction(name, function, integratedfunction, integratedfunctionOverT):
    if name in fitFunctions:
        raise Exception("This function name is already in use!")
    fitFunctions[name] = function
    fitFunctions[name+"Integrated"] = integratedfunction
    fitFunctions[name+"IntegratedOverT"] = integratedfunctionOverT

#When registering a new fit function, you need to add the function, its integral, and the integral of (function/T)
registerFitFunction("Poly",
                    #The function (f)
                    lambda T, C : sum([constant * (T ** order) for constant, order in C]),
                    #The integrated function (f) (without the integration constant)
                    lambda T, C : sum([constant * math.log(T) if order == -1 else constant * (T ** (order+1)) / (order+1) for constant, order in C]),
                    #The integrated function over temperature (f/T)
                    lambda T, C : sum([constant * math.log(T) if order == 0 else constant * (T ** (order)) / (order) for constant, order in C]))

####################################################################
# Thermodynamic Data Structures
####################################################################
#SpeciesData is a dictionary of a list of named tuples of of Cp/H/S
#coefficient data.  Tmin and Tmax are the valid temperature limits for
#each data entry. fitFunction is the name of the function used to
#correlate the data (see next section). constants are the coefficients
#for the fit function. HConst is the integration constant in Hf for
#enthalpy calculations and SConst is the integration constant for the
#entropy.
from collections import namedtuple
ThermoDataType = namedtuple('ThermoData', ['Tmin', 'Tmax', 'fitFunction', 'constants', 'HConst', 'SConst'])
SpeciesDataType = namedtuple('SpeciesData', ['elementalComposition','dataset'])

def inDataRange(component, T):
    for Tmin, Tmax, func, C, Hconst, Sconst in speciesData[component].dataset:
        if T >= Tmin and T <= Tmax:
            return True
    return False

def Cp(component, T):#J / mol K
    for Tmin, Tmax, func, C, Hconst, Sconst in speciesData[component].dataset:
        if T >= Tmin and T <= Tmax:
            return R * fitFunctions[func](T, C)
    raise Exception("Cannot find valid Cp expression for "+str(T)+"K")

def Hf(component, T):#J / mol
    for Tmin, Tmax, func, C, Hconst, Sconst in speciesData[component].dataset:
        if T >= Tmin and T <= Tmax:
            return R * (fitFunctions[func+"Integrated"](T, C) + Hconst)
    raise Exception("Cannot find valid Hf expression for "+str(T)+"K")

def S(component, T):#J / K
    for Tmin, Tmax, func, C, Hconst, Sconst in speciesData[component].dataset:
        if T >= Tmin and T <= Tmax:
            return R * (fitFunctions[func+"IntegratedOverT"](T, C) + Sconst)
    raise Exception("Cannot find valid expression for the temperature T"+str(T))

class PrexistingComponentError(Exception):
    pass

def registerComponent(component, dataset, elementalComposition):
    if component in speciesData:
        raise PrexistingComponentError("Species \""+component+"\" already exists in the database")
    speciesData[component] = SpeciesDataType(elementalComposition, dataset)
    #Check the elemental composition has no unknown elements
    for element in elementalComposition:
        if element not in elements:
            raise Exception("Component "+component+" with composition "+str(elementalComposition)+" has an unknown element, "+element)

####################################################################
# Test functions
####################################################################
HfMaxError=0.01 #1% error
MWMaxError=0.0003 #0.03% error

def relativeError(val, ref):
    return math.fabs((val - ref) / (ref + (ref==0)))

####################################################################
# NASA Glenn Thermodynamic Database
####################################################################

def parseFortanFloat(string):
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

def parseNASADataFile(filename, quiet=True):
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
        component=line.split()[0]
        comments=line.split()[1:]
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
                coeffs.append(ThermoDataType(Tmin, Tmax, fitFunction, constants, HConst, Sconst))
            
            #Skip non-gas phases for now
            if phase != 0:
                continue
            
            try:
                registerComponent(component, coeffs, Components(MolecularFormula))
            except PrexistingComponentError:
                if not quiet:
                    print "Failed to add component",component," to database as component already exists:File:",filename,"at line",linecount
                    
            error = relativeError(Components({component:1}).avgMolarMass(), MW)
            if error > MWMaxError:
                raise Exception("Component \""+component+"\" has a relative MW error of "+str(error))

            if inDataRange(component, 298.15):
                error = relativeError(Hf(component, 298.15), HfRef)
                if error > HfMaxError:
                    raise Exception("Component \""+component+"\" has a relative error of "+str(error)+" for Hf at 298.15")

            print "\rLoaded",len(speciesData),"thermodynamic components",

        except Exception as e:
            if not quiet:
                import traceback
                print traceback.print_exc()
                print "Error parsing record for",component,"in file:",filename,"at line",linecount
                print "   ",e.message
            while linecount < startline + 1 + intervals * 3:
                linecount += 1
                lineit.next()

import os
###thermo.inp database takes priority over the Burcat database
parseNASADataFile(os.path.join(os.path.dirname(__file__), 'thermo.inp'), quiet=False)
parseNASADataFile(os.path.join(os.path.dirname(__file__), 'NEWNASA.inp'), quiet=True)
print ""
