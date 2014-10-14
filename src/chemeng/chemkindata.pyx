#!/usr/bin/env python
import math
#Ensure that the elemental database has been loaded first
import chemeng.elementdata
from chemeng.components import Components
from chemeng.speciesdata import speciesData, registerSpecies, relativeError
from chemeng.speciesdata cimport SpeciesDataType,ThermoConstantsType

####################################################################
# Physical constants
####################################################################
cdef public double R = 8.31451
T0 = 273.15 + 25.0
P0 = 1.0e5

cdef class ChemKinPolynomial(ThermoConstantsType):
    cdef double a[7]
    def __init__(ChemKinPolynomial self, double Tmin, double Tmax, a, comments):
        ThermoConstantsType.__init__(self, Tmin, Tmax, comments)
        cdef int i
        for i in range(7):
            self.a[i] = a[i]

    def __str__(self):
        retval = "ChemKinPolynomial{Tmin="+str(self.Tmin)+", Tmax="+str(self.Tmax)+", "
        for i in range(7):
            retval+=str(self.a[i])+", "
        return retval[:-2]+", notes='"+self.comments+"'}"
    
    cpdef double Cp0(self, double T):
        return R * (self.self.a[0] + self.a[1] * T + self.a[2] * T**2 + self.a[3] * T**3 + self.a[4] * T**4)

    cpdef double Hf0(self, double T):
        return R * (self.a[5] + T * (self.a[0] + self.a[1] * T / 2 + self.a[2] * T**2 / 3 + self.a[3] * T**3 / 4 + self.a[4] * T**4 / 5))

    cpdef double S0(self, double T):
        return R * (self.self.a[0] * math.log(T) + self.a[1] * T + self.a[2] * T**2 / 2 + self.a[3] * T**3 / 3 + self.a[4] * T**4 / 4 + self.a[6])

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
    commonT = float(line[10:20])
    #Begin parsing records
    linecount=0
    while True:
        try:
            #Keep reading from the file until we run out of lines
            try:
                linecount += 1
                line = lineit.next()
            except StopIteration:
                break
            
            if line.strip() == "END":
                break
            
            species = line[0:18].strip()
            comments = line[18:24].strip()
            if not quiet:
                print ">>>Parsing:'"+species+"' '"+comments+"'"
                        
            if not quiet:
                print "  Chemical Formula:",line,
            MolecularFormula = Components({})
            print line,
            for entry in range(4):
                shift = 5 * entry
                print line[24+shift:26+shift],line[26+shift:29+shift]
                atom = line[24+shift:26+shift].strip()
                if atom == "":
                    continue
                if atom == "E":
                    atom = 'e-'
                if len(atom) == 2:
                    atom = atom[0] + atom[1].lower()
                num=line[26+shift:29+shift]
                if num[-1] == '.':
                    num = num[:-1]
                natom=int(num)
                MolecularFormula[atom] = natom
            print MolecularFormula
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
                C.append(parseFortanFloat(line[i*15:(i+1)*15]))
            
            try:
                linecount += 1
                line = lineit.next()
            except StopIteration:
                break
            
            for i in range(5):
                C.append(parseFortanFloat(line[i*15:(i+1)*15]))
            
            try:
                linecount += 1
                line = lineit.next()
            except StopIteration:
                break
            
            for i in range(4):
                C.append(parseFortanFloat(line[i*15:(i+1)*15]))
            
            registerSpecies(species, MolecularFormula)
            sp = speciesData[species]
            sp.registerPhase(phase)
            sp.registerPhaseCoeffs(ChemKinPolynomial(lowT, midT, C[:7], comments), phase)
            sp.registerPhaseCoeffs(ChemKinPolynomial(midT, highT, C[7:14], comments), phase)
        except Exception as e:
            import traceback
            print traceback.print_exc()
            print "Error parsing record for",species,"in file:",filename,"at line",linecount
            print "   ",e.message
            raise

def initDataDir(directory):
    import os
    parseCHEMKINDataFile(os.path.join(directory, 'BurcatCHEMKIN.DAT'), quiet=False)
    print "Loaded Burcat ChemKin dataset (database at",len(speciesData),"species)"

import chemeng.config
initDataDir(chemeng.config.datadir)
