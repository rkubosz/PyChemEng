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
#Constants used in the NASA data set! (If changed, will need to
#rescale the data)
cdef public double R = 8.31451
T0 = 273.15 + 25.0
P0 = 1.0e5

cdef class NASAPolynomial(ThermoConstantsType):
    cdef double a[7]
    cdef double b[2]
    def __init__(NASAPolynomial self, double Tmin, double Tmax, a, b, comments):
        ThermoConstantsType.__init__(self, Tmin, Tmax, comments)
        cdef int i
        for i in range(7):
            self.a[i] = a[i]
        for i in range(2):
            self.b[i] = b[i]

    def __str__(self):
        retval = "NASAPolynomial{Tmin="+str(self.Tmin)+", Tmax="+str(self.Tmax)+", notes='"+self.comments+"', a=["
        for i in range(7):
            retval+=str(self.a[i])+", "
        retval = retval[:-2] + "], b=["
        for i in range(2):
            retval+=str(self.b[i])+", "
        return retval[:-2]+"]}"

    def __repr__(self):
        return self.__str__()
    
    cpdef double Cp0(self, double T):
        return R * (self.a[0] * T ** (-2) + self.a[1] / T + self.a[2] + self.a[3] * T + self.a[4] * T**2 + self.a[5] * T**3 + self.a[6] * T**4)

    cpdef double Hf0(self, double T):
        return R * (-self.a[0] / T + self.a[1] * math.log(T) + self.a[2] * T + self.a[3] * T**2 / 2 + self.a[4] * T**3 / 3 + self.a[5] * T**4 / 4 + self.a[6] * T**5 / 5 + self.b[0])

    cpdef double S0(self, double T):
        return R * (-self.a[0] * T ** (-2) / 2 - self.a[1] / T + self.a[2] * math.log(T) + self.a[3] * T + self.a[4] * T**2 / 2 + self.a[5] * T**3 / 3 + self.a[6] * T**4 / 4 + self.b[1])


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
    linecount=0
    while True:
        linecount += 1
        #Keep reading from the file until we run out of lines
        try:
            line = lineit.next()
        except StopIteration:
            break

        #Skip blank lines, comments, section markers
        if (line.strip() == "") or (line[0] == "!") or (line.strip() == "END PRODUCTS") or (line.strip() == "END REACTANTS"):
            continue

        ######NEED A SMART COMPONENT PARSER#####
        species = line.split()[0]
        comments = line[len(species):-1].strip()
        species = species.strip()
        if not quiet:
            print ">>>Parsing:'"+species+"' '"+comments+"'"
        ########################################
        linecount += 1
        try:
            line = lineit.next()
        except StopIteration:
            break
        try:
            if not quiet:
                print "  Chemical Formula:",line
            refdatecode=line[3:9].strip()
            ####Parse Chemical formula?
            phase=int(line[50:52])
            MW=parseFortanFloat(line[52:65])#g/mol
            HfRef=parseFortanFloat(line[65:80])# @298.15K , in J/mol
            coeffs = []
            MolecularFormula = Components({})
            for offset in range(5):
                shift = 8 * offset
                atom = line[10+shift:12+shift].strip()
                natom = float(line[12+shift:18+shift].strip())
                if atom == "E ":
                    atom = 'e-'
                if len(atom) == 1 and atom[0] == "E":
                    atom = atom[0].lower()+"-"
                if len(atom)==2:
                    atom = atom[0]+atom[1].lower()
                if len(atom)==3:
                    atom = atom[0]+atom[1].lower()+atom[2].lower()
                if natom != 0:
                    MolecularFormula[atom] = natom
            
            intervals=int(line[1:2])
            if intervals == 0:
                #If there are no intervals, there is still the first
                #line of the record, but we skip this data
                if not quiet:
                    print "Skipping Record for",species
                linecount += 1
                line = lineit.next()
                continue

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

                linecount += 1
                line = lineit.next()

                C=[]
                for offset in range(5):
                    C.append(parseFortanFloat(line[16 * offset: 16 * (offset + 1)]))

                linecount += 1
                line = lineit.next()
                for offset in range(Ncoeffs - 5):
                    C.append(parseFortanFloat(line[16 * offset: 16 * (offset + 1)]))

                for i in range(7 - Ncoeffs):
                    C.append(0.0)

                HConst = parseFortanFloat(line[48:48 + 16])
                SConst = parseFortanFloat(line[64:64 + 16])
                coeffs.append(NASAPolynomial(Tmin, Tmax, C, [HConst, SConst], comments))
            
            ##Apply fixes for NASA naming Cl->CL and Al->AL
            if "Cl" in MolecularFormula:
                species = species.replace('CL', 'Cl')

            if "Al" in MolecularFormula:
                species = species.replace('AL', 'Al')

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
                    
            if "(" in species and ")" in species:
                indcomma =species.find(")")+1 
                if len(species) > indcomma+1 and species[indcomma] == ",":
                    phasename = species[species.find("(") + 1: species.find(")")]
                    species = species.replace('('+phasename+')','')
                    if phasename == "L":
                        phasename = "Liquid"
                    else:
                        phasename = str(phase)+phasename
                        
            registerSpecies(species, MolecularFormula)
            sp = speciesData[species]
            sp.registerPhase(phasename)
            
            for C in coeffs:
                sp.registerPhaseCoeffs(C, phasename)
                    
            if sp.inDataRange(T0, phasename):
                HfCalc = sp.Hf0(T0, phasename)
                error = relativeError(HfCalc, HfRef)
                if error > HfMaxError:
                    print "Warning: Species \""+species+"\" and phase "+phasename+" in file:"+filename+" at line "+str(linecount)+" has a Hf of "+str(HfRef)+" but a calculated value of "+str(HfCalc)+" a relative error of "+str(error)+" for Hf at "+str(T0)

        except Exception as e:
            import traceback
            print traceback.print_exc()
            print "Error parsing record for",species,"in file:",filename,"at line",linecount
            print "   ",e.message
            raise

def initDataDir(directory):
    import os
    parseNASADataFile(os.path.join(directory, 'NASA_CEA.inp'))
    print "Loaded NASA dataset (database at",len(speciesData),"species)"
    #parseNASADataFile(os.path.join(directory, 'NEWNASA.TXT'), quiet=False)
    #print "Loaded NEW_NASA dataset (database at",len(speciesData),"species)"

import chemeng.config
initDataDir(chemeng.config.datadir)
