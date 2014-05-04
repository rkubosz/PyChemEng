#!/usr/bin/env python

cdef class IsotopeData:
    """
    This structure holds information on an isotope, including
    its mass, and abundance when available.
    """
    cdef public int Z
    cdef public int N
    cdef public double mass
    cdef public double mass_uncertainty
    cdef public double abundance
    cdef public str name
    def __init__(self, name, Z, N, mass, mass_uncertainty, abundance = 0):
        self.name = name
        self.Z = Z
        self.N = N
        self.mass = mass #Mass of the isotope in U
        self.mass_uncertainty = mass_uncertainty
        self.abundance = abundance #Natural occurance of the isotope (fraction) if available

    def __str__(self):
        return "Isotope{"+self.name+", Z="+str(self.Z)+", N="+str(self.N)+", M="+str(self.mass)+"("+str(self.mass_uncertainty)+"), P="+str(self.abundance)+"}"
    
    def __repr__(self):
        return self.__str__()

cdef class ElementData:
    """
    This class holds information on an element, including its
    average atomic mass, and its known isotopes.
    """

    cdef public int Z
    cdef public int N
    cdef public double mass
    cdef public str name
    cdef public dict isotopes

    
    def __init__(self, name, Z, mass=0):
        self.name = name
        self.Z = Z
        self.mass = mass
        self.isotopes = {}

    def __str__(self):
        output="Element{"+self.name+", Z="+str(self.Z)+", AW="+str(self.mass)
        if len(self.isotopes) != 0:
            output+=", "+str(len(self.isotopes))+" isotopes"
        return output+"}"
    
    def __repr__(self):
        return self.name

####################################################################
# Element data structures
####################################################################
cdef class ElementDatabaseType:
    """
    This class is a dictionary of information on the elements. It can
    be accessed using four types of array access, e.g.

    #Returns the element with a proton count of 6 (carbon) and all isotopes
    ElementDatabaseType[6] 

    #Returns the "C" (carbon) element data with all isotopes
    ElementDatabaseType["C"]

    #Three ways to access the Deuterium isotope of Hydrogen (second argument is neutron number) (no isotopes will be returned for these accesses)
    ElementDatabaseType['D']
    ElementDatabaseType[(1,1)]
    ElementDatabaseType[('H',1)]

    The data returned is of sub type ElementData.
    """
    
    cdef public dict nameIndex #A way to look up the proton number (Z) from an element name
    cdef public dict data

    def __getitem__(self, key):
        if type(key) is tuple:
            element, N = key
            elementdata = self.data[element]
            isotopedata = elementdata.isotopes[N]
            return ElementData(name=isotopedata.name, Z=isotopedata.Z, mass=isotopedata.mass)
        if type(key) is str:
            return self[self.nameIndex[key]]
        return self.data[key]

    def __setitem__(self, key, value):
        if type(key) is str:
            self.data[self.nameIndex[key]] = value
        else:
            self.data[key] = value

    def __contains__(self, key):
        if type(key) is str:
            return key in self.nameIndex
        if type(key) is tuple:
            element, N = key
            if element not in self:
                return False
            elementdata = self.data[element]
            return N in elementdata.isotopes
        else:
            return key in self.data
        
    def ParseTableData(self, filename):
        """This function initialises the data on all isotopes (the
        data and isotopes member variables). This file contains no
        information on abundances which must be obtained from another
        file."""
        self.data = {}
        self.nameIndex = {}
        file = open(filename, "r")
        lineit = iter(file)
        try:
            #First, skip to the third page
            pagecount = 0
            while pagecount < 2:
                if lineit.next()[0] == '1':
                    pagecount += 1
            #Skip the second header line
            lineit.next()

        #Begin record by record parsing of the elements
            A=0
            while True:
                line = lineit.next()
                cc = line[0:1].strip()
                NZ = line[1:4].strip()
                N = int(line[4:9].strip()) #neutron number
                Z = int(line[9:14].strip()) #proton count
                newA = line[14:19].strip()
                if newA != "":#The A values are carried over in the data
                    A = int(newA)
                #discard = line[19:20].strip()
                element = line[20:23].strip()
                #origin = line[23:27].strip()
                #discard = line[27:28].strip()
                #mass = line[28:41].strip()
                #massunc = line[41:52].strip()
                #binding = line[52:63].strip()
                #bindingunc = line[63:72].strip()
                #discard = line[72:73].strip()
                #B = line[73:75].strip()
                #beta = line[75:86].strip()
                #betaunc = line[86:95].strip()
                #discard = line[95:96].strip()
                mass_int = int(line[96:99].strip())
                #discard = line[99:100].strip()
                mass_frac = float(line[100:112].strip().split('#')[0])
                mass_unc = float(line[112:123].strip().split('#')[0])
                #discard = line[123:124].strip()
                if Z not in self.data:
                    self.data[Z] = ElementData(name=element, Z=Z)
                    self.nameIndex[element] = Z

                self[Z].isotopes[N] = IsotopeData(name=element, N=N, Z=Z, mass=mass_int + mass_frac * 1e-6, mass_uncertainty=mass_unc * 1e-6)
        except StopIteration:
            pass

    def ParseIsotopicCompositionsFile(self, filename):
        def ParseLine(expectedtext, line):
            data = line.split("=")
            if len(data) != 2:
                raise Exception("Could not split line '"+line+"'")
            if data[0] != expectedtext:
                raise Exception("Expected '"+expectedtext+"' but got '"+str(data)+"'")
            return data[1].strip()

        file = open(filename, "r")
        lineit = iter(file)
        try:
            while True:
                Z = int(ParseLine("Atomic Number ", lineit.next()))
                Symbol = ParseLine("Atomic Symbol ", lineit.next())
                A = int(ParseLine("Mass Number ", lineit.next()))
                RelMass = ParseLine("Relative Atomic Mass ", lineit.next())
                abundance = ParseLine("Isotopic Composition ", lineit.next())
                AW = ParseLine("Standard Atomic Weight ", lineit.next())
                Notes = ParseLine("Notes ", lineit.next())
                if lineit.next().strip() != "":
                    raise Exception("Expected newline missing!")
                N = A - Z

                #We use information in this file to set abundance information
                if abundance != "":
                    self[Z].isotopes[N].abundance = float(abundance.split('(')[0])

                if AW[0] == '[':
                    AW = AW[1:-1]
                AW = AW.split('(')[0]
                self[Z].mass = float(AW)
                #Figure out which element is representative
        except StopIteration:
            pass

####################################################################
# Periodic Table of elements data
####################################################################
#elements is a dictionary of elements and their masses and proton counts (including electrons)

elements = ElementDatabaseType()
def initDataDir(directory):
    import os
    elements.ParseTableData(os.path.join(directory, 'mass.mas03round.txt'))
    elements.ParseIsotopicCompositionsFile(os.path.join(directory, 'isotopicCompositions.inp'))

    ####Add some special cases
    elements['H'].isotopes[1].name = "D"
    elements['H'].isotopes[2].name = "T"
    elements.nameIndex["D"] = (1,1)
    elements.nameIndex["T"] = (1,2)
    ####Also add a special entry for electrons (using Z=-1 as Z=0 is taken
    ####by neutrons). This is needed for reactions involving ions.
    elements.nameIndex["e-"] = -1
    elements[-1] = ElementData(name="e-", Z=-1, mass=5.4857990943e-4)


import sys
import os.path
initDataDir('/usr/local/PyChemEng/data')
