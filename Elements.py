#!/usr/bin/env python
import os

####################################################################
# Element data structures
####################################################################
#elements is a dictionary of elements and their masses and proton counts (including electrons)
from collections import namedtuple

class ElementDatabaseType:
    """
    This class is a (read-only) database of information on the
    elements. It can be accessed using three types of array access, e.g.
    ElementDatabaseType[6] #Returns the stable isotope of the element with a proton count of 6 (carbon)
    ElementDatabaseType["C"] #Returns the stable isotope of the "C" element data (carbon)
    ElementDatabaseType[(1,2)] #Returns the element with one proton and two neutrons (heavy isotop of Hydrogen known as Tritium).

    The data returned is of type ElementData.
    """
    data = {} #Data on each isotope, using keys of (Z, N)
    NameIndex = {} #A way to look up the stable isotope key from an element name
    ZIndex = {} #A dictionary of the stable isotope key from a proton count (Z)
    isotopes = {} #A dictionary of isotopes from a proton count (Z)
    avgmass = {} #A dictionary of the avg mass for a proton count
    
    class ElementData:
        """This structure holds information on an element/isotope"""
        def __init__(self, name, mass, mass_uncertainty, Z, N, abundance = None):
            self.name = name
            self.mass = mass #Mass of the isotope in U
            self.mass_uncertainty = mass_uncertainty
            self.abundance = abundance #Natural occurance of the isotope (fraction) if available
            self.Z = Z
            self.N = N
            self.Nrep = None
            
    def __getitem__(self, key):
        if type(key) is str:
            return data[NameIndex[key]]
        elif type(key) is int:
            return data[key]
        elif type(key) is tuple:
            return data[key]
        else:
            raise Exception("Unknown key '"+str(key)+"'")
        
    def ParseTableData(self, filename):
        """This function initialises the data on all isotopes (the
        data and isotopes member variables). This file contains no
        information on abundances which must be obtained from another
        file."""
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
                self.data[(Z, N)] = self.ElementData(name = element, 
                                                     mass = mass_int + mass_frac * 1e-6,
                                                     mass_uncertainty = mass_unc * 1e-6,
                                                     Z = Z,
                                                     N = N)
                if Z not in self.isotopes:
                    self.isotopes[Z] = []
                self.isotopes[Z].append(N)
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
                    self.data[(Z,N)].abundance = float(abundance.split('(')[0])
                if Z not in self.avgmass:
                    if AW[0] == '[':
                        AW = AW[1:-1]
                    AW = AW.split('(')[0]
                    self.avgmass[Z] = float(AW)
                #Figure out which element is representative
        except StopIteration:
            pass

####################################################################
# Periodic Table of elements data
####################################################################

elements = ElementDatabaseType()
elements.ParseTableData(os.path.join(os.path.dirname(__file__), 'datafiles/mass.mas03round.txt'))
elements.ParseIsotopicCompositionsFile(os.path.join(os.path.dirname(__file__), 'datafiles/isotopicCompositions.inp'))
print elements.data
