#!/usr/bin/env python
import os

####################################################################
# Element data structures
####################################################################
#elements is a dictionary of elements and their masses and proton counts (including electrons)
from collections import namedtuple

Element = namedtuple('Element', ['mass', 'protoncount'])
elementToProtonCount = {}

elements={}

def registerElement(name, mass, protoncount):
    """Registers an element (or an isotope) with a given proton count, two-letter name, and mass"""
    if name in elementToProtonCount:
        raise Exception("This element already exists ("+name+")")
    elements[name] = Element(mass, protoncount)

####################################################################
# Periodic Table of elements data
####################################################################
def ParseTableData(filename):
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
        print lineit.next()
    except StopIteration:
        pass

ParseTableData(os.path.join(os.path.dirname(__file__), 'datafiles/mass.mas03round.txt'))

####################################################################
# Isotopic compositions
####################################################################
def ParseIsotopicCompositionsFile(filename):
    def ParseLine(expectedtext, line):
        data = line.split("=")
        if len(data) != 2:
            raise Exception("Could not split line '"+line+"'")
        if data[0] != expectedtext:
            raise Exception("Expected '"+expectedtext+"' but got '"+expectedtext)
        return data[1].strip()

    file = open(filename, "r")
    lineit = iter(file)
    try:
        while True:
            Z = int(ParseLine("Atomic Number ", lineit.next()))
            Symbol = ParseLine("Atomic Symbol ", lineit.next())
            Mass = ParseLine("Mass Number ", lineit.next())
            RelMass = ParseLine("Relative Atomic Mass ", lineit.next())
            Composition = ParseLine("Isotopic Composition ", lineit.next())
            AW = ParseLine("Standard Atomic Weight ", lineit.next())
            Notes = ParseLine("Notes ", lineit.next())
            if lineit.next().strip() != "":
                raise Exception("Expected newline missing!")
            print Z, Symbol, Mass, RelMass, Composition, AW, Notes
    except StopIteration:
        pass

ParseIsotopicCompositionsFile(os.path.join(os.path.dirname(__file__), 'datafiles/isotopicCompositions.inp'))
