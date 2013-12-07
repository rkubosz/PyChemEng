#!/usr/bin/env python
import math
from Elements import elements
from Data import speciesData

####################################################################
# Components class
####################################################################
class Components(dict):
    """A class which represents a set of molar flows OR mole fractions OR elemental composition"""
    def __add__(self, other):
        """A + operator to allow mixing of components"""
        output=Components(self)
        for component in other:
            if component in output:
                output[component] += other[component]
            else:
                output[component] = other[component]
        return output

    def __mul__(self, scale):
        """A * operator to allow scaling of components (e.g. Components * 2)"""
        scale=scale+0.0 #make sure the scale is a float!
        output=Components()
        for component, flow in self.iteritems():
            output[component] = scale * flow
        return output

    __rmul__ = __mul__ #to allow right multiplications (e.g. 2 * Components)

    def __div__(self, scale):
        """A / operator to allow scaling of components (e.g. Components / 2)"""
        return self * (1.0 / scale)

    def total(self):#mol
        """Returns the total molar flow of all the components"""
        return sum(flow for component, flow in self.iteritems())

    def normalised(self):
        """Creates a copy of the component stream which is normalised (total molar flow of 1)"""
        return self / self.total()

    def totalMass(self):#g
        """Returns the total mass flow"""
        return sum(flow * elements[element].mass for element, flow in self.elementalComposition().iteritems())

    def avgMolarMass(self): #g / mol
        return  self.totalMass() / self.total()

    def elementalComposition(self): #mol
        output = Components()
        for component, flow in self.iteritems():
            output += flow * speciesData[component].elementalComposition
        return output

    def __str__(self):
        output="C{"
        for component, flow in self.iteritems():
            output += "\'%s\':%g, " % (component, flow)
        return output[:-2]+"}"
