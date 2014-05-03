#!/usr/bin/env python
# distutils: language = c++
# cython: profile=True

from libcpp.string cimport string
from libcpp.map cimport map

from chemeng.elementdata import elements
from chemeng.speciesdata import speciesData
   
cdef class Components:
   def __init__(self, dict data):
       for key, entry in data.iteritems():
           self._list[key] = entry
       
   cpdef Components copy(Components self):
       cdef Components retval = Components.__new__(Components)
       retval._list = self._list.copy()
       return retval

   cpdef Components mix(self, Components other):
       for entry in other._list:
           it = self._list.find(entry.first)
           if it == self._list.end():
               self._list.insert(entry)
           else:
               self._list[entry.first] += entry.second

   def __add__(self, Components other):
       cdef Components copy = self.copy()
       copy.mix(other)
       return copy

   def __mul__(self, double factor):
       cdef Components copy = self.copy()
       copy.scale(factor)
       return copy

   cpdef double totalMass(self):
       cdef double sum = 0.0
       for entry in self._list:
           sum += entry.second * speciesData[entry.first].mass
       return sum

   cpdef double avgMolarMass(self): #g / mol
       return self.totalMass() / self.total()

   cpdef Components elementalComposition(Components self): #mol
       cdef Components retval = Components.__new__(Components)
       for entry in self._list:
           retval.mix(speciesData[entry.first].elementalComposition * entry.second)
       return retval

   cpdef Components scale(Components self, double factor):
       """A * operator to allow scaling of components (e.g. Components * 2)"""
       for entry in self._list:
           self._list[entry.first] *= factor
       return self

   cpdef double total(self):#mol
       """Returns the total molar flow of all the components"""
       cdef double total = 0
       for entry in self._list:
           total += entry.second
       return total

   cpdef Components normalised(self):
       """Creates a copy of the component stream which is normalised (total molar flow of 1)"""
       cdef Components copy = self.copy()
       copy.scale(1.0 / self.total())
       return copy

   def __str__(self):
       output="C{"
       for entry in self._list:
           output += "\'%s\':%g, " % (entry.first, entry.second)
       return output[:-2]+"}"
   
   def __len__(self):
       return self._list.size()

   def __getitem__(self, string key):
       cdef map[string, double].iterator it
       it = self._list.find(key)
       if it == self._list.end():
           return 0.0
       cimport cython.operator.dereference
       return cython.operator.dereference(it).second

   def __setitem__(self, string key, double value):
       self._list[key] = value

   def __contains__(self, string key):
       cdef map[string, double].iterator it
       it = self._list.find(key)
       return it != self._list.end()
