#!/usr/bin/env python
#distutils: language = c++

from libcpp.map cimport map
from libcpp.string cimport string

cdef class Components:
   cdef map[string, double] _list
   cpdef Components copy(Components)
   cpdef Components mix(Components, Components)
   cpdef totalMass(Components)
   cpdef avgMolarMass(Components)
   cpdef Components elementalComposition(Components)
   cpdef Components scale(Components, double)
   cpdef double total(Components)
   cpdef Components normalised(Components)
