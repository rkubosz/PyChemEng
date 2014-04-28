#!/usr/bin/env python
# distutils: language = c++

from libcpp.string cimport string
from chemeng.components cimport Components

cdef class SpeciesDataType:
    cdef public string name
    cdef public double mass
    cdef public Components elementalComposition
    cdef dict phases
    cdef public list antoineData

    cpdef bint inDataRange(SpeciesDataType self, double T, string phase)
    cpdef double Cp0(SpeciesDataType self, double T, string phase)
    cpdef double Hf0(SpeciesDataType self, double T, string phase)
    cpdef double S0(SpeciesDataType self, double T, string phase)
    cpdef double Psat(SpeciesDataType self, double T)
    cpdef PsatTRange(SpeciesDataType self)
    cpdef double Gibbs0(SpeciesDataType self, double T, string phase)
