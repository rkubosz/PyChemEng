#!/usr/bin/env python
# distutils: language = c++

from libcpp.string cimport string
from chemeng.components cimport Components

cdef class ThermoConstantsType:
    cdef public double Tmin
    cdef public double Tmax    
    cpdef double Cp0(ThermoConstantsType, double T)
    cpdef double S0(ThermoConstantsType, double T)
    cpdef double Hf0(ThermoConstantsType, double T)

cdef class SpeciesDataType:
    cdef public string name
    cdef public double mass
    cdef public Components elementalComposition
    cdef public dict phases
    cdef public list antoineData

    cpdef bint inDataRange(SpeciesDataType self, double T, string phase)
    cpdef double Cp0(SpeciesDataType self, double T, string phase)
    cpdef double Hf0(SpeciesDataType self, double T, string phase)
    cpdef double S0(SpeciesDataType self, double T, string phase)
    cpdef double Gibbs0(SpeciesDataType self, double T, string phase)
    cpdef str validRanges(self, phase)
