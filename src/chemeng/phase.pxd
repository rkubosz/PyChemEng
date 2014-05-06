#!/usr/bin/env python
#distutils: language = c++

from chemeng.components cimport Components

cdef class Phase:
    """A base class which holds fundamental methods and members of a single phase which may contain multiple components"""
    cdef public double T
    cdef public double P
    cdef public str phase
    cdef public Components components

    cpdef double Cp(Phase) except +
    cpdef double enthalpy(Phase) except +
    cpdef double entropy(Phase) except +
    cpdef double gibbsFreeEnergy(Phase)
    cpdef double internalEnergy(Phase)
    cpdef double volume(Phase)
    cpdef double helmholtzFreeEnergy(Phase)

cdef class IdealGasPhase(Phase):
    cpdef IdealGasPhase copy(IdealGasPhase)
    cpdef double entropy(IdealGasPhase)
    cpdef double volume(IdealGasPhase)
    cpdef double Cv(IdealGasPhase)

cdef class IncompressiblePhase(Phase):
    cdef public double molarvolume
    cpdef IncompressiblePhase copy(IncompressiblePhase)
    cpdef double enthalpy(IncompressiblePhase)
    cpdef double volume(IncompressiblePhase)
    cpdef double Cv(IncompressiblePhase)

