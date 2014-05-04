#!/usr/bin/env python
# distutils: language = c++
# cython: profile=True

from chemeng.phase cimport Phase
from chemeng.phase import IdealGasPhase
from chemeng.components cimport Components

cpdef Components getTotalSpecies(list phaselist):
    cdef Components retval
    for phase in phaselist:
        retval += phase.components


def minimiseGibbs(phaselist, elemental = True, T=None, P=None):
    #We need to provide bounds on the optimizer parameters, which may
    #include the number of species, the temperature, and the pressure.

    #Determine what the maximum number of moles of any species is. We
    #estimate this by looking for the maximum number of any single
    #element is (assuming that there is no fractional element species,
    #and that there is not a significant production of electrons/ions).
    cdef Components totalSpecies = getTotalSpecies(phaselist)
    cdef Components totalElements = totalSpecies.elementalComposition()
    cdef double maxmoles = max(totalElements.values())
    
    ####### Numerical optimisation
    ### Step size
    #The step size should be some small fraction of the available
    #moles in the system. The same value is used for temperature, so
    #we limit it to a max of 1 mol/K.
    cdef double totalInputMoles = totalSpecies.total()
    cdef double stepsize = min(1.0, totalInputMoles * 0.01)

    def getState():
        variables =[]
        for phase in phaselist:
            variables += phase.getStateVector()
        return variables
        
    def restoreState(variables):
        cdef int var_offset = 0
        cdef int var_count
        for phase in phaselist:
            var_count = len(phase.components)+2
            phase.restoreStateVector(variables[var_offset:var_offset+var_count])

    cdef double initialOptVal = 0
    def optimisation_func(variables):
        restoreState(variables)
        cdef double Gtotal = 0
        for phase in phaselist:
            Gtotal += phase.gibbsFreeEnergy()
        return Gtotal / totalInputMoles - initialOptVal

    initialOptVal = optimisation_func
    #Create the constraints
    def constraintFunc(variables):
        retval=[]
        cdef Components currentSpecies = getTotalSpecies(phaselist)
        cdef Components currentElements = totalSpecies.elementalComposition()
        if not elemental:
            for key in totalSpecies.keys():
                retval.append(currentSpecies[key] - totalSpecies[key])
        else:
            for key in totalElements.keys():
                retval.append(currentElements[key] - totalElements[key])
        return retval
