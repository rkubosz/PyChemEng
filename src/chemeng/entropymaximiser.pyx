#!/usr/bin/env python
# distutils: language = c++
# cython: profile=True

from chemeng.phase cimport Phase
from chemeng.phase import IdealGasPhase
from chemeng.components cimport Components
from itertools import izip

cpdef Components getTotalSpecies(list phaselist):
    cdef Components retval = Components({})
    for phase in phaselist:
        retval += phase.components
    return retval

cdef class EquilibriumFinder:
    cdef public list inputPhases
    cdef public list outputPhases
    cdef public bint constT
    cdef public bint constP
    cdef public bint constH
    cdef public bint elemental
    cdef public double initialS
    cdef public double inputMoles
    cdef public Components inputSpecies
    cdef public Components inputElements
    cdef public list constraintTargets

    cpdef list getStateVector(EquilibriumFinder self):
        cdef list variables = []
        cdef Phase phase
        cdef double moles
        for phase in self.outputPhases:
            variables += [moles for moles in phase.components.values()]
            if not self.constT:
                variables.append(phase.T)
            if not self.constP:
                variables.append(phase.P)
        return variables

    cpdef restoreStateVector(EquilibriumFinder self, statevec):
        cdef int var_count = 0
        cdef Phase phase
        for phase in self.outputPhases:
            for key in phase.components.keys():
                phase.components[key] = statevec[var_count]
                var_count +=1
            if not self.constT:
                phase.T = statevec[var_count]
                var_count +=1
            if not self.constP:
                phase.P = statevec[var_count]
                var_count +=1            

    #Create the constraints
    cpdef list constraintCalc(EquilibriumFinder self, variables):
        self.restoreStateVector(variables)
        cdef Components currentSpecies = getTotalSpecies(self.outputPhases)
        cdef Components currentElements = currentSpecies.elementalComposition()
        constraintvals=[]
        if self.elemental:
            for key in self.inputElements.keys():
                constraintvals.append(currentElements[key])
        else: #not elemental
            for key in self.inputSpecies.keys():
                constraintvals.append(currentSpecies[key])

        cdef double sumH = 0.0
        cdef Phase phase
        if self.constH:
            for phase in self.outputPhases:
                sumH += phase.enthalpy()
            constraintvals.append(sumH)
        return constraintvals

    cpdef list constraint_func(EquilibriumFinder self, variables):
        cdef double current
        cdef double target
        cdef list retval = []
        for current, target in izip(self.constraintCalc(variables), self.constraintTargets):
            retval.append(current - target)
        return retval

    cpdef double optimisation_func(EquilibriumFinder self, variables):
        self.restoreStateVector(variables)
        cdef double Stotal = 0
        for phase in self.outputPhases:
            Stotal += phase.entropy()
        return - Stotal / self.inputMoles - self.initialS

    def __init__(EquilibriumFinder self, list inputPhases, bint constT = False, bint constP = False, bint constH = False, bint elemental = False, double Tmax = 20000, double Pmax = 500):
        self.inputPhases = inputPhases
        self.outputPhases = []
        cdef Phase phase
        for phase in self.inputPhases:
            self.outputPhases.append(phase.copy())
        self.constT = constT
        self.constP = constP
        self.constH = constH
        self.elemental = elemental
        self.inputSpecies = getTotalSpecies(self.inputPhases)
        self.inputElements = self.inputSpecies.elementalComposition()

        #We need to provide bounds on the optimizer parameters, which may
        #include the number of species, the temperature, and the pressure.
        #Determine what the maximum number of moles of any species is. We
        #estimate this by looking for the maximum number of any single
        #element is (assuming that there is no fractional element species,
        #and that there is not a significant production of electrons/ions).
        cdef double maxmoles = max(self.inputElements.values())

        ### Step size
        #The step size should be some small fraction of the available
        #moles in the system. The same value is used for temperature, so
        #we limit it to a max of 0.1 mol/K.
        self.inputMoles = self.inputSpecies.total()
        cdef double stepsize = min(0.1, self.inputMoles * 0.01)

        cdef list initialState = self.getStateVector()
        self.initialS = 0
        self.initialS = self.optimisation_func(initialState)

        cdef list variableBounds = []
        for phase in self.inputPhases:
            for i in range(len(phase.components)):
                variableBounds += [[0, maxmoles]]
            if not self.constT:
                variableBounds.append([0, Tmax])
            if not self.constP:
                variableBounds.append([0, Pmax])
        
        self.constraintTargets = self.constraintCalc(initialState)

        ### Optimisation
        from scipy.optimize import fmin_slsqp
        optimisedState, fx, its, imode, smode = fmin_slsqp(self.optimisation_func,
                                                           initialState,
                                                           f_eqcons=self.constraint_func,
                                                           bounds=variableBounds,
                                                           iprint=5,
                                                           full_output=True)
        if imode != 0:
            error_string = "Minimisation error:("+str(imode)+")-'"+smode+"' after "+str(its)+" iterations\n"
            
            raise Exception(error_string + {
                    2:'Try adding more species to the initial phases.'
                    }.get(imode, "Not sure what to advise here, must exit!"))
        self.restoreStateVector(optimisedState)

def findEquilibrium(*args, **kwrdargs):
    solver = EquilibriumFinder(*args, **kwrdargs)
    return solver.outputPhases
