#!/usr/bin/env python
# distutils: language = c++
# cython: profile=True

from chemeng.phase cimport Phase
from chemeng.phase import IdealGasPhase
from chemeng.components cimport Components
from itertools import izip

cpdef double R = 8.3144621

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
    cdef public bint constV
    cdef public bint constH
    cdef public bint constU
    cdef public bint constS
    cdef public bint elemental
    cdef public double initialObjVal
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

        #Add the temperature and pressure of the first phase
        if not self.constT:
            variables.append(self.outputPhases[0].T)
        if not self.constP:
            variables.append(self.outputPhases[0].P)

        return variables

    cpdef restoreStateVector(EquilibriumFinder self, statevec):
        cdef int var_count = 0
        cdef Phase phase
        for phase in self.outputPhases:
            for key in phase.components.keys():
                phase.components[key] = statevec[var_count]
                var_count +=1
                
        if not self.constT:
            for phase in self.outputPhases:
                phase.T = statevec[var_count]
            var_count +=1

        if not self.constP:
            for phase in self.outputPhases:
                phase.P = statevec[var_count]
            var_count +=1            

    #Create the constraints
    cpdef list constraint_func(EquilibriumFinder self):
        cdef Components currentSpecies = getTotalSpecies(self.outputPhases)
        cdef Components currentElements = currentSpecies.elementalComposition()
        constraintvals=[]
        if self.elemental:
            for key in self.inputElements.keys():
                constraintvals.append(currentElements[key])
        else: #not elemental
            for key in self.inputSpecies.keys():
                constraintvals.append(currentSpecies[key])

        cdef double sumX
        cdef Phase phase
        if self.constH:
            sumX = 0
            for phase in self.outputPhases:
                sumX += phase.enthalpy()
            constraintvals.append(sumX)
        if self.constU:
            sumX = 0
            for phase in self.outputPhases:
                sumX += phase.internalEnergy()
            constraintvals.append(sumX)
        if self.constV:
            sumX = 0
            for phase in self.outputPhases:
                sumX += phase.volume()
            constraintvals.append(sumX)
        if self.constS:
            sumX = 0
            for phase in self.outputPhases:
                sumX += phase.entropy()
            constraintvals.append(sumX)
        return constraintvals

    cpdef double optimisation_func(EquilibriumFinder self) except +:
        cdef double total = 0
        cdef Phase phase
        if (self.constT and self.constP):
            for phase in self.outputPhases:
                total += phase.gibbsFreeEnergy()
            #Scale the problem, G / (R * T)
            total /= R * self.inputPhases[0].T
        elif (self.constT and self.constV):
            for phase in self.outputPhases:
                total += phase.helmholtzFreeEnergy()
            #Scale the problem, A / (R * T)
            total /= R * self.inputPhases[0].T
        elif (self.constH and self.constP) or (self.constU and self.constV):
            for phase in self.outputPhases:
                total -= phase.entropy()
            #Scale the problem, S / R
            total /= R
        elif (self.constS and self.constV):
            for phase in self.outputPhases:
                total += phase.internalEnergy()
        elif (self.constS and self.constP):
            for phase in self.outputPhases:
                total += phase.enthalpy()
        else:
            raise Exception("Error, unanticipated combination of constant state variables")
        return total

    cpdef obj_func(EquilibriumFinder self, variables):
        self.restoreStateVector(variables)
        #Remove the large constant offset of the objective function (if present)
        cdef double f = self.optimisation_func() - self.initialObjVal
        cdef list constraints = self.constraint_func()
        cdef int i
        for i in range(len(constraints)):
            constraints[i] -= self.constraintTargets[i]
        return f, constraints, 0

    def __init__(EquilibriumFinder self, list inputPhases, bint constT = False, bint constP = False, bint constH = False, bint constV = False, bint constU = False, bint constS = False, bint elemental = False, double Tmax = 20000, double Pmax = 500, double xtol=1e-8):
        self.inputPhases = inputPhases
        self.outputPhases = []
        cdef Phase phase
        for phase in self.inputPhases:
            newphase = phase.copy()
            #Ensure all phases have the same temperature and pressure
            newphase.T = self.inputPhases[0].T
            newphase.P = self.inputPhases[0].P
            self.outputPhases.append(newphase)
            
        if (constT + constP + constH + constV + constU + constS) != 2:
            raise Exception("EquilibriumFinder requires exactly 2 constant state variables from T, P, H, V, U, and S.")
        self.constT = constT
        self.constP = constP
        self.constH = constH
        self.constV = constV
        self.constU = constU
        self.constS = constS
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
        cdef double stepsize = min(0.1, self.inputMoles * xtol)

        cdef list initialState = self.getStateVector()
        self.initialObjVal = self.optimisation_func()
        self.constraintTargets = self.constraint_func()

        cdef list variableBounds = []
        for phase in self.inputPhases:
            for i in range(len(phase.components)):
                variableBounds += [[0, maxmoles]]
            if not self.constT:
                variableBounds.append([0, Tmax])
            if not self.constP:
                variableBounds.append([0, Pmax])
        
        import pyOpt
        opt_prob = pyOpt.Optimization('Entropy maximisation',lambda x : self.obj_func(x))
        opt_prob.addObj('-S', value=0.0)
        for i,var in enumerate(initialState):
            opt_prob.addVar('x'+str(i), type='c', value=var, lower=0)
        for i,cons in enumerate(self.constraintTargets):
            opt_prob.addCon('cons'+str(i), type='e')
            
            
        opt = pyOpt.pySLSQP.SLSQP()
        opt.setOption('ACC', stepsize)

        cdef bint debug = False
        if debug:
            print opt_prob
            opt.setOption('IPRINT', 0)

        [fstr, xstr, inform] = opt(opt_prob, sens_type='FD', disp_opts=debug)

        if debug:
            print fstr
            print xstr
            print inform
            print opt_prob._solutions[0]

        self.restoreStateVector(xstr)

def findEquilibrium(*args, **kwrdargs):
    solver = EquilibriumFinder(*args, **kwrdargs)
    return solver.outputPhases
