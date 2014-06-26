#!/usr/bin/env python
# distutils: language = c++
# cython: profile=True

from chemeng.phase cimport Phase
from chemeng.phase import IdealGasPhase
from chemeng.components cimport Components
import pyOpt
import math

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
    cdef public bint logMolar
    cdef public double initialObjVal
    cdef public double inputMoles
    cdef public Components inputSpecies
    cdef public Components inputElements
    cdef public list constraintTargets

    cpdef restoreStateVector(EquilibriumFinder self, statevec):
        '''This loads the optimisation variables (stored in a
        list/array for the optimizer) from statevec into
        self.outputPhases so that we might calculate thermodynamic
        properties used for the constraints and objective function.'''
        cdef int var_count = 0
        cdef Phase phase
        for phase in self.outputPhases:
            for key in phase.components.keys():
                if self.logMolar:
                    phase.components[key] = math.exp(statevec[var_count])
                else:
                    phase.components[key] = statevec[var_count]
                var_count +=1
        
        #Here we can "constrain" the temperature and/or pressure by
        #deciding whether to use it as a variable of the minimisation
        if not self.constT:
            for phase in self.outputPhases:
                phase.T = statevec[var_count] * self.inputPhases[0].T
            var_count +=1

        if not self.constP:
            for phase in self.outputPhases:
                phase.P = statevec[var_count] * self.inputPhases[0].P
            var_count +=1

    #Create the constraints
    cpdef list constraint_func(EquilibriumFinder self):
        '''This function returns a list of the constrained values of
        the system, depending on what is held constant'''
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
        #Constant temperature and pressure are handled in
        #restoreStateVector, where we decided wheter to use them as
        #optimisation variables or not
        return constraintvals

    cpdef double optimisation_func(EquilibriumFinder self) except +:
        '''This function returns the required thermodynamic potential
        to be minimised depending on the process studied (what is held constant).'''
        cdef double total = 0
        cdef Phase phase
        if (self.constT and self.constP):
            for phase in self.outputPhases:
                total += phase.gibbsFreeEnergy()
            #Using the dimensionless form, G / (R * T)
            total /= R * self.inputPhases[0].T
        elif (self.constT and self.constV):
            for phase in self.outputPhases:
                total += phase.helmholtzFreeEnergy()
            #Using the dimensionless form, A / (R * T)
            total /= R * self.inputPhases[0].T
        elif (self.constH and self.constP) or (self.constU and self.constV):
            for phase in self.outputPhases:
                total -= phase.entropy()
            #Using the dimensionless form, S / R
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
        '''Returns the objective function to be minimised, a list of
        constraint values which must be zero, and a
        success(0)/fail(-1) parameter.'''
        cdef double f
        cdef list constraints
        cdef int i
        try:
            self.restoreStateVector(variables)
            #Remove any constant offset of the objective function
            f = self.optimisation_func() - self.initialObjVal
            #Take the current constraint values and subtract the
            #target values, as the optimizer expects constraints to
            #have a value of zero when they are satisfied.
            constraints = self.constraint_func()
            for i in range(len(constraints)):
                constraints[i] -= self.constraintTargets[i]
            return f, constraints, 0
        except:
            return 0, self.constraintTargets, -1

    def __init__(EquilibriumFinder self, list inputPhases, bint constT = False, bint constP = False, bint constH = False, bint constV = False, bint constU = False, bint constS = False, bint elemental = False, double Tmax = 20000, double Pmax = 500, double xtol=1e-5, bint logMolar = False, bint debug=False, double Tmin = 200.00001, iterations=50):
        '''Sets up and runs the optimisation process to find equilibrium'''

        #Sanity checks
        if (constT + constP + constH + constV + constU + constS) != 2:
            raise Exception("EquilibriumFinder requires exactly 2 constant state variables from T, P, H, V, U, and S.")

        #Data member initialisation
        self.constT = constT
        self.constP = constP
        self.constH = constH
        self.constV = constV
        self.constU = constU
        self.constS = constS
        self.logMolar = logMolar
        self.elemental = elemental
        self.inputPhases = inputPhases
        self.outputPhases = []
        self.inputSpecies = getTotalSpecies(self.inputPhases)
        self.inputElements = self.inputSpecies.elementalComposition()
        self.inputMoles = self.inputSpecies.total()

        #In the minimisation, we try to scale all variables to be of
        #order unity. This means we divide all molar amounts by the
        #number of input moles, and rescale back
        #afterwards. Temperatures and pressures are scaled by the
        #input values when being passed back and forth from the
        #optimizer. Objective functions, such as the Gibb's free
        #energy, are scaled to their dimensionless forms where
        #possible.

        #Setup of the "output" phases. These are working copies of the
        #"input" phases to calculate thermodynamic properties, and
        #will hold the final solution of the problem.
        cdef Phase phase
        for phase in self.inputPhases:
            #We must copy the phases, as python tends to take things
            #by reference. We cannot assume that we can change
            #inputPhases (and we need the properties of the input
            #phases anyway).
            newphase = phase.copy()
            #Scale the phases by the total number of input moles
            newphase.components /= self.inputMoles
            #Ensure all phases have the same temperature and pressure
            #(taken from the first phase)
            newphase.T = self.inputPhases[0].T
            newphase.P = self.inputPhases[0].P
            self.outputPhases.append(newphase)
        
        #Determine the initial value of the optimisation function, so
        #that we can remove the constant offset of its value.
        self.initialObjVal = self.optimisation_func()

        #Grab the current values of the constrained variables. These
        #are the target values the optimisation must satisfy.
        self.constraintTargets = self.constraint_func()
        
        #########  Create Optimisation System
        opt_prob = pyOpt.Optimization('Entropy maximisation', self.obj_func)
        opt_prob.addObj('ThermoPotential')
        opt = pyOpt.pySLSQP.SLSQP()

        ######### Step size
        cdef double stepsize = xtol
        opt.setOption('ACC', stepsize)
        opt.setOption('MAXIT', iterations)

        #########  Load optimisation variables and their bounds        
        #Here, we also name the variables within the optimiser to
        #enable debugging later on.
        cdef int phase_count = 0
        for phase in self.outputPhases:
            phase_count += 1
            for key, moles in phase.components.iteritems():
                if self.logMolar:
                    if moles == 0.0:
                        opt_prob.addVar(str(phase_count)+':ln('+key+')', type='c', value=math.log(xtol))
                    else:
                        opt_prob.addVar(str(phase_count)+':ln('+key+')', type='c', value=math.log(moles))
                else:
                    opt_prob.addVar(str(phase_count)+':'+key, type='c', value=moles, lower=0.0)
             
        if not self.constT:
            opt_prob.addVar('T/Tinitial', type='c', value=1.0, lower= Tmin / self.inputPhases[0].T)
            
        if not self.constP:
            opt_prob.addVar('P/Pinitial', type='c', value=1.0, lower=0)
        

        #########  Load optimisation constraints
        if self.elemental:
            for key in self.inputElements.keys():
                opt_prob.addCon('Element:'+key, type='e')
        else: #not elemental
            for key in self.inputSpecies.keys():
                opt_prob.addCon('Species:'+key, type='e')

        if self.constH:
            opt_prob.addCon('Enthalpy', type='e')
        if self.constU:
            opt_prob.addCon('InternalE', type='e')
        if self.constV:
            opt_prob.addCon('Volume', type='e')
        if self.constS:
            opt_prob.addCon('Entropy', type='e')

        ######### Run optimisation
        if debug:
            print opt_prob
            opt.setOption('IPRINT', 0)

        [fstr, xstr, inform] = opt(opt_prob, sens_type='FD', disp_opts=debug)
    
        if inform['value'] != 0:
            raise Exception("findEquilibrium error:"+str(inform['value'])+": "+inform['text'])

        if debug:
            print fstr #The optimial minimum value of the objective function
            print xstr #The final state vector
            print inform #dictionary of 
            print opt_prob._solutions[0]

        #recreate the output system, ready for output
        self.restoreStateVector(xstr)
        for phase in self.outputPhases:
            phase.components *= self.inputMoles

def findEquilibrium(*args, **kwrdargs):
    solver = EquilibriumFinder(*args, **kwrdargs)
    return solver.outputPhases
