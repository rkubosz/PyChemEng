#!/usr/bin/env python
#distutils: language = c++

from chemeng.phase cimport Phase
from chemeng.phase import R, IdealGasStream
from chemeng.components cimport Components

cpdef tuple mapToIndex(Components components):
    cdef dict keyToIndex = {}
    cdef list indexToKey = []
    for entry in components._list:
        keyToIndex[entry.first] = len(indexToKey)
        indexToKey.append(entry.first)
    return keyToIndex, indexToKey

#Describe how to turn from the "variables" back into a stream for
#the gibbs calculations.
cpdef Phase variablesToStream(variables, Phase InputStream, bint constP, list indexToComponent):
    """Turns the optimised variables into a Stream class for
    calculations."""
    output = InputStream.copy()
    for i in range(len(indexToComponent)):
        output.components[indexToComponent[i]] = variables[i]
    if not constP:
        output.P = InputStream.P * output.components.total() * output.T / (InputStream.components.total() * InputStream.T)
    return output

#Generate a vector of the elemental composition
cpdef constraintCalc(variables, InputStream, bint constP, indexToComponent, indexToElement):
    stream = variablesToStream(variables, InputStream, constP, indexToComponent)
    cdef Components elementalcomposition = stream.components.elementalComposition()
    cdef list constraints = []
    for i in range(len(indexToElement)):
        constraints.append(elementalcomposition[indexToElement[i]])
    return constraints

cpdef react(Phase InputStream, set extraComponents, bint constT=False, bint constP=False, int outputlevel=0):
    """Performs gas-phase reaction equilibrium calculations. If the
    temperature is constant, this simply calls the IsothermalReact
    function. If non-isothermal, this function iterates an energy
    balance to converge the temperature. """

    #Create an output stream with the input flow rates and all of the
    #required products. We have to copy to stop the input stream being
    #modified
    OutputStream = InputStream.copy()
    for component in extraComponents:
        if component not in InputStream.components:
            OutputStream.components[component] = 0

    if constT:
        return IsothermalReact(OutputStream, constP, outputlevel)

    import math
    for iteration in range(100):
        Told = OutputStream.T
        oldVol = OutputStream.volume()
        OutputStream = IsothermalReact(OutputStream, constP, outputlevel)
        if constP:
            OutputStream.setEnthalpy(InputStream.enthalpy())
        else:
            OutputStream.setInternalEnergy(InputStream.internalEnergy())
            OutputStream.P = OutputStream.components.total() * R * OutputStream.T / oldVol
        if math.fabs((OutputStream.T - Told) / Told) < 0.0001:
            return OutputStream

    raise Exception("Failed to converge to temperature in reaction calculation")
    
def IsothermalReact(InputStream, bint constP=False, int outputlevel=0):
    """Performs a Gibbs or Helmholz free energy minimisation to
    calculate the equilibrium concentrations of a set of species"""

    #Define a mapping from component keys to sequential indices (and
    #back) as we need to be able to turn everything into a matrix and
    #vector representation.
    componentToIndex, indexToComponent = mapToIndex(InputStream.components)

    #Perform the same operation as above, but for elements as these
    #are used to constrain the optimisation to ensure elements (and
    #charge) are conserved.
    elementToIndex, indexToElement = mapToIndex(InputStream.components.elementalComposition())

    #Set up the initial "state" of the system. The initial state
    #"variables" are all of the product (and reactant) species
    #concentrations at the start.
    initialStateVariables=[InputStream.components[component] for component in indexToComponent]

    #Store the initial elemental concentrations to use this as a
    #constraint (it must be preserved)
    constraintTargets = constraintCalc(initialStateVariables, InputStream, constP, indexToComponent, indexToElement)

    ####### Numerical optimisation
    ### Step size
    #The step size should be some small fraction of the available
    #moles in the system. The same value is used for temperature, so
    #we limit it to a max of 1 mol/K.
    stepsize = min(1.0, InputStream.components.total() * 0.01)

    ### Bounds on optimised variables
    #We assume that there cannot be more moles of a single component
    #than there are total moles of elements in the input. This may be
    #broken if free electrons are generated in significant amounts?
    maxmoles = InputStream.components.elementalComposition().total()
    variableBounds = [(0.0, maxmoles)  for i in range(len(indexToComponent))]

    ### Problem scaling and optimisation function
    #The SLSQP algorithm seems to assume the function/gradients are of
    #order "1", so we also need to scale our problem carefully. We
    #divide by the number of moles to make the problem intensive. We
    #also subtract the initial Gibbs free energy and divide by 1000 to bring the energy down to O(100).
    totalInputMoles = 1000.0 * InputStream.components.total()
    optfunc = None
    if constP:
        initialval = InputStream.gibbsFreeEnergy()
        optfunc = lambda variables: (variablesToStream(variables, InputStream, constP, indexToComponent).gibbsFreeEnergy() - initialval) / totalInputMoles
    else:
        initialval = InputStream.gibbsFreeEnergy()
        optfunc = lambda variables: (variablesToStream(variables, InputStream, constP, indexToComponent).helmholtzFreeEnergy() - initialval) / totalInputMoles

    ### Constraints
    # This function returns an array of the deviations of the constraints
    constraintFunc = lambda variables : [(val1 - val2) for val1, val2 in zip(constraintCalc(variables, InputStream, constP, indexToComponent, indexToElement), constraintTargets)]
    
    ### Optimisation
    from scipy.optimize import fmin_slsqp
    optimisedState = fmin_slsqp(optfunc,
                                initialStateVariables,
                                f_eqcons=constraintFunc,
                                bounds=variableBounds,
                                iprint=outputlevel)

    final = variablesToStream(optimisedState, InputStream, constP, indexToComponent)

    ###Return a Stream with the optimised composition (and temperature)
    return final
