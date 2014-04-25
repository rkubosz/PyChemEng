#!/usr/bin/env python
from Stream import  R, IdealGasStream
from Components import Components

def dictToIndex(dictionary):
    keyToIndex={}
    indexToKey=[]
    for key in dictionary:
        keyToIndex[key] = len(indexToKey)
        indexToKey.append(key)
    return keyToIndex, indexToKey

def react(InputStream, extraComponents, constT=False, constP=False, outputlevel=0):
    """Performs gas-phase reaction equilibrium calculations. If the
    temperature is constant, this simply calls the IsothermalReact
    function. If non-isothermal, this function iterates an energy
    balance to converge the temperature. """

    #Create an output stream with the input flow rates and all of the
    #required products.
    for component in extraComponents:
        if component not in InputStream.components:
            InputStream.components[component] = 0

    if constT:
        return IsothermalReact(InputStream, constP, outputlevel)

    import math
    Output=InputStream
    for iteration in range(100):
        Told = Output.T
        oldVol = Output.volume()
        Output = IsothermalReact(Output, constP, outputlevel)
        if constP:
            Output.setEnthalpy(InputStream.enthalpy())
        else:
            Output.setInternalEnergy(InputStream.internalEnergy())
            Output.P = Output.components.total() * R * Output.T / oldVol
        if math.fabs((Output.T - Told) / Told) < 0.0001:
            return Output

    raise Exception("Failed to converge to temperature in reaction calculation")
    
def IsothermalReact(InputStream, constP=False, outputlevel=0):
    """Performs a Gibbs or Helmholz free energy minimisation to
    calculate the equilibrium concentrations of a set of species"""

    #Define a mapping from component keys to sequential indices (and
    #back) as we need to be able to turn everything into a matrix and
    #vector representation.
    componentToIndex, indexToComponent = dictToIndex(InputStream.components)

    #Perform the same operation as above, but for elements as these
    #are used to constrain the optimisation to ensure elements (and
    #charge) are conserved.
    elementToIndex, indexToElement = dictToIndex(InputStream.components.elementalComposition())

    #Set up the initial "state" of the system. The initial state
    #"variables" are all of the product (and reactant) species
    #concentrations at the start.
    initialStateVariables=[InputStream.components[component] for component in indexToComponent]

    #Describe how to turn from the "variables" back into a stream for
    #the gibbs calculations.
    def variablesToStream(variables):
        """Turns the optimised variables into a Stream class for
        calculations."""
        from copy import copy 
        output = copy(InputStream)
        output.components = Components({component: max(variables[index], 0.0) for index, component in enumerate(indexToComponent)})
        if not constP:
            output.P = InputStream.P * output.components.total() * output.T / (InputStream.components.total() * InputStream.T)
        return output

    #Generate a vector of the elemental composition
    def constraintCalc(variables):
        stream = variablesToStream(variables)
        elementalcomposition = stream.components.elementalComposition()
        return [elementalcomposition[element] for index,element in enumerate(indexToElement)]

    #Store the initial elemental concentrations to use this as a
    #constraint (it must be preserved)
    constraintTargets = constraintCalc(initialStateVariables)

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
        initialval = variablesToStream(initialStateVariables).gibbsFreeEnergy()
        optfunc = lambda variables: (variablesToStream(variables).gibbsFreeEnergy() - initialval) / totalInputMoles
    else:
        initialval = variablesToStream(initialStateVariables).gibbsFreeEnergy()
        optfunc = lambda variables: (variablesToStream(variables).helmholtzFreeEnergy() - initialval) / totalInputMoles

    ### Constraints
    # This function returns an array of the deviations of the constraints
    constraintFunc = lambda variables : [(val1 - val2) for val1, val2 in zip(constraintCalc(variables), constraintTargets)]
    
    ### Optimisation
    from scipy.optimize import fmin_slsqp
    optimisedState = fmin_slsqp(optfunc,
                                initialStateVariables,
                                f_eqcons=constraintFunc,
                                bounds=variableBounds,
                                iprint=outputlevel,
                                )
    
    ###Return a Stream with the optimised composition (and temperature)
    return variablesToStream(optimisedState)
