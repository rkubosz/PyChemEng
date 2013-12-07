#!/usr/bin/env python
from Stream import  R, IdealGasStream
from Components import Components

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
    for iteration in range(10):
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
    componentIndex={}
    componentKey={}
    nextindex=0
    for component in InputStream.components:
        if component not in componentIndex:
            componentIndex[component] = nextindex
            componentKey[nextindex] = component
            nextindex += 1

    #Perform the same operation as above, but for elements as these
    #are used to constrain the optimisation to ensure elements (and
    #charge) are conserved.
    elementIndex={}
    elementKey={}
    nextindex=0
    for element in InputStream.components.elementalComposition():
        if element not in elementIndex:
            elementIndex[element] = nextindex
            elementKey[nextindex] = element
            nextindex += 1

    #Set up the initial "state" of the system. The initial state
    #"variables" are all of the product (and reactant) species
    #concentrations at the start.
    initialStateVariables=[]
    for i in range(len(componentIndex)):
        initialStateVariables.append(InputStream.components[componentKey[i]])

    #Describe how to turn from the "variables" back into a stream for
    #the gibbs calculations.
    def variablesToStream(variables):
        """Turns the optimised variables into a Stream class for
        calculations."""
        output = IdealGasStream(InputStream.T, {componentKey[i]: max(variables[i], 0.0) for i in range(len(variables))})
        if constP:
            output.P = InputStream.P
        else:
            output.P = InputStream.P * output.components.total() * output.T / (InputStream.components.total() * InputStream.T)
        return output

    #Describe how to determine the elemental concentrations and, if
    #required, the enthalpy of a set of state variables.
    def constraintCalc(variables):
        stream = variablesToStream(variables)
        elementalcomposition = stream.components.elementalComposition()
        constraints=[]
        for i in range(len(elementIndex)):
            constraints.append(elementalcomposition[elementKey[i]])
        return constraints

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
    variableBounds = [(0.0, maxmoles)  for i in range(len(componentIndex))]

    ### Problem scaling
    #The SLSQP algorithm might assume the function/gradients is of
    #order "1", so we also need to scale our problem. We choose to
    #scale by the initial Gibbs free energy.
    import math
    scale = 1.0 / math.fabs(InputStream.gibbsFreeEnergy())

    ### Constraints
    # This function returns an array of the deviations of the constraints
    constraintFunc = lambda variables : [(val1 - val2) for val1, val2 in zip(constraintCalc(variables), constraintTargets)]

    func = None
    if constP:
        func = lambda variables: variablesToStream(variables).gibbsFreeEnergy() * scale
    else:
        func = lambda variables: variablesToStream(variables).helmholtzFreeEnergy() * scale
        
    ### Optimisation
    from scipy.optimize import fmin_slsqp
    optimisedState = fmin_slsqp(func,
                                initialStateVariables,
                                f_eqcons=constraintFunc,
                                bounds=variableBounds,
                                iprint=outputlevel)
    
    ###Return a Stream with the optimised composition (and temperature)
    return variablesToStream(optimisedState)
