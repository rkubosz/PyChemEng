#!/usr/bin/env python
from Stream import *
import math
from collections import namedtuple
from Elements import elements
from Components import Components
from Data import speciesData, R, T0, P0
import copy

def relativeError(val, ref):
    import math
    return math.fabs((val - ref) / (ref + (ref==0)))
  
def validate(output, expected, error=0.05): ######Maximum allowable Error at 5%######
    if relativeError(output, expected) > error:
        raise Exception("Failed test with error of "+str(relativeError(output, expected)))

class PrexistingComponentError(Exception):
    pass

ViscosityDataType = namedtuple('speciesViscosityData', ['Tmin','Tmax','dataset'])
ThermcondDataType = namedtuple('speciesThermCondData', ['Tmin','Tmax','dataset'])
speciesViscosityData = {}
speciesThermCondData = {}

def registerviscositydata(component, dataset):
    if component in speciesViscosityData:
        raise PrexistingComponentError("Species \""+component+"\" already exists in the database")
    speciesViscosityData[component] = dataset
 
def registerthermconddata(component, dataset):
    if component in speciesThermCondData:
        raise PrexistingComponentError("Species \""+component+"\" already exists in the database")
    speciesThermCondData[component] = dataset

file = open("trans.inp", "r")
lineit = iter(file)
#Skip first line
lineit.next()
while True:
    try:
        line = lineit.next()
    except StopIteration:
        break
    if line[0:3] == "end":
        break
    component1=line[0:15].strip() #Main component
    component2=line[16:31].strip() #Second component
    comments=line[40:80].strip()
    entries=int(line[35])+int(line[37])
    visc_coeffs=[]
    therm_coeffs=[]
    for entry in range(entries):
        try:
            line = lineit.next()
        except StopIteration:
            raise Exception("Ran out of data while reading entry "+component1+" "+component2)
        Tlow=float(line[2:11])
        Thigh=float(line[11:20])
        C=[]
        coeffs=[]
        for offset in range(4):
            start=20 + 15 * offset
            end = start + 15
            data=line[start:end]
            mantissa=line[start:end][:-3].strip()
            exponent=line[start:end][-3:].strip()
            C.append(float(mantissa+exponent))
        if line[1] == 'V':
            visc_coeffs.append(ViscosityDataType(Tlow, Thigh, C))
        elif line[1] == 'C':
            therm_coeffs.append(ThermcondDataType(Tlow, Thigh, C))
        else:
            raise Exception("Unexpected data type")

    if component2 != "":#Skip the mixture data for now
        continue
    registerviscositydata(component1,visc_coeffs)
    registerthermconddata(component1,therm_coeffs) 

def Viscosity(component,T): # Viscosity in Pascal.s
    for Tmin, Tmax, C in speciesViscosityData[component]:
        if T >= Tmin and T <= Tmax:
            return math.exp((C[0] * math.log(T) + C[1]/T + C[2]/(T**2) + C[3]))/1e7 # Divided by 1 e7 to convert to Pascal*S
    raise Exception("Cannot find valid Viscosity expression for "+str(T)+"K")

def ThermalConductivity(component,T) : # Thermal Conductivity in W/m*K 
    for Tmin, Tmax, C in speciesThermCondData [component]:
        if T >= Tmin and T <= Tmax:
	    return (math.exp(C[0] * math.log(T) + C[1]/T + C[2]/(T**2) + C[3]))/1e4 # Divided by 1 e4 to convert to watts per meter kelvin
    raise Exception("Cannot find valid Viscosity expression for "+str(T)+"K")

def IdealGasDensity(component,T): # kg/m3
    return (101325*(Components({component:1}).avgMolarMass())/(R*T))/1000 ###### NOT SURE ABOUT THIS ######

def KinematicViscosity(component,T):# in m2/s
    return Viscosity(component,T)/(IdealGasDensity(component,T))

def Phi (component1, component2, T):
    return 0.25*(1 + ((Viscosity(component1,T)/Viscosity(component2,T))**0.5) * (((Components({component2:1}).avgMolarMass()) / (Components({component1:1}).avgMolarMass()))**0.25))**2 * (2*(Components({component2:1}).avgMolarMass()) / ((Components({component1:1}).avgMolarMass()) + (Components({component2:1}).avgMolarMass())))**0.25

def Psi (component1, component2,T):
    return Phi(component1,component2,T) * (1 + (2.41*((Components({component1:1}).avgMolarMass()) - (Components({component2:1}).avgMolarMass())) * ((Components({component1:1}).avgMolarMass()) - 0.142*(Components({component2:1}).avgMolarMass()))) / (((Components({component1:1}).avgMolarMass()) + (Components({component2:1}).avgMolarMass()))**2)) 

def ViscosityofMixtureGordon(Mixture, T):
    Mixture = Mixture.normalised()
    components = []
    molefractions = []
    runningsum=0
    for component1, value1 in Mixture.iteritems():
        visc1 = Viscosity(component1, T)
        numerator = value1 * visc1
        denominator = value1
        for component2, value2 in Mixture.iteritems():
            if component1 == component2:
                continue
            denominator += value2 * Phi(component1, component2, T)
        runningsum += numerator / denominator
    return runningsum

def ThermCondofMixtureGordon(Mixture, T):
    Mixture = Mixture.normalised()
    components = []
    molefractions = []
    runningsum=0
    for component1, value1 in Mixture.iteritems():
        visc1 = ThermalConductivity(component1, T)
        numerator = value1 * visc1
        denominator = value1
        for component2, value2 in Mixture.iteritems():
            if component1 == component2:
                continue
            denominator += value2 * Psi(component1, component2, T)
        runningsum += numerator / denominator
    return runningsum

def GasMixtureKinematicViscosityGordon(Mixture,T): 
    Mixture = Mixture.normalised()
    components = []
    molefractions = []
    runningsum=0
    for component1, value1 in Mixture.iteritems():
        visc1 = KinematicViscosity(component1, T)
        numerator = value1 * visc1
        denominator = value1
        for component2, value2 in Mixture.iteritems():
            if component1 == component2:
                continue
            denominator += value2 * Phi(component1, component2, T)
        runningsum += numerator / denominator
    return runningsum

def ViscosityofmixtureHZE(Mixture,T): # Herning Zipperer Equation
    Mixture = Mixture.normalised()
    components = Mixture.keys()
    molefractions = Mixture.values()
    NM = len(Mixture)
    return sum(molefractions[i] * Viscosity(components[i],T) * (Components({components[i]:1}).avgMolarMass())**0.5 for i in range (NM))/sum(molefractions[i] * (Components({components[i]:1}).avgMolarMass())**0.5 for i in range (NM))

def ThermCondofmixtureHZE(Mixture,T): # Herning Zipperer Equation
    Mixture = Mixture.normalised()
    components = Mixture.keys()
    molefractions = Mixture.values()
    NM = len(Mixture)
    return sum(molefractions[i] * ThermalConductivity(components[i],T) * (Components({components[i]:1}).avgMolarMass())**0.5 for i in range (NM))/sum(molefractions[i] * (Components({components[i]:1}).avgMolarMass())**0.5 for i in range (NM))

def GasMixtureKinematicViscosityHZE(Mixture,T): # Herning Zipperer Equation
    Mixture = Mixture.normalised()
    components = Mixture.keys()
    molefractions = Mixture.values()
    NM = len(Mixture)
    return sum(molefractions[i] * KinematicViscosity(components[i],T) * (Components({components[i]:1}).avgMolarMass())**0.5 for i in range (NM))/sum(molefractions[i] * (Components({components[i]:1}).avgMolarMass())**0.5 for i in range (NM))

GasMixtureViscosity = ViscosityofMixtureGordon #Selected this equation for Viscosity of gaseous mixtures
GasMixtureThermalConductivity = ThermCondofMixtureGordon #Selected this equation for Thermal Conductivity of gaseous mixtures
GasMixtureKinematicViscosity = GasMixtureKinematicViscosityGordon #Selected this equation for Kinematic Viscosity of gaseous mixtures

def PrandtlFr(Mixture,T):
    return IdealGasStream(T, Mixture).Cp() / (Mixture.totalMass()/1000) * GasMixtureViscosity(Mixture,T)/GasMixtureThermalConductivity(Mixture,T)

DryAir = Components({"N2":78.084, "O2":20.946,  "Ar":0.934, "CO2":0.0397, "Ne":0.001818, "He": 0.000524, "CH4":0.000179, "Kr":0.000114})
DryAir = DryAir.normalised()

validate(PrandtlFr(DryAir,293.15),0.713)

validate(GasMixtureKinematicViscosity(DryAir,300.0),1.5918E-05)# Chemkin @ http://navier.engr.colostate.edu/~dandy/code/code-2/
validate(GasMixtureKinematicViscosity(DryAir,700.0),6.6405E-5)# Chemkin @ http://navier.engr.colostate.edu/~dandy/code/code-2/
validate(GasMixtureKinematicViscosity(DryAir,1500.0),2.3489E-4)# Chemkin @ http://navier.engr.colostate.edu/~dandy/code/code-2/
validate(GasMixtureKinematicViscosity(DryAir,1000.0),1.2026E-4)# Chemkin @ http://navier.engr.colostate.edu/~dandy/code/code-2/

validate(GasMixtureViscosity(DryAir,300.0),18.57e-6) #Kadoya Matsunga and Nagashima Viscosity and Thermal Conductivity of Dryair 
validate(GasMixtureThermalConductivity(DryAir,300.0),26.23e-3)#Kadoya Matsunga and Nagashima Viscosity and Thermal Conductivity of Dryair 
validate(GasMixtureViscosity(DryAir,400.0),23.10e-6) #Kadoya Matsunga and Nagashima Viscosity and Thermal Conductivity of Dryair 
validate(GasMixtureThermalConductivity(DryAir,400.0),33.28e-3)#Kadoya Matsunga and Nagashima Viscosity and Thermal Conductivity of Dryair 
validate(GasMixtureViscosity(DryAir,800.0),37.47e-6) #Kadoya Matsunga and Nagashima Viscosity and Thermal Conductivity of Dryair 
validate(GasMixtureThermalConductivity(DryAir,800.0),56.99e-3)#Kadoya Matsunga and Nagashima Viscosity and Thermal Conductivity of Dryair 
validate(GasMixtureViscosity(DryAir,1500.0),56.48e-6) #Kadoya Matsunga and Nagashima Viscosity and Thermal Conductivity of Dryair 
validate(GasMixtureThermalConductivity(DryAir,1500.0),92.96e-3)#Kadoya Matsunga and Nagashima Viscosity and Thermal Conductivity of Dryair 
validate(GasMixtureViscosity(DryAir,2000.0),67.91e-6) #Kadoya Matsunga and Nagashima Viscosity and Thermal Conductivity of Dryair 
validate(GasMixtureThermalConductivity(DryAir,2000.0),117.5e-3)#Kadoya Matsunga and Nagashima Viscosity and Thermal Conductivity of Dryair 

