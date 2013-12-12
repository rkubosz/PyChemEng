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
  
def validate(output, expected, error=0.025):
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

def IdealDensity(component,T): # kg/m3
    return (101325*(Components({component:1}).avgMolarMass())/(R*T))/1000

def KinematicViscosity(component,T):# in m2/s
    return Viscosity(component,T)/(IdealDensity(component,T))
  

def Phi (component1, component2, T):
    return 0.25*(1 + ((Viscosity(component1,T)/Viscosity(component2,T))**0.5) * (((Components({component2:1}).avgMolarMass()) / (Components({component1:1}).avgMolarMass()))**0.25))**2 * (2*(Components({component2:1}).avgMolarMass()) / ((Components({component1:1}).avgMolarMass()) + (Components({component2:1}).avgMolarMass())))**0.25

def Psi (component1, component2,T):
    return Phi(component1,component2,T) * (1 + (2.41*((Components({component1:1}).avgMolarMass()) - (Components({component2:1}).avgMolarMass())) * ((Components({component1:1}).avgMolarMass()) - 0.142*(Components({component2:1}).avgMolarMass()))) / (((Components({component1:1}).avgMolarMass()) + (Components({component2:1}).avgMolarMass()))**2)) 

def ViscosityofMixture(Mixture, T):
    Mixture = Mixture.normalised()
    components = []
    molefractions = []
    for key, value in Mixture.iteritems():
        components.append(key)
        molefractions.append(value)
    components2 = copy.deepcopy(components)
    molefractions2 = copy.deepcopy(molefractions)
    NM = len(Mixture)
    firstitemC = components2[0]
    firstitemM = molefractions2[0]
    components2.pop(0)
    components2.append(firstitemC)
    molefractions2.pop(0)
    molefractions2.append(firstitemM)
    return sum((molefractions[i] * Viscosity(components[i],T)) / (molefractions[i] + sum(molefractions2[j] * Phi(components[i],components2[j],T) for j in range (NM) )) for i in range (NM)) 


def ThermCondofMixture(Mixture, T):
    Mixture = Mixture.normalised()
    components = []
    molefractions = []
    for key, value in Mixture.iteritems():
        components.append(key)
        molefractions.append(value)
    components2 = copy.deepcopy(components)
    molefractions2 = copy.deepcopy(molefractions)
    NM = len(Mixture)
    firstitemC = components2[0]
    firstitemM = molefractions2[0]
    components2.pop(0)
    components2.append(firstitemC)
    molefractions2.pop(0)
    molefractions2.append(firstitemM)
    return sum(molefractions[i] * ThermalConductivity(components[i],T) / (molefractions[i] + sum(molefractions2[j] * Psi(components[i],components2[j],T) for j in range (NM) )) for i in range (NM)) 
  
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
  
def KinematicViscosityofmixtureHZE(Mixture,T): # Herning Zipperer Equation
    Mixture = Mixture.normalised()
    components = Mixture.keys()
    molefractions = Mixture.values()
    NM = len(Mixture)
    return sum(molefractions[i] * KinematicViscosity(components[i],T) * (Components({components[i]:1}).avgMolarMass())**0.5 for i in range (NM))/sum(molefractions[i] * (Components({components[i]:1}).avgMolarMass())**0.5 for i in range (NM))
  
def PrandtlFr(Mixture,T):
    return IdealGasStream(T, Mixture).Cp() / (Mixture.totalMass()/1000) * ViscosityofmixtureHZE(Mixture,T)/ThermCondofmixtureHZE(Mixture,T)



DryAir = Components({"N2":78.084, "O2":20.946,  "Ar":0.934, "CO2":0.0397, "Ne":0.001818, "He": 0.000524, "CH4":0.000179, "Kr":0.000114})
DryAir = DryAir.normalised()
#DryAir = Components({"N2":78.0, "O2":21.0,"Ar":1.0})



validate(PrandtlFr(DryAir,293.15),0.713)

validate(KinematicViscosityofmixtureHZE(DryAir,300.0),1.5918E-05)# Chemkin @ http://navier.engr.colostate.edu/~dandy/code/code-2/
validate(ThermCondofmixtureHZE(DryAir,300.0),2.6351E-02)# Chemkin @ http://navier.engr.colostate.edu/~dandy/code/code-2/
validate(KinematicViscosityofmixtureHZE(DryAir,1100.0),1.4270E-4)# Chemkin @ http://navier.engr.colostate.edu/~dandy/code/code-2/

#validate(ThermCondofmixtureHZE(DryAir,1100.0),7.5246E-02 )# Chemkin @ http://navier.engr.colostate.edu/~dandy/code/code-2/
print Phi("N2","O2",300)
print Phi("O2","N2",300)
print ViscosityofMixture(DryAir,300.0)
print ViscosityofmixtureHZE(DryAir,300.0)
print ThermCondofMixture(DryAir,300.0)
print ThermCondofmixtureHZE(DryAir,300.0)
print KinematicViscosityofmixtureHZE(DryAir,300.0)



