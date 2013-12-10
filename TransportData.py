#!/usr/bin/env python
from Stream import *
import math
from collections import namedtuple
from Elements import elements
from Components import Components
from Data import speciesData, R, T0, P0
import copy

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

def IdealKinematicViscosity(component,T):# in m2/s
    return Viscosity(component,T)/(IdealDensity(component,T))
  

def Phi (component1, component2, T):
    return 0.25*(1 + ((Viscosity(component1,T)/Viscosity(component2,T))**0.5) * (((Components({component2:1}).avgMolarMass()) / (Components({component1:1}).avgMolarMass()))**0.25))**2 * (2*(Components({component2:1}).avgMolarMass()) / ((Components({component1:1}).avgMolarMass()) + (Components({component2:1}).avgMolarMass())))**0.25

def Psi (component1, component2,T):
    return Phi(component1,component2,T) * (1 + (2.41*((Components({component1:1}).avgMolarMass()) - (Components({component2:1}).avgMolarMass())) * ((Components({component1:1}).avgMolarMass()) - 0.142*(Components({component2:1}).avgMolarMass()))) / (((Components({component1:1}).avgMolarMass()) + (Components({component2:1}).avgMolarMass()))**2)) 

def ViscosityofMixture(Components, T):
    Components = Components.normalised()
    components = []
    molefractions = []
    for key, value in Components.iteritems():
        components.append(key)
        molefractions.append(value)
    components2 = copy.deepcopy(components)
    molefractions2 = copy.deepcopy(molefractions)
    NM = len(Components)
    firstitemC = components2[0]
    firstitemM = molefractions2[0]
    components2.pop(0)
    components2.append(firstitemC)
    molefractions2.pop(0)
    molefractions2.append(firstitemM)
    return sum((molefractions[i] * Viscosity(components[i],T)) / (molefractions[i] + sum(molefractions2[i] * Phi(components[i],components2[i],T) for i in range (NM) if components[i] != components2[i])) for i in range (NM)) 


def ThermCondofMixture(Components, T):
    Components = Components.normalised()
    components = []
    molefractions = []
    for key, value in Components.iteritems():
        components.append(key)
        molefractions.append(value)
    components2 = copy.deepcopy(components)
    molefractions2 = copy.deepcopy(molefractions)
    NM = len(Components)
    firstitemC = components2[0]
    firstitemM = molefractions2[0]
    components2.pop(0)
    components2.append(firstitemC)
    molefractions2.pop(0)
    molefractions2.append(firstitemM)
    return sum(molefractions[i] * ThermalConductivity(components[i],T) / (molefractions[i] + sum(molefractions2[i] * Psi(components[i],components2[i],T) for i in range (NM) if components[i] != components2[i])) for i in range (NM)) 
  
def ViscosityofmixtureHerningandZippererequation(Components,T):
    Components = Components.normalised()
    components = Components.keys()
    molefractions = Components.values()
    NM = len(Components)
    print components, molefractions, NM
    return sum(molefractions[i] * Viscosity(components[i],T) * (Components({components[i]:1}).avgMolarMass())**0.5 for i in range (NM))/sum(molefractions[i] * (Components({components[i]:1}).avgMolarMass())**0.5 for i in range (NM))




DryAir = Components({"N2":78.084, "O2":20.946,  "Ar":0.934, "CO2":0.0397, "Ne":0.001818, "He": 0.000524, "CH4":0.000179, "Kr":0.000114})
#DryAir = Components({"N2":78.0, "O2":21.0,  "Ar":1.0})


print ViscosityofMixture(DryAir,300.0)
print ViscosityofmixtureHerningandZippererequation(DryAir,300.0)




