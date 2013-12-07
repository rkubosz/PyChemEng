#!/usr/bin/env python
from Streams import *

speciesTransportData = {}

def registerTransportProperty(component, dataset):
    if component in speciesTransportData:
        raise PrexistingComponentError("Species \""+component+"\" already exists in the database")
    speciesTransportData[component] = dataset

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
            visc_coeffs.append([[Tlow, Thigh], C])
        elif line[1] == 'C':
            therm_coeffs.append([[Tlow, Thigh], C])
        else:
            raise Exception("Unexpected data type")

    if component2 != "":#Skip the mixture data for now
        continue
    coeffs.append([therm_coeffs,visc_coeffs])
    registerTransportProperty(component1,coeffs)   
     
print speciesTransportData["O2"]
