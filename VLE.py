#!/usr/bin/env python
from scipy import optimize
import numpy as np
import AntoineData
from Data import speciesData

A="C6H6"
B="C7H8"
psatA = speciesData[A].Psat
psatB = speciesData[B].Psat
TRangeA = speciesData[A].PsatTRange()
TRangeB = speciesData[B].PsatTRange()
PRangeA = map(psatA, TRangeA)
PRangeB = map(psatB, TRangeB)

TRange = [max(TRangeA[0],TRangeB[0]), min(TRangeA[1],TRangeB[1])]
PRange = [max(psatA(TRange[0]), psatB(TRange[0])), min(psatA(TRange[1]), psatB(TRange[1]))]

print "Tranges are",TRangeA, TRangeB, TRange
print "PRanges are",PRangeA, PRangeB, PRange

def solveForT(targetP, TRange, P, xA):
    Tlow, Thigh = TRange
    if P(Tlow, xA) > targetP or P(Thigh, xA) < targetP:
        raise Exception("Out of data range")
    while (Thigh - Tlow) > 0.00001:
        Tmid = 0.5 * (Thigh + Tlow)
        if P(Tmid, xA) > targetP:
            Thigh=Tmid
        else:
            Tlow=Tmid
    return 0.5 * (Thigh + Tlow)


pressure = 1.013*1e5

def Pmix(T, xA):
    return xA * psatA(T) + (1.0 - xA) * psatB(T)

def PA(T, xA):
    return psatA(T)

def PB(T, xA):
    return psatB(T)

xs=[0.0]
ys=[0.0]
Ts=[]
hVs=[]
hLs=[]
T = solveForT(pressure, TRangeB, PB, 0)
Ts.append(T)
hVs.append(speciesData[B].Hf0(T, 0))
hLs.append(speciesData[B].Hf0(T, 1))
for xA in np.arange(0, 1, 0.05):
    try:
        T = solveForT(pressure, TRange, Pmix, xA)
        y = xA * psatA(T) / pressure
        hV = xA * speciesData[A].Hf0(T, 0) + (1.0 - xA) * speciesData[B].Hf0(T, 0)
        hL = xA * speciesData[A].Hf0(T, 1) + (1.0 - xA) * speciesData[B].Hf0(T, 1)
        xs.append(xA)
        ys.append(y)
        Ts.append(T)
        hVs.append(hV)
        hLs.append(hL)
    except:
        pass
T = solveForT(pressure, TRangeA, PA, 0)
xs.append(1.0)
ys.append(1.0)
Ts.append(T)
hVs.append(speciesData[A].Hf0(T, 0))
hLs.append(speciesData[A].Hf0(T, 1))

import matplotlib.pyplot as plt
plt.plot(xs,ys, 'x-')
plt.show()

plt.plot(xs,Ts, 'x-')
plt.plot(ys,Ts, 'x-')
plt.show()

plt.plot(xs,hLs, 'x-')
plt.plot(ys,hVs, 'x-')
plt.show()
