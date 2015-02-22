#!/usr/bin/env python

from chemeng import *
import chemeng.NASAdata
import scipy.optimize
import math

#Taken from     
def pVAP2(T):
    g=[-2.9912729e3,
       -6.0170128e3,
       +1.887643854e1,
       -2.8354721e-2,
       +1.7838301e-5,
       -8.4150417e-10,
       +4.4412543e-13,
       2.858487]
    sumval=g[7] * math.log(T)
    for i in range(7):
        sumval += g[i] * math.pow(T, i-2)
    return math.exp(sumval)

#Table by T
f=open("table.tex", "w")
print >>f,"\\documentclass{report}"
print >>f,"\\usepackage{longtable}"
print >>f,"\\begin{document}"
print >>f,"\\begin{longtable}{|c|c|c|c|c|c|c|}"
print >>f,"\\caption{Thermodynamic properties of saturated steam, calculated
using the NASA CEA database and the vapour pressure data of
Wexler}\\\\\\hline"
print >>f,"$T$ & $P$ & $C_{p,l}$ & $C_{p,v}$ & $h_l$ & $h_{lv}$ & $h_v$\\\\"
print >>f,"(K) & (bar) & (kJ/kg~K) & (kJ/kg) & (kJ/kg) & (kJ/kg) & (kJ/kg) \\\\\\hline\hline"
print >>f,"\\endhead"
Tref=0.01
Pref = pVAP2(Tref+273.15)
href = IncompressiblePhase({"H2O":1}, "Liquid", T=Tref+273.15, P=Pref, molarvolume=18.01528/999.97e3).enthalpy()


import numpy
#[Tref]+[274.15+1*i for i in range(30)]+[305.15+2*i for i in range(10)]+[328.15+5*i for i in range(10)] +[373.15+10*i for i in range(11)]
TCs = map(float, [Tref]+list(numpy.arange(1,21,1))+list(numpy.arange(25,85,5))+list(numpy.arange(90,210,10)))
for TC in TCs:
    T=TC+273.15
    P = pVAP2(T)
    Liq = IncompressiblePhase({"H2O":1}, "Liquid", T=T, P=P, molarvolume=18.01528/999.97e3)
    Gas = IdealGasPhase({"H2O":1}, T=T, P=P)
    print >>f," & ".join(("%4g"% TC,"%4g"%(P/1e5), "%4g"%(Liq.Cp()/18.01528), "%4g"%(Gas.Cp()/18.01528), "%4g"%((Liq.enthalpy()-href)/18.01528), "%4g"%((Gas.enthalpy()- Liq.enthalpy())/18.01528), "%4g"%((Gas.enthalpy()-href)/18.01528)))+"\\\\\\hline"
print >>f,"\\end{longtable}"
print >>f,"\\end{document}"
