#!/usr/bin/env python

from chemeng import *
import chemeng.NASAdata
import scipy.optimize
import math

#Taken from     
def pVAP2(T):
    if T < 379:
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
    elif T <= 360+273.15:
        Tc = 647.096
        Pc = 22064.0 #kpa
        a = [-7.85951783, 1.84408259, -11.7866497, 22.6807411, -15.9618719, 1.80122502]
        t = 1-T/Tc
        sumval = 0
        for i,power in zip(range(6), [1,1.5,3,3.5,4,7.5]):
            sumval += a[i]*math.pow(t, power)
        return 1e3*Pc*math.exp(sumval*(Tc/T))
    else:
        raise Exception("Not valid!")

#Table by T
f=open("table.tex", "w")
complete=False
if complete:
    print >>f,"\\documentclass{report}"
    print >>f,"\\usepackage{longtable}"
    print >>f,"\\usepackage[a4paper,margin=20mm]{geometry}"
    print >>f,"\\begin{document}"
print >>f,"\\LTcapwidth=\\textwidth"
print >>f,"\\begin{longtable}{|c|c|c|c|c|c|c|c|c|}"
print >>f,"\\caption{\label{Tab:steambyT}Thermodynamic properties of saturated steam by temperature, calculated using the NASA CEA database and the vapour pressure data of Wexler or Wagner and Pruss (1990). The reference state is the triple point of saturated liquid water.}\\\\\\hline"
print >>f,"$T$ & $P$ & $C_{p,l}$ & $C_{p,v}$ & $h_l$ & $h_{lv}$ & $h_v$ & $s_l$ & $s_v$\\\\"
print >>f,"(${}^\circ$C) & (bar) & (kJ/kg~K) & (kJ/kg) & (kJ/kg) & (kJ/kg) & (kJ/kg) & (kJ/kg~K) & (kJ/kg~K)\\\\\\hline\hline"
print >>f,"\\endfirsthead"
print >>f,"\\caption*{Table~\\ref{Tab:steambyp} continued: Thermodynamic properties of saturated steam by temperature.}\\\\\\hline"
print >>f,"$T$ & $P$ & $C_{p,l}$ & $C_{p,v}$ & $h_l$ & $h_{lv}$ & $h_v$ & $s_l$ & $s_v$\\\\"
print >>f,"(${}^\circ$C) & (bar) & (kJ/kg~K) & (kJ/kg) & (kJ/kg) & (kJ/kg) & (kJ/kg) & (kJ/kg~K) & (kJ/kg~K)\\\\\\hline\hline"
print >>f,"\\endhead"
Tref=0.01
Pref = pVAP2(Tref+273.15)
Lref = IncompressiblePhase({"H2O":1}, "Liquid", T=Tref+273.15, P=Pref, molarvolume=18.01528/999.97e3)
href = Lref.enthalpy()
sref = Lref.entropy()

import numpy
#[Tref]+[274.15+1*i for i in range(30)]+[305.15+2*i for i in range(10)]+[328.15+5*i for i in range(10)] +[373.15+10*i for i in range(11)]
TCs = map(float, [Tref]+list(numpy.arange(1,11,1))+list(numpy.arange(12,22,2))+list(numpy.arange(25,105,5))+list(numpy.arange(110,210,10))+list(numpy.arange(250,350,50)))
for TC in TCs:
    T=TC+273.15
    P = pVAP2(T)
    Liq = IncompressiblePhase({"H2O":1}, "Liquid", T=T, P=P, molarvolume=18.01528/999.97e3)
    Gas = IdealGasPhase({"H2O":1}, T=T, P=P)
    print >>f," & ".join(("%4g"% TC,"%.4g"%(P/1e5), "%4g"%(Liq.Cp()/18.01528), "%4g"%(Gas.Cp()/18.01528), "%4g"%((Liq.enthalpy()-href)/18.01528), "%4g"%((Gas.enthalpy()- Liq.enthalpy())/18.01528), "%4g"%((Gas.enthalpy()-href)/18.01528), "%4g"%((Liq.entropy() - sref)/18.01528),"%4g"%((Gas.entropy()-sref)/18.01528)))+"\\\\\\hline"
print >>f,"\\end{longtable}"

import scipy.optimize as opt
Ps = map(float, [0.006117e5]+list(numpy.arange(0.010e5,0.10e5,0.005e5))+list(numpy.arange(0.12e5, 0.50e5, 0.02e5))+list(numpy.arange(0.50e5, 1.0e5, 0.05e5))+list(numpy.arange(1.0e5, 2.0e5, 0.1e5))+list(numpy.arange(2.0e5, 5.0e5, 0.5e5))+list(numpy.arange(5.0e5, 10.0e5, 1.0e5))+list(numpy.arange(10.0e5, 50.0e5, 5.0e5))+list(numpy.arange(50.0e5, 100.0e5, 10.0e5))+list(numpy.arange(100.0e5, 130.0e5, 20.0e5)))
print >>f,"\\clearpage\\begin{longtable}{|c|c|c|c|c|c|c|c|c|}"
print >>f,"\\caption{\label{Tab:steambyp}Thermodynamic properties of saturated steam by pressure, calculated using the NASA CEA database and the vapour pressure data of Wexler or Wagner and Pruss (1990). The reference state is the triple point of saturated liquid water.}\\\\\\hline"
print >>f,"$P$ & $T$ & $C_{p,l}$ & $C_{p,v}$ & $h_l$ & $h_{lv}$ & $h_v$ & $s_l$ & $s_v$\\\\"
print >>f,"(bar) & (${}^\circ$C) & (kJ/kg~K) & (kJ/kg) & (kJ/kg) & (kJ/kg) & (kJ/kg) & (kJ/kg~K) & (kJ/kg~K)\\\\\\hline\hline"
print >>f,"\\endfirsthead"
print >>f,"\\caption*{Table~\\ref{Tab:steambyp} continued: Thermodynamic properties of saturated steam by pressure.}\\\\\\hline"
print >>f,"$P$ & $T$ & $C_{p,l}$ & $C_{p,v}$ & $h_l$ & $h_{lv}$ & $h_v$ & $s_l$ & $s_v$\\\\"
print >>f,"(bar) & (${}^\circ$C) & (kJ/kg~K) & (kJ/kg) & (kJ/kg) & (kJ/kg) & (kJ/kg) & (kJ/kg~K) & (kJ/kg~K)\\\\\\hline\hline"
print >>f,"\\endhead"
for P in Ps:
    T = opt.brentq(lambda T,P: pVAP2(T)-P, 273.15, 360.0+273.15, args=(P))
    TC = T-273.15
    Liq = IncompressiblePhase({"H2O":1}, "Liquid", T=T, P=P, molarvolume=18.01528/999.97e3)
    Gas = IdealGasPhase({"H2O":1}, T=T, P=P)
    print >>f," & ".join(("%.4g"%(P/1e5),"%.4g"% TC, "%4g"%(Liq.Cp()/18.01528), "%4g"%(Gas.Cp()/18.01528), "%4g"%((Liq.enthalpy()-href)/18.01528), "%4g"%((Gas.enthalpy()- Liq.enthalpy())/18.01528), "%4g"%((Gas.enthalpy()-href)/18.01528), "%4g"%((Liq.entropy() - sref)/18.01528),"%4g"%((Gas.entropy()-sref)/18.01528)))+"\\\\\\hline"
print >>f,"\\end{longtable}"
if complete:
    print >>f,"\\end{document}"
