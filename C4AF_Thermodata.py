#!/usr/bin/env python
from chemeng import *
import chemeng.NASAdata
import csv
import math
import chemeng.cementdata

################ Thermodynamics of reactions among ..... by Zhu Jiang Guo Yang And............Energy Requirements of using oil shale in the production, .....Allaboun, Al-Otoom

registerSpecies("Ca4Fe2Al2O10", Components({'Ca':4, 'Al':2, 'Fe':2, 'O':10}))

# Before Rubbish [1574.509104, 0.07308832384, 365441.6192,-5553914.12218,-9217.080078125]
# Zhu 
class C4AFThermoData(ThermoConstantsType):
    def __init__(self, Tmin, Tmax,a, notes=""):
        ThermoConstantsType.__init__(self, Tmin, Tmax, notes) #Required 
        self.a = a
    
    def Cp0(self, T):
        return self.a[0] + self.a[1]*T + self.a[2]/(T**2)

    def Hf0(self, T):
        return self.a[0] * T + (self.a[1] * T**2)/2.0 - self.a[2]/T + self.a[3] 

    def S0(self, T):
        return self.a[0]*math.log(T) + self.a[1]*T - self.a[2]/(2.0*T**2) + self.a[4]

    def __str__(self):
        retval = "CementPolynomial{Tmin="+str(self.Tmin)+", Tmax="+str(self.Tmax)+", notes='"+self.comments+"', a=["
        for i in range(5):
            retval+=str(self.a[i])+", "
        retval = retval[:-2] + "]}"
        return retval

    def __repr__(self):
        return self.__str__()

[374.42616,0.0728016,0.0,-5.19088701e+06,  -1.72775463e+03]

      
speciesData["Ca4Fe2Al2O10"].registerPhase("Crystal")  
speciesData["Ca4Fe2Al2O10"].registerPhaseCoeffs(C4AFThermoData(298.0, 1863.0,[374.42616,0.0728016,0.0,-5.19088701e+06,  -1.72775463e+03], "G_zhu_H_thorvaldsen_Cp_bab"), "Crystal") 
speciesData["Ca4Fe2Al2O10"].registerPhaseCoeffs(C4AFThermoData(298.0, 1863.0,[  6.75447047e+02,  -1.26071150e-01,  -5.38417399e+07, -5.45638014e+06,  -3.68426527e+03], "Zhu"), "Crystal") 
speciesData["Ca4Fe2Al2O10"].registerPhaseCoeffs(C4AFThermoData(298.0, 1863.0,[374.42616,0.0728016,0.0,-5195083.7513,-1828.68305648], "Babushkin"), "Crystal") # WARNING, TF = T melt ???
 

