#!/usr/bin/env python
from chemeng import *
import chemeng.NASAdata
import csv
import math
import chemeng.cementdata


registerSpecies("Ca4Al6SO16", Components({'Ca':4, 'Al':6, 'S':1, 'O':16}))
registerSpecies("Ca5Si2SO12", Components({'Ca':5, 'Si':2, 'S':1, 'O':12}))

class PredThermoData(ThermoConstantsType):
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
      
speciesData["Ca4Al6SO16"].registerPhase("Yeelemite")            
speciesData["Ca4Al6SO16"].registerPhaseCoeffs(PredThermoData(298.0, 1800.0,[554.05, 0.14334, -11340000,-8602785.545,-2794.056201], ""), "Yeelemite")
speciesData["Ca4Al6SO16"].registerPhase("Yeelemite_VP")            
speciesData["Ca4Al6SO16"].registerPhaseCoeffs(PredThermoData(298.0, 1673.15,[  5.64804789e+02,   8.45578516e-02,  -8.89463497e+06,  -8.77732852e+06,  -2.89997330e+03], ""), "Yeelemite_VP")
speciesData["Ca5Si2SO12"].registerPhase("Ternesite")            
speciesData["Ca5Si2SO12"].registerPhaseCoeffs(PredThermoData(298.0, 1800.0,[418.684116,0.1065058693,-8400000.0,-6224806.013,-2116.497219], ""), "Ternesite")
speciesData["Ca5Si2SO12"].registerPhase("Ternesite_VP")            
speciesData["Ca5Si2SO12"].registerPhaseCoeffs(PredThermoData(298.0, 1800.0,[  3.96207710e+02,   1.07140634e-01,  -4.77030171e+07,-6.38968628e+06,  -2.06827921e+03], ""), "Ternesite_VP")