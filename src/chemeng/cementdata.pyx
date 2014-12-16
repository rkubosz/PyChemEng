#!/usr/bin/env python
from chemeng import *
import chemeng.NASAdata
import csv
import math

################ NIST data source

registerSpecies("CaAl2SiO6", Components({'Ca':1, 'Al':2, 'Si':1, 'O':6}))
registerSpecies("CaAl2Si2O8", Components({'Ca':1, 'Al':2, 'Si':2, 'O':8}))
registerSpecies("Ca2Al2Si3O10(OH)2", Components({'Ca':2, 'Al':2, 'Si':3, 'O':12, 'H':2}))
registerSpecies("Ca2Al2SiO7", Components({'Ca':2, 'Al':2, 'Si':1, 'O':7}))
registerSpecies("Ca3Al2Si3O12", Components({'Ca':3, 'Al':2, 'Si':3, 'O':12}))
registerSpecies("Al2Si4O10(OH)2", Components({'Al':2, 'Si':4, 'O':12, 'H':2}))
registerSpecies("Al2Si2O5(OH)4", Components({'Al':2, 'Si':2, 'O':9, 'H':4}))
registerSpecies("Ca2Al3Si3O12(OH)", Components({'Ca':2, 'Al':3, 'Si':3, 'O':13, 'H':1}))
registerSpecies("CaAl4Si2O10(OH)2", Components({'Ca':1, 'Al':4, 'Si':2, 'O':12, 'H':2}))
registerSpecies("CaSiO3", Components({'Ca':1, 'Si':1, 'O':3}))
registerSpecies("Ca2SiO4", Components({'Ca':2, 'Si':1, 'O':4}))
registerSpecies("Ca3SiO5", Components({'Ca':3, 'Si':1, 'O':5}))
registerSpecies("Ca3Si2O7", Components({'Ca':3, 'Si':2, 'O':7}))


class CementThermoData(ThermoConstantsType):
    def __init__(self, Tmin, Tmax, a, notes=""):
        ThermoConstantsType.__init__(self, Tmin, Tmax, notes) #Required 
        self.a = a
    
    def Cp0(self, T):
        return self.a[0] * T**(-2) + self.a[2] * T**(-0.5) + self.a[4] +2 * self.a[5] * T + self.a[6] * T ** 2

    def Hf0(self, T):
        return -self.a[0] / T + self.a[1] + 2 * self.a[2] * T**(0.5) + self.a[4] * T + self.a[5] * T**2 + self.a[6] * T ** 3 / 3

    def S0(self, T):
        return -self.a[0] / (2 * T**2) - 2 * self.a[2] * T**(-0.5) + self.a[3] + self.a[4] * math.log(T)+ 2 *  self.a[5] * T + self.a[6] * T**(2) / 2

    def __str__(self):
        retval = "CementPolynomial{Tmin="+str(self.Tmin)+", Tmax="+str(self.Tmax)+", notes='"+self.comments+"', a=["
        for i in range(7):
            retval+=str(self.a[i])+", "
        retval = retval[:-2] + "]}"
        return retval

    def __repr__(self):
        return self.__str__()


import chemeng.config
import os
with open(os.path.join(chemeng.config.datadir, 'Cement.csv'), 'rb') as datafile:
    reader = csv.reader(filter(lambda row: row[0]!='!', datafile), delimiter=',', quotechar='"')
    reader.next() #Skip the header
    for row in reader:
        #Load data
        species = row[0]
        phase = row[1]
        Tmin = float(row[6])
        Tmax = float(row[7])
        a = map(float, row[8:15])
        if len(row[3].strip()) != 0:
            Hf0 = float(row[3]) * 1000.0
            a[1] += Hf0
        notes = row[15]
        speciesData[species].registerPhase(phase)
        speciesData[species].registerPhaseCoeffs(CementThermoData(Tmin, Tmax, a, notes), phase)
        if len(row[4].strip()) != 0:
            V0 = float(row[4])



def Validate_NIST_Data(species,phase,T,Cp,S,HmHr,GmHoT,Htr,percent_error):
    if math.fabs((speciesData[species].phases[phase].constants[-1].Cp0(T) - Cp)*100.0 / Cp) > percent_error:
        print "\nError in Cp "+species+" "+phase+" @ "+str(T)+"K of " + str(abs((speciesData[species].phases[phase].constants[-1].Cp0(T) - Cp)*100.0 / Cp ))+" percent"
    if math.fabs((speciesData[species].phases[phase].constants[-1].S0(T) - S)*100.0 / S) > percent_error:
        print "\nError in S "+species+" "+phase+" @ "+str(T)+"K of " + str(abs((speciesData[species].phases[phase].constants[-1].S0(T) - S)*100.0 / S) )+" percent"
    if T!=298.15: 
        if math.fabs(((speciesData[species].phases[phase].constants[-1].Hf0(T) - Htr)  - HmHr)*100.0 / HmHr) > percent_error:
            print "\nError in H-Hr "+species+" "+phase+" @ "+str(T)+"K of " + str(abs(((speciesData[species].phases[phase].constants[-1].Hf0(T) - Htr)  - HmHr)*100.0 / HmHr))+" percent"    
    if math.fabs(((speciesData[species].phases[phase].constants[-1].Hf0(T) - T*speciesData[species].phases[phase].constants[-1].S0(T) - Htr)/T  - GmHoT)*100.0 / GmHoT) > percent_error:
        print "\nError in G-Htr/T "+species+" "+phase+" @ "+str(T)+"K of " + str(abs(((speciesData[species].Gibbs0(T, phase) - Htr)/T  - GmHoT)*100.0 / GmHoT))+" percent"    


import chemeng.config
import os
with open(os.path.join(chemeng.config.datadir, 'NistData.csv'), 'rb') as datafile:
    reader = csv.reader(filter(lambda row: row[0]!='!', datafile), delimiter=',', quotechar='"')
    reader.next() #Skip 1st line
    for row in reader:
        if row[0] == "Si" and row[1] == "Crystal":
            continue
        if row[0] == "H2O" and row[1] == "Liquid":
            continue  
        species = row[0]
        phase = row[1]
        T = float(row[2])
        Cp = float(row[3])
        S = float(row[4])
        GmHoT = float(row[5])
        HmHr = float(row[6])
        Hfe = float(row[7])
        Htr = float(row[13])
        if row[8] != "NA":
            Gfe = float(row[8])
        else:
            Gfe = row[8]
        if row[9] != "NA":
            logkfe = float(row[9])
        else:
            logkfe = row[9] 
        if row[10] != "NA":
            Hfox = float(row[10])
        else:
            Hfox = row[10]   
        if row[11] != "NA":
            Gfox = float(row[11])
        else:
            Gfox = row[11]
        if row[12] != "NA":
            logkfox = float(row[12])
        else:
            logkfox = row[12]   
        Validate_NIST_Data(species,phase,T,Cp,S,HmHr,GmHoT,Htr,1.0)

################################################ Mainly ACS
registerSpecies("Li2TiO3", Components({'Li':2, 'Ti':1,  'O':3}))
registerSpecies("TiS2", Components({'Ti':1, 'S':2}))
registerSpecies("FeCl3", Components({'Fe':1, 'Cl':3}))
registerSpecies("Ca3B2O6", Components({'Ca':3, 'B':2 , 'O':6}))
registerSpecies("Ca2B2O5", Components({'Ca':2, 'B':2 , 'O':5}))
registerSpecies("CaB2O4", Components({'Ca':1, 'B':2 , 'O':4}))
registerSpecies("CaB4O7", Components({'Ca':1, 'B':4 , 'O':7}))
registerSpecies("Na2TiO3", Components({'Na':2, 'Ti':1 , 'O':3}))
registerSpecies("Fe2SiO4", Components({'Fe':2, 'Si':1, 'O':4}))
registerSpecies("CaTiSiO5", Components({'Ca':1, 'Ti':1, 'Si':1,'O':5}))
registerSpecies("Ca3Al2O6", Components({'Ca':3, 'Al':2, 'O':6}))
registerSpecies("Ca12Al14O33", Components({'Ca':12, 'Al':14, 'O':33}))
registerSpecies("CaAl2O4", Components({'Ca':1, 'Al':2, 'O':4}))
registerSpecies("CaAl4O7", Components({'Ca':1, 'Al':4, 'O':7}))
registerSpecies("CaFe2O4", Components({'Ca':1, 'Fe':2, 'O':4}))
registerSpecies("Ca2Fe2O5", Components({'Ca':2, 'Fe':2, 'O':5}))
registerSpecies("MgFe2O4", Components({'Mg':1, 'Fe':2, 'O':4}))
registerSpecies("Al2TiO5", Components({'Al':2, 'Ti':1, 'O':5}))
registerSpecies("Fe2TiO4", Components({'Fe':2, 'Ti':1, 'O':4}))
registerSpecies("Zn2TiO4", Components({'Zn':2, 'Ti':1, 'O':4}))
registerSpecies("CaTiO3", Components({'Ca':1, 'Ti':1,  'O':3}))
registerSpecies("FeTiO3", Components({'Fe':1, 'Ti':1,  'O':3}))

class CemThermoData(ThermoConstantsType):
    def __init__(self, Tmin, Tmax, a, notes=""):
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


def Validate_Cem(species,phase,a,H298,S298,T,error,HTmH298,STmS298):
    H = a[0] * T + (a[1] * T**2)/2.0 - a[2]/T + a[3] 
    S = a[0]*math.log(T) + a[1]*T - a[2]/(2.0*T**2) + a[4]
    if math.fabs(((H - H298) - HTmH298)/HTmH298 ) > error:
        print "Error for H at ",T, species,phase,math.fabs(((H - H298) - HTmH298)/HTmH298 )*100.0, " percent"  
    if math.fabs(((S - S298) - STmS298)/STmS298 ) > error:
        print "Error for H at ",T, species,phase,  math.fabs(((S - S298) - STmS298)/STmS298 )*100.0, " percent"  
  
    
with open(os.path.join(chemeng.config.datadir,'Cement_Therm_New2.csv'), 'rb') as datafile:
    reader = csv.reader(filter(lambda row: row[0]!='!', datafile), delimiter=',', quotechar='"')
    reader.next() #Skip 1st line
    reader.next() #Skip 2nd line
    for row in reader:
        if row[0] == "Exit":
            break
        else:
            if row[2] and row[4] and row[6] and row[7] and row[8] and row[9] and row[10] and row[11]:
                #Load data
                species = row[0]
                phase = row[1]
                S0 = float(row[2])*4.184
                Hf0 = float(row[4])*4.184
                a = float(row[6])*4.184
                b = float(row[7])*2.0*4.184e-3
                c = float(row[8])*-4.184e5
                d = float(row[9])*4.184 + Hf0
                Tmin = float(row[10])
                Tmax = float(row[11])
                notes = "S0:"+str(row[3])+" Hf0:"+str(row[5])+" HT-H298:"+str(row[12])
                if Tmin<=298.15<=Tmax:
                    e = S0 - (a*math.log(298.15) + b*298.15 - c/(2.0*298.15**2))
                else:
                    e = float(row[13])*4.184 - ( a*math.log(Tmin) + b*Tmin - c/(2.0*Tmin**2)) + S0                   
                with open(os.path.join(chemeng.config.datadir,'Cement_New_Tests.csv'), 'rb') as Testfile:
                    Testreader = csv.reader(filter(lambda row: row[0]!='!', Testfile), delimiter=',', quotechar='"')
                    Testreader.next() #Skip 1st line 
                    for Testrow in Testreader:
                        if Testrow[0] == species and Testrow[1] == phase:
                            Validate_Cem(species,phase,[a,b,c,d,e],Hf0,S0,float(Testrow[2]),0.03,float(Testrow[3])*4.184,float(Testrow[4])*4.184)      
                speciesData[species].registerPhase(phase)
                speciesData[species].registerPhaseCoeffs(CemThermoData(Tmin, Tmax, [a,b,c,d,e], notes), phase)


############# Holland, T. J. B., and R. Powell. "An improved and extended internally consistent thermodynamic dataset for phases of petrological interest, involving a new equation of state for solids." Journal of Metamorphic Geology 29.3 (2011): 333-383.
registerSpecies("Ca5Si2CO11", Components({'Ca':5, 'Si':2, 'C':1, 'O':11}))
registerSpecies("CaMgC2O6", Components({'Ca':1, 'Mg':1, 'C':2, 'O':6}))
registerSpecies("Ca5Si2C2O13", Components({'Ca':5, 'Si':2, 'C':2, 'O':13}))
registerSpecies("NaAlSi3O8", Components({'Na':1, 'Al':1, 'Si':3, 'O':8}))
registerSpecies("CaMgSi2O6", Components({'Ca':1, 'Mg':1, 'Si':2, 'O':6}))
registerSpecies("Mg2Si2O6", Components({ 'Mg':2, 'Si':2, 'O':6}))
registerSpecies("Fe2SiO4", Components({ 'Fe':2, 'Si':1, 'O':4}))
registerSpecies("K3Fe0.5Al4Si19.5O47", Components({ 'K':3.0,'Fe':0.5, 'Si':19.5, 'O':47,'Al':4.0}))
registerSpecies("KAlSi3O8", Components({ 'K':1, 'Si':3, 'O':8,'Al':1.0}))
registerSpecies("K3Mg0.5Al4Si19.5O47", Components({ 'K':3.0,'Mg':0.5, 'Si':19.5, 'O':47,'Al':4.0}))
registerSpecies("NaFeSi2O6", Components({ 'Na':1, 'Fe':1, 'Si':2,'O':6}))
registerSpecies("NaAlSi2O6", Components({'Na':1, 'Al':1, 'Si':2, 'O':6}))
registerSpecies("Fe2Si2O6", Components({'Fe':2, 'Si':2, 'O':6}))

class HPThermoData(ThermoConstantsType):
    def __init__(self, Tmin, Tmax, a, notes=""):
        ThermoConstantsType.__init__(self, Tmin, Tmax, notes) #Required 
        self.a = a
    
    def Cp0(self, T):
        return self.a[0] + self.a[1] * T + self.a[2]/(T**2) + self.a[3]/(T**0.5)

    def Hf0(self, T):
        return self.a[0]*T + 0.5*self.a[1]*T**2.0 - self.a[2]/T + 2.0*self.a[3]*T**0.5 + self.a[4]

    def S0(self, T):
        return self.a[0]*math.log(T) + self.a[1]*T - 0.5*self.a[2]/(T**2.0) - 2.0*self.a[3]/(T**0.5) + self.a[5] 

    def __str__(self):
        retval = "CementPolynomial{Tmin="+str(self.Tmin)+", Tmax="+str(self.Tmax)+", notes='"+self.comments+"', a=["
        for i in range(6):
            retval+=str(self.a[i])+", "
        retval = retval[:-2] + "]}"
        return retval

    def __repr__(self):
        return self.__str__()
      
speciesData["Ca5Si2CO11"].registerPhase("Spurrite")            
speciesData["Ca5Si2CO11"].registerPhaseCoeffs(HPThermoData(298.0, 2000.0,[614.1,-0.003508,-2493100.0,-4168.0,-5897053.0,-3664.64123237], "37"), "Spurrite")

speciesData["CaMgC2O6"].registerPhase("Dolomite")            
speciesData["CaMgC2O6"].registerPhaseCoeffs(HPThermoData(298.0, 2000.0,[358.9,-0.004905,0.0,-3456.2,-2311991.47077,-2287.7288762], "37"), "Dolomite")

speciesData["Ca5Si2C2O13"].registerPhase("Tilleyite")            
speciesData["Ca5Si2C2O13"].registerPhaseCoeffs(HPThermoData(298.0, 2000.0,[741.7,-0.005345,-1434600,-5878.5,-6391033.598,-4523.276465], "37"), "Tilleyite")

speciesData["NaAlSi3O8"].registerPhase("Liquid")            
speciesData["NaAlSi3O8"].registerPhaseCoeffs(HPThermoData(298.0, 2000.0,[358.5,0.0,0.0,0.0,-4031406.775,-1890.58842253], "37"), "Liquid")

speciesData["NaAlSi3O8"].registerPhase("Albite")            
speciesData["NaAlSi3O8"].registerPhaseCoeffs(HPThermoData(298.0, 2000.0,[452.0,-0.013364,-1275900.0,-3953.6,-3936515.41852,-2826.44236975], "37"), "Albite")

speciesData["NaAlSi3O8"].registerPhase("highAlbite")            
speciesData["NaAlSi3O8"].registerPhaseCoeffs(HPThermoData(298.0, 2000.0,[452.0,-0.013364,-1275900.0,-3953.6,-3926755.41852,-2813.04236975], "37"), "highAlbite")

speciesData["CaAl2Si2O8"].registerPhase("Liquid")            
speciesData["CaAl2Si2O8"].registerPhaseCoeffs(HPThermoData(298.0, 2000.0,[417.5,0.0,0.0,0.0,-4406117.625,-2339.74662875], "37"), "Liquid")

speciesData["CaMgSi2O6"].registerPhase("Liquid")            
speciesData["CaMgSi2O6"].registerPhaseCoeffs(HPThermoData(298.0, 2000.0,[345.3,0.0,0.0,0.0,-3312281.195,-1944.38014589], "37"), "Liquid")

speciesData["CaMgSi2O6"].registerPhase("Diopside")            
speciesData["CaMgSi2O6"].registerPhaseCoeffs(HPThermoData(298.0, 2000.0,[314.5,0.000041,-2745900,-2020.1,-3235757.57446,-1898.63491189], "37"), "Diopside")

speciesData["Mg2Si2O6"].registerPhase("Liquid")            
speciesData["Mg2Si2O6"].registerPhaseCoeffs(HPThermoData(298.0, 2000.0,[354.9,0.0,0.0,0.0,-3191463.435,-2021.07707436], "37"), "Liquid")

speciesData["Mg2Si2O6"].registerPhase("Enstatite")            
speciesData["Mg2Si2O6"].registerPhaseCoeffs(HPThermoData(298.0, 2000.0,[356.2,-0.002990,-596900.0,-3185.3,-3088328.86538,-2268.39597837], "37"), "Enstatite")

speciesData["Fe2SiO4"].registerPhase("Liquid")            
speciesData["Fe2SiO4"].registerPhaseCoeffs(HPThermoData(298.0, 2000.0,[239.7,0.0,0.0,0.0,-1506736.555,-1247.71393272], "37"), "Liquid")

speciesData["Mg2SiO4"].registerPhase("Liquid")            
speciesData["Mg2SiO4"].registerPhaseCoeffs(HPThermoData(298.0, 2000.0,[267.9,0.0,0.0,0.0,-2242104.385,-1552.3861601], "37"), "Liquid")

speciesData["K3Fe0.5Al4Si19.5O47"].registerPhase("Liquid")            
speciesData["K3Fe0.5Al4Si19.5O47"].registerPhaseCoeffs(HPThermoData(298.0, 2000.0,[2375.0,0.0,0.0,0.0,-23632876.25,-13311.7921995], "37"), "Liquid")

speciesData["K3Mg0.5Al4Si19.5O47"].registerPhase("Liquid")            
speciesData["K3Mg0.5Al4Si19.5O47"].registerPhaseCoeffs(HPThermoData(298.0, 2000.0,[2386.0,0.0,0.0,0.0,-23819995.9,-13410.4657633], "37"), "Liquid")

speciesData["KAlSi3O8"].registerPhase("Liquid")            
speciesData["KAlSi3O8"].registerPhaseCoeffs(HPThermoData(298.0, 2000.0,[367.3,0.0,0.0,0.0,-4081390.495,-1949.72727363], "37"), "Liquid")

speciesData["KAlSi3O8"].registerPhase("Sanidine")            
speciesData["KAlSi3O8"].registerPhaseCoeffs(HPThermoData(298.0, 2000.0,[448.8,-0.010075,-1007300.0,-3973.1,-3964433.22114,-2789.93851863], "37"), "Sanidine")

speciesData["Al2SiO5"].registerPhase("Liquid")            
speciesData["Al2SiO5"].registerPhaseCoeffs(HPThermoData(298.0, 2000.0,[237.6,0.0,0.0,0.0,-2576200.44,-1286.74897962], "37"), "Liquid")

speciesData["NaAlSi2O6"].registerPhase("Jadeite")            
speciesData["NaAlSi2O6"].registerPhaseCoeffs(HPThermoData(298.0, 2000.0,[301.1,0.010143,-2239300.0,-2055.1,-3054593.52897,-1835.70351849], "37"), "Jadeite")

speciesData["NaFeSi2O6"].registerPhase("Acmite")            
speciesData["NaFeSi2O6"].registerPhaseCoeffs(HPThermoData(298.0, 2000.0,[307.1,0.016758,-1685500.0,-2125.8,-2611237.43295,-1839.83541625], "37"), "Acmite")

speciesData["Fe2SiO4"].registerPhase("Fayalite")            
speciesData["Fe2SiO4"].registerPhaseCoeffs(HPThermoData(298.0, 2000.0,[201.1,0.01733,-1960600.0,-900.9,-1514412.39343,-1115.33066401], "37"), "Fayalite")

speciesData["Fe2Si2O6"].registerPhase("Ferrosilite")            
speciesData["Fe2Si2O6"].registerPhaseCoeffs(HPThermoData(298.0, 2000.0,[398.7, -0.006579, 1290100.0, -4058.0,-2362863.83439,-2541.84281492], "37"), "Ferrosilite")

############## Mullite Thermo Data #Mullite data taken from Thermodynamic properties of mullite.... (Waldbaum) Harvard 1965
registerSpecies("Al6Si2O13", Components({'Al':6, 'Si':2, 'O':13}))

class MulliteThermoData(ThermoConstantsType):
    def __init__(self, Tmin, Tmax, notes=""):
        ThermoConstantsType.__init__(self, Tmin, Tmax, notes) #Required 
        self.a = [ -6.34719844e+06,  -6.96694385e+06,  -6.89076313e+02, -2.08446637e+03,   3.77248082e+02,   1.30438133e-01,-8.70793049e-05]
    
    def Cp0(self, T):
        return self.a[0] * T**(-2) + self.a[2] * T**(-0.5) + self.a[4] +2 * self.a[5] * T + self.a[6] * T ** 2

    def Hf0(self, T):
        return -self.a[0] / T + self.a[1] + 2 * self.a[2] * T**(0.5) + self.a[4] * T + self.a[5] * T**2 + self.a[6] * T ** 3 / 3

    def S0(self, T):
        return -self.a[0] / (2 * T**2) - 2 * self.a[2] * T**(-0.5) + self.a[3] + self.a[4] * math.log(T)+ 2 *  self.a[5] * T + self.a[6] * T**(2) / 2

    def __str__(self):
        retval = "CementPolynomial{Tmin="+str(self.Tmin)+", Tmax="+str(self.Tmax)+", notes='"+self.comments+"', a=["
        for i in range(7):
            retval+=str(self.a[i])+", "
        retval = retval[:-2] + "]}"
        return retval

    def __repr__(self):
        return self.__str__()
      
speciesData["Al6Si2O13"].registerPhase("Mullite")            
speciesData["Al6Si2O13"].registerPhaseCoeffs(MulliteThermoData(298.0, 2000.0, ""), "Mullite")

