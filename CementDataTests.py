#!/usr/bin/env python
from chemeng import *
#import chemeng.NASAdata
import chemeng.cementdata
import math
from cementdata2 import *

def Standard_States_Cement_enthalpy(element,T):
    if element == "Al" and T <= 933:
        return speciesData[element].phases["Crystal"].constants[-1].Hf0(T)
    if element == "Al" and T > 933:
        return speciesData[element].phases["Liquid"].constants[-1].Hf0(T)
    if element == "Ca" and T <= 720:
        return speciesData[element].phases["alpha"].constants[-1].Hf0(T)
    if element == "Ca" and  720 < T <= 1112:
        return speciesData[element].phases["beta"].constants[-1].Hf0(T)
    if element == "Ca" and  1112 < T <= 1755:
        return speciesData[element].phases["Liquid"].constants[-1].Hf0(T)
    if element == "Ca" and  1755 < T <= 1800:
        return speciesData[element].phases["Gas"].constants[-1].Hf0(T)  
    if element == "Si" and T <= 1685:
        return speciesData[element].phases["1cr"].constants[-1].Hf0(T) # Using Nasa here
    if element == "Si" and T > 1685:
        return speciesData[element].phases["Liquid"].constants[-1].Hf0(T)
    if element == "O":
        return speciesData["O2"].phases["Gas"].constants[-1].Hf0(T)/2.0
    if element == "H":
        return speciesData["H2"].phases["Gas"].constants[-1].Hf0(T)/2.0

def Standard_States_enthalpy(element,T, P = 1e5): #WARNING some dont reach reaquired T range, e.g. Standard_States_enthalpy("Al",1000.0)
    if element in ["O","N","H","F","Cl"] :
        return IdealGasPhase({element+"2":0.5}, T=T, P=P).enthalpy()
    elif element == "Br":
        return IncompressiblePhase({"Br2":0.5}, T=T, P=P, phaseID="Liquid").enthalpy()
    elif element == "I":
        return IncompressiblePhase({"I2":0.5}, T=T, P=P, phaseID="1cr").enthalpy()
    elif element in ["He","Ne","Ar","Kr","Xe","Rn"]:
        return IdealGasPhase({element:1.0}, T=T, P=P).enthalpy()
    elif element == "Hg":
        return IncompressiblePhase({"Hg":1.0}, T=T, P=P, phaseID="Liquid").enthalpy()
    elif element == "B":
        return IncompressiblePhase({"B":1.0}, T=T, P=P, phaseID="1b").enthalpy()
    elif element == "C":
        return IncompressiblePhase({"C":1.0}, T=T, P=P, phaseID="1gr").enthalpy() 
    elif element in ["Be","S","Ca","Sc","Ti","Mn","Fe","Co","Sr","Zr"]:
        return IncompressiblePhase({element:1.0}, T=T, P=P, phaseID="1a").enthalpy()
    elif element in ["Li","Na","Al","Mg","Si","P","K","V","Cr","Ni","Cu","Zn","Ga","Rb","Nb","Mo","Ag","Cd","In","Sn","Cs","Ba","Ta","W","Pb"]:
        return IncompressiblePhase({element:1.0}, T=T, P=P, phaseID="1cr").enthalpy()
    else:
        print "Sorry, no data for element: ", element
      

def Validate_Cp(Formula,phase,T,expval,error = 0.01):
    Cp_sim = speciesData[Formula].phases[phase].constants[-1].Cp0(T)
    if math.fabs(1.0 - Cp_sim/expval) > error:
        print "WARNING: Cp sim ",Formula," , ",phase," at , ", T ," = " , Cp_sim, " compared to ", expval
    
def Validate_S(Formula,phase,T,expval,error = 0.01):
    S_sim =  speciesData[Formula].phases[phase].constants[-1].S0(T)
    if math.fabs(1.0 - S_sim/expval) > error:
        print "WARNING: S sim ",Formula,", ",phase," at , ", T ," = " , S_sim, " compared to ", expval
        
def Validate_Hf(Formula,phase,T,expval,error = 100,verbose = False): #0.01%
    Hf_sim = speciesData[Formula].phases[phase].constants[-1].Hf0(T) + sum(-Standard_States_Cement_enthalpy(keys,T)*vals for keys,vals in Components({Formula:1.0}).elementalComposition().iteritems())
    if verbose:
        print "Hf at ",T, " for ", Formula," ",phase, "sim =",Hf_sim ,"exp = ",expval, "diff = ", Hf_sim - expval
    if math.fabs(Hf_sim- expval) > error:
        print "WARNING: Hf sim ",Formula,",",phase,"at ,", T ," =" , Hf_sim, " compared to ", expval, ",diff = ", Hf_sim - expval

def Validate_dH(Formula,phase,T,expval,error = 0.01,verbose = False):
    Hf_sim =  speciesData[Formula].phases[phase].constants[-1].Hf0(T) - speciesData[Formula].phases[phase].constants[-1].Hf0(298.15)
    if verbose:
        print Hf_sim,expval,Hf_sim-expval
    if math.fabs(1.0 - Hf_sim/expval) > error:
        print "WARNING: dH sim ",Formula,", ",phase," at , ", T ," = " , Hf_sim, " compared to ", expval, ",diff = ", Hf_sim - expval

############################################################
Validate_Hf("Al","Crystal",298.15,0.0000000000000000000001)  

Validate_Cp("Al","Crystal",300.0,24.338)
Validate_S("Al","Crystal",300.0,28.500)
Validate_dH("Al","Crystal",300.0,45.0)
Validate_Hf("Al","Crystal",300,0.0000000000000000000001)

Validate_Cp("Al","Crystal",550.0,27.255)
Validate_S("Al","Crystal",550.0,44.094)
Validate_dH("Al","Crystal",550.0,6512)

Validate_Cp("Al","Liquid",1300.0,31.756)
Validate_S("Al","Liquid",1300.0,81.939)
#Validate_dH("Al","Liquid",1300.0,40472)
Validate_Hf("Al","Liquid",1300,0.0000000000000000000001)

############################################################

Validate_Hf("Ca","alpha",298.15,0.00000000000000001)

Validate_Cp("Ca","alpha",300.0,25.354)
Validate_S("Ca","alpha",300.0,41.773)
Validate_dH("Ca","alpha",300.0,47.0)
Validate_Hf("Ca","alpha",300.0,0.00000000000000001)

Validate_Cp("Ca","alpha",550.0,28.560)
Validate_S("Ca","alpha",550.0,57.858)
Validate_dH("Ca","alpha",550.0,6721)
Validate_Hf("Ca","alpha",550.0,0.00000000000000001)

Validate_Cp("Ca","beta",750.0,30.581)
Validate_S("Ca","beta",750.0,68.512)
#Validate_dH("Ca","beta",750.0,13700)
Validate_Hf("Ca","beta",750.0,0.00000000000000001)

Validate_Cp("Ca","beta",1100.0,45.029)
Validate_S("Ca","beta",1100.0,82.816)
#Validate_dH("Ca","beta",1100.0,26933)
Validate_Hf("Ca","beta",1100.0,0.00000000000000001)

Validate_Cp("Ca","Liquid",1300.0,29.275)
Validate_S("Ca","Liquid",1300.0,95.541)
#Validate_dH("Ca","Liquid",1300.0,41499)
Validate_Hf("Ca","Liquid",1300.0,0.00000000000000001)

Validate_Cp("Ca","Gas",1755.0,20.851)
Validate_S("Ca","Gas",1755.0,191.628)
#Validate_dH("Ca","Gas",1755.0,208032)
Validate_Hf("Ca","Gas",1755.0001,0.00000000000000001)

############################################################
# Continue Using Nasa it matches cement experimental

Validate_Cp("Si","1cr",550.0,23.970)
Validate_S("Si","1cr",550.0,32.365)
Validate_dH("Si","1cr",550.0,5626)
Validate_Hf("Si","1cr",550.0,0.00000000000000001)
Validate_Cp("Si","Liquid",1700.0,25.522)
Validate_S("Si","Liquid",1700.0,92.031)
Validate_Hf("Si","Liquid",1755.0001,0.00000000000000001)

############################################################

Validate_Cp("O2","Gas",550.0,31.722)
Validate_S("O2","Gas",550.0,223.664)
Validate_dH("O2","Gas",550.0,7694)
Validate_Hf("O2","Gas",550.0,0.00000000000000001)

Validate_Cp("O2","Gas",1700.0,36.826)
Validate_S("O2","Gas",1700.0,262.581)
Validate_dH("O2","Gas",1700.0,47914)
Validate_Hf("O2","Gas",1700,0.00000000000000001)

############################################################

Validate_Cp("H2","Gas",550.0,29.263)
Validate_S("H2","Gas",550.0,148.429)
Validate_dH("H2","Gas",550.0,7350)
Validate_Hf("H2","Gas",550.0,0.00000000000000001)

Validate_Cp("H2","Gas",1700.0,33.099)
Validate_S("H2","Gas",1700.0,182.820)
Validate_dH("H2","Gas",1700.0,42812)
Validate_Hf("H2","Gas",1700,0.00000000000000001)

############################################################
Validate_Cp("HAlO2","Diaspore",300.0,53.374)
Validate_S("HAlO2","Diaspore",300.0,35.668)
Validate_dH("HAlO2","Diaspore",300.0,98)
Validate_Hf("HAlO2","Diaspore",300.0,-999484) 

Validate_Cp("HAlO2","Diaspore",550.0,77.591)
Validate_S("HAlO2","Diaspore",550.0,75.645)
Validate_dH("HAlO2","Diaspore",550.0,16890)
Validate_Hf("HAlO2","Diaspore",550.0,-1000447) 

Validate_Cp("HAlO2","Diaspore",750.0,87.817)
Validate_S("HAlO2","Diaspore",750.0,101.333)
Validate_dH("HAlO2","Diaspore",750.0,33507)
Validate_Hf("HAlO2","Diaspore",750.0,-998984) 
############################################################
Validate_Hf("HAlO2","Boehmite",298.15,-990424)

Validate_Cp("HAlO2","Boehmite",300.0,65.846)
Validate_S("HAlO2","Boehmite",300.0,48.84)
Validate_dH("HAlO2","Boehmite",300.0,121)
Validate_Hf("HAlO2","Boehmite",300.0,-990428) 

Validate_Cp("HAlO2","Boehmite",550.0,98.917)
Validate_S("HAlO2","Boehmite",550.0,99.01)
Validate_dH("HAlO2","Boehmite",550.0,21230)
Validate_Hf("HAlO2","Boehmite",550.0,-987075) 

Validate_Cp("HAlO2","Boehmite",750.0,113.611)
Validate_S("HAlO2","Boehmite",750.0,132.014)
Validate_dH("HAlO2","Boehmite",750.0,42587)
Validate_Hf("HAlO2","Boehmite",750.0,-980872) 
############################################################
Validate_Hf("Al(OH)3","Gibbsite",298.15,-1293334)

Validate_Cp("Al(OH)3","Gibbsite",300.0,92.226)
Validate_S("Al(OH)3","Gibbsite",300.0,69.009)
Validate_dH("Al(OH)3","Gibbsite",300.0,170)
Validate_Hf("Al(OH)3","Gibbsite",300.0,-1293370)

Validate_Cp("Al(OH)3","Gibbsite",550.0,142.361)
Validate_S("Al(OH)3","Gibbsite",550.0,139.954)
Validate_dH("Al(OH)3","Gibbsite",550.0,30061)
Validate_Hf("Al(OH)3","Gibbsite",550.0,-1292350)

Validate_Cp("Al(OH)3","Gibbsite",750.0,169.716)
Validate_S("Al(OH)3","Gibbsite",750.0,188.303)
Validate_dH("Al(OH)3","Gibbsite",750.0,61379)
Validate_Hf("Al(OH)3","Gibbsite",750.0,-1285312)
############################################################
Validate_Hf("CaAl2SiO6","Ca-Al Clinopyroxene",298.15,-3298956)

Validate_Cp("CaAl2SiO6","Ca-Al Clinopyroxene",550.0,220.652)
Validate_S("CaAl2SiO6","Ca-Al Clinopyroxene",550.0,261.723)
Validate_dH("CaAl2SiO6","Ca-Al Clinopyroxene",550.0,50093)
Validate_Hf("CaAl2SiO6","Ca-Al Clinopyroxene",550.0,-3297315)

Validate_Cp("CaAl2SiO6","Ca-Al Clinopyroxene",1050.0,252.925)
Validate_S("CaAl2SiO6","Ca-Al Clinopyroxene",1050.0,415.727)
Validate_dH("CaAl2SiO6","Ca-Al Clinopyroxene",1050.0,170028)
Validate_Hf("CaAl2SiO6","Ca-Al Clinopyroxene",1050,-3310493)

Validate_Cp("CaAl2SiO6","Ca-Al Clinopyroxene",1650.0,268.038)
Validate_S("CaAl2SiO6","Ca-Al Clinopyroxene",1650.0,533.626)
Validate_dH("CaAl2SiO6","Ca-Al Clinopyroxene",1650.0,326876)
Validate_Hf("CaAl2SiO6","Ca-Al Clinopyroxene",1650.0,-3300141)

############################################################
Validate_Hf("CaAl2Si2O8","Anorthite",298.15,-4227833)

Validate_Cp("CaAl2Si2O8","Anorthite",550.0,279.430)
Validate_S("CaAl2Si2O8","Anorthite",550.0,351.246)
Validate_dH("CaAl2Si2O8","Anorthite",550.0,63366)
Validate_Hf("CaAl2Si2O8","Anorthite",550.0,-4226238)

Validate_Cp("CaAl2Si2O8","Anorthite",1050.0,321.182)
Validate_S("CaAl2Si2O8","Anorthite",1050.0,546.185)
Validate_dH("CaAl2Si2O8","Anorthite",1050.0,215181)
Validate_Hf("CaAl2Si2O8","Anorthite",1050.0,-4237076)

Validate_Cp("CaAl2Si2O8","Anorthite",1550.0,362.132)
Validate_S("CaAl2Si2O8","Anorthite",1550.0,677.944)
Validate_dH("CaAl2Si2O8","Anorthite",1550.0,384985)
Validate_Hf("CaAl2Si2O8","Anorthite",1550.0,-4222328)

############################################################
Validate_Hf("Ca2Al2Si3O10(OH)2","Prehnite",298.15,-6193631)

Validate_Cp("Ca2Al2Si3O10(OH)2","Prehnite",550.0,441.543)
Validate_S("Ca2Al2Si3O10(OH)2","Prehnite",550.0,531.382)
Validate_dH("Ca2Al2Si3O10(OH)2","Prehnite",550.0,99563)
Validate_Hf("Ca2Al2Si3O10(OH)2","Prehnite",550.0,-6190923)

Validate_Cp("Ca2Al2Si3O10(OH)2","Prehnite",1050.0,501.798)
Validate_S("Ca2Al2Si3O10(OH)2","Prehnite",1050.0,840.083)
Validate_dH("Ca2Al2Si3O10(OH)2","Prehnite",1050.0,339866)
Validate_Hf("Ca2Al2Si3O10(OH)2","Prehnite",1050.0,-6192430)

############################################################
Validate_Hf("Ca2Al2SiO7","Gehlenite",298.15,-3981707)

Validate_Cp("Ca2Al2SiO7","Gehlenite",550.0,263.707)
Validate_S("Ca2Al2SiO7","Gehlenite",550.0,354.815)
Validate_dH("Ca2Al2SiO7","Gehlenite",550.0,60351)
Validate_Hf("Ca2Al2SiO7","Gehlenite",550.0,-3980375)

Validate_Cp("Ca2Al2SiO7","Gehlenite",1050.0,297.965)
Validate_S("Ca2Al2SiO7","Gehlenite",1050.0,537.664)
Validate_dH("Ca2Al2SiO7","Gehlenite",1050.0,202637)
Validate_Hf("Ca2Al2SiO7","Gehlenite",1050.0,-3997596)

############################################################
Validate_Hf("Ca3Al2Si3O12","Grossular",298.15,-6636338)

Validate_Cp("Ca3Al2Si3O12","Grossular",550.0,438.472)
Validate_S("Ca3Al2Si3O12","Grossular",550.0,494.171)
Validate_dH("Ca3Al2Si3O12","Grossular",550.0,99347)
Validate_Hf("Ca3Al2Si3O12","Grossular",550.0,-6633217)

Validate_Cp("Ca3Al2Si3O12","Grossular",1050.0,490.538)
Validate_S("Ca3Al2Si3O12","Grossular",1050.0,797.963)
Validate_dH("Ca3Al2Si3O12","Grossular",1050.0,335622)
Validate_Hf("Ca3Al2Si3O12","Grossular",1050.0,-6641917)

############################################################
Validate_Hf("Al2Si4O10(OH)2","Pyrophyllite",298.15,-5642023)

Validate_Cp("Al2Si4O10(OH)2","Pyrophyllite",550.0,402.465)
Validate_S("Al2Si4O10(OH)2","Pyrophyllite",550.0,455.069)
Validate_dH("Al2Si4O10(OH)2","Pyrophyllite",550.0,90070)
Validate_Hf("Al2Si4O10(OH)2","Pyrophyllite",550.0,-5640993)

Validate_Cp("Al2Si4O10(OH)2","Pyrophyllite",950.0,491.469)
Validate_S("Al2Si4O10(OH)2","Pyrophyllite",950.0,697.168)
Validate_dH("Al2Si4O10(OH)2","Pyrophyllite",950.0,268821)
Validate_Hf("Al2Si4O10(OH)2","Pyrophyllite",950.0,-5639980)

############################################################
Validate_Hf("Al2Si2O5(OH)4","Dickite",298.15,-4118475)

Validate_Cp("Al2Si2O5(OH)4","Dickite",550.0,327.337)
Validate_S("Al2Si2O5(OH)4","Dickite",550.0,373.166)
Validate_dH("Al2Si2O5(OH)4","Dickite",550.0,73556)
Validate_Hf("Al2Si2O5(OH)4","Dickite",550.0,-4118516)

Validate_Cp("Al2Si2O5(OH)4","Dickite",950.0,348.591)
Validate_S("Al2Si2O5(OH)4","Dickite",950.0,560.983)
Validate_dH("Al2Si2O5(OH)4","Dickite",950.0,211396)
Validate_Hf("Al2Si2O5(OH)4","Dickite",950.0,-4130058)

############################################################
Validate_Hf("Al2Si2O5(OH)4","Halloysite",298.15,-4101028)

Validate_Cp("Al2Si2O5(OH)4","Halloysite",550.0,326.629)
Validate_S("Al2Si2O5(OH)4","Halloysite",550.0,380.696)
Validate_dH("Al2Si2O5(OH)4","Halloysite",550.0,73988)
Validate_Hf("Al2Si2O5(OH)4","Halloysite",550.0,-4100636)

Validate_Cp("Al2Si2O5(OH)4","Halloysite",950.0,353.307)
Validate_S("Al2Si2O5(OH)4","Halloysite",950.0,568.841)
Validate_dH("Al2Si2O5(OH)4","Halloysite",950.0,212166)
Validate_Hf("Al2Si2O5(OH)4","Halloysite",950.0,-4111840)


############################################################
Validate_Hf("Al2Si2O5(OH)4","Kaolinite",298.15,-4119780)

Validate_Cp("Al2Si2O5(OH)4","Kaolinite",550.0,326.624)
Validate_S("Al2Si2O5(OH)4","Kaolinite",550.0,382.631)
Validate_dH("Al2Si2O5(OH)4","Kaolinite",550.0,74100)
Validate_Hf("Al2Si2O5(OH)4","Kaolinite",550.0,-4119277)

Validate_Cp("Al2Si2O5(OH)4","Kaolinite",950.0,353.585)
Validate_S("Al2Si2O5(OH)4","Kaolinite",950.0,570.744)
Validate_dH("Al2Si2O5(OH)4","Kaolinite",950.0,212259)
Validate_Hf("Al2Si2O5(OH)4","Kaolinite",950.0,-4130500)

############################################################
Validate_Hf("Al2O3","Corundum",298.15,-1675711)

Validate_Cp("Al2O3","Corundum",550.0,109.195)
Validate_S("Al2O3","Corundum",550.0,109.425)
Validate_dH("Al2O3","Corundum",550.0,24443)
Validate_Hf("Al2O3","Corundum",550.0,-1675833)

Validate_Cp("Al2O3","Corundum",950.0,124.323)
Validate_S("Al2O3","Corundum",950.0,173.663)
Validate_dH("Al2O3","Corundum",950.0,71728)
Validate_Hf("Al2O3","Corundum",950.0,-1694168)

############################################################
Validate_Hf("Al2SiO5","Kyanite",298.15,-2594269)

Validate_Cp("Al2SiO5","Kyanite",300.0,122.936)
Validate_S("Al2SiO5","Kyanite",300.0,85.224)
Validate_dH("Al2SiO5","Kyanite",300.0,227)
Validate_Hf("Al2SiO5","Kyanite",300.0,-2594305)


Validate_Cp("Al2SiO5","Kyanite",550.0,170.145)
Validate_S("Al2SiO5","Kyanite",550.0,175.066)
Validate_dH("Al2SiO5","Kyanite",550.0,37871)
Validate_Hf("Al2SiO5","Kyanite",550.0,-2594282)

Validate_Cp("Al2SiO5","Kyanite",950.0,196.033)
Validate_S("Al2SiO5","Kyanite",950.0,275.840)
Validate_dH("Al2SiO5","Kyanite",950.0,112092)
Validate_Hf("Al2SiO5","Kyanite",950.0,-2609095)

Validate_Hf("Al2SiO5","Kyanite",430.46,-2595232)

############################################################
Validate_Hf("Al2SiO5","Andalusite",298.15,-2590270)
Validate_Hf("Al2SiO5","Andalusite",430.46,-2591268)

Validate_Cp("Al2SiO5","Andalusite",550.0,169.117)
Validate_S("Al2SiO5","Andalusite",550.0,184.150)
Validate_dH("Al2SiO5","Andalusite",550.0,37774)
Validate_Hf("Al2SiO5","Andalusite",550.0,-2590381)

Validate_Cp("Al2SiO5","Andalusite",950.0,190.256)
Validate_S("Al2SiO5","Andalusite",950.0,283.059)
Validate_dH("Al2SiO5","Andalusite",950.0,110537)
Validate_Hf("Al2SiO5","Andalusite",950.0,-2606651)

############################################################
Validate_Hf("Al2SiO5","Sillimanite",298.15,-2587774)

Validate_Cp("Al2SiO5","Sillimanite",550.0,168.102)
Validate_S("Al2SiO5","Sillimanite",550.0,186.585)
Validate_dH("Al2SiO5","Sillimanite",550.0,37775)
Validate_Hf("Al2SiO5","Sillimanite",550.0,-2587883)

Validate_Cp("Al2SiO5","Sillimanite",1600.0,204.042)
Validate_S("Al2SiO5","Sillimanite",1600.0,389.776)
Validate_dH("Al2SiO5","Sillimanite",1600.0,240794)
Validate_Hf("Al2SiO5","Sillimanite",1600.0,-2591083)

############################################################
Validate_Hf("Ca2Al3Si3O12(OH)","Zoisite",298.15,-6891117)

Validate_Cp("Ca2Al3Si3O12(OH)","Zoisite",550.0,465.329)
Validate_S("Ca2Al3Si3O12(OH)","Zoisite",550.0,548.022)
Validate_dH("Ca2Al3Si3O12(OH)","Zoisite",550.0,105160)
Validate_Hf("Ca2Al3Si3O12(OH)","Zoisite",550.0,-6889496)

Validate_Cp("Ca2Al3Si3O12(OH)","Zoisite",950.0,532.537)
Validate_S("Ca2Al3Si3O12(OH)","Zoisite",950.0,822.107)
Validate_dH("Ca2Al3Si3O12(OH)","Zoisite",950.0,306989)
Validate_Hf("Ca2Al3Si3O12(OH)","Zoisite",950.0,-6906700)


############################################################
Validate_Hf("CaAl4Si2O10(OH)2","Margarite",298.15,-6239610)

Validate_Cp("CaAl4Si2O10(OH)2","Margarite",550.0,439.497)
Validate_S("CaAl4Si2O10(OH)2","Margarite",550.0,499.626)
Validate_dH("CaAl4Si2O10(OH)2","Margarite",550.0,98535)
Validate_Hf("CaAl4Si2O10(OH)2","Margarite",550.0,-6238608)

Validate_Cp("CaAl4Si2O10(OH)2","Margarite",950.0,505.305)
Validate_S("CaAl4Si2O10(OH)2","Margarite",950.0,759.380)
Validate_dH("CaAl4Si2O10(OH)2","Margarite",950.0,289829)
Validate_Hf("CaAl4Si2O10(OH)2","Margarite",950.0,-6264413)

############################################################
Validate_Hf("CaO","Lime",298.15,-635094)

Validate_Cp("CaO","Lime",550.0,49.707)
Validate_S("CaO","Lime",550.0,66.470)
Validate_dH("CaO","Lime",550.0,11764)
Validate_Hf("CaO","Lime",550.0,-633898)

Validate_Cp("CaO","Lime",950.0,53.550)
Validate_S("CaO","Lime",950.0,94.764)
Validate_dH("CaO","Lime",950.0,32541)
Validate_Hf("CaO","Lime",950.0,-633686)

############################################################
Validate_Hf("CaSiO3","Wollastonite",298.15,-1634766)

Validate_Cp("CaSiO3","Wollastonite",550.0,110.331)
Validate_S("CaSiO3","Wollastonite",550.0,142.113)
Validate_dH("CaSiO3","Wollastonite",550.0,25413)
Validate_Hf("CaSiO3","Wollastonite",550.0,-1633240)

Validate_Cp("CaSiO3","Wollastonite",950.0,123.322)
Validate_S("CaSiO3","Wollastonite",950.0,206.234)
Validate_dH("CaSiO3","Wollastonite",950.0,72576)
Validate_Hf("CaSiO3","Wollastonite",950.0,-1630055)

############################################################
Validate_Hf("CaSiO3","Cyclowollastonite",298.15,-1627614)

Validate_Cp("CaSiO3","Cyclowollastonite",950.0,122.030)
Validate_S("CaSiO3","Cyclowollastonite",950.0,210.641)
Validate_dH("CaSiO3","Cyclowollastonite",950.0,71503)
Validate_Hf("CaSiO3","Cyclowollastonite",950.0,-1623977)

Validate_Cp("CaSiO3","Cyclowollastonite",1550.0,131.270)
Validate_S("CaSiO3","Cyclowollastonite",1550.0,272.747)
Validate_dH("CaSiO3","Cyclowollastonite",1550.0,147847)
Validate_Hf("CaSiO3","Cyclowollastonite",1550.0,-1624302)

############################################################

Validate_Cp("Ca2SiO4","alpha",1650.0,199.600)
Validate_S("Ca2SiO4","alpha",1650.0,426.418)
Validate_Hf("Ca2SiO4","alpha",1650.0,-2278342)

Validate_Cp("Ca2SiO4","alpha",1750.0,199.600)
Validate_S("Ca2SiO4","alpha",1750.0,438.162)
Validate_Hf("Ca2SiO4","alpha",1750.0,-2324801)

############################################################

Validate_Cp("Ca2SiO4","alphaprime",1050.0,182.454)
Validate_S("Ca2SiO4","alphaprime",1050.0,329.098)
Validate_Hf("Ca2SiO4","alphaprime",1050.0,-2297850)

Validate_Cp("Ca2SiO4","alphaprime",1650.0,213.067)
Validate_S("Ca2SiO4","alphaprime",1650.0,417.454)
Validate_Hf("Ca2SiO4","alphaprime",1650.0,-2293654)

############################################################
Validate_Cp("Ca2SiO4","gamma",550.0,155.417)
Validate_S("Ca2SiO4","gamma",550.0,206.706)
Validate_dH("Ca2SiO4","gamma",550.0,35833)
Validate_Hf("Ca2SiO4","gamma",550.0,-2315156)


Validate_Cp("Ca2SiO4","gamma",1050.0,182.952)
Validate_S("Ca2SiO4","gamma",1050.0,316.884)
Validate_dH("Ca2SiO4","gamma",1050.0,121783)
Validate_Hf("Ca2SiO4","gamma",1050.0,-2311529)

############################################################
Validate_Cp("Ca2SiO4","beta",550.0,160.388)
Validate_S("Ca2SiO4","beta",550.0,215.634)
Validate_dH("Ca2SiO4","beta",550.0,36978)
Validate_Hf("Ca2SiO4","beta",550.0,-2304174)

Validate_Cp("Ca2SiO4","beta",950.0,181.741)
Validate_S("Ca2SiO4","beta",950.0,309.393)
Validate_dH("Ca2SiO4","beta",950.0,105983)
Validate_Hf("Ca2SiO4","beta",950.0,-2299712)


############################################################
Validate_Cp("SiO2","alpha",550.0,62.160)
Validate_S("SiO2","alpha",550.0,74.204)
Validate_dH("SiO2","alpha",550.0,13690)
Validate_Hf("SiO2","alpha",550.0,-910329)

Validate_Cp("SiO2","beta",950.0,68.450)
Validate_S("SiO2","beta",950.0,112.464)
Validate_Hf("SiO2","beta",950.0,-905514)
#################################################################

Validate_Cp("Ca3SiO5","Crystal",550.0,213.200)
Validate_S("Ca3SiO5","Crystal",550.0,287.139)
Validate_dH("Ca3SiO5","Crystal",550.0,49284)
Validate_Hf("Ca3SiO5","Crystal",550.0,-2928875)

Validate_Cp("Ca3SiO5","Crystal",950.0,239.686)
Validate_S("Ca3SiO5","Crystal",950.0,411.303)
Validate_dH("Ca3SiO5","Crystal",950.0,140635)
Validate_Hf("Ca3SiO5","Crystal",950.0,-2922633)
#################################################################

Validate_Cp("Ca3Si2O7","Rankinite",550.0,267.039)
Validate_S("Ca3Si2O7","Rankinite",550.0,359.181)
Validate_dH("Ca3Si2O7","Rankinite",550.0,61782)
Validate_Hf("Ca3Si2O7","Rankinite",550.0,-3969712)

Validate_Cp("Ca3Si2O7","Rankinite",950.0,293.498)
Validate_S("Ca3Si2O7","Rankinite",950.0,513.288)
Validate_dH("Ca3Si2O7","Rankinite",950.0,175045)
Validate_Hf("Ca3Si2O7","Rankinite",950.0,-3965020)

#####################################################################

Validate_Cp("H2O","Gas",550.0,35.686)
Validate_S("H2O","Gas",550.0,209.801)
Validate_dH("H2O","Gas",550.0,8692)
Validate_Hf("H2O","Gas",550.0,-244340)

Validate_Cp("H2O","Gas",950.0,40.624)
Validate_S("H2O","Gas",950.0,230.504)
Validate_dH("H2O","Gas",950.0,23934)
Validate_Hf("H2O","Gas",950.0,-247566)

#####################################################################

Validate_S("Ca3Al2O6","Crystal",500.0,(28.41+49.1)*4.1840)
Validate_dH("Ca3Al2O6","Crystal",500.0,11180*4.1840)

Validate_S("Ca3Al2O6","Crystal",1000.0,(72.42+49.1)*4.1840)
Validate_dH("Ca3Al2O6","Crystal",1000.0,43110*4.1840)

Validate_S("Ca3Al2O6","Crystal",1500.0,(99.61+49.1)*4.1840)
Validate_dH("Ca3Al2O6","Crystal",1500.0,76660*4.1840)

#####################################################################

Validate_S("Ca12Al14O33","alpha",500.0,(147.8+249.7)*4.1840)
Validate_dH("Ca12Al14O33","alpha",500.0,58250*4.1840)

Validate_S("Ca12Al14O33","alpha",1000.0,(382.9+249.7)*4.1840)
Validate_dH("Ca12Al14O33","alpha",1000.0,228750*4.1840)

Validate_S("Ca12Al14O33","beta",1500.0,(533.2+249.7)*4.1840)
print math.fabs(1.0 - ((speciesData["Ca12Al14O33"].phases["beta"].constants[-1].Hf0(1500)  - speciesData["Ca12Al14O33"].phases["alpha"].constants[-1].Hf0(298.15))/(414400*4.1840))) < 0.01

#####################################################################

Validate_S("CaAl2O4","Crystal",500.0,(16.76+27.3)*4.1840)
Validate_dH("CaAl2O4","Crystal",500.0,6610*4.1840)

Validate_S("CaAl2O4","Crystal",1000.0,(43.59+27.3)*4.1840)
Validate_dH("CaAl2O4","Crystal",1000.0,26100*4.1840)

Validate_S("CaAl2O4","Crystal",1500.0,(60.99+27.3)*4.1840)
Validate_dH("CaAl2O4","Crystal",1500.0,47610*4.1840)

#####################################################################

Validate_S("MgAl2O4","Crystal",500.0,(16.85+19.26)*4.1840)
Validate_dH("MgAl2O4","Crystal",500.0,6650*4.1840)

Validate_S("MgAl2O4","Crystal",1000.0,(43.98+19.26)*4.1840)
Validate_dH("MgAl2O4","Crystal",1000.0,26390*4.1840)

Validate_S("MgAl2O4","Crystal",1500.0,(61.95+19.26)*4.1840)
Validate_dH("MgAl2O4","Crystal",1500.0,48620*4.1840)

#####################################################################

Validate_S("CaFe2O4","Crystal",500.0,(20.18+34.7)*4.1840)
Validate_dH("CaFe2O4","Crystal",500.0,7900*4.1840)

Validate_S("CaFe2O4","Crystal",1000.0,(49.22+34.7)*4.1840)
Validate_dH("CaFe2O4","Crystal",1000.0,28950*4.1840)

Validate_S("CaFe2O4","Liquid",1600.0,(88.11+34.7)*4.1840)
print math.fabs(1.0 - ((speciesData["CaFe2O4"].phases["Liquid"].constants[-1].Hf0(1600)  - speciesData["CaFe2O4"].phases["Crystal"].constants[-1].Hf0(298.15))/(82810*4.1840))) < 0.01

#####################################################################

Validate_S("Ca2Fe2O5","Crystal",500.0,(25.92+45.1)*4.1840)
Validate_dH("Ca2Fe2O5","Crystal",500.0,10210*4.1840)

Validate_S("Ca2Fe2O5","Crystal",1000.0,(65.55+45.1)*4.1840)
Validate_dH("Ca2Fe2O5","Crystal",1000.0,38890*4.1840)

Validate_S("Ca2Fe2O5","Liquid",1800.0,(121.00+45.1)*4.1840)
print math.fabs(1.0 - ((speciesData["Ca2Fe2O5"].phases["Liquid"].constants[-1].Hf0(1800)  - speciesData["Ca2Fe2O5"].phases["Crystal"].constants[-1].Hf0(298.15))/(122580*4.1840))) < 0.01

#####################################################################

Validate_S("MgFe2O4","alpha",500.0,(20.03+28.3)*4.1840)
Validate_dH("MgFe2O4","alpha",500.0,7870*4.1840)

Validate_S("MgFe2O4","beta",1000.0,(51.83+28.3)*4.1840)
print math.fabs(1.0 - ((speciesData["MgFe2O4"].phases["beta"].constants[-1].Hf0(1000)  - speciesData["MgFe2O4"].phases["alpha"].constants[-1].Hf0(298.15))/(30810*4.1840))) < 0.01

Validate_S("MgFe2O4","gamma",1500.0,(70.27+28.3)*4.1840)
print math.fabs(1.0 - ((speciesData["MgFe2O4"].phases["gamma"].constants[-1].Hf0(1500)  - speciesData["MgFe2O4"].phases["alpha"].constants[-1].Hf0(298.15))/(53520*4.1840))) < 0.01


#####################################################################

Validate_S("Al2TiO5","Crystal",500.0,(19.3+26.2)*4.1840)
Validate_dH("Al2TiO5","Crystal",500.0,7620*4.1840)

Validate_S("Al2TiO5","Crystal",1000.0,(50.73+26.2)*4.1840)
Validate_dH("Al2TiO5","Crystal",1000.0,30450*4.1840)

Validate_S("Al2TiO5","Crystal",1500.0,(70.80+26.2)*4.1840)
Validate_dH("Al2TiO5","Crystal",1500.0,55260*4.1840)

#####################################################################

Validate_S("Fe2TiO4","Crystal",500.0,(19.40+169.034)*4.1840)
Validate_dH("Fe2TiO4","Crystal",500.0,7610*4.1840)

Validate_S("Fe2TiO4","Crystal",1000.0,(49.29+169.034)*4.1840)
Validate_dH("Fe2TiO4","Crystal",1000.0,29400*4.1840)

Validate_S("Fe2TiO4","Crystal",1500.0,(70.27+169.034)*4.1840)
Validate_dH("Fe2TiO4","Crystal",1500.0,55450*4.1840)
#####################################################################

Validate_S("Zn2TiO4","Crystal",500.0,(18.55+32.8)*4.1840)
Validate_dH("Zn2TiO4","Crystal",500.0,7290*4.1840)

Validate_S("Zn2TiO4","Crystal",1000.0,(47.89+32.8)*4.1840)
Validate_dH("Zn2TiO4","Crystal",1000.0,28660*4.1840)

Validate_S("Zn2TiO4","Crystal",1500.0,(66.76+32.8)*4.1840)
Validate_dH("Zn2TiO4","Crystal",1500.0,51950*4.1840)
