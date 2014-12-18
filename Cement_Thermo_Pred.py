#!/usr/bin/env python
from chemeng.standarddefinitions import DryAir
from chemeng import *
import chemeng.cementdata
from C4AF_Thermodata import *
from Yeelemite_and_Ternesite_Data import*
import pickle
from operator import add
from matplotlib import pyplot as plt
from matplotlib import interactive
import colorsys
import sys
from matplotlib.widgets import Button
from datetime import datetime

# Auto refine plot , True for ON or False for off
AutoRefinement = True
# Default SO2 Partial pressure of 10-6 atm ,True for ON or False for off
Default_SO2_pp = False

T_start = 1050
T_Finish = 1300
T_step = 30

###Factor to multiply by to keep partial pressure constant
Factor = 100000

Air_mass_flow_rate = 2.624* Factor 
SO2_mass_flow_rate = 0.105 * Factor

####Clinker 1
mass_Al2O3 = 32.172
mass_CaO =  46.205
mass_Fe2O3 =3.286
mass_SiO2 = 10.465
mass_MgO = 0.0
mass_K2O =0.0
mass_Na2O = 0.0
mass_TiO2 = 0.0
mass_CaSO4 = 0.0

#Clinker 2
#mass_Al2O3 = 14.629 #g
#mass_CaO =  43.107  #g
#mass_Fe2O3 = 3.286  #g
#mass_SiO2 = 20.698  #g

## INPUT AREA ENDED 
##########################################################################################################################################################################################################################################################################################################################################################################################################################################
print "######################################################"
print 
print
print "CEMENT THERMO PRED SOFTWARE STARTED: ", str(datetime.now())
print

def RawInput(ref,Statement):
    val = raw_input(Statement)
    if val == "":
        val = ref
    val = float(val)
    return val

maxinteger = float(sys.maxint)
plotstyles=['xr-', 'xg-', 'xb-', 'xc-', 'xm-', 'xy-', 'xk-', 'or-', 'og-', 'ob-', 'oc-', 'om-', 'oy-', 'ok-', '+r-', '+g-', '+b-', '+c-', '+m-', '+y-', '+k-']

unwantedlist = ["Gas"]

### Dryair Components normalised to 1 mole air
DryAir = DryAir.normalised()

Air_molar_flow_rate = Air_mass_flow_rate/ DryAir.totalMass()
SO2_molar_flow_rate = SO2_mass_flow_rate/ Components({'SO2':1.0}).totalMass()



if Default_SO2_pp and SO2_molar_flow_rate/(SO2_molar_flow_rate + Air_molar_flow_rate) < 1e-6:
    print 
    print "!!! Default SO2 partial pressure activated !!!"
    print 
    print
    
    SO2_molar_flow_rate = (Air_molar_flow_rate*0.000001)/0.999999

    

moles_Al2O3 = mass_Al2O3 / Components({'Al2O3':1.0}).totalMass()  
moles_CaO = mass_CaO / Components({'CaO':1.0}).totalMass()
moles_Fe2O3 = mass_Fe2O3 / Components({'Fe2O3':1.0}).totalMass() 
moles_SiO2 = mass_SiO2 / Components({'SiO2':1.0}).totalMass() 
moles_MgO = mass_MgO / Components({'MgO':1.0}).totalMass()
moles_K2O = mass_K2O / Components({'K2O':1.0}).totalMass() 
moles_Na2O = mass_Na2O / Components({'Na2O':1.0}).totalMass() 
moles_TiO2 = mass_TiO2/ Components({'TiO2':1.0}).totalMass() 
moles_CaSO4 = mass_CaSO4/Components({'CaSO4':1.0}).totalMass()

###Gaseous phases considered , Missing F and Cl containg Gases
phases=[IdealGasPhase(DryAir*Air_molar_flow_rate, T=298.15, P=1e5) + IdealGasPhase({'CO2':0.0}, T=298.15, P=1e5) +IdealGasPhase({'H2O':0.0}, T=298.15, P=1e5) + IdealGasPhase({'CO':0.0}, T=298.15, P=1e5) + IdealGasPhase({'SO3':0.0}, T=298.15, P=1e5) + IdealGasPhase({'SO2':SO2_molar_flow_rate }, T=298.15, P=1e5) +  IdealGasPhase({'OH':0.0}, T=298.15, P=1e5) +IdealGasPhase({'NaOH':0.0}, T=298.15, P=1e5) +IdealGasPhase({'Na':0.0}, T=298.15, P=1e5)+IdealGasPhase({'KOH':0.0}, T=298.15, P=1e5)+IdealGasPhase({'K':0.0}, T=298.15, P=1e5) + IdealGasPhase({'Na2SO4':0.0}, T=298.15, P=1e5)   + IdealGasPhase({'K2SO4':0.0}, T=298.15, P=1e5)  ]

names=["Gas"]

#### Possible Incompressible phases that can be formed WARNING KEEP OR REMOVE LIQUID PHASES
Incompressiblephases = [
    ('H2O', "Liquid", 0),
    ('NaAlSi3O8', "Liquid", 0),
    ('NaAlSi3O8', "Albite", 0),
    ('NaAlSi3O8', "highAlbite", 0),
    ('NaFeSi2O6', "Acmite", 0),
    ('NaAlSi2O6', "Jadeite", 0),
    ('CaAl2Si2O8', "Liquid", 0),
    ('CaMgSi2O6', "Liquid", 0),
    ('CaMgSi2O6', "Diopside", 0),
    ('CaMgC2O6', "Dolomite", 0),
    ('Mg2Si2O6', "Liquid", 0),
    ('Mg2Si2O6', "Enstatite", 0),
    ('Fe2SiO4', "Liquid", 0),
    ('Fe2SiO4', "Fayalite", 0),
    ('Fe2Si2O6', "Ferrosilite", 0),
    ('K3Fe0.5Al4Si19.5O47', "Liquid", 0),
    ('K3Mg0.5Al4Si19.5O47', "Liquid", 0),
    ('KAlSi3O8', "Liquid", 0),
    ('KAlSi3O8', "Sanidine", 0),
    ('Al2SiO5', "Liquid", 0),
    ('Fe2O3', "1cr", moles_Fe2O3),
    ('SiO2', "alpha",moles_SiO2), 
    ('SiO2', "beta", moles_SiO2),
    ('Fe2O3', "2cr", moles_Fe2O3),
    ('MgCO3', "1cr",0.0 ),
    ('KAlO2', "1II",0.0 ),
    ('KAlO2', "2I",0.0 ),
    ('Na2CO3', "1a",0.0 ),
    ('Na2CO3', "2b",0.0 ),
    ('Na2CO3', "3c",0.0 ),
    ('Na2CO3', "Liquid",0.0 ),
    ('NaAlO2', "1a",0.0 ),
    ('NaAlO2', "2b",0.0 ),
    ('K2O', "1c",moles_K2O ),
    ('K2O', "2b",moles_K2O ),
    ('K2O', "3a",moles_K2O ),
    ('K2O', "Liquid",moles_K2O ),
    ('Na2O', "1c",moles_Na2O),
    ('Na2O', "2b",moles_Na2O ),
    ('Na2O', "3a",moles_Na2O ),
    ('Na2O', "Liquid",moles_Na2O ),
    ('Na2SO4', "1V",0.0 ),
    ('Na2SO4', "2IV",0.0 ),
    ('Na2SO4', "3I",0.0 ),
    ('Na2SO4', "Liquid",0.0 ),
    ('K2CO3', "1a",0.0),
    ('K2CO3', "2b",0.0 ),
    ('K2CO3', "Liquid",0.0 ),
    ('K2SO4', "1II",0.0 ),
    ('K2SO4', "2I",0.0 ),
    ('K2SO4', "Liquid",0.0 ),
    ('K2Si2O5', "1a",0 ),
    ('K2Si2O5', "2b",0 ),
    ('K2Si2O5', "3c",0 ),
    ('K2Si2O5', "Liquid",0 ),
    ('K2SiO3', "1cr",0 ),
    ('K2SiO3', "Liquid",0 ),  
    ('MgCO3', "Liquid",0.0 ),
    ('MgSiO3', "1I",0.0 ),
    ('MgSiO3', "2II",0.0 ),
    ('MgSiO3', "3III",0.0 ),
    ('MgSiO3', "Liquid",0.0 ),
    ('MgFe2O4', "alpha",0.0 ),
    ('MgFe2O4', "beta",0.0 ),
    ('MgFe2O4', "gamma",0.0 ),
    ('MgAl2O4', "1cr",0.0 ),
    ('Mg2SiO4', "1cr",0.0 ),
    ('MgO', "1cr",moles_MgO ),
    ('CaAl2SiO6', "Ca-Al Clinopyroxene", 0),
    ('CaFe2O4', "Liquid", 0),
    ('CaFe2O4', "Crystal", 0),
    ('Ca4Fe2Al2O10', "Crystal", 0),
    ('Al2O3', "Corundum", moles_Al2O3),
    ('Ca3Al2Si3O12', "Grossular", 0),
    ('CaSiO3', "Cyclowollastonite", 0),
    ('CaSiO3', "Wollastonite", 0), 
    ('Al2SiO5', "Sillimanite", 0),
    ('Al2SiO5', "Andalusite", 0),
    ('Al2SiO5', "Kyanite", 0),
    ('CaAl2Si2O8', "Anorthite", 0),
    ('CaAl4Si2O10(OH)2', "Margarite", 0),
    ('Ca2Al3Si3O12(OH)', "Zoisite", 0),
    ('Ca2Al2Si3O10(OH)2', "Prehnite", 0),
    ('Ca2SiO4', "beta", 0),              
    ('Al2Si4O10(OH)2', "Pyrophyllite", 0),
    ('HAlO2', "Boehmite", 0),
    ('HAlO2', "Diaspore", 0),
    ('Ca2Al2SiO7', "Gehlenite", 0),
    ('Ca3Si2O7', "Rankinite", 0), 
    #('Al2TiO5', "Crystal", 0), From Barin Before
    #('Fe2TiO4', "Crystal", 0),From Barin Before
    ('MgTiO3', "1cr", 0),
    ('Mg2TiO4', "1cr", 0),
    ('MgTi2O5', "1cr", 0),
    ('Ti3O5', "1a", 0),
    ('Ti3O5', "2b", 0),
    ('TiO2', "1cr", moles_TiO2),
    ('TiO', "1a", 0),
    ('TiO', "2b", 0),
    ('TiO', "3c", 0),
    ('Ti2O3', "1I", 0),
    ('Ti2O3', "2I'", 0),
    ('Ti4O7', "1cr", 0),
    #('Zn2TiO4', "Crystal", 0),
    ('CaTiO3', "alpha", 0),
    ('CaTiO3', "beta", 0),
    ('FeTiO3', "Crystal", 0),
    ('FeTiO3', "Liquid", 0),
    ('Ca2Fe2O5', "Liquid", 0), 
    ('Ca2Fe2O5', "Crystal", 0),
    ('Al6Si2O13', "Mullite", 0),
    ('Ca5Si2C2O13', "Tilleyite", 0),   
    ('Ca12Al14O33', "alpha", 0),
    ('Ca12Al14O33', "beta", 0),
    ('CaAl2O4', "Crystal", 0),
    ('CaAl4O7', "Crystal", 0),  
    ('Ca3Al2O6', "Crystal", 0),
    ('Ca5Si2CO11', "Spurrite", 0),
    ('Ca2SiO4', "alphaprime", 0.0),
    ('Ca2SiO4', "gamma", 0),
    ('Ca2SiO4', "alpha", 0),
    ('Ca3SiO5', "Crystal", 0),
    ('Ca5Si2SO12', "Ternesite_VP", 0),
    ('Ca4Al6SO16', "Yeelemite_VP", 0),
    ('CaO', "Lime", moles_CaO),
    ('CaCO3', "1cr",0.0 ),
    ('CaCO3', "Liquid",0.0 ), 
    ('CaSO4', "1II", moles_CaSO4),
    ('CaSO4', "3I", moles_CaSO4),
    ('CaSO4', "Liquid", moles_CaSO4),
    ]


def checkEqual2(iterator):
   return len(set(iterator)) <= 1

### Get Elemental Composition of Everything in System
Elemental_Composition =  phases[0].components.elementalComposition() + Components({'Al2O3':moles_Al2O3, 'CaO':moles_CaO,'Fe2O3':moles_Fe2O3,'SiO2':moles_SiO2,'MgO':moles_MgO,'K2O':moles_K2O,'Na2O':moles_Na2O,'TiO2':moles_TiO2}).elementalComposition()
Unavailable_Elements = []
### Find Elements not available in System
for comp, amnt in Elemental_Composition.iteritems():
    if amnt == 0:
        Unavailable_Elements.append(comp)
indexes = []
### Get indexes of Species that cannot form in System as Elements not available
for species,phase,amnt in Incompressiblephases:
    for element,amount in Components({species:amnt}).elementalComposition().iteritems():
        if element in Unavailable_Elements:
	    if Incompressiblephases.index((species,phase,amnt)) not in indexes:
	        indexes.append(Incompressiblephases.index((species,phase,amnt)))
### Ammend Incompressible Phase List	    
for index in sorted(indexes, reverse=True):
    del Incompressiblephases[index]

### Add Incompressible Phases to Phase List	    
for species,phase,amount in Incompressiblephases:
    phases.append(IncompressiblePhase({species:amount}, T=298.15, P=1e5, phaseID=phase))
    names.append(species+':'+phase)

### Make Dictionaries to store Data
Ts_gas = {}
pp_gas = {}
Ts = {}
Moles = {}
Mass = {}
Mass_percent = {}
Total_Solid_Mass = []
T_solved=[]


for name in names:
    Ts[name] = []
    Moles[name] = []
    Mass[name] = []
    Mass_percent[name] = []
    
    
for comp,amnt in phases[0].components.iteritems():
    Ts_gas[comp] = []
    pp_gas[comp] = []
    
def updatePlot():
    legendnames = []
    legendProxies = []
    i=0
    for name in names:
	if name not in unwantedlist and sum(Moles[name])>0.01:
	    legendnames.append(name)
	    legendProxies.append(plt.Rectangle((0, 0), 1, 1, fc=colorsys.hls_to_rgb(abs(hash(name))*360/maxinteger, 0.5, abs(hash(name))/maxinteger)))
	    i+=1

    ax1.clear();
    ax1.legend(legendProxies, legendnames,loc='center left', bbox_to_anchor=(1, 0.5),fancybox=True, shadow=True)
    ax1.set_xlim([T_start,T_Finish])
    ax1.set_xlabel("Temperature $^\circ$C", fontsize = 20)
    ax1.set_ylabel("Proportions by weight", fontsize = 20)

    Saved_mass = [0]*len(Mass[names[0]])   	    
    for name in names: 
	if name not in unwantedlist and sum(Moles[name])>0.001:
	    try:
		ax1.fill_between(Ts[name], Saved_mass, [sum(x) for x in zip(Saved_mass, Mass[name])], facecolor = colorsys.hls_to_rgb(abs(hash(name))*360/maxinteger, 0.5, abs(hash(name))/maxinteger) ,label = name, picker=True) 
		Saved_mass = [sum(x) for x in zip(Saved_mass, Mass[name])]
		
	    except:
		print "ERRROROROR"
    for T in T_solved:
	ax1.plot([T-273.15,T-273.15],[0,Total_Solid_Mass[0]],'--k')
    fig.canvas.mpl_connect('pick_event', pick_handler)	
    plt.draw()
    plt.pause(0.0001)

def solveAtT(T):
    T_solved.append(T)
    currentphases=[]
    currentnames=[]
    for phase,name in zip(phases,names):
        if speciesData[phase.components.keys()[0]].inDataRange(T, phase.phase) == False:
            continue
        currentphases.append(phase)
        currentnames.append(name)
    currentphases[0].T = T
    try:
        result = findEquilibrium(currentphases, constP=True, constT=True, elemental = True,iterations=300, xtol=1e-6) 
    except Exception as e:
        print "Did not converge for T=",T,'\n', e
        return
    try:
        result = findEquilibrium(result, constP=True, constT=True, elemental = True,iterations=300, xtol=1e-8) 
    except Exception as e: 
        print "Did not converge for T=",T,'\n', e
        return
    for comp,amnt in result[0].components.iteritems():
        if amnt/result[0].components.total() > 0.0000499:
            Ts_gas[comp].append(T-273.15)
            pp_gas[comp].append(math.log10(amnt/result[0].components.total()))	

    for phase,name,datum in zip(currentphases, currentnames, result):
        Ts[name].append(T - 273.15)
        Moles[name].append(datum.components.total())
        Mass[name].append(datum.components.totalMass())
        Mass_percent[name].append(100*datum.components.totalMass()/(sum (result[i].components.totalMass() for i in range(2,len(result)))))
        
    Total_Solid_Mass.append(sum (result[i].components.totalMass() for i in range(1,len(result))))# Skip gas phase # WARNING Incorrect, interpretation, need if statement    

    for name in names: 
        for vals in Ts["Gas"]:
            if vals not in Ts[name]:
                Ts[name].append(vals)
                Moles[name].append(0.0)
                Mass[name].append(0.0)
                Mass_percent[name].append(0.0)
                    
    for name in names:
        sorted_lists = sorted(zip(Mass[name], Ts[name],Moles[name],Mass_percent[name]), reverse=False, key=lambda x: x[1])
        Mass[name],Ts[name],Moles[name],Mass_percent[name] = [[x[i] for x in sorted_lists] for i in range(4)]            
    updatePlot()
    

plt.ion()
fig = plt.figure()
ax1 = fig.add_subplot(111)
box = ax1.get_position()
ax1.set_position([box.x0, box.y0, box.width * 0.8, box.height])

    
axRef = plt.axes([0.81, 0.05, 0.1, 0.075])
bRef = Button(axRef, 'Show PP')
    
def on_button_clicked(event):
    
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    gasplotstyles = ["^","v",">","<","D","x","8","s","p","*","+","H","o","^","v",">","<","D","x","8","s","p","*","+","H","o"]
    sasplotline = ["-","-","-","-","-","-","-","-","-","-","-","-","-","--","--","--","--","--","--","--","--","--","--","--","--","--"] 

    for name in pp_gas: 
	sorted_lists2 = sorted(zip(pp_gas[name], Ts_gas[name]), reverse=False, key=lambda x: x[1])
	pp_gas[name], Ts_gas[name] = [[x[i] for x in sorted_lists2] for i in range(2)]     
    ii = 0
    for name in pp_gas:
	if len(pp_gas[name]) > 2:
	    plt.plot(Ts_gas[name], pp_gas[name], linestyle=sasplotline[ii], marker=gasplotstyles[ii], color='k', markevery = 5,label = name)
	    ii+=1
    plt.tick_params(axis='x', labelsize=16)    
    plt.tick_params(axis='y', labelsize=16) 
    plt.title("Partial Pressures",fontsize = 28)

    plt.yticks([-4.0,-3.0,-2.0,-1.0,0.0])

    plt.xlabel("Temperature $^\circ$C", fontsize = 18)
    plt.ylabel("$log_{10}$(P/atm) ", fontsize = 18)
    box = ax1.get_position()
    ax1.set_position([box.x0, box.y0, box.width * 0.83, box.height])    

    ax1.legend(loc='center left', bbox_to_anchor=(1, 0.5),fancybox=True, shadow=True)
    plt.show()

bRef.on_clicked(on_button_clicked)

#(xm,ym),(xM,yM)=bnext.label.clipbox.get_points()

#def on_press(event):
    #if xm<event.x<xM and ym<event.y<yM:
	#print "Button clicked, do nothing. This triggered event is useless."
#fig.canvas.mpl_connect('button_press_event', on_press)


txt = fig.suptitle("Phase Diagram",fontsize = 25)
Handled_data = []    
def pick_handler(event):
    if event.mouseevent.xdata not in Handled_data:
	txt.set_text(event.artist.get_label())
	T = (event.mouseevent.xdata)
	print "\nTemperature =",T , "C"
	solveAtT(event.mouseevent.xdata+273.15)
	for name in names:
	    indexT = min(range(len(Ts[name])), key=lambda i: abs(Ts[name][i]-T))  
	    if Mass_percent[name][indexT]> 0.01 and name != "Gas":
		print name,Mass_percent[name][indexT], "Percent"
	print
	Handled_data.append(event.mouseevent.xdata)

plt.pause(0.0001)
#Generate temperatures
T_vals = [T_start,T_Finish]
deltaT = T_Finish-T_start
denom = 2
while ((deltaT / denom) > T_step):
    for i in range(1,denom,2):
	T_vals.append(T_Finish - i * deltaT / denom)
    denom = denom * 2

for T in map(lambda T : T+273.15, T_vals):
    solveAtT(T)


if AutoRefinement:
    i=0
    while i< len( Ts["Gas"])-1:
	names1 = []
	names2 = []
	for name in names:
	    if Moles[name][i]>0.001:
		names1.append(name)
	    if Moles[name][i+1]>0.001:
		names2.append(name)
	if names1 != names2:
	    solveAtT((Ts["Gas"][i] + Ts["Gas"][i+1])/2.0 + 273.15)   
	i+=1    


if AutoRefinement:
    i=0
    while i< len( Ts["Gas"])-1:
	names1 = []
	names2 = []
	for name in names:
	    if Moles[name][i]>0.001:
		names1.append(name)
	    if Moles[name][i+1]>0.001:
		names2.append(name)
	if names1 != names2:
	    solveAtT((Ts["Gas"][i] + Ts["Gas"][i+1])/2.0 + 273.15)   
	i+=1  
  
data1 = Ts,Moles,Mass,Total_Solid_Mass,Ts_gas,pp_gas,Mass_percent
output = open(("CEMENTSOFTWARE") +'.pkl', 'wb')
pickle.dump(data1, output)
output.close()   

  
pkl_file = open("CEMENTSOFTWARE" +'.pkl', 'rb')
data1 = pickle.load(pkl_file)
pkl_file.close()
Ts,Moles,Mass,Total_Solid_Mass,Ts_gas,pp_gas,Mass_percent = data1    

plt.ioff()
plt.show()

