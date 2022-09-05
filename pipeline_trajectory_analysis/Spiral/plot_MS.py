"Plot MSD and plain displacement."

import matplotlib.pyplot as plt
import glob
import numpy as np
import pandas as pd

fmonoS59o = '/home_d/malin/Sim/1osl/3_results_mono_slide/Model12_Lys59off/'
fmonoS = '/home_d/malin/Sim/1osl/3_results_mono_slide/final/'
fmonoR = '/home_d/malin/Sim/1efa/3_results_mono_slide/final/'
ffullS = '/home_d/malin/Sim/full-length/2-MD/S-model/'
ffullR = '/home_d/malin/Sim/full-length/2-MD/R-model/'
fNS = '/home_d/malin/Sim/N-core/2-MD/S-model/'
fNR = '/home_d/malin/Sim/N-core/2-MD/R-model/'
transitionmodel = '/trajectories/malin/Sim/TransitionModel/onBDNA/'

fR = '/home_d/malin/Sim/dimer_model/2-MD/R-helix/AprilkBT04/6contacts/04kBT/'
fR59off = '/home_d/malin/Sim/dimer_model/2-MD/R-helix/AprilkBT04/6contacts/Lys59off_6contacts/'
#fRN = '/home_d/malin/Sim/dimer_model/2-MD/R_Ncore/salt04/'
#fS = '/home_d/malin/Sim/dimer_model/2-MD/Model12/kBT04/'
#fS = '/home_d/malin/Sim/dimer_model/2-MD/Model12/symmetric/sym/'
fS = '/home_d/malin/Sim/dimer_model/2-MD/Model12/final/'
fS59off = '/home_d/malin/Sim/dimer_model/2-MD/Model12/nochargeLys59unimore/'

#c = ['01','02','03','04']
c = ['01','02','03','04','05','06']

StoR = [transitionmodel + x + '/output/Traj/' for x in c]

fullS = [ffullS + x + '/output/Analisis/' for x in c]
fullR = [ffullR + x + '/output/Analisis/' for x in c]

NS = [fNS + x + '/output/Analisis/' for x in c]
NR = [fNR + x + '/output/Analisis/' for x in c]

monoS = [fmonoS + x + '/output/Analisis/' for x in c]
monoR = [fmonoR + x + '/output/Analisis/' for x in c]
monoS59o = [fmonoS59o + x + '/output/Traj/' for x in c]

R = [fR + x + '/output/Analisis/' for x in c]
R59o = [fR59off + x + '/output/Analisis/' for x in c]
#RN = [fRN + x + '/output/Traj/' for x in c]
S = [fS + x + '/output/Analisis/' for x in c]  
S59o = [fS59off + x + '/output/Analisis/' for x in c]  

def plot_MS(folder):
    M = [] 
    S = [] 
    for Conc in folder:
        "Plot Displacement for each trajectory."
        filesA = glob.glob(Conc+'*Displacement_recA.txt')
        filesB = glob.glob(Conc+'*Displacement_recB.txt')
        files = filesA + filesB
        CZ = []
        CXY = []
        ang = []
        for file in files:
            df = pd.read_csv(file, sep='\s+',header=None)
            CZ.append(df[0])
            CXY.append(df[1])
            ang.append(df[2])
        CXYnp = np.array([x for y in CXY for x in y])
        Mean = CXYnp.mean()
        Std = CXYnp.std()
        M.append(Mean)
        #This is the standard error of the mean
        S.append(Std/(np.sqrt(len(M))))
    return(M,S)


#monorS = plot_MS(monoS)
#monorR = plot_MS(monoR)
#monor59o = plot_MS(monoS59o)

fullrS = plot_MS(fullS)
fullrR = plot_MS(fullR)

NrS = plot_MS(NS)
NrR = plot_MS(NR)

#Rr = plot_MS(R)
#Rr59o = plot_MS(R59o)
#RrN = plot_MS(RN)
#Sr = plot_MS(S)
#Sr59o = plot_MS(S59o)

StoRr = plot_MS(StoR)

#ci = [0.01,0.02,0.03,0.04] 
ci = [0.01,0.02,0.03,0.04,0.05,0.06] 

fig, ax = plt.subplots(figsize = (4,2.5))

#plt.errorbar(ci,monorS[0],yerr=monorS[1], label = 'Monomer, S-state')
#plt.errorbar(ci,monorR[0],yerr=monorR[1], label = 'Monomer, R-state')
#plt.errorbar(ci,monor59o[0],yerr=monor59o[1], label = 'Monomer_K59o, S-state')
#plt.errorbar(ci,Rr[0],yerr=Rr[1], label = 'R-state')
#plt.errorbar(ci,Rr59o[0],yerr=Rr59o[1], label = 'R-state_K59-')
#plt.errorbar(ci,RrN[0],yerr=RrN[1], label = 'R-state_288')
#plt.errorbar(ci,Sr[0],yerr=Sr[1], label = 'S-state')
#plt.errorbar(ci,Sr59o[0],yerr=Sr59o[1], label = 'S-state_K59-')

#plt.errorbar(ci,monorS[0],yerr=monorS[1], linestyle='-.',label = 'Monomer, S-state',linewidth=0.5,capsize=3, markeredgewidth=0.5, color='k')
#plt.errorbar(ci,monorR[0],yerr=monorR[1], linestyle='-.',label = 'Monomer, R-state',linewidth=0.5,capsize=3, markeredgewidth=0.5, color = 'grey')
#plt.errorbar(ci,monor59o[0],yerr=monor59o[1],linestyle='-', label = 'Monomer_K59o, S-state',linewidth=0.5,capsize=3, markeredgewidth=0.5)
#plt.errorbar(ci,Rr59o[0],yerr=Rr59o[1],linestyle='-',linewidth=0.5,label = 'R-state_K59-',capsize=3, markeredgewidth=0.5)
#plt.errorbar(ci,RrN[0],yerr=RrN[1], label = 'R-state_288',,capsize=3, markeredgewidth=0.5,linewidth=0.5)
#plt.errorbar(ci,Sr[0],yerr=Sr[1],linestyle='-',linewidth=0.5, label = 'S-state',capsize=3,markeredgewidth=0.5, color = 'b')
#plt.errorbar(ci,Rr[0],yerr=Rr[1], linestyle='-',linewidth=0.5,label = 'R-state',capsize=3, markeredgewidth=0.5, color ='g')
#plt.errorbar(ci,Sr59o[0],yerr=Sr59o[1],linestyle='-',linewidth=0.5, label = 'S-state_K59-',capsize=3,markeredgewidth=0.5)

plt.errorbar(ci,fullrS[0],yerr=fullrS[1],linestyle='-',linewidth=0.5, label = 'S-state',capsize=3,markeredgewidth=0.5, color = 'b')
plt.errorbar(ci,fullrR[0],yerr=fullrS[1], linestyle='-',linewidth=0.5,label = 'R-state',capsize=3, markeredgewidth=0.5, color ='g')

#plt.errorbar(ci,NrS[0],yerr=NrS[1],linestyle='-',linewidth=0.5, label = 'S-state',capsize=3,markeredgewidth=0.5, color = 'b')
#plt.errorbar(ci,NrR[0],yerr=NrS[1], linestyle='-',linewidth=0.5,label = 'R-state',capsize=3, markeredgewidth=0.5, color ='g')

plt.errorbar(ci,StoRr[0],yerr=StoRr[1], linestyle='-',linewidth=0.5,label = 'Mixed model',capsize=3, markeredgewidth=0.5, color ='r')

axes = plt.gca()
axes.set_xlim(0,0.08)
axes.set_ylim(0,200)
plt.xlabel('Ionic strength (M)')
plt.ylabel('Distance (A) from DNA axis')

ax.tick_params(
axis='both',          # changes apply to the x-axis
which='both',      # both major and minor ticks are affected
bottom=True,      # ticks along the bottom edge are off
top=True,         # ticks along the top edge are off
labelbottom=True,
labeltop = False,
right=True, left=True, labelleft=True,labelsize='small',
width = 0.5,
length = 3, color='k',direction ='in')

plt.legend()
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.savefig('Cdep.png',dpi=400,bbox_inches='tight')

#fig, ax = plt.subplots(figsize = (4,2.5))
#plt.errorbar(ci[0:3],RMN[0:3],yerr=RSN[0:3], label = 'R model 288')
#plt.errorbar(ci[0:3],SM[0:3],yerr=SS[0:3], label = 'S model 118')
#plt.errorbar(ci[0:3],RM[0:3],yerr=RS[0:3], label = 'R model 118')
#axes = plt.gca()
#plt.xlabel('Ionic strength (M)')
#plt.ylabel('Distance (A) from DNA axis')

#ax.tick_params(
#axis='both',          # changes apply to the x-axis
#which='both',      # both major and minor ticks are affected
#bottom=True,      # ticks along the bottom edge are off
#top=True,         # ticks along the top edge are off
#labelbottom=True,
#labeltop = False,
#right=True, left=True, labelleft=True,labelsize='small',
#width = 0.5,
#length = 3, color='k',direction ='in')

#plt.legend()
#plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
#plt.savefig('CdepNear.png',dpi=400,bbox_inches='tight')

