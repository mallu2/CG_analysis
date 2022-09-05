import numpy as np
import pandas as pd
import glob
import matplotlib.pyplot as plt


fR = '/home_d/malin/Sim/dimer_model/2-MD/R-helix/AprilkBT04/'
fRN = '/home_d/malin/Sim/dimer_model/2-MD/R_Ncore/salt04/'
fS = '/home_d/malin/Sim/dimer_model/2-MD/Model12/'

c = ['01','02','03','04']

R = [fR + x + '/output/EbyType/' for x in c]
RN = [fRN + x + '/output/EbyType/' for x in c]
S = [fS + x + '/output/EbyType/' for x in c]  


def read_energy(folder):
    files = glob.glob(folder+'EbyType*.dat')
    Energies = []
    bonds= []
    angles= []
    dih = []
    lj = []
    rep = []
    el= []
    for file in files:
        df = pd.read_csv(file, sep='\s+',header=None,engine='python')
        bonds.append(df[0])
        angles.append(df[1])
        dih.append(df[2])
        lj.append(df[3])
        rep.append(df[4])
        el.append(df[6])
        Energies.append(df)
    #combine all 5 replica
    Bonds = np.array([x for y in bonds for x in y])
    BM, BS = Bonds.mean(), Bonds.std()
    Angles = np.array([x for y in angles for x in y])
    AM, AS = Angles.mean(), Angles.std()
    Dih = np.array([x for y in dih for x in y])
    DM, DS = Dih.mean(), Dih.std()
    LJ = np.array([x for y in lj for x in y])
    LJM, LJS = LJ.mean(), LJ.std()
    Rep = np.array([x for y in rep for x in y])
    RM, RS = Rep.mean(), Rep.std()
    El = np.array([x for y in el for x in y])
    EM, ES = El.mean(), El.std()

    #return([(BM,BS),()
    return(Energies)

def plot_energy_vs_conc(folderR,folderS,type):
    dataR = read_energy(folderR)
    dataS = read_energy(folderS)
    
    fig, ax = plt.subplots(figsize = (4,2.5))
    ci = [0.01,0.02,0.03,0.04] 
    plt.errorbar(ci[2],yerr=dataS[1],label= 'electrostatics S')
    plt.errorbar(ci[2],dataR[0],yerr=dataR[1],label= 'electrostatics R')
    axes = plt.gca()
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
    plt.xlabel("N Frames")
    plt.ylabel("Energy ()")
    plt.legend()
    fig.savefig('El.png',dpi=400,bbox_inches='tight')

def plot_energy(folderR,folderS,name): #for whole dataframe
    dataR = read_energy(folderR)
    dataS = read_energy(folderS)
    
    fig, ax = plt.subplots(figsize = (4,2.5))
  #  plt.plot(dataS[0][0],label= 'bonds S'i,linewidth = 0.5)
  #  plt.plot(dataR[0][0],label= 'bonds R',linewidth = 0.5)
  #  plt.plot(dataS[0][1],label= 'angles S',linewidth = 0.5)
  #  plt.plot(dataR[0][1],label= 'angles R',linewidth = 0.5)
  #  plt.plot(dataS[0][2],label= 'dihedrals S',linewidth = 0.5)
  #  plt.plot(dataR[0][2],label= 'dihedrals R',linewidth = 0.5)
  #  plt.plot(dataS[0][3],label= 'LJ S',linewidth = 0.5)
  #  plt.plot(dataR[0][3],label= 'LJ R',linewidth = 0.5)
  #  plt.plot(dataS[0][4],label= 'repulsion S',linewidth = 0.5)
  #  plt.plot(dataR[0][4],label= 'repulsion R',linewidth = 0.5)
    plt.plot(dataS[0][6],label= 'electrostatics S',linewidth = 0.5)
    plt.plot(dataR[0][6],label= 'electrostatics R',linewidth = 0.5)
    axes = plt.gca()
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
    plt.xlabel("N Frames")
    plt.ylabel("Energy (kcal/mol)")
    plt.title("{}".format(name))
    plt.legend()
    fig.savefig('{}'.format(name)+'El.png',dpi=400,bbox_inches='tight')

#plot_energy(R[3],S[3],2)
