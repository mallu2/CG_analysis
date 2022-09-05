import glob
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np


fR = '/home_d/malin/Sim/dimer_model/2-MD/R-helix/AprilkBT04/6contacts/04kBT/'
fR59off = '/home_d/malin/Sim/dimer_model/2-MD/R-helix/AprilkBT04/6contacts/Lys59off_6contacts/'
#fRN = '/home_d/malin/Sim/dimer_model/2-MD/R_Ncore/salt04/'
fS = '/home_d/malin/Sim/dimer_model/2-MD/Model12/kBT04/'
fS59off = '/home_d/malin/Sim/dimer_model/2-MD/Model12/nochargeLys59unimore/'

c = ['01','02','03','04']

R = [fR + x + '/output/Analisis/' for x in c]
R59o = [fR59off + x + '/output/Analisis/' for x in c]
#RN = [fRN + x + '/output/Traj/' for x in c]
S = [fS + x + '/output/Traj/' for x in c]  
S59o = [fS59off + x + '/output/Traj/' for x in c]  

def read_dis(folder):
    folderDist = folder
    files = glob.glob(folderDist+'*DistanceRec.txt')
    Dist = []
    for file in files:
        df = pd.read_csv(file, sep='\s+',header=None,skipfooter=3,engine='python')
        dist = df[0]
        Dist.append(np.array(dist))
    Distf = [x for y in Dist for x in y]
    return(pd.DataFrame(Distf))

def plot_dis(folderR,folderS):
    dataR = read_dis(folderR)
    dataS = read_dis(folderS)
    df = pd.concat([dataR,dataS], ignore_index=True, axis=1)
    df.rename(columns={0:'R',1:'S'}, inplace=True)
    colorS = ['brown','indianred','lightcoral','darksalmon','lightsalmon']
    colorR = ['darkcyan','darkturquoise','turquoise','cyan','paleturquoise']    
    colors = [colorR[2],colorS[1]]
    
    fig, ax = plt.subplots(figsize = (4,2.5))
    df.plot.hist(ax=ax,bins=range(21,45,2), color=colors, density=True, alpha = 0.7)
    df.plot.kde(ax=ax,color=colors, alpha = 1, legend=False)
    #dataR.plot.hist( bins=range(21,45,2), color=colors[0], density=True, label=['R'])
    #dataS.plot.hist( bins=range(21,45,2), color=colors[1], density=True, label=['S'])
    axes = plt.gca()
    axes.set_ylim([0,0.25])
    axes.set_xlim([5,45])
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
    plt.xlabel("Distance between recognition helixes (A)")
    plt.ylabel("Frequency")
    plt.title("K59-, ci = 0.04 M")
    fig.savefig('DistcancerecArecB.png',dpi=400,bbox_inches='tight')


#plot_dis(R[3],S[3])
plot_dis(R59o[3],S59o[3])

