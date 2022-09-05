"Plot MSD and plain displacement."

import matplotlib.pyplot as plt
import glob
import pandas as pd

def readlinesfromfile(file):
    lines = []
    with open(file, "r") as f:
        lines.append(f.readlines())
    return(lines)

class plot:
    def plot(folder):
        files = glob.glob(folder+'*MSD.txt')
        fig, ax = plt.subplots(figsize = (4,2.5))
        for i in range(len(files)):
            data = readlinesfromfile(files[i])
            MSD = [float(x) for x in data[0]]
            plt.plot(MSD,label='{}'.format(i+1))
    
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
    
        plt.xlabel('tau value')
        plt.ylabel('MSD (A)')
        plt.legend()
        plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        plt.savefig('{}'.format(folder)+'MSD.png',dpi=400,bbox_inches='tight')
        
        "Plot Displacement for each trajectory."
        CZ = []
        CXY = []
        ang = []

        files = glob.glob(folder+'*Displacement.txt')
        for file in files:
            df = pd.read_csv(file, sep='\s+',header=None)
            CZ.append(df[0])
            CXY.append(df[1])
            ang.append(df[2])

        for i in range(len(files)):
            fig, ax = plt.subplots(figsize = (4,2.5))
            plt.plot(CZ[i], label = 'Z trajectory')
            plt.plot(CXY[i], label= 'XY trjectory')
            axes = plt.gca()
            axes.set_ylim([-150,150])
            plt.xlabel('Displacement (A)')
            plt.ylabel('Rotation Angle (deg)')
    
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
            plt.savefig('{}'.format(files[i])+'Disp.png',dpi=400,bbox_inches='tight')
        
        "Plot angle vx rotation"
        import numpy as np
        for i in range(5):
            colors = np.arange(len(CZ[0]))
            fig, ax = plt.subplots(figsize = (4.5,2.5))
            plt.scatter(CZ[i],ang[i]*(180/np.pi),s=0.5,c=colors,edgecolors='none')
            cb = plt.colorbar()
          #  plot(ax, i)
            axes = plt.gca()
            axes.set_xlim([-150,150])
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
            plt.xlabel('Z Displacement (A)')
            plt.ylabel('Rotation Angle (deg)')
            plt.savefig('{}'.format(files[i])+'Scatter.png',dpi=400,bbox_inches='tight')
    
