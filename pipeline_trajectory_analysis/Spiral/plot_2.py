"Plot MSD and plain displacement."

import matplotlib.pyplot as plt
import glob
import pandas as pd
import numpy as np

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
        plt.ylabel('MSD (A^2)')
        plt.legend()
        plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        plt.savefig('{}'.format(folder)+'MSD.png',dpi=400,bbox_inches='tight')
        
        "Plot Displacement for each trajectory."
        CZ = []
        CXY = []
        ang = []
        CZA = []
        CXYA = []
        angA = []
        CZB = []
        CXYB = []
        angB = []

        files = glob.glob(folder+'*Displacement.txt')
        files.sort()
        filesA = glob.glob(folder+'*Displacement_H2A.txt')
        filesA.sort()
        filesB = glob.glob(folder+'*Displacement_H2B.txt')
        filesB.sort()
        print(files,filesA) 
        for file in files:
            df = pd.read_csv(file, sep='\s+',header=None)
            CZ.append(df[0])
            CXY.append(df[1])
            ang.append(df[2])
        
        for fileA in filesA:
            dfA = pd.read_csv(fileA, sep='\s+',header=None)
            CZA.append(dfA[0])
            CXYA.append(dfA[1])
            angA.append(dfA[2])
        
        for fileB in filesB:
            dfB = pd.read_csv(fileB, sep='\s+',header=None)
            CZB.append(dfB[0])
            CXYB.append(dfB[1])
            angB.append(dfB[2])
        
        "Plot xy radius of H2A and H2B against the Z-displacement of H2A and H2B" 
        for i in range(len(files)):
            fig, ax = plt.subplots(figsize = (4,2.5))
           # plt.plot(CZ[i], label = 'Z trajectory')
           # plt.plot(CXY[i], label= 'XY trjectory')
            plt.plot(CZA[i], label = 'Z trajectory A',linewidth = 0.5,color='k')
            plt.plot(CXYA[i], label= 'XY trjectory A',linewidth = 0.5,color='r')
            plt.plot(CZB[i], label = 'Z trajectory B',linewidth = 0.5,color='grey')
            plt.plot(CXYB[i], label= 'XY trjectory B',linewidth = 0.5,color='yellow')
            axes = plt.gca()
            axes.set_ylim([-150,150])
            plt.xlabel('Simulation Time step')
            plt.ylabel('Displacement (A)') 
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
        
        "Plot xy radius of H2A against H2B." 
        for i in range(len(files)):
            fig, ax = plt.subplots(figsize = (4,2.5))
            color = np.arange(len(CXYA[0]))
            plt.scatter(CXYB[i], CXYA[i],s = 0.5,c=color)
            rr = [10,20,30,40,50,60]
            ax.set_xticks(rr)
            axes = plt.gca()
            axes.set_aspect('equal', adjustable='box')
            #axes.set_ylim([-150,150])
            ax.set_yticks(rr)
            # ax.set_xticklabels(['90:A','180:A','270:A','90:B','180:B','270:B'])
            #ax.axis('equal')
            #ax.set_aspect('equal')
            
            plt.xlabel('Radius H2B (A)')
            plt.ylabel('Radius H2A (A)') 
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
            cb = plt.colorbar() 
            #plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
            plt.savefig('{}'.format(files[i])+'H2AH2B.png',dpi=400,bbox_inches='tight')
        
        "Plot angle vs rotation"
        for i in range(len(files)):
            colors = np.arange(len(CZ[0]))
            fig, ax = plt.subplots(figsize = (4.5,2.5))
            #plotting in degrees
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
    
        "Plot angle vx rotation"
        for i in range(len(files)):
            colors = np.arange(len(CZ[0]))
            fig, ax = plt.subplots(figsize = (4.5,2.5))
            plt.scatter(CZA[i],angA[i]*(180/np.pi),s=0.5,c=colors,edgecolors='none')
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
            plt.savefig('{}'.format(files[i])+'ScatterA.png',dpi=400,bbox_inches='tight')
    
        "Plot angle vx rotation"
        for i in range(len(files)):
            colors = np.arange(len(CZ[0]))
            fig, ax = plt.subplots(figsize = (4.5,2.5))
            plt.scatter(CZB[i],angB[i]*(180/np.pi),s=0.5,c=colors,edgecolors='none')
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
            plt.savefig('{}'.format(files[i])+'ScatterB.png',dpi=400,bbox_inches='tight')
    
