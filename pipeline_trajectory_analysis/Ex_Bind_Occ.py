#!/usr/bin/python
import sys
sys.path.append('/home_d/malin/scripts/AnalysisPy/')

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import math
import glob

from Spiral import COoperations
from Spiral import geom
from Spiral import MSD
from Spiral import Coordinates 
from Spiral import plot 
from Spiral import Region
from Spiral import Iterator_full as Iterator

topo_N = '/home_d/malin/scripts/AnalysisPy/Binding/ADNAex.dat'
targetC = '/home_d/malin/scripts/AnalysisPy/Binding/fullLength.txt'

class P:
    def readlinesfromfile(file):
        """Returns list of lines in the file."""
        lines = []
        with open(file, "r") as f:
            lines.append(f.readlines())
        return(lines)
    
    def P_coordinates(DNA):
        """Get the coordinates for the P in the DNA in float form."""
        DNA_lines = P.readlinesfromfile(DNA)
        P_lines_C = [line for line in DNA_lines[0] if 'P    C' in line]
        P_lines_G = [line for line in DNA_lines[0] if 'P    G' in line]
        P_lines_A = [line for line in DNA_lines[0] if 'P    A' in line]
        P_lines_T = [line for line in DNA_lines[0] if 'P    T' in line]
        P_lines = P_lines_C + P_lines_G + P_lines_T + P_lines_A
        P_x = []
        P_y = []
        P_z = []
        P_index = []
        for line in P_lines:
            P_index.append(float(line[5:9]))
            P_x.append(float(line[16:24]))
            P_y.append(float(line[24:32]))
            P_z.append(float(line[32:40]))
        dfP1 = pd.DataFrame(P_index)
        dfP = dfP1.rename({0:'ID'}, axis='columns')
        dfP['X'] = P_x
        dfP['Y'] = P_y
        dfP['Z'] = P_z
        
        return(dfP)

class analyse_traj:
    
    def COM(data):
        """Get the Coordinates and do the an operation on it."""
        COMl = []
        for frame in data:
            COMl.append(COoperations.COM(Coordinates.get_coord_traj(frame))) 
        return(COMl)
    
    def CZ(data):
        """Get the Coordinates and do the an operation on it."""
        CZl= []
        for frame in data:
             CZl.append(COoperations.CZ(Coordinates.get_coord_traj(frame))) 
        return(CZl)
    
    def CXY(data):
        """Get the Coordinates and do the an operation on it."""
        CXYl = []
        for frame in data:
            CXYl.append(COoperations.CXY(Coordinates.get_coord_traj(frame))) 
        return(CXYl)
    
    def COXY(data):
        """Get the Coordinates and do the an operation on it."""
        COXYl = []
        for frame in data:
            COXYl.append(COoperations.COXY(Coordinates.get_coord_traj(frame))) 
        return(COXYl)

    def CXYCA(data):
        """Get the Coordinates and do the an operation on it."""
        CXYCAl = []
        for frame in data:
            CXYCAl.append(COoperations.CXYCA(Coordinates.get_coord_traj(frame))) 
        return(CXYCAl)
    
           
    def Displacement(folder,framecount, step):
        """Get path of all trajectory files in folder and return the Coordinates for all atoms and different regions."""

        Trajs = glob.glob(folder+'Traj*.dat')
        Trajs.sort()
        
        iterator = Iterator.traj
        
        for Traj in Trajs:
            CZ = analyse_traj.CZ(iterator(Traj,framecount,step))
            CXY = analyse_traj.CXY(iterator(Traj,framecount,step))
            COXY = analyse_traj.COXY(iterator(Traj,framecount,step))
            ang = [] 
            v0 = COXY[0]
            for i in range(len(COXY)):
                 ang.append(geom.Angle(COXY[i],v0))
                 uang = np.unwrap(ang) 
            f3 = open('{}'.format(Traj)+'Displacement.txt','w')
            for i in range(len(COXY)):
                f3.write('{l:8.3f}  {m:8.3f}  {n:8.3f} \n'.format(l = CZ[i],m=CXY[i],n=uang[i] ))
            f3.close()

    def Footprint(folder,framecount, step):
        """Get path of all trajectory files in folder and return the Coordinates for all atoms and different regions."""

        Trajs = glob.glob(folder+'Traj*.dat')
        Trajs.sort()
        iterator = Iterator.traj
        print(iterator)
        CXI_frame = []
        for Traj in Trajs:
            """Take XY coordinate for each CA bead."""
            CXI_frame.append(analyse_traj.CXYCA(iterator(Traj,framecount,step)))
        Mean_all = []
        STD_all = []
        for frameData in CXI_frame:
            dfi = pd.DataFrame(frameData)
            Mean_all.append(dfi.mean())
            STD_all.append(dfi.std())
        print(Mean_all)
        Mean1=pd.DataFrame(Mean_all)
        Mean = Mean1.mean()
        print(Mean)
        Stdev1 = pd.DataFrame(STD_all)
        Stdev = Stdev1.mean()
        f4 = open('Footprint.txt','w')
        for i in list(range(len(Mean))):
            f4.write('{l:8.3f}  {m:8.3f} \n'.format(l = Mean[i],m=Stdev[i]))
        f4.close()

        return(Mean,Stdev)        


    def CountCont(folder, framecount, step, targetC, topo):
        iterator = Iterator.traj
        """Analyse the formation of contacts to the target side."""    
        "Get the list of taget side contacts"
        Contacts = pd.read_csv(targetC, sep='\s+', header=None) 
        """Get the DNA-Positions for each relevant atom form the topology file."""
        P_coord = P.P_coordinates(topo)
        """Get the relavent coordinates for the contacts you defines."""
        iPs = list(Contacts[2])
        P_contacts = P_coord.loc[P_coord['ID'].isin(iPs)]
        Trajs = glob.glob(folder+'Traj*.dat')
        Trajs.sort()
        dfs = []
        for Traj in Trajs:
            "Analyse each Trajectory individually"
            CoordinatesTraj= []
            for frame in iterator(Traj,framecount,step):
                CoordinatesTraj.append(Coordinates.get_coord_traj(frame)) 
            """Analyse the distances from the dat file."""
            fractionformed = []
            Contactformedall = []
            for Frame in CoordinatesTraj:
                """Iterate over frame and find the atom iR of the repressor,  get the position of the partner P and 
                calculate the distance. Next evaluate of it is a contact or not."""
                Contactformed = []
                for R in list(Contacts.index):
                    iR = int(Contacts[1][R])
                    iP = int(Contacts[2][R])
                    XYZR = (Frame[0][iR-1],Frame[1][iR-1],Frame[2][iR-1])
                    P_contact = P_contacts.loc[P_contacts['ID'] == iP]
                    XYZP = (float(P_contact['X']),float(P_contact['Y']),float(P_contact['Z']))
                    DistanceRP = geom.dist(XYZR,XYZP)
                    iMax =int( Contacts[0][R]+3 )
                    if DistanceRP <= iMax:
                        Contactformed.append(1)
                    else:
                        Contactformed.append(0)
                    Contactformedarray = np.array(Contactformed)
                fractionformed.append(Contactformedarray.sum()/len(Contacts.index))
                Contactformedall.append(Contactformedarray)
            f5 = open('{}'.format(Traj)+'Contacts.txt','w')
            for i in list(range(len(fractionformed))):
                f5.write('{l:8.3f} \n'.format(l = float(fractionformed[i])))
            f5.close()
             
            CDF = pd.DataFrame(Contactformedall)
            CDF.rename(columns={ 0:'Leu6', 1:'Ser16',2:'Tyr17',3:'Thr19',4:'Arg22',5:'Asn25',6:'Ser31',7:'Thr34',8:'Tyr47',9:'Asn50', 10:'Gln54', 11:'Leu6', 12:'Ser16', 13:'Tyr17',14:'Thr19',15:'Arg22',16:'Asn25',17:'Ser31',18:'Thr34',19:'Tyr47',20:'Asn50', 21:'Gln54'}, inplace=True)
            allcont = len(CDF.index)
            residueoccupancies1 = np.array(CDF.sum())
            residueoccupancies = residueoccupancies1/allcont
            print(residueoccupancies)	    
            #df=pd.DataFrame(residueoccupancies)
            CDF.to_csv('{}'.format(Traj)+'Occ.txt',sep=' ',float_format='%.3f',header=('Leu6','Ser16','Tyr17','Thr19','Arg22','Asn25','Ser31','Thr34','Tyr47','Asn50','Gln54','Leu6','Ser16','Tyr17','Thr19','Arg22','Asn25','Ser31','Thr34','Tyr47','Asn50','Gln54'))
            dfSum = pd.DataFrame(CDF.sum())
            dfSum['Occ'] = residueoccupancies
            dfSum.drop(columns=0,inplace=True)
            dfSum.rename_axis('residues', axis='index', inplace=True)            
            dfs.append(dfSum)
        dfs1 = pd.DataFrame(dfs[0])
        for i in range(1,len(dfs)):
            dfsc = pd.concat((dfs1,pd.DataFrame(dfs[i])),axis=1)
            dfs1 = dfsc 
        dfs1M = dfs1.mean(axis=1, numeric_only=True)  
        dfs1S = dfs1.std(axis=1, numeric_only=True)  
        MeanStd = pd.concat((dfs1M,dfs1S),axis=1)
        MeanStd.to_csv('Contacts_split.txt',sep=' ',float_format='%.3f',header=('Mean','Std'))
            

        print('{}'.format(len(Trajs))+'trajectories analysed')                                           

from pathlib import Path
mypath = Path().absolute()
folder='{}/'.format(mypath)

analyse_traj.CountCont(folder, 50000, 10, targetC, topo_N)
#analyse_traj.Footprint(folder,50000,10)
#analyse_traj.Displacement(folder,50000,10)
#plot.plot(folder)
