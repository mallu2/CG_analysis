#!/usr/bin/python
import sys
sys.path.append('/home_d/malin/scripts/AnalysisPy/')

import numpy as np
import matplotlib.pyplot as plt
import math
import glob

from Spiral import COoperations
from Spiral import geom
from Spiral import MSD
from Spiral import Coordinates 
from Spiral import plot 
from Spiral import Region
from Spiral import Iterator_Nb as Iterator

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

    def writer_dis(com,num):
        COMcoord = []
        """Number of frames to analyse"""
        Nof = num
        """List of frame indexes to analyse"""
        numL = list(range(0,num,1)) 
        com = com[0:Nof]
        for Coordinates, n in zip(com,numL):
            coord = 'ATOM  {n:4.0f}  CA  AAA A   1    {l[0]:8.3f} {l[1]:8.3f} {l[2]:8.3f}  0.00  0.00\n'.format(l = Coordinates, n=n )
            COMcoord.append(coord)
        return(COMcoord)
    
    def writer_traj(com,num):
        COMcoord = []
        """Number of frames to analyse"""
        Nof = num
        """List of frame indexes to analyse"""
        numL = list(range(0,num,1)) 
        com = com[0:Nof]
        for Coordinates, n in zip(com,numL):
            coord = 'ATOM  1     CA  AAA A   1    {l[0]:8.3f} {l[1]:8.3f} {l[2]:8.3f}  0.00  0.00\n'.format(l = Coordinates, n=n )
            COMcoord.append(coord)
        return(COMcoord)
    
    def glob_COM(folder,framecount, step):

        """Get path of all trajectory files in folder and return the Coordinates for all atoms and different regions."""
        Trajs = glob.glob(folder+'Traj*.dat')
        Trajs.sort()
        print(Trajs)
        
        All = []
        H2A = []
        H2B = []
        hiA = []
        hiB = []
        
        iterator = Iterator.traj
        
        Dist_all = []
        Ang_all = []
        for Traj in Trajs:

            Alllines = iterator(Traj,framecount,step)
            All = (analyse_traj.COM(Alllines))
          
            filename='{}'.format(Traj)+'COM.pdb'
            filenameT='{}'.format(Traj)+'COM_Traj.pdb'
            f1 = open(filename,'w')
            for model in range(len(All)):
                f1.write(analyse_traj.writer_dis(All,framecount)[model])
            f1.close()
            fT1 = open(filenameT,'w')
            for model in range(len(All)):
                fT1.write('MODEL {:8d}\n'.format(int(model)+1))
                fT1.write(analyse_traj.writer_traj(All,framecount)[model])
                fT1.write('ENDMDL\n')
            fT1.close()
             
            H2A = (analyse_traj.COM(Region.H2A(Alllines))) #get only lines of H2A
            H2AXY = (analyse_traj.COXY(Region.H2A(Alllines))) #get only lines of H2A
            filename='{}'.format(Traj)+'COM_H2A.pdb'
            filenameT='{}'.format(Traj)+'COM_Traj_H2A.pdb'
            f1 = open(filename,'w')
            for model in range(len(H2A)):
                f1.write(analyse_traj.writer_dis(H2A,framecount)[model])
            f1.close()
            fT1 = open(filenameT,'w')
            for model in range(len(H2A)):
                fT1.write('MODEL {:8d}\n'.format(int(model)+1))
                fT1.write(analyse_traj.writer_traj(H2A,framecount)[model])
                fT1.write('ENDMDL\n')
            fT1.close()
            
            H2B = (analyse_traj.COM(Region.H2B(Alllines))) #get only lines of H2B
            H2BXY = (analyse_traj.COXY(Region.H2B(Alllines))) #get only lines of H2B
            filename='{}'.format(Traj)+'COM_H2B.pdb'
            filenameT='{}'.format(Traj)+'COM_Traj_H2B.pdb'
            f1 = open(filename,'w')
            for model in range(len(H2B)):
                f1.write(analyse_traj.writer_dis(H2B,framecount)[model])
            f1.close()
            fT1 = open(filenameT,'w')
            for model in range(len(H2B)):
                fT1.write('MODEL {:8d}\n'.format(int(model)+1))
                fT1.write(analyse_traj.writer_traj(H2B,framecount)[model])
                fT1.write('ENDMDL\n')
            fT1.close()
            
            f2 = open('{}'.format(Traj)+'DistanceRec.txt','w')
           
            Dist=[]
            Ang=[]
            for i in list(range(0,len(H2A))):
                dist = geom.dist(H2A[i],H2B[i])
                ang = geom.AngleN(H2AXY[i],H2BXY[i])
                Dist.append( dist )
                Dist_all.append( dist )
                Ang.append( ang *(180/np.pi) )
                Ang_all.append( ang *(180/np.pi) )
                f2.write('{l:7.3f}  {a:7.3f}\n'.format(l = Dist[i], a = Ang[i] ))
            f2.write('Avergae over all frames:\n')
            f2.write('{m:7.3f} {s:7.3f}\n'.format(m=np.mean(Dist),s=np.std(Dist)))
            f2.write('{a:7.3f} {t:7.3f}\n'.format(a=np.mean(Ang),t=np.std(Ang)))
            f2.close()
        
        print(np.mean(Dist_all),np.std(Dist_all))
        print(np.mean(Ang_all),np.std(Ang_all))
    
    def Displacement(folder,framecount, step):
        """Get path of all trajectory files in folder and return the Coordinates for all atoms and different regions."""

        Trajs = glob.glob(folder+'Traj*.dat')
        Trajs.sort()
        
        iterator = Iterator.traj
        
        for Traj in Trajs:
            CZ = analyse_traj.CZ(iterator(Traj,framecount,step))
            CXY = analyse_traj.CXY(iterator(Traj,framecount,step))
            COXY = analyse_traj.COXY(iterator(Traj,framecount,step))
            MSDt = MSD.MSDtau(CZ)
            print(MSDt)
            ang = [] 
            v0 = COXY[0]
            for i in range(len(COXY)):
                 ang.append(geom.Angle(COXY[i],v0))
                 uang = np.unwrap(ang) 
                 #rang = [x*(180/np.pi) for x in uang] 
            f3 = open('{}'.format(Traj)+'Displacement.txt','w')
            for i in range(len(COXY)):
                f3.write('{l:8.3f}  {m:8.3f}  {n:8.3f} \n'.format(l = CZ[i],m=CXY[i],n=uang[i] ))
            f3.close()
            f4 = open('{}'.format(Traj)+'MSD.txt','w')
            for i in range(len(MSDt)):
                f4.write('{l:8.3f} \n'.format(l = MSDt[i] ))
            f4.close()

        for Traj in Trajs:
            CZA = analyse_traj.CZ(Region.H2A(iterator(Traj,framecount,step)))
            CXYA = analyse_traj.CXY(Region.H2A(iterator(Traj,framecount,step)))
            COXYA = analyse_traj.COXY(Region.H2A(iterator(Traj,framecount,step)))
            angA = [] 
            v0A = COXY[0]
            for i in range(len(COXY)):
                 angA.append(geom.Angle(COXYA[i],v0A))
                 uangA = np.unwrap(angA) 
                 #rang = [x*(180/np.pi) for x in uang] 
            f3 = open('{}'.format(Traj)+'Displacement_H2A.txt','w')
            for i in range(len(COXYA)):
                f3.write('{l:8.3f}  {m:8.3f}  {n:8.3f} \n'.format(l = CZA[i],m=CXYA[i],n=uangA[i] ))
            f3.close()

        for Traj in Trajs:
            CZB = analyse_traj.CZ(Region.H2B(iterator(Traj,framecount,step)))
            CXYB = analyse_traj.CXY(Region.H2B(iterator(Traj,framecount,step)))
            COXYB = analyse_traj.COXY(Region.H2B(iterator(Traj,framecount,step)))
            angB = [] 
            v0B = COXYB[0]
            for i in range(len(COXYB)):
                 angB.append(geom.Angle(COXYB[i],v0B))
                 uangB = np.unwrap(angB) 
                 #rang = [x*(180/np.pi) for x in uang] 
            f3 = open('{}'.format(Traj)+'Displacement_H2B.txt','w')
            for i in range(len(COXYB)):
                f3.write('{l:8.3f}  {m:8.3f}  {n:8.3f} \n'.format(l = CZB[i],m=CXYB[i],n=uangB[i] ))
            f3.close()

from pathlib import Path
mypath = Path().absolute()
folder='{}/'.format(mypath)

#analyse_traj.glob_COM(folder,100000,10)
#analyse_traj.Displacement(folder,10000,1)
#plot.plot(folder)

from sklearn.linear_model import LinearRegression

from sklearn.metrics import mean_squared_error, r2_score

def readlinesfromfile(file):
    lines = []
    with open(file, "r") as f:
        lines.append(f.readlines())
    return(lines)

def calc_D(folder,start,stop):
    "Calculate the diffusion coefficient from MSD."
    
    files = glob.glob(folder+'*MSD.txt')
    slope, D = [],[]
    for i in range(len(files)):
        data = readlinesfromfile(files[i])
        MSD = [float(x) for x in data[0]]
    
        y = MSD[start:stop+1]
        x = np.arange(start+1,stop,1).reshape(-1,1)
        
        # sckit-learn implementation
        # Model initialization
        regression_model = LinearRegression()

        # Fit the data(train the model)
        regression_model.fit(x, y)

        # Predict
        y_predicted = regression_model.predict(x)

        # model evaluation
        rmse = mean_squared_error(y, y_predicted)
        r2 = r2_score(y, y_predicted)

        # printing values
        print('Slope:' ,regression_model.coef_)
        print('Intercept:', regression_model.intercept_)
        print('Root mean squared error: ', rmse)
        print('R2 score: ', r2)

        slope.append(regression_model.coef_)
    slope = np.array(slope)
    D = slope/2  #in A^2 tu^-1
    print('Diffusion constant: ', D)

calc_D(folder,50,200)
