#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
import math
import glob
class COoperations:
    
    def COXYZ(coord):
        """Returns x,y,z coordinates as floats."""
        x = []
        y = []
        z = []
        CO_coord = []
        for xi in coord[0]:
            x.append(float(xi))
        for yi in coord[1]:
            y.append(float(yi))
        for zi in coord[2]:
            z.append(float(zi))
        CO_coord = (x,y,z)
        return(CO_coord)
    
    def COM(coord):
        """Returns the center of mass of the defined region."""
        x = []
        y = []
        z = []
        for xi in coord[0]:
            x.append(float(xi))
        for yi in coord[1]:
            y.append(float(yi))
        for zi in coord[2]:
            z.append(float(zi))
        atoms = len(x)
        x_com = sum(x)/atoms
        y_com = sum(y)/atoms
        z_com = sum(z)/atoms
        COM_coord = (x_com,y_com,z_com)
        return(COM_coord)
    
    def CXY(coord):
        """Calculates the displacement with respect to the x and y-axis for in a given simulation frame. 
        Returns the root mean square of both points."""
        x = []
        y = []
        for xi in coord[0]:
            x.append(float(xi))
        for yi in coord[1]:
            y.append(float(yi))
        atoms = len(x)
        x_com = sum(x)/atoms
        y_com = sum(y)/atoms
        xy_coord = np.sqrt(x_com**2+y_com**2)
        CXY_coord = xy_coord
        return(CXY_coord)
    
    def COXY(coord):
        """Calculates the displacement with respect to the x and y-axis for in a given simulation frame.
        Returns x and y coordinate."""
        x = []
        y = []
        for xi in coord[0]:
            x.append(float(xi))
        for yi in coord[1]:
            y.append(float(yi))
        atoms = len(x)
        x_com = sum(x)/atoms
        y_com = sum(y)/atoms
        xy_coord =( x_com,y_com)
        return(xy_coord)
    
    def CZ(coord):
        """Returns z-coordinate."""
        z = []
        for zi in coord[2]:
            z.append(float(zi))
        atoms = len(z)
        z_com = sum(z)/atoms
        CZ_coord = z_com
        return(CZ_coord)

class MSD:
    def MSD(tau,coordinates):
        """Calculated Mean square displacement at a certain tau value."""
        step = (lambda x,x0 : (x-x0)**2)
        iterations = len(coordinates)//tau
        steps = []
        for i in range(0,len(coordinates)-tau,tau):
            steps.append(step(coordinates[i+tau],coordinates[i]))
        MSD = sum(steps)/iterations
        #propably also evaluate the variance of this data
        return(MSD)
    
    def MSDtau(coordinates):
        """Calculates the mans sqaure displacement for tau values from 1 to 200 in steps of 1."""
        taus = np.arange(1,200,1)
        MSDtau=[]
        for i in taus:
            msd = MSD.MSD(i,coordinates)
            MSDtau.append(msd)
        return(MSDtau)

class geom:
    def determinant(v1,v2):
        "Returns the determinant of two vectors."
        return(v1[0]*v2[1]-v1[1]*v2[0])
    
    def Angle(v1,v2):
        "Returns the angle between two vectors."
        cosine_ang = np.dot(v1,v2)/(math.sqrt(v1[0]**2 + v1[1]**2)*math.sqrt(v2[0]**2 + v2[1]**2))
        "Make sure the values are between 1 and -1."
        cosine_ang = np.clip(cosine_ang, -1, 1)

        angleA = np.arccos(cosine_ang)

        if geom.determinant(v1,v2)<0: #this is a property of the det. If the det < 0 then v1 is clockwise of v2
            return(angleA)
        else: # if the det > 0 then v1 is immediately clockwise of v2
            return(2*(np.pi)-angleA)
    
    def dist(P1,P2):
        "Retruns the distance between P1 and P2."
        P1x = P1[0]
        P1y = P1[1]
        P1z = P1[2]
        P2x = P2[0]
        P2y = P2[1]
        P2z = P2[2]
        D = np.sqrt( (P2x - P1x)**2 + (P2y - P1y)**2 + (P2z - P1z)**2 )
        return(D)

class Iterator:
    """Iterate over coordinates--> adjust for each model: monomer: 59
                                                    dimer: 118
                                                    N-core: 288
                                                    full-lenth: 
    Define the number of models and set step.
    Retruns the coordinates fo the beads in the trajectory file.
    lp = 288 length of the Ncore sequence,
    lp = 59 length of the monomer sequence"""
    
    def readlinesfromfile(file):
        """Returns list of lines in the file."""
        lines = []
        with open(file, "r") as f:
            lines.append(f.readlines())
        return(lines)
    
    def traj(file,framecount,step):
        lines = Iterator.readlinesfromfile(file)
        lp = 118
        mods = []
        for i in range(0,framecount,step):
            mods.append(lines[0][(lp+606)+i*(lp+2):(lp+606)+(i+1)*(lp+2)][1:(lp+1)])
        return(mods)
    

    def pdb(file):
        """Retruns coordinates of CA atoms in pdb."""
        lines = Iterator.readlinesfromfile(file)
        mods = []
        for i in range(len(lines)):
            if 'CA' in lines[i]:
                mods.append(lines[i])
        return(mods)
    
class Region():
    """Retruns coordinates of distinct regions in lacI protein."""

    def recA(data): 
        recA = []
        for frame in data:
            recA.append(frame[16:24])
        return(recA)
    
    def recB(data):
        lp = len(data[0])/2
        mono=lp/2 
        recB = []
        B1 = int(16+mono)
        B2 = int(24+mono)
        for frame in data:
            recB.append(frame[B1:B2])
        return(recB)
    
    def hiA(data):
        hiA = []
        for frame in data:
            hiA.append(frame[50:56])
        return(hiA)
    
    def hiB(data):
        lp = len(data[0])/2
        mono=lp/2 
        Bh1 = int(50+mono)
        Bh2 = int(56+mono)
        hiB = []
        for frame in data:
            hiB.append(frame[Bh1:Bh2])
        return(hiB)

#class coordinates_pdb():
#    def get_coord(self):
#        """Reruns the coordinates in float form for calculation"""
#        Xcord = []
#        Ycord = []
#        Zcord = []
#        for line in self:
#            Xcord.append(line[30:38])
#            Ycord.append(line[38:46])
#            Zcord.append(line[46:54])
#        return(Xcord,Ycord,Zcord)
#    def coord_F(self):
#        X,Y,Z = coordinates_pdb.get_coord(self)
#        x=[]
#        y=[]
#        z=[]
#        for i in range(len(X)):
#            x.append(float(X[i]))
#            y.append(float(Y[i]))
#            z.append(float(Z[i]))
#        return(x,y,z)
#    
#class coordinates():
#    def get_coord(self):
#        """Reruns the coordinates in float form for calculation"""
#        Xcord = []
#        Ycord = []
#        Zcord = []
#        for line in self:
#            Xcord.append(line[0:8])
#            Ycord.append(line[8:16])
#            Zcord.append(line[16:24])
#        x = []
#        y = []
#        z = []
#        for xi in Xcord:
#            x.append(float(xi))
#        for yi in Ycord:
#            y.append(float(yi))
#        for zi in Zcord:
#            z.append(float(zi))
#        return(x, y, z)
#    
#    def All(iterator,filepath,framecount,step):
#        """Combines the iterator and the folat converter to get the coordinates for each frame."""
#        frames = iterator(filepath,framecount,step=step)
#        Coordinates = []
#        for frame in frames:
#            Coordinates.append(coordinates.get_coord(frame))
#        return(Coordinates)
#
#    def Region(region,iterator,filepath,framecount,step):
#        """Combines the iterator and the float converter to get the coordinates for each frame."""
#        frames = region(iterator(filepath,framecount,step=step))
#        Coordinates = []
#        for frame in frames:
#            Coordinates.append(coordinates.get_coord(frame))
#        return(Coordinates)
#
#    
class Coordinates:   
    def get_coord_traj(lines):
        """Rertuns the coordinates in float form for calculation"""
        Xcord = []
        Ycord = []
        Zcord = []
        for line in lines:
            Xcord.append(line[0:8])
            Ycord.append(line[8:16])
            Zcord.append(line[16:24])
        x = [float(x) for x in Xcord]
        y = [float(y) for y in Ycord]
        z = [float(z) for z in Zcord]
        return(x,y,z)
    
    def get_coord_pdb(lines):
        """Rertuns the coordinates in float form for calculation"""
        Xcord = []
        Ycord = []
        Zcord = []
        for line in lines:
            Xcord.append(line[30:38])
            Ycord.append(line[38:46])
            Zcord.append(line[46:54])
        x = [float(x) for x in Xcord]
        y = [float(y) for y in Ycord]
        z = [float(z) for z in Zcord]
        return(x,y,z)

    def coord_traj(frames):
        """Combines the iterator and the float converter to get the coordinates for each frame."""
        Co = []
        for frame in frames:
            Co.append(Coordinates.get_coord_traj(frame))
        return(Co)
    
    def coord_pdb(frames):
        """Combines the iterator and the float converter to get the coordinates for each frame."""
        Co = []
        for frame in frames:
            Co.append(Coordinates.get_coord_pdb(frame))
        return(Co)


class analyse_traj:
    
    def COM(data):
        """Get the coordinates and do the an operation on it."""
        COMl = []
        COM = COoperations.COM
        for frame in data:
            COMl.append(COM(Coordinates.get_coord_traj(frame))) 
        return(COMl)
    
    def CZ(data):
        """Get the coordinates and do the an operation on it."""
        CZl= []
        CZ = COoperations.CZ
        for frame in data:
             CZl.append(CZ(Coordinates.get_coord_traj(frame))) 
        return(CZl)
    
    def CXY(data):
        """Get the coordinates and do the an operation on it."""
        CXYl = []
        CXY = COoperations.CXY
        for frame in data:
            CXYl.append(CXY(Coordinates.get_coord_traj(frame))) 
        return(CXYl)
    
    def COXY(data):
        """Get the coordinates and do the an operation on it."""
        COXYl = []
        COXY = COoperations.COXY
        for frame in data:
            COXYl.append(COXY(Coordinates.get_coord_traj(frame))) 
        return(COXYl)

    def writer_dis(com,num):
        COMcoord = []
        """Number of frames to analyse"""
        Nof = num
        """List of frame indexes to analyse"""
        numL = list(range(0,num,1)) 
        com = com[0:Nof]
        for coordinates, n in zip(com,numL):
            coord = 'ATOM  {n:4.0f}  CA  AAA A   1    {l[0]:8.3f} {l[1]:8.3f} {l[2]:8.3f}  0.00  0.00\n'.format(l = coordinates, n=n )
            COMcoord.append(coord)
        return(COMcoord)
    
    def writer_traj(com,num):
        COMcoord = []
        """Number of frames to analyse"""
        Nof = num
        """List of frame indexes to analyse"""
        numL = list(range(0,num,1)) 
        com = com[0:Nof]
        for coordinates, n in zip(com,numL):
            coord = 'ATOM  1     CA  AAA A   1    {l[0]:8.3f} {l[1]:8.3f} {l[2]:8.3f}  0.00  0.00\n'.format(l = coordinates, n=n )
            COMcoord.append(coord)
        return(COMcoord)
    
    def glob_COM(folder,framecount, step):


        """Get path of all trajectory files in folder and return the c
oordinates for all atoms and different regions."""
        Trajs = glob.glob(folder+'Traj*.dat')
        Trajs.sort()
        print(Trajs)
        
        All = []
        recA = []
        recB = []
        hiA = []
        hiB = []
        
        iterator = Iterator.traj
        
        for Traj in Trajs:
            Alllines = iterator(Traj,framecount,step)
            All = (analyse_traj.COM(Alllines))
          
            filename='{}'.format(Traj)+'COM.pdb'
            filenameT='{}'.format(Traj)+'COM_Traj.pdb'
            f1 = open(filename,'w')
            for model in range(len(All)):
                f1.write(analyse_traj.writer_dis(All,framecount)[model])
            f1.close()
            f1 = open(filenameT,'w')
            for model in range(len(All)):
                f1.write('MODEL {:8d}\n'.format(int(model)+1))
                f1.write(analyse_traj.writer_traj(All,framecount)[model])
                f1.write('ENDMDL\n')
            f1.close()
             
            recA = (analyse_traj.COM(Region.recA(Alllines))) #get only lines of recA
            filename='{}'.format(Traj)+'COM_recA.pdb'
            filenameT='{}'.format(Traj)+'COM_Traj_recA.pdb'
            f1 = open(filename,'w')
            for model in range(len(recA)):
                f1.write(analyse_traj.writer_dis(recA,framecount)[model])
            f1.close()
            f1 = open(filenameT,'w')
            for model in range(len(recA)):
                f1.write('MODEL {:8d}\n'.format(int(model)+1))
                f1.write(analyse_traj.writer_traj(recA,framecount)[model])
                f1.write('ENDMDL\n')
            f1.close()
            
            recB = (analyse_traj.COM(Region.recB(Alllines))) #get only lines of recB
            filename='{}'.format(Traj)+'COM_recB.pdb'
            filenameT='{}'.format(Traj)+'COM_Traj_recB.pdb'
            f1 = open(filename,'w')
            for model in range(len(recB)):
                f1.write(analyse_traj.writer_dis(recB,framecount)[model])
            f1.close()
            f1 = open(filenameT,'w')
            for model in range(len(recB)):
                f1.write('MODEL {:8d}\n'.format(int(model)+1))
                f1.write(analyse_traj.writer_traj(recB,framecount)[model])
                f1.write('ENDMDL\n')
            f1.close()
            
            f2 = open('{}'.format(Traj)+'DistanceRec.txt','w')
            Dist=[]
            for i in list(range(0,len(recA))):
                Dist.append( geom.dist(recA[i],recB[i]))
                f2.write('{l:7.3f}\n'.format(l = Dist[i] ))
            print(np.mean(Dist),np.std(Dist))
    
    def Displacement(folder,framecount, step):
        """Get path of all trajectory files in folder and return the coordinates for all atoms and different regions."""
        Trajs = glob.glob(folder+'Traj*.dat')
        Trajs.sort()
        print(Trajs)
        
        iterator = Iterator.traj
        
        for Traj in Trajs:
            CZ = analyse_traj.CZ(iterator(Traj,framecount,step))
            CXY = analyse_traj.CXY(iterator(Traj,framecount,step))
            COXY = analyse_traj.COXY(iterator(Traj,framecount,step))
            MSDt = MSD.MSDtau(CZ)
            ang = [] 
            v0 = COXY[0]
            for i in range(len(COXY)):
                 ang.append(geom.Angle(v0,COXY[i])*(180/np.pi))
            f3 = open('{}'.format(Traj)+'Displacement.txt','w')
            for i in range(len(COXY)):
                f3.write('{l:8.3f}  {m:8.3f}  {n:8.3f} \n'.format(l = CZ[i],m=CXY[i],n=ang[i] ))
            f3.close()
            f4 = open('{}'.format(Traj)+'MSD.txt','w')
            for i in range(len(MSDt)):
                f4.write('{l:8.3f} \n'.format(l = MSDt[i] ))
            f4.close()
            print(ang)             

foldS='/Users/mallu899/Desktop/Traj_Slong/'
analyse_traj.Displacement(foldS,100000,1)
