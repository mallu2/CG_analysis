#!/Users/mallu899/anaconda3/bin/python
import numpy as np
import math

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
    
    def CXYCA(coord):
        """Calculates the displacement with respect to the x and y-axis for in a given simulation frame.
        Returns x and y coordinate for each CA."""
        x = []
        y = []
        for xi in coord[0]:
            x.append(float(xi))
        for yi in coord[1]:
            y.append(float(yi))
        xarray = np.array(x)
        yarray = np.array(y)

        xy_coord = np.sqrt(xarray**2+yarray**2)
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
        for i in range(0,len(coordinates)-tau,1):
            steps.append(step(coordinates[i+tau],coordinates[i]))
        iterations = len(steps)
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
    
    def AngleN(v1,v2):
        "Returns the angle between two vectors."
        cosine_ang = np.dot(v1,v2)/(math.sqrt(v1[0]**2 + v1[1]**2)*math.sqrt(v2[0]**2 + v2[1]**2))
        "Make sure the values are between 1 and -1."
        cosine_ang = np.clip(cosine_ang, -1, 1)

        angleA = np.arccos(cosine_ang)
        return(angleA)
    
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
