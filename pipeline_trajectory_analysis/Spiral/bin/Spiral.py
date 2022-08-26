#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
import math
import glob

from classes import COoperations
from classes import geom
from classes import MSD
from classes import coordinates 

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


