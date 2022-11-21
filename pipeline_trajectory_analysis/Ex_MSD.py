#!/usr/bin/python
import numpy as np
import matplotlib.pyplot as plt
import math
import glob
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
