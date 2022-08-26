hhkk
def calc_D(CZ,start,stop):

    "Calculate the diffusion coefficient from MSD."

    y = MSD.MSDtau(CZ1)[start:stop+1]

    x = np.arange(start,stop,1).reshape(-1,1)


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

    slope = regression_model.coef_

    D = slope/2

    D #in A^2 tu^-1

    print('Diffusion constant: ', D)

