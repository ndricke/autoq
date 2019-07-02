import pandas as pd
import numpy as np
import os
import re
import sys
import getMullikan
import rdkitTest
import doesitbind
import matplotlib.pyplot as plt
import seaborn
import rdkit
from rdkit import Chem
from sklearn.model_selection import train_test_split
from sklearn import linear_model
from sklearn.metrics import mean_squared_error, r2_score
from sklearn import svm
from math import sqrt

def getData(cat, O2, plusone, mol, xyz):
    returndata = doesitbind.collectData(cat, O2, plusone, mol, xyz)
    return(returndata)
    
def getDataFromCSV(csv):
    df = pd.read_csv(csv)
    return (df)

def linearRegression(alldata, f1, f2, f3, f4, f5, f6, f7, f8, f9, f10):

    print(alldata.size)
    alldata = alldata[(alldata["Doesitbind"] == True) & (alldata["BindingEnergy"]>-1.5)]
    print(alldata.size)
        
    cols = []
    features = [f1, f2, f3, f4, f5, f6, f7, f8, f9, f10]
    for feature in features:
        if feature != None:
            cols.append(feature)
    
    print("For the features ", cols)

    X=alldata[cols]
    X = X.values
    y=alldata['BindingEnergy'].astype('float')
    y = y.values
    
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3, random_state=0)

    lin_model = linear_model.LinearRegression()
    svr = svm.SVR(kernel = 'rbf', gamma = 'auto')
    # Train the model using the training sets
    lin_model.fit(X_train, y_train)
    svr.fit(X_train, y_train)
    # Make predictions using the testing set
    y_pred_lin = lin_model.predict(X_test)
    y_pred_svr = svr.predict(X_test)
    # The coefficients
    print('Coefficients: \n', lin_model.coef_)
    # The mean squared error
    print("Root mean squared error, linear: %.2f"
          % sqrt(mean_squared_error(y_test, y_pred_lin)))
    # Explained variance score: 1 is perfect prediction
    print('Variance score: %.2f' % r2_score(y_test, y_pred_lin))
    print("Root mean squared error, svr: %.2f"
          % sqrt(mean_squared_error(y_test, y_pred_svr)))
    # Explained variance score: 1 is perfect prediction
    print('Variance score: %.2f' % r2_score(y_test, y_pred_svr))
    # Plot outputs
    #print(X_test, y_test)
    if len(cols)==1:
        plt.scatter(X, y,  color='black')
        plt.plot(X_test, y_pred_lin, color='blue', linewidth=3)
        #plt.plot(X_test, y_pred_svr, color = 'red', linewidth = 3)

        #plt.xticks(())
        #plt.yticks(())
        plt.xlabel(cols[0])
        plt.ylabel("Binding Energy")

        plt.show()
    else:
        plt.scatter(y_test, y_pred_lin, color='blue')
        #plt.scatter(y_pred_svr, y_test, color = 'red')

        #plt.xticks(())
        #plt.yticks(())
        
        plt.xlabel("Binding energy (actual) (eV)")
        plt.ylabel("Binding energy (predicted) (eV)")

        plt.show()

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser("Linear regression of data")
    parser.add_argument('-csv', help = 'doesitbind output csv', type=str, default = None)
    parser.add_argument('-cat', help='input file/directory (bare catalyst .out)', type=str)
    parser.add_argument('-O2', help='O2 bound directory', type=str)
    parser.add_argument('-plusone', help='Plus one data for bare catalyst to grab charge difference data', default = None, type=str)
    parser.add_argument('-mol', help='Mol files for bare catalyst to grab neighbor data', default = None, type=str)
    parser.add_argument('-xyz', help='XYZ files for bare catalyst to grab bond data', default = None, type=str)
    parser.add_argument('-f1', help = 'First feature to compare against, default to spin density', default = 'SpinDensity', type = str)
    parser.add_argument('-f2', help = "Second feature to compare against, like ChElPGNeutralCharge, Doesitbind, BindingEnergy, NeutralFreeEnergy, NeighborSpinDensity, NeighborChElPGCharge, NeighborChargeDifference", type = str)
    parser.add_argument('-f3', help = 'Third feature', default = None)
    parser.add_argument('-f4', help = 'Fourth feature', default = None)
    parser.add_argument('-f5', help = 'Fifth feature', default = None)
    parser.add_argument('-f6', help = 'Sixth feature', default = None)
    parser.add_argument('-f7', help = 'Seventh feature', default = None)
    parser.add_argument('-f8', help = 'Eighth feature', default = None)
    parser.add_argument('-f9', help = 'Ninth feature', default = None)
    parser.add_argument('-f10', help = 'Tenth feature', default = None)


    args = parser.parse_args()

    if (args.csv != None):
        alldata = getDataFromCSV(args.csv)
    else:
        alldata = getData(args.cat, args.O2, args.plusone, args.mol, args.xyz)
        
    linearRegression(alldata, args.f1, args.f2, args.f3, args.f4, args.f5, args.f6, args.f7, args.f8, args.f9, args.f10)
