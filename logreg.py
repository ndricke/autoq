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
from sklearn.linear_model import LogisticRegression
from sklearn import metrics
from sklearn.model_selection import train_test_split
from sklearn.neural_network import MLPClassifier
from sklearn import tree
from sklearn.svm import SVC
from sklearn.svm import LinearSVC
from sklearn.preprocessing import PolynomialFeatures
from sklearn.preprocessing import OneHotEncoder
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import confusion_matrix
from sklearn.utils.multiclass import unique_labels

def plot_coefficients(classifier, feature_names, top_features=5):
    coef = classifier.coef_.ravel()
    top_positive_coefficients = np.argsort(coef)[-top_features:]
    top_negative_coefficients = np.argsort(coef)[:top_features]
    top_coefficients = np.hstack([top_negative_coefficients, top_positive_coefficients])
    # create plot
    plt.figure(figsize=(15, 8))
    colors = ['red' if c < 0 else 'blue' for c in coef[top_coefficients]]
    plt.bar(np.arange(1, (len(feature_names)+1)), coef[top_coefficients], color=colors)
    feature_names = np.array(feature_names)
    plt.xticks(np.arange(1, 1 + 2*top_features), feature_names[top_coefficients], rotation=20, ha='right', fontsize = 10)
    plt.show()

def plot_confusion_matrix(y_true, y_pred,
                          normalize=False,
                          title=None,
                          cmap=plt.cm.Blues):

    if not title:
        if normalize:
            title = 'Normalized confusion matrix'
        else:
            title = 'Confusion matrix, without normalization'

    # Compute confusion matrix
    cm = confusion_matrix(y_true, y_pred)
    # Only use the labels that appear in the data
    #classes = classes[unique_labels(y_true, y_pred)]
    if normalize:
        cm = cm.astype('float') / cm.sum(axis=1)[:, np.newaxis]
        print("Normalized confusion matrix")
    else:
        print('Confusion matrix, without normalization')

    print(cm)

    fig, ax = plt.subplots()
    im = ax.imshow(cm, interpolation='nearest', cmap=cmap)
    ax.figure.colorbar(im, ax=ax)
    # We want to show all ticks...
    ax.set(xticks=np.arange(cm.shape[1]),
           yticks=np.arange(cm.shape[0]),
           # ... and label them with the respective list entries
           #xticklabels=classes, yticklabels=classes,
           title=title,
           ylabel='True label',
           xlabel='Predicted label')

    # Rotate the tick labels and set their alignment.
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right",
             rotation_mode="anchor")

    # Loop over data dimensions and create text annotations.
    fmt = '.2f' if normalize else 'd'
    thresh = cm.max() / 2.
    for i in range(cm.shape[0]):
        for j in range(cm.shape[1]):
            ax.text(j, i, format(cm[i, j], fmt),
                    ha="center", va="center",
                    color="white" if cm[i, j] > thresh else "black")
    fig.tight_layout()
    np.set_printoptions(precision=2)
    plt.show()


def plotData(X, y, dtc, f1, f2):
    x_min, x_max = X[:, 0].min() - .5, X[:, 0].max() + .5
    y_min, y_max = X[:, 1].min() - .5, X[:, 1].max() + .5
    h = .01  # step size in the mesh
    xx, yy = np.meshgrid(np.arange(x_min, x_max, h), np.arange(y_min, y_max, h))
    Z = dtc.predict(np.c_[xx.ravel(), yy.ravel()])

    # Put the result into a color plot
    Z = Z.reshape(xx.shape)
    plt.figure(1, figsize=(4, 3))
    plt.pcolormesh(xx, yy, Z, cmap=plt.cm.Paired)

   # Plot also the training points
    plt.scatter(X[:, 0], X[:, 1], c=y, edgecolors='k', cmap=plt.cm.Paired, s = 80)
    plt.xlabel(f1)
    plt.ylabel(f2)

    plt.xlim(xx.min(), xx.max())
    plt.ylim(yy.min(), yy.max())
    plt.xticks(())
    plt.yticks(())

    plt.show()


def getData(cat, O2, plusone, mol, xyz):
    returndata = doesitbind.collectData(cat, O2, plusone, mol, xyz)
    return(returndata)
    
def getDataFromCSV(csv):
    df = pd.read_csv(csv)
    return (df)

def processData(alldata, features):
    cols = []
    scaledCols = []
    oneHotCols = []
    alreadyProcessedCols = []
    for feature in features:
        if feature != None:
            if feature == "NumberOfHydrogens":
                oneHotCols.append(feature)
                cols.append(feature)
            elif feature in ["OrthoOrPara", "Meta", "FartherThanPara", "IsInRingSize6", "IsInRingSize5", "SpinDensity",
            "Substructure1", "Substructure2", "Substructure3", "Substructure4", "Substructure5", "Substructure6", "Substructure7", "Substructure8", "Substructure9", "Substructure10"]:
                alreadyProcessedCols.append(feature)
                cols.append(feature)
            else:
                scaledCols.append(feature)
                cols.append(feature)
    
    print(scaledCols)
    print(oneHotCols)
    print(alreadyProcessedCols)
    
    scaler = StandardScaler()
    oneHotEncoder = OneHotEncoder(categories = "auto", sparse = False)
    
    scaled_columns = scaler.fit_transform(alldata[scaledCols])
    #print(scaled_columns)
    print(scaled_columns.shape)
    encoded_columns = oneHotEncoder.fit_transform(alldata[oneHotCols])
    #print(encoded_columns)
    print(encoded_columns.shape)
    already_processed_columns = alldata[alreadyProcessedCols]
    already_processed_columns = already_processed_columns.values
    processed_data = np.concatenate((already_processed_columns, scaled_columns, encoded_columns), axis = 1)
    print(processed_data)
    return(processed_data)
    
def classify(alldata, which, f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, c, poly, onehot, confusion, plotcoef, checkmissed, testdata):
    
    #print(alldata.loc[(alldata["SpinDensity"]>0.27)  & (alldata["Doesitbind"]==False)])
    #print(alldata)
    testcolumn = alldata['Doesitbind']
    trues = 0
    falses = 0
    for value in testcolumn:
        if value == True:
            trues +=1
        else:
            falses +=1
    
    print("The number of binding active sites is ", trues)
    print("The number of nonbinding active sites is ", falses)
    
    cols = []
    scaledCols = []
    oneHotCols = []
    alreadyProcessedCols = []
    features = [f1, f2, f3, f4, f5, f6, f7, f8, f9, f10]
    processed_data = processData(alldata, features)
    if testdata != None:
        processed_test_data = processData(testdata, features)
    for feature in features:
        if feature != None:
            cols.append(feature)    
    
    print("For the features ", cols)
    # if poly != 0 :
        # print("With polynomial basis set")
        # poly = PolynomialFeatures()
        # X = alldata[cols]
        # X = poly.fit_transform(X)
        # y = alldata['Doesitbind'].astype('int')
        # y = y.values
    # if onehot!=0 :
        # print("With one-hot encoding")
        # onehot = OneHotEncoder()
        # X = alldata[cols]
        # X = onehot.fit_transform(X)
        # y = alldata['Doesitbind'].astype('int')
        # y = y.values
    # else:
    X=processed_data
    y=alldata['Doesitbind'].astype('int')
    y = y.values

    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3, random_state=0)
    logregression = LogisticRegression(solver = 'lbfgs')
    logregression.fit(X_train, y_train)
    print('Accuracy of logistic regression classifier on test set: {:.2f}'.format(logregression.score(X_test, y_test)))

    #plotData(X, y, logregression)
    
    dtc = tree.DecisionTreeClassifier(max_depth = 2)
    dtc.fit(X_train, y_train)
    print('Accuracy of decision tree classifier on test set: {:.2f}'.format(dtc.score(X_test, y_test)))
    
    mlp = MLPClassifier(max_iter=2500, hidden_layer_sizes = (2048,15))
    mlp.fit(X_train, y_train)
    print('Accuracy of MLP classifier on test set: {:.2f}'.format(mlp.score(X_test, y_test)))
    
    linsvc = LinearSVC(C = float(c))
    linsvc.fit(X_train, y_train)
    print('Accuracy of LinearSVC classifier on test set: {:.2f}'.format(linsvc.score(X_test, y_test)))
    
    svc = SVC(C = float(c), kernel = 'rbf', gamma = 'scale')
    svc.fit(X_train, y_train)
    print('Accuracy of SVC on test set: {:.2f}'.format(svc.score(X_test, y_test)))
    
    
    y_pred = logregression.predict(X_test)
    print(y_pred)

    print(X_train)
    if testdata != None:
        test_pred = linsvc.predict(processed_test_data)
        #for i in range(len(test_pred)):
        success = []
        successes = np.where(test_pred == 1)
        for index in successes:
            for value in index:
                test_spin_density = (processed_test_data[value][0])
                for index, row in testdata.iterrows():
                    if row['SpinDensity'] == test_spin_density:
                        success.append((index, row['Catalyst Name'], row['DistanceToN'], row['SpinDensity']))
        for tuple in success:
            index = tuple[0]
            name = tuple[1]
            distance = tuple[2]
            spin = tuple[3]
            print(name, " at index ", index, "\nDistance to N is ", distance, "\nSpin density is", spin, "\n")
    
    misclassified = np.where(y_test != y_pred)
    #print(misclassified)
    #print(X_test)
    #(misclassified)
    #(len(y_pred))
    missed = []
    if checkmissed!= 0:
        for index in misclassified:
            #print(index)
            for value in index:
                #print(value)
                test_spin_density = (X_test[value][0])
                print(test_spin_density)
                #print(test_spin_density)
                for index, row in alldata.iterrows():
                    if row['SpinDensity'] == test_spin_density:
                        missed.append((index, row['Catalyst Name'], row['Doesitbind'], row['DistanceToN'], row['SpinDensity']))
        for tuple in missed:
            index = tuple[0]
            name = tuple[1]
            bind = tuple[2]
            distance = tuple[3]
            spin = tuple[4]
            print(name, " at index ", index, "\nDoes it bind? ", bind, "\nDistance to N is ", distance, "\nSpin density is", spin, "\n")
    #print(svc.dual_coef_)
    #print(svc.get_params())
    #for index in misclassified:
     #   for index, row in alldata.iterrows():
      #      if row['SpinDensity'
    
    classifier_dict = {"dtc" : dtc, "log" : logregression, "svc" : svc, "mlp" : mlp, "linsvc" : linsvc}


    if len(cols) == 2 and which != None:
        plotData(X, y, classifier_dict[which], f1, f2)
    
    if confusion!=0:
        plot_confusion_matrix(y_test, y_pred) 
    
    if plotcoef!=0:
        plot_coefficients(linsvc, cols)
    
    return(svc.score(X_test, y_test))
    


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser("Grab data about bare catalysts and store in csv")
    parser.add_argument('-csv', help = 'doesitbind output csv', type=str, default = None)
    parser.add_argument('-cat', help='input file/directory (bare catalyst .out)', type=str, default = None)
    parser.add_argument('-O2', help='O2 bound directory', type=str, default = None)
    parser.add_argument('-plusone', help='Plus one data for bare catalyst to grab charge difference data', default = None, type=str)
    parser.add_argument('-mol', help='Mol files for bare catalyst to grab neighbor data', default = None, type=str)
    parser.add_argument('-xyz', help='xyz files for bare catalyst to grab bond length data', default = None, type=str)
    parser.add_argument('-f1', help = 'First feature to compare against, default to spin density', default = 'SpinDensity', type = str)
    parser.add_argument('-f2', help = "Second feature to compare against, like ChElPGNeutralCharge, Doesitbind, BindingEnergy, NeutralFreeEnergy, NeighborSpinDensity, NeighborChElPGCharge, NeighborChargeDifference", type = str)
    parser.add_argument('-plot', help = 'Choose classifier to plot', default = None, type = str)
    parser.add_argument('-f3', help = 'Third feature', default = None)
    parser.add_argument('-f4', help = 'Fourth feature', default = None)
    parser.add_argument('-f5', help = 'Fifth feature', default = None)
    parser.add_argument('-f6', help = 'Sixth feature', default = None)
    parser.add_argument('-f7', help = 'Seventh feature', default = None)
    parser.add_argument('-f8', help = 'Eighth feature', default = None)
    parser.add_argument('-f9', help = 'Ninth feature', default = None)
    parser.add_argument('-f10', help = 'Tenth feature', default = None)

    parser.add_argument('-c', help = 'Penalty for svc', default = 1)
    parser.add_argument('-poly', help = 'Enable polynomial basis', default = 0)
    parser.add_argument('-onehot', help = "Enable one-hot encoder", default = 0)
    parser.add_argument('-confusion', help = "Enable confusion matrix", default = 0)
    parser.add_argument('-plotcoef', help = "Enable coefficient plot", default = 0)
    parser.add_argument('-misclassified', help = "Enable misclassification analysis", default = 0)
    
    parser.add_argument('-test', help = 'csv for unclassified data', default = None)

    args = parser.parse_args()
    
    if (args.csv != None):
        alldata = getDataFromCSV(args.csv)
    else:
        alldata = getData(args.cat, args.O2, args.plusone, args.mol, args.xyz)
    
    classify(alldata, args.plot, args.f1, args.f2, args.f3, args.f4, args.f5, args.f6, args.f7, args.f8, args.f9, args.f10, args.c, args.poly, args.onehot, args.confusion, args.plotcoef, args.misclassified, args.test)