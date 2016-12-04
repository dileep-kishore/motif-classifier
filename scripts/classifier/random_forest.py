#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Dec  3 19:48:21 2016

@author: Manu
"""

# Random Forest Classification
import pandas as pd
from sklearn.cross_val_score import cross_val_score
from sklearn.ensemble import RandomForestClassifier
import math

def run_classifier(feature_df, labels_df):
#    url = "https://archive.ics.uci.edu/ml/machine-learning-databases/pima-indians-diabetes/pima-indians-diabetes.data"
#    names = ['preg', 'plas', 'pres', 'skin', 'test', 'mass', 'pedi', 'age', 'class']
#    dataframe = pandas.read_csv(url, names=names)
#    array = dataframe.values
    
#    num_folds = 10
#    num_instances = len(feature_df)
#    seed = 7
    num_trees = 100
    max_features = int(math.sqrt(feature_df.shape[1]))
        
#    kfold = cross_validation.KFold(n=num_instances, n_folds=num_folds, random_state=seed)
    model = RandomForestClassifier(n_estimators=num_trees, max_features=max_features)
#    results = cross_validation.cross_val_score(model, X, Y, cv=kfold)
#    print(results.mean())
    
    scores = cross_val_score(model, feature_df , labels_df, cv=5, scoring='f1')
    print(scores)
#MAX FEATURES parameter
# max_features=sqrt(n_features) for classification tasks (where n_features is the number of features in the data). 