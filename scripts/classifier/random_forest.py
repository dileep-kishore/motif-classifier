#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Dec  3 19:48:21 2016

@author: Manu
"""

# Random Forest Classification
import pandas as pd
from sklearn.model_selection import cross_val_score
from sklearn.ensemble import RandomForestClassifier
from sklearn.dummy import DummyClassifier
from sklearn.model_selection import ShuffleSplit
import math
import numpy as np

def run_classifier(feature_path, labelled_path):
    training_data = pd.read_csv(feature_path)
    labelled_data = pd.read_csv(labelled_path)
    feature_df = training_data
    headers = feature_df.columns
    labels_df = np.array(labelled_data['label'], dtype=int)
    feature_df = feature_df.as_matrix()
    num_trees = 100
    max_features = int(math.sqrt(feature_df.shape[1]))
    model = RandomForestClassifier(n_estimators=num_trees, max_features=max_features)
    dummy_model = DummyClassifier(constant=None, random_state=0, strategy='most_frequent')
    scores = cross_val_score(model, feature_df , labels_df, cv=ShuffleSplit(n_splits=3, train_size=0.7, random_state=0))
    # scores = cross_val_score(model, feature_df , labels_df, cv=10, verbose=1)
    dummy_scores = cross_val_score(dummy_model, feature_df, labels_df, cv=ShuffleSplit(n_splits=3, train_size=0.7, random_state=0))
    print('randomforest_scores=',np.mean(scores), np.std(scores))
    print('dummy_scores=',np.mean(dummy_scores), np.std(dummy_scores))
    np.random.shuffle(feature_df)
    dummy_model = dummy_model.fit(feature_df[:400,:], labels_df[:400])
    model = model.fit(feature_df[:400,:], labels_df[:400])
    # print(model.predict(feature_df[:100,:]))
    print('randomforest',model.score(feature_df[400:,:], labels_df[400:]))
    print('dummy', dummy_model.score(feature_df[400:,:], labels_df[400:]))
    # return scores, model
    return None, None

#MAX FEATURES parameter
# max_features=sqrt(n_features) for classification tasks (where n_features is the number of features in the data). 
