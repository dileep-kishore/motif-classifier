#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Dec  3 19:48:21 2016

@author: Manu
"""

# Random Forest Classification
import pandas as pd
from sklearn.model_selection import cross_val_score
#from sklearn.model_selection import cross_val_predict
from sklearn.ensemble import RandomForestClassifier
from sklearn.dummy import DummyClassifier
#from sklearn.model_selection import ShuffleSplit
from sklearn.model_selection import KFold
import math
import numpy as np
#from sklearn.metrics import classification_report
from sklearn.metrics import matthews_corrcoef
from sklearn.metrics import make_scorer
from sklearn.metrics import confusion_matrix
import matplotlib.pyplot as plt


def prepare_data_for_classifier(feature_path, labelled_path, randomize = False, only_columns = None):
    # We read the files
    features_data = pd.read_csv(feature_path)
    labelled_data = pd.read_csv(labelled_path)
    # We save the features names
    headers = features_data.columns

    # Filter unwanted features
    if only_columns is not None:
        features_data = features_data[only_columns]
    
    #Randomize dataset    
    if randomize :
        rd_index = np.random.permutation(features_data.index)
        features_data = features_data.reindex(rd_index)
        labelled_data = labelled_data.reindex(rd_index)

    # Scikit uses arrays, so we convert the dataframes 
    labels_array = np.array(labelled_data['label'], dtype=int)
    feature_array = features_data.as_matrix()
    return (feature_array, labels_array, headers)    
    
def cross_validation(the_feature_df, the_labels_df):
    # Build a Dummy estimator (classifier)
    dummy_type = 'most_frequent'
#    dummy_type = 'stratified'
    dummy_model = DummyClassifier(constant=None, random_state=0, strategy=dummy_type)
    
    #Build a Random Forest Classifier
    num_trees = 100
    # Number of features: use sqrt, except in case of high-dimensional dataset
    number_of_features = 'sqrt'
#    number_of_features = 'log2'
    forest_model = RandomForestClassifier(n_estimators=num_trees, max_features=number_of_features)
    
    # We create an scorer from the metric "Mathew correlation"
    # The Matthews correlation coefficient (+1 represents a perfect
    #    prediction, 0 an average random prediction and -1 and inverse
    #    prediction).
    #
    # More info https://en.wikipedia.org/wiki/Matthews_correlation_coefficient
    # Sci-kit documentation -> http://scikit-learn.org/stable/modules/generated/sklearn.metrics.matthews_corrcoef.html#sklearn.metrics.matthews_corrcoef
    # Paper that suggest use (according to Wikipedia)
    # http://www.flinders.edu.au/science_engineering/fms/School-CSEM/publications/tech_reps-research_artfcts/TRRA_2007.pdf
    # Didactic explanation of MCC
    # https://lettier.github.io/posts/2016-08-05-matthews-correlation-coefficient.html
    matthews_scorer = make_scorer(matthews_corrcoef, greater_is_better=True, needs_proba=False, needs_threshold=False)

    # We define the cross_validators
    cross_validator_forest = KFold(n_splits=10, shuffle=True, random_state=None)
    cross_validator_dummy = KFold(n_splits=10, shuffle=True, random_state=None)
    # Note: StratifiedKFold could be used also
    #
    # Note: ShuffleSplit could be used, BUT the documentation says:
    #   Note: contrary to other cross-validation strategies, random splits
    #   do not guarantee that all folds will be different, although this is
    #   still very likely for sizeable datasets.
    #
    # Eg: ShuffleSplit(n_splits=10, test_size=0.1, train_size=None, random_state=None)
    #
    # Other cross-validators, check:
    # http://scikit-learn.org/stable/modules/classes.html#module-sklearn.model_selection
    
    # Cross validators will train and test models, returning their score
    score_dummy = cross_val_score(dummy_model, the_feature_df, the_labels_df, cv=cross_validator_dummy, scoring = matthews_scorer)
    score_forest = cross_val_score(forest_model, the_feature_df, the_labels_df, cv=cross_validator_forest, scoring = matthews_scorer)
    
    return (score_forest, score_dummy)
    
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
<<<<<<< HEAD
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
=======
>>>>>>> f05746cd2d0bc30be3396fbfad1ee86b4788b4b9
    # return scores, model
    
    # MANU ADDED THIS
    run_once(feature_path, labelled_path, features_to_use=['tfs_D_fw', 'tfs_D_rv', 'tfs_U_fw', 'tfs_U_fw.1'])
    print(run_cross_validation(feature_path, labelled_path, features_to_use=['tfs_D_fw', 'tfs_D_rv', 'tfs_U_fw', 'tfs_U_fw.1']))
    
    return None, None

<<<<<<< HEAD
=======
def run_once(feature_path,labelled_path,features_to_use = None):
#    features_to_use = ['intergenetic', 'compA_d', 'compT_d', 'compG_d', 'compC_d', 'compA_u', 'compT_u', 'compG_u', 'compC_u', 'tfs_D_fw', 'tfs_D_rv', 'tfs_U_fw', 'tfs_U_fw.1']    
#    features_to_use = ['tfs_D_fw', 'tfs_D_rv', 'tfs_U_fw', 'tfs_U_fw.1']
    features, labels, headers = prepare_data_for_classifier(feature_path,labelled_path, randomize = True, only_columns=features_to_use)

    # We manually split the randomized dataset into Training and Testing sets
    train_features = features[:500,:]
    train_labels = labels[:500]
    test_features = features[500:,:]
    test_labels = labels[500:]
    
    # We create some estimators and we trained them
    forest = RandomForestClassifier(n_estimators=100, max_features='sqrt')
    forest = forest.fit(train_features, train_labels)
    dummy = DummyClassifier(constant=None, random_state=0, strategy='most_frequent')
    dummy = dummy.fit(train_features,train_labels)
    
    # We use them to clasify 
    predicted_labels_forest = forest.predict(test_features)
    predicted_labels_dummy = dummy.predict(test_features)

    # We get the MCC scores (1 is perfect classification, 0 is random, -1 is inverse prediction)
    forest_score = matthews_corrcoef(test_labels, predicted_labels_forest, sample_weight=None)
    dummy_score = matthews_corrcoef(test_labels, predicted_labels_dummy, sample_weight=None)
    
    # We generate the Confusion Matrix 
    # True negatives is C_{0,0}
    # False negatives is C_{1,0}
    # True positives is C_{1,1}
    # False positives is C_{0,1}
    forest_matrix = confusion_matrix(test_labels, predicted_labels_forest, labels=None, sample_weight=None)
    dummy_matrix = confusion_matrix(test_labels, predicted_labels_dummy, labels=None, sample_weight=None)
    
    # We print everything
    print("## RUN ONCE ##\n")    
    print("Forest score",forest_score)
    print("Dummy score",dummy_score)
    
    print("\nForest Matrix")
    print(forest_matrix)
    print("\nDummy Matrix")
    print(dummy_matrix)
    
    # And we print the Feature Importance (with plot)
    forest_feature_importance(features, headers, forest, True)    

#    forest_s , dummy_s = cross_validation(feature_df, labels_df)
#    print(forest_s , "\n", dummy_s)
    
#    predictions = cross_val_predict(model, feature_df, labels_df, cv=KFold(n_splits=10, shuffle=True, random_state=None))
#    np.set_printoptions(threshold=np.nan)
#    print(predictions)
#    print(np.sum(predictions))
#    print(len(predictions))
    ##
    
def forest_feature_importance(features_array, features_names, the_trained_forest, plot=False):
    importances = the_trained_forest.feature_importances_
    std = np.std([tree.feature_importances_ for tree in the_trained_forest.estimators_],
             axis=0)
    indices = np.argsort(importances)[::-1]

    # Print the feature ranking
    print("\nFeature ranking:")
    
    for f in range(features_array.shape[1]):
        print("%d. feature %s (%f)" % (f + 1, features_names[indices[f]], importances[indices[f]]))

    features_sorted_by_importance = [ features_names[i] for i in indices]
        
    # Plot the feature importances of the forest
    if plot:
        plt.figure()
        plt.title("Feature importances")
        plt.bar(range(features_array.shape[1]), importances[indices],
               color="r", yerr=std[indices], align="center")
        plt.xticks(range(features_array.shape[1]), tuple(features_sorted_by_importance))
        plt.xlim([-1, features_array.shape[1]])
        plt.show()

def summurize_cross_score(scores_array):
    best = scores_array.max()
    worst = scores_array.min()
    mean = np.mean(scores_array)
    assert (scores_array < 0).sum > 0, "Danger: we have a negative MCC score. Our way of summurizing MCC scores can not handle this"
    return (best,worst,mean)
    
def run_cross_validation(feature_path, labelled_path, features_to_use = None):
    # We process the features and labels so as to Scikit like it
    features_array, labels_array, features_names = prepare_data_for_classifier(feature_path, labelled_path, randomize=True, only_columns = features_to_use)
    # We run the cross-validation
    forest_scores , dummy_scores = cross_validation(features_array, labels_array)
    # We return a summury of the MCC scores
    print("\n ## Running CROSS validation ##")
    return summurize_cross_score(forest_scores)
>>>>>>> f05746cd2d0bc30be3396fbfad1ee86b4788b4b9
