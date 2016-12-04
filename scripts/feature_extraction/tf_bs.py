#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  2 01:59:51 2016

@author: Manuel Gimenez
"""

#import paths
import pandas as pd
from functools import partial

def get_distance_to_other_tfs(the_data, a_series):
    start = a_series[0]
    end = a_series[1]
    position = (start + end)/2
    mask = (the_data['TFBS_center_position']-position).abs().argsort()
    tfs_around_position = the_data.iloc[mask]
    tfs_upstream_position = tfs_around_position[tfs_around_position['TFBS_center_position'] >= position]
    tfs_downstream_position = tfs_around_position[tfs_around_position['TFBS_center_position'] < position]
    
    tfs_upstream_position_fw = tfs_upstream_position[tfs_upstream_position['DNA_strand']=="forward"]
    tfs_upstream_position_rv = tfs_upstream_position[tfs_upstream_position['DNA_strand']=="reverse"]
    
    tfs_downstream_position_fw = tfs_downstream_position[tfs_downstream_position['DNA_strand']=="forward"]
    tfs_downstream_position_rv = tfs_downstream_position[tfs_downstream_position['DNA_strand']=="reverse"]
        
    #Some consistency checks
    assert len(tfs_downstream_position_fw.index) > 0, "DANGER: No downstream-forward other-TFs binding site!"
    assert len(tfs_downstream_position_rv.index) > 0, "DANGER: No downstream-reverse other-TFs binding site!"

    ## Manual MANU debug
    #    print("The start is: "+str(start)+"\n")
#    print("The end is: "+str(end)+"\n")

    #    if start == 4636959:
#        print(tfs_upstream_position_rv)
#        print(tfs_upstream_position_rv.shape[0])
#        assert False
    
    ## VERY UGLY FIX in case there are NO tf upstream
    if not len(tfs_upstream_position_fw.index) > 0 :
        tfs_upstream_position_fw = the_data
    if not len(tfs_upstream_position_rv.index) > 0 :
        tfs_upstream_position_rv = the_data
#    assert len(tfs_upstream_position_fw.index) > 0, "DANGER: No upstream-forward other-TFs binding site!"
#    assert len(tfs_upstream_position_rv.index) > 0, "DANGER: No upstream-reverse other-TFs binding site!"
    
    the_list = [tfs_downstream_position_fw, tfs_downstream_position_rv, tfs_upstream_position_fw, tfs_upstream_position_rv]
    new_list = []
    for df in the_list:
        df = df[:1]['TFBS_center_position']-position
    #    print(df)
        new_list.append(df)
    #tfs_downstream_position_fw[:1]['TFBS_center_position']-position
    res = pd.concat(new_list, axis=0)
    res = res.abs()
#    res = data_filtered.iloc[:,1]
    return pd.Series(res.tolist(), index=['tfs_D_fw', 'tfs_D_rv', 'tfs_U_fw', 'tfs_U_fw'])

    
def get_df_distances_to_other_tfs(labelled_data, feature_path):
    # FIX THIS! It should be read from the package
    features_data_path = feature_path
    positions_df = pd.read_csv(labelled_data,comment='#')
    
    # The name of the file with the RegulonDB data of TFs and their binding site
    regulondb_tfbs_filename = features_data_path + "regulondb-tfbs.txt"
    
    # Im not using this at the moment
#    ecoli_genome_highest_position = 4700000.0
    
    # The name of the columns of the data (check the datafile for more detail)
    column_names =['TF_id', 'TF_name', 'TFBS_id', 'TFBS_start', 'TFBS_end', 'DNA_strand', 'interaction_id', 'unit_regulated','effect','p_name', 'relative_center', 'TFBS_sequence', 'evidence', 'evidence_level']
    
    # We read the data file into a Pandas Dataframe
    regulondb_df = pd.read_csv(regulondb_tfbs_filename, sep='\t', comment='#', skip_blank_lines=True, names=column_names)
    regulondb_df = regulondb_df[regulondb_df['evidence_level'] == "Strong"]
    #regulondb_df = regulondb_df.sort_values(by='TFBS_start', ascending=True)
    regulondb_df['TFBS_center_position'] = (regulondb_df['TFBS_start'] + regulondb_df['TFBS_end'])/2
    
    distances_getter = partial(get_distance_to_other_tfs, regulondb_df)
    res = positions_df.apply(distances_getter, axis=1)
    return res

# df = pd.DataFrame({'x' : [5, 32, 50, 64], 'y' : [13, 40, 58, 72]})
# temp = get_df_distances_to_other_tfs(df)

#path_to_data = "../../results/motif_training/labelled_data.csv"
##input_data = pd.read_csv("../../results/motif_training/labelled_data.csv")
#path_to_features = "../../data/features/"
##
#temp = get_df_distances_to_other_tfs(path_to_data, path_to_features)