#!/usr/bin/env python3
# @Author: Dileep Kishore <dileep>
# @Date:   2016-12-04T16:35:17-05:00
# @Filename: extract_features.py
# @Last modified by:   dileep
# @Last modified time: 2016-12-06T15:49:17-05:00

"""Script to extract features for classification"""

from subprocess import call
import os
import paths
import pandas as pd
from feature_extraction import intergenetic_spaces
from feature_extraction import motif_scores
from feature_extraction.Base_Comp import get_df_base_composition
from feature_extraction.tf_bs import get_df_distances_to_other_tfs
from feature_extraction.dna_topology import get_df_topologies_of_bindingsites

def clean_up(out_path):
    command = 'rm -rf ' + out_path
    call(command, shell=True)
    return None

def main(feature_data_dir, out_dir, labelled_data):
    # Intergenetic space feature extraction
    intergene_file = feature_data_dir + 'Intergenetic_Spaces.txt'
    intergene_data = intergenetic_spaces.extract_features(intergene_file, labelled_data)
    # Motif score feature extraction
    n_file = paths.motif_training_path + 'best_n.txt'
    fimo_dir = paths.motif_training_path
    motif_data = motif_scores.extract_features(n_file, fimo_dir, labelled_data)
    # Base composition feature extraction
    base_comp = get_df_base_composition(labelled_data, paths.genome_path)
    # Distance to other features
    tf_dist = get_df_distances_to_other_tfs(labelled_data, feature_data_dir)
    # Topologies of binding sites
    tf_topologies = get_df_topologies_of_bindingsites(labelled_data, feature_data_dir)
    # for col in tf_topologies.columns:
    #     if 'mgw' in col:
    #         if 'mgw_tops' in locals():
    #             mgw_tops = pd.concat([mgw_tops, tf_topologies[col]], axis=1)
    #         else:
    #             mgw_tops = tf_topologies[col]
    # final_dataframe = pd.concat([intergene_data, base_comp, tf_dist, mgw_tops], axis=1)
    # Joining all the data_frames
    final_dataframe = pd.concat([intergene_data, motif_data, base_comp, tf_dist, tf_topologies], axis=1)
    op_file = out_dir + 'final_dataframe.csv'
    final_dataframe.to_csv(op_file, index=False)
    return intergene_data, base_comp, tf_topologies

if __name__ == '__main__':
    ans = input('Do you want to clean the results directory? (Y/N)\n')
    out_dir = paths.out_path + 'feature_extraction/'
    if ans == 'Y':
        clean_up(out_dir)
    feature_data_dir = paths.features_data_path
    if not os.path.exists(out_dir):
        call('mkdir -p ' + out_dir, shell=True)
    labelled_data = paths.labelled_data
    a, b, c = main(feature_data_dir, out_dir, labelled_data)
