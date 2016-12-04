#!/usr/bin/env python3
"""Script to extract features for classification"""

from subprocess import call
import paths
from feature_extraction import intergenetic_spaces
from feature_extraction.Base_Comp import get_df_base_composition
from feature_extraction.tf_bs import get_df_distances_to_other_tfs
from feature_extraction.dna_topology import get_df_topologies_of_bindingsites

def clean_up(out_path):
    results_dir = out_path + 'feature_extraction/'
    command = 'rm -rf ' + results_dir
    call(command, shell=True)
    return None

def main(feature_data_dir, out_dir, labelled_data):
    # Intergenetic space feature extraction
    intergene_file = feature_data_dir + 'Intergenetic_Spaces.txt'
    intergene_data = intergenetic_spaces.extract_features(intergene_file, labelled_data)
    # Base composition feature extraction
    base_comp = get_df_base_composition(labelled_data, paths.genome_path)
    # Distance to other features
    tf_dist = get_df_distances_to_other_tfs(labelled_data, feature_data_dir)
    # Topologies of binding sites
    tf_topologies = get_df_topologies_of_bindingsites(labelled_data, feature_data_dir)
    # Joining all the data_frames
    final_dataframe = intergene_data.join(base_comp, tf_dist, tf_topologies)
    op_file = out_dir + 'final_dataframe.csv'
    final_dataframe.to_csv(op_file)
    return final_dataframe

if __name__ == '__main__':
    ans = input('Do you want to clean the results directory? (Y/N)\n')
    if ans == 'Y':
        clean_up(paths.out_path)
    feature_data_dir = paths.features_data_path
    out_dir = paths.out_path + 'feature_extraction/'
    labelled_data = paths.labelled_data
    main(feature_data_dir, out_dir, labelled_data)
