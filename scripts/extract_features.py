#!/usr/bin/env python3
"""Script to extract features for classification"""
import paths
from feature_extraction.intergenetic_spaces

def clean_up(out_path):
    results_dir = out_path + 'feature_extraction'
    command = 'rm -rf ' + results_dir
    call(command, shell=True)
    return None

def main(in_dir, out_dir):
    # Intergenetic space feature extraction
    intergene_file = in_dir + 'Intergenetic_Spaces.txt'
    intergenetic_spaces.extract_features(intergene_file, labelled_data)

if __name__ == '__main__':
    ans = input('Do you want to clean the results directory? (Y/N)\n')
    if ans == 'Y':
        clean_up(paths.out_path)
    in_dir = paths.features_data_path
    out_dir = paths.out_path + 'feature_extraction/'
    labelled_data = paths.labelled_data
    main(in_dir, out_dir, labelled_data)
