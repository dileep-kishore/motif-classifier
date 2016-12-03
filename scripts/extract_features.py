#!/usr/bin/env python3
"""Script to extract features for classification"""
import paths
from feature_extraction.intergenetic_spaces import intergenetic_space

def clean_up(out_path):
    results_dir = out_path + 'feature_extraction'
    command = 'rm -rf ' + results_dir
    call(command, shell=True)
    return None

def main():
    features_data_path = paths.features_data_path
    intergenetic_space(ip_file, features_data_path)

if __name__ == '__main__':
    ans = input('Do you want to clean the results directory? (Y/N)\n')
    if ans == 'Y':
        clean_up(paths.out_path)
    main()
