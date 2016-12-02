#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  2 16:53:04 2016

@author: Manu
"""

import pandas as pd

# FIX THIS! It should be read from the package
features_data_path = "../../data/features/"

# The name of the file with the RegulonDB data of TFs and their binding site
regulondb_tfbs_filename = features_data_path + "regulondb-tfbs.txt"

column_names = ['position', 'value']