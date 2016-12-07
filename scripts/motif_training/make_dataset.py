# @Author: Dileep Kishore <dileep>
# @Date:   2016-12-03T03:28:55-05:00
# @Filename: make_dataset.py
# @Last modified by:   dileep
# @Last modified time: 2016-12-06T22:06:03-05:00



#!/usr/bin/env python3
"""Make the dataset of tf binding sites using fitness data"""

import pandas as pd
from Bio import motifs

def get_chipseq_ranges(chip_fimo, chip_all, motif_file):
    chip_all_data = pd.read_csv(chip_all)
    chip_fimo_data = pd.read_table(chip_fimo)
    motif_id = set(chip_fimo_data['#pattern name']).pop()
    with open(motif_file, 'r') as fid:
        motif_record = motifs.parse(fid, 'MEME')
    mot_len = int(motif_record[motif_id-1].length)
    mot_range = []
    for pos in list(chip_all_data['Position']):
        mot_range.append((pos-mot_len-10, pos+mot_len+10))
    return mot_range

def label_genome_data(genome_fimo, motif_range):
    genome_fimo_data = pd.read_table(genome_fimo)
    mot_start = list(genome_fimo_data["start"])
    mot_stop = list(genome_fimo_data["stop"])
    label = ['' for _ in range(len(mot_stop))]
    for i in range(len(mot_start)):
        if any(start <= mot_start[i] <= stop for (start, stop) in motif_range):
            label[i] = 'true'
            continue
        if any(start <= mot_stop[i] <= stop for (start, stop) in motif_range):
            label[i] = 'true'
            continue
        label[i] = 'false'
    return label

def make_dataframe(results_dir, genome_fimo, labels):
    genome_fimo_data = pd.read_table(genome_fimo)
    mot_start = list(genome_fimo_data["start"])
    mot_stop = list(genome_fimo_data["stop"])
    d = {'start': mot_start, 'stop': mot_stop, 'label': labels}
    label_df = pd.DataFrame(data=d, index=None)
    cols = ['start', 'stop', 'label']
    label_df = label_df[cols]
    out_file = results_dir + 'labelled_data.csv'
    label_df.to_csv(out_file, index=False)

def make_dataset(n, results_dir):
    curr_dir = results_dir + str(n) + '/'
    fimo_dir = curr_dir + 'fimo/'
    chip_all = curr_dir + 'chip_out_all.csv'
    motif_file = curr_dir + 'meme/meme.txt'
    chip_fimo = fimo_dir + 'chip/fimo.txt'
    genome_fimo = fimo_dir + 'genome/fimo.txt'
    motif_range = get_chipseq_ranges(chip_fimo, chip_all, motif_file)
    labels = label_genome_data(genome_fimo, motif_range)
    make_dataframe(results_dir, genome_fimo, labels)
