"""Script to id which motif in the motif file is the actual motif"""

from Bio import motifs
import sys

def motif_identifier(tf, meme_folder):
    meme_file = meme_folder + '/meme.txt'
    with open(meme_file, 'r') as fid:
        try:
            records = motifs.parse(fid, 'MEME')
        except ValueError:
            return '1'
    motif_list = []
    evalue_list = []
    for record in records:
        curr_evalue = record.evalue
        for curr_motif in record.instances:
            # if curr_motif.sequence_name == tf.lower():
            motif_list.append(curr_motif.motif_name)
            evalue_list.append(curr_evalue)
    motif_list = zip(motif_list, evalue_list)
    sorted_list = sorted(motif_list, key=lambda x: x[1])
    motif_name = sorted_list[0][0]
    return motif_name[-1]

if __name__ == '__main__':
    op_folder = sys.argv[1]
    print(motif_identifier(op_folder))
