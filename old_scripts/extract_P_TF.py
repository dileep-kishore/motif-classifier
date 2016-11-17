"""Script to extract promoter seqs and TF_promoter seqs. Input to meme"""

import sys
import os
import subprocess
from Bio import SeqIO
from Bio import pairwise2
import pandas as pd

def get_promoter_locs(promoter_file, f_type='text'):
    """Read promoter locations from the file
       :params promoter_file: File contains locs and coverage info
    """
    promoter_locs = []
    if f_type == 'text':
        with open(promoter_file, 'r') as fid:
            header = 1
            for line in fid.readlines():
                if header:
                    header = 0
                    continue
                line_info = line.strip().split()
                promoter_locs.append(int(line_info[0]))
    elif f_type == 'excel':
        xl_file = pd.read_excel(promoter_file)
        promoter_locs = list(xl_file['Position'])
    elif f_type == 'csv':
        csv_file = pd.read_csv(promoter_file)
        promoter_locs = list(csv_file['Locations'])
    return promoter_locs

def check_if_TF(feature):
    """Check whether the given CDS is a promoter or note
       :param feature: Feature set of a given CDS
       :return Whether the given CDS is a TF or not
    """
    TF_criteria = ['regulator', 'transcription factor', 'DNA binding', 'two-component', 'repressor', 'activator']
    TF_anti_criteria = 'phage'
    flag = False
    try:
        function_info = feature.qualifiers['function']
    except:
        function_info = []
    try:
        product_info = feature.qualifiers['product']
    except:
        product_info = []
    try:
        note_info = feature.qualifiers['note']
    except:
        note_info  = []
    info_dict = {'function': function_info, 'product': product_info, 'note': note_info}
    for qualifier in info_dict:
        for criteria in TF_criteria:
            for list_item in info_dict[qualifier]:
                if TF_anti_criteria in list_item.lower():
                    return False
                if criteria in list_item.lower():
                    flag = True
    if flag:
        return True
    else:
        return False

def get_TF_locs(gbk_file, GO_term):
    """Read TF locations from file using the GO_term
       :params gbk_file: Genbank file that contains info about all the genes
       :params GO_term: Gene Ontology term for TF
    """
    genome = SeqIO.read(gbk_file, "genbank")
    TF_data = []
    for ind, feature in enumerate(genome.features):
        if feature.type == 'CDS':
            if check_if_TF(feature):
                start = feature.location.start.position
                end = feature.location.end.position
                strand = feature.strand
                gene = feature.qualifiers['gene']
                TF_data.append([gene[0], start, end, strand])
    return TF_data

def get_complement(seq):
    """Generate the reverse complement of the given sequence
       :params seq: Nucleotide sequence to be complemented
    """
    rev_dict = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
    seq = seq[-1::-1]
    rev_comp = [rev_dict[nucleotide] for nucleotide in seq]
    return ''.join(rev_comp)

def get_promoter_seqs(seq_file, promoter_locs, promoter_size, complement=0):
    """Read promoter sequences from sequence file using promoter locations
       :params promoter_locs: Promoter coordinates obtained through ChIP-seq
       :params promoter_size: Size of the promoter on either side of the coord
    """
    with open(seq_file, 'r') as fid:
        seq = [line.strip() for line in fid.readlines()]
        sequence = ''.join(seq[1:])
        promoters = []
        for loc in promoter_locs:
            start = loc - promoter_size//2
            end = loc + promoter_size//2
            if complement:
                promoters.append(get_complement(sequence[start:end+1]))
            else:
                promoters.append(sequence[start:end+1])
    return promoters

def get_TF_pseqs(seq_file, TF_data, TF_psize):
    """Read TF promoter sequences from sequence file using locs in TF_data
       :params TF_data: Contains gene_name, start, end, strand_direction
       :params TF_psize: The size of the TF promoter to be used
    """
    with open(seq_file, 'r') as fid:
        seq = [line.strip() for line in fid.readlines()]
        sequence = ''.join(seq[1:])
        TF_ps = []
        TF_names = []
        for TF in TF_data:
            TF_names.append(TF[0])
            if TF[3] == -1:
                start = TF[2]
                end = start + TF_psize
                TF_ps.append(get_complement(sequence[start:end+1]))
            elif TF[3] == 1:
                end = TF[1]
                start = TF[1] - TF_psize
                TF_ps.append(sequence[start:end+1])
            else:
                print(TF)
                sys.exit("Error")
    return TF_names, TF_ps

def write2fasta(folname, promoters, tf_names, tf_ps):
    """Write promoter sequences and TF promoter sequences in fasta files for meme i/p
       :params promoters: Promoter sequences to be input to meme
       :params tf_ps: Promoter sequences of TFs to be input to meme
       :params TF_names: Gene names of the TFs
    """
    folder_name = 'meme_in/' + folname + '/'
    subprocess.call(["mkdir", folder_name])
    for tf_ind, ps in enumerate(tf_ps):
        fname = folder_name + tf_names[tf_ind] + '.fa'
        with open(fname, 'w') as fid:
            for i, promoter in enumerate(promoters):
                # if pairwise2.align.localxx(promoter, ps)[0][2]/len(ps) > 0.75:
                #     continue
                header = '>promoter' + str(i) + '\n'
                fid.write(header)
                fid.write(promoter)
                fid.write('\n')
            header = '>TF_' + tf_names[tf_ind] + '\n'
            fid.write(header)
            fid.write(ps)
            fid.write('\n')
    return None

def write_TF_list(fname, TF_data):
    with open(fname, 'w') as fid:
        fid.write('Transcription factor,Start,End,Strand')
        fid.write('\n')
        for tf_info in TF_data:
            fid.write(tf_info[0])
            fid.write(',')
            fid.write(str(tf_info[1]))
            fid.write(',')
            fid.write(str(tf_info[2]))
            fid.write(',')
            fid.write(str(tf_info[3]))
            fid.write('\n')
    return None

if __name__ == '__main__':
    # Get locations of promoters
    if len(sys.argv) < 2:
        print("Run as extract_P_TF.py <text/excel/csv>")
    f_type = sys.argv[1]
    os.chdir('../')
    # qseq_path = 'original_data/nac_query/'
    qseq_path = 'rnaseq_query/'
    qseq_files = os.listdir(qseq_path)
    qseq_files.sort()
    answer = input("Do you want to delete input folder Y/N? ")
    if answer=='Y':
        try:
            subprocess.call(["rm", "-rf", "meme_in"])
            subprocess.call(["mkdir", "meme_in"])
        except:
            pass
    for ind, curr_qseq in enumerate(qseq_files):
        promoter_file = qseq_path + curr_qseq
        promoter_locs = get_promoter_locs(promoter_file, f_type)
        # Get locations of TFs
        gbk_file = 'original_data/U00096.2.gbk'
        GO_term = 'GO:0016563'
        TF_data = get_TF_locs(gbk_file, GO_term)
        # Get sequences
        seq_file = 'original_data/U00096.2.fa'
        # Get promoter sequences
        promoter_size = 200
        promoter_seqs = get_promoter_seqs(seq_file, promoter_locs, promoter_size)
        rev_promoter_seqs = get_promoter_seqs(seq_file, promoter_locs, promoter_size, 1)
        # Get TF promoter sequences
        TF_psize = 200
        TF_names, TF_pseqs = get_TF_pseqs(seq_file, TF_data, TF_psize)
        # Write files for input into meme
        if answer == 'Y':
            write2fasta('forward'+str(ind), promoter_seqs, TF_names, TF_pseqs)
            # write2fasta('reverse', rev_promoter_seqs, TF_names, TF_pseqs)
        tf_file = 'meme_in/TF_list'+str(ind)+'.csv'
        # write_TF_list(tf_file, TF_data)
