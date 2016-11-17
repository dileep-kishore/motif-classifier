"""Script to find motifs in rna_seq data"""

import os
import sys
from subprocess import call
from exp_corr import read_expdata, parse_expdata
from Bio import SeqIO
import csv

def select_genes(exp_dict, fold_change, num_genes):
    """Obtain $num_genes genes with highest fold_change
        :param fold_change: increase/decrease/no_Change
        :param exp_dict: Dict with genes under respective fold_change types
        :param num_genes: Number of genes to be selected
        :return selected_genes: List of top genes with corr fold_change
    """
    if fold_change == 'both':
        genes = exp_dict['increase'] + exp_dict['decrease']
    else:
        genes = exp_dict[fold_change]
    # Sort the data
    genes = sorted(genes, key=lambda x: abs(x[-1]), reverse=True)
    #FIXME: Check for tf in list?
    if num_genes == 'all':
        selected_genes = [gene[1] for gene in genes]
        selected_ids = [gene[0] for gene in genes]
    else:
        selected_genes = [gene[1] for gene in genes[:num_genes]]
        selected_ids = [gene[0] for gene in genes[:num_genes]]
    return selected_genes, selected_ids

def check_if_gene(feature, gene_ids):
    curr_id = feature.qualifiers['db_xref'][-1].split(':')[-1]
    return bool(curr_id in gene_ids)

def get_gene_locs(gbk_file, genes, gene_ids, prom_length):
    genome = SeqIO.read(gbk_file, "genbank")
    gene_data = dict()
    for ind, feature in enumerate(genome.features):
        # NOTE:Some of the most regulated genes code for tRNA
        if feature.type == 'CDS':
            if check_if_gene(feature, gene_ids):
                gene = feature.qualifiers['gene'][0]
                start = feature.location.start.position
                end = feature.location.end.position
                strand = feature.strand
                if strand == 1:
                    pos = start - prom_length//2
                elif strand == -1:
                    pos = end + prom_length//2
                if pos <= 100:
                    pos = 101
                gene_data[gene] = str(pos)
                # gene_data.append([str(pos), gene, 'intergenic'])
    return gene_data

def write_query_file(fname, gene_data):
    # locations, coverage, type
    header = 'Locations,Gene,Geneid\n'
    with open(fname, 'w') as fid:
        fid.write(header)
        for gene in gene_data:
            fid.write(','.join(gene))
            fid.write('\n')
    return None

def get_expdata(fname):
    exp_folder = 'new_data/testbindingdata/'
    # Expression file needs to be changed for different TFs
    exp_file = exp_folder + fname
    exp_sheet = 2
    exp_data = read_expdata(exp_file, exp_sheet)
    exp_dict, not_present = parse_expdata(exp_data, exp_data['Symbol'])
    return exp_dict

def get_codinggenes(genes, gene_ids, gene_dict, n):
    promoter_data = []
    count = 0
    for i, gene in enumerate(genes):
        if gene in gene_dict.keys():
            promoter_data.append([gene_dict[gene], gene, gene_ids[i]])
            count += 1
        if count >= n:
            break
    return promoter_data

def main(n, gene_dict, genes, gene_ids):
    gbk_file = 'original_data/U00096.2.gbk'
    prom_length = 200
    query_file = 'rnaseq_query/'+str(n)+'.csv'
    exp_dict = get_expdata('Nac.xlsx')
    # genes, gene_ids = select_genes(exp_dict, 'both', 'all')
    promoter_data = get_codinggenes(genes, gene_ids, gene_dict, n)
    prom_file = 'rnaseq_prom/promoter_seqs_'+str(n)+'.csv'
    pfid = open(prom_file, 'w')
    wr = csv.writer(pfid)
    wr.writerows(list(zip(genes, gene_ids)))
    pfid.close()
    write_query_file(query_file, promoter_data)

if __name__ == '__main__':
    os.chdir('../')
    try:
        call('rm -rf rnaseq_query rnaseq_prom', shell=True)
        call('mkdir rnaseq_query rnaseq_prom', shell=True)
    except:
        pass
    gbk_file = 'original_data/U00096.2.gbk'
    prom_length = 200
    exp_dict = get_expdata('Nac.xlsx')
    genes, gene_ids = select_genes(exp_dict, 'both', 'all')
    gene_file = 'rnaseq_prom/genes.csv'
    gene_dict = get_gene_locs(gbk_file, genes, gene_ids, prom_length)
    for n_seq in [10, 25, 50, 100]:
        main(n_seq, gene_dict, genes, gene_ids)
