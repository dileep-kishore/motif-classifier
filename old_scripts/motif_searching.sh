#!/bin/bash

# Getting correct input and output parameters
if [[ $# < 2 ]]; then
    echo "Run as motif_searching.sh <input=text/excel> <output=text/html>"
fi

# Generating the promoter location list from the rna-seq file
echo "Generating promoter location list from rna-seq data..."
python3 rnaseq_motifs.py

# Extracting promoters and transcription factors
echo "Extracting promoter and TF promoter sequences..."
python3 extract_P_TF.py $1

# Use meme to obtain common motifs
echo "Running meme..."
python3 get_meme_output.py $2

# Generate summary of the motifs obtained from meme
echo "Analyzing meme output..."
python3 meme_output_analysis.py $2

# Generating a html file with clustered motifs
echo "Generating html file..."
python3 write_html.py

# Run meme on the generated fasta files
## Each fasta file has one TF promoter and the rest are promoters given in the files
