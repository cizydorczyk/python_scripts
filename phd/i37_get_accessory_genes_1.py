import argparse
import pandas as pd
from Bio import SeqIO
import os.path

parser = argparse.ArgumentParser()

parser.add_argument("--input_roary", help="input roary csv file")
parser.add_argument("--st_list", help="st list to work with")
# parser.add_argument("--rast_fna_dir", help="directory with rast dna features files")
parser.add_argument("--output_dir", help="directory where to output individual isolate accessory gene fastas")
parser.add_argument("--group_alignments_dir", help="panaroo folder with group alignments")
parser.add_argument("--st", help="sequence type/identifier")

args = parser.parse_args()

with open(args.st_list, 'r') as infile1:
    st_list = [line.strip()+"_formatted" for line in infile1]

roary_data = pd.read_csv(filepath_or_buffer = args.input_roary, sep=",", index_col = 0)

roary_data = roary_data[st_list]

roary_data = roary_data[roary_data.isnull().any(axis=1)]
roary_data = roary_data[roary_data.isnull().sum(axis=1) < len(st_list)]

roary_data_filepath = args.output_dir + args.st + "_accessory_presence_absence.csv"
roary_data.to_csv(roary_data_filepath, sep=",")

for isolate in st_list:
    isolate_accessory_fasta_seqs = []
    isolate_data = roary_data[isolate]
    isolate_data = isolate_data[isolate_data.notna()].index

    print("Working on isolate " + isolate + "...")

    for gene in isolate_data:
        gene_alignment = args.group_alignments_dir + gene + ".aln.fas"

        # If only 1 sample has a gene in an alignment file, that file
        # is named "group_XX.fasta" instead of "group_XX.aln.fasta"
        # elif catches this exception
        if os.path.isfile(gene_alignment):
            seqs = list(SeqIO.parse(gene_alignment, "fasta"))
            to_append = [i for i in seqs if isolate in i.id]
            if len(to_append) > 1:
                print("\tIsolate " + isolate + " had a split gene: " + gene)

        elif not os.path.isfile(gene_alignment):
            gene_alignment = args.group_alignments_dir + gene + ".fasta"
            seqs = list(SeqIO.parse(gene_alignment, "fasta"))
            to_append = [i for i in seqs if isolate in i.id]
            if len(to_append) > 1:
                print("\tIsolate " + isolate + " had a split gene: " + gene)

        isolate_accessory_fasta_seqs += to_append

    print("\n\tNote: for each split gene, there is 1 additional sequence added to the total written to file.\n")

    print("\n\tNumber of accessory genes: " + str(len(isolate_data)))
    print("\tNumber of individual sequences written to file: " + str(len(isolate_accessory_fasta_seqs)) + "\n")
    isolate_output_fasta = args.output_dir + isolate + "_accessory_genes.fasta"
    SeqIO.write(isolate_accessory_fasta_seqs, isolate_output_fasta, "fasta")
