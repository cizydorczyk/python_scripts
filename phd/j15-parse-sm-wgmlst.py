import argparse
from Bio import SeqIO

parser = argparse.ArgumentParser()

parser.add_argument("--input_fasta", help="input fasta to split")
parser.add_argument("--output_dir", help="output directory for split fasta")

args = parser.parse_args()

genes_dict = {}

for record in SeqIO.parse(args.input_fasta, "fasta"):
    gene_id = record.id.split("_")[1]

    try:
        genes_dict[gene_id].append(record)
    except:
        genes_dict[gene_id] = [record]

for gene in genes_dict:
    output_file = args.output_dir + str(gene) + ".fasta"

    with open(output_file, "w") as outfile1:
        SeqIO.write(genes_dict[gene], outfile1, "fasta")
