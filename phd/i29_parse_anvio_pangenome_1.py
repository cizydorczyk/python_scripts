import argparse
import pandas as pd

parser = argparse.ArgumentParser()

parser.add_argument("--input_file", help="anvio summary file")
parser.add_argument("--isolate_list", help="isolate list to analyze")
parser.add_argument("--output_dir", help="dir for output files")

args = parser.parse_args()

isolate_list = []
with open(args.isolate_list, 'r') as infile0:
    for line in infile0:
        isolate_list.append(line.strip())

input_data = pd.read_csv(filepath_or_buffer = args.input_file, sep="\t", header=0, index_col=0)

# Create pangenome fasta alignment:
print("Creating pangenome alignment...")
gene_clusters2 = sorted(list(set(input_data["gene_cluster_id"])))

output_dict2 = {}
for gene in gene_clusters2:
    unique_genomes = sorted(list(set(input_data[input_data["gene_cluster_id"] == gene]["genome_name"])))

    output_dict2[gene] = []
    for isolate in isolate_list:
        if isolate in unique_genomes:
            output_dict2[gene].append("A")
        elif isolate not in unique_genomes:
            output_dict2[gene].append("C")

output_df2 = pd.DataFrame.from_dict(data=output_dict2, orient="columns")
output_df2.index = isolate_list

fasta_dict2 = output_df2.transpose().to_dict(orient="list")

fasta_list2 = []
for key in fasta_dict2:
    header = ">" + key
    seq = "".join(fasta_dict2[key])
    fasta_list2.append(header + "\n" + seq)

## write to file
output_fasta_file2 = args.output_dir + "/pangenome_pseudo_full_aln.fasta"
with open(output_fasta_file2, 'w') as outfile2:
    outfile2.write("\n".join(fasta_list2))

################################################
print("Creating pangenome distance matrix...")
# Get distance counts
input_data = input_data[input_data["num_genomes_gene_cluster_has_hits"] < len(isolate_list)]

gene_clusters = sorted(list(set(input_data["gene_cluster_id"])))

output_dict = {}
for gene in gene_clusters:
    unique_genomes = sorted(list(set(input_data[input_data["gene_cluster_id"] == gene]["genome_name"])))

    output_dict[gene] = []
    for isolate in isolate_list:
        if isolate in unique_genomes:
            output_dict[gene].append("A")
        elif isolate not in unique_genomes:
            output_dict[gene].append("C")

output_df = pd.DataFrame.from_dict(data=output_dict, orient="columns")
output_df.index = isolate_list

# create hamming distance matrix
# code for function from https://stackoverflow.com/questions/42752610/python-how-to-generate-the-pairwise-hamming-distance-matrix
def compute_HammingDistance(X):
    return (X[:, None, :] != X).sum(2)

dist_matrix = compute_HammingDistance(pd.DataFrame.as_matrix(output_df))
to_write_dist_matrix = pd.DataFrame(dist_matrix, index=isolate_list, columns=isolate_list)
#print(to_write_dist_matrix)

## write to file
output_dist_matrix = args.output_dir + "/pangenome_dist_matrix.txt"
to_write_dist_matrix.to_csv(path_or_buf = output_dist_matrix, sep="\t")


# create pseudo fasta file in case required
print("Creating accessory genome alignment...")
fasta_dict = output_df.transpose().to_dict(orient="list")

fasta_list = []
for key in fasta_dict:
    header = ">" + key
    seq = "".join(fasta_dict[key])
    fasta_list.append(header + "\n" + seq)

## write to file
output_fasta_file = args.output_dir + "/pangenome_pseudo_fasta.fasta"
with open(output_fasta_file, 'w') as outfile1:
    outfile1.write("\n".join(fasta_list))
