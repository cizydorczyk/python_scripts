import argparse

parser = argparse.ArgumentParser()

parser.add_argument("--output_file", help="output file")
parser.add_argument("--isolate_list", help="isolate list; one isolat per line")

parser.add_argument("--input_amrfinderplus_path", help="path to dir with amrfinderplus input files")
parser.add_argument("--input_cardrgi_path", help="path to dir with card input files")
parser.add_argument("--input_resfinder_path", help="path to dir with resfinder input files")

parser.add_argument("--input_amrfinderplus_suffix", help="everything after isolate name in amrfinderplus output file")
parser.add_argument("--input_cardrgi_suffix")
parser.add_argument("--input_resfinder_suffix")

args = parser.parse_args()

# Read isolate list:
with open(args.isolate_list, "r") as infile1:
    isolate_list = [line.strip() for line in infile1]

# Create genes list:
genes_list = []
for isolate in isolate_list:
    amrfinderplus_file = args.input_amrfinderplus_path + isolate + "-" + args.input_amrfinderplus_suffix
    cardrgi_file = args.input_cardrgi_path + isolate + "-" + args.input_cardrgi_suffix
    resfinder_file = args.input_resfinder_path + isolate + "/" + args.input_resfinder_suffix

    with open(amrfinderplus_file, "r") as infile1:
        for line in infile1:
            if not line.startswith("Name"):
                gene = line.strip().split("\t")[2].lower()
                if gene not in genes_list:
                    genes_list.append(gene)

    with open(cardrgi_file, "r") as infile2:
        for line in infile2:
            if not line.startswith("ORF_ID"):
                gene = line.strip().split("\t")[8].lower()
                if gene not in genes_list:
                    genes_list.append(gene)

    with open(resfinder_file, "r") as infile3:
        for line in infile3:
            if not line.startswith("Resistance"):
                gene = line.strip().split("\t")[0].lower()
                if gene not in genes_list:
                    genes_list.append(gene)

# Parse isolate files:
genes_list = sorted(genes_list)
output_dict = {}
for gene in genes_list:
    output_dict[gene] = []

for isolate in isolate_list:
    amrfinderplus_file = args.input_amrfinderplus_path + isolate + "-" + args.input_amrfinderplus_suffix
    cardrgi_file = args.input_cardrgi_path + isolate + "-" + args.input_cardrgi_suffix
    resfinder_file = args.input_resfinder_path + isolate + "/" + args.input_resfinder_suffix

    # Get list of isolate-specific genes:
    amrfinder_genes = []
    with open(amrfinderplus_file, "r") as infile1:
        for line in infile1:
            if not line.startswith("Name"):
                gene = line.strip().split("\t")[2].lower()
                if gene not in amrfinder_genes:
                    amrfinder_genes.append(gene)

    cardrgi_genes = []
    with open(cardrgi_file, "r") as infile2:
        for line in infile2:
            if not line.startswith("ORF_ID"):
                gene = line.strip().split("\t")[8].lower()
                if gene not in cardrgi_genes:
                    cardrgi_genes.append(gene)

    resfinder_genes = []
    with open(resfinder_file, "r") as infile3:
        for line in infile3:
            if not line.startswith("Resistance"):
                gene = line.strip().split("\t")[0].lower()
                if gene not in resfinder_genes:
                    resfinder_genes.append(gene)

    for gene in output_dict:
        key = ""
        if gene in amrfinder_genes:
            key += "A"
        if gene in cardrgi_genes:
            key += "C"
        if gene in resfinder_genes:
            key += "R"

        output_dict[gene].append(key)

# Write output:
header = "Gene\t" + "\t".join(isolate_list)
to_write = []
for gene in output_dict:
    line = gene + "\t" + "\t".join(output_dict[gene])
    to_write.append(line)
with open(args.output_file, "w") as outfile1:
    outfile1.write(header + "\n" + "\n".join(to_write))












#     isolate_genes = []
#
#     with open(amrfinderplus_file, "r") as infile4:
#         for line in infile4:
#             if not line.startswith("Name"):
#                 gene = line.strip().split("\t")[2].lower()
#                 if gene not in isolate_genes:
#                     isolate_genes.append(gene)
#
#     with open(cardrgi_file, "r") as infile5:
#         for line in infile5:
#             if not line.startswith("ORF_ID"):
#                 gene = line.strip().split("\t")[8].lower()
#                 if gene not in isolate_genes:
#                     isolate_genes.append(gene)
#
#     with open(resfinder_file, "r") as infile6:
#         for line in infile6:
#             if not line.startswith("Resistance"):
#                 gene = line.strip().split("\t")[0].lower()
#                 if gene not in isolate_genes:
#                     isolate_genes.append(gene)
#
#     # Record presence/absence for each gene per isolate:
#     for gene in output_dict:
#         if gene in isolate_genes:
#             output_dict[gene].append("1")
#         elif gene not in isolate_genes:
#             output_dict[gene].append("0")
#
# # Write output:
# header = "Gene\t" + "\t".join(isolate_list)
# to_write =
