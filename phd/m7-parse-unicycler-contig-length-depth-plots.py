import argparse

parser = argparse.ArgumentParser()

parser.add_argument("--input_assembly", help="input fasta assembly")
parser.add_argument("--input_contig_depths", help="input contig depths summary from m7 bash script")
parser.add_argument("--isolate", help="current isolate being worked on")
parser.add_argument("--output_dir", help="output directory")

args = parser.parse_args()


contigs_dict = {}
with open(args.input_assembly, "r") as infile1:
    for line in infile1:
        if line.startswith(">"):
            line_elements = line.strip().split(" ")
            contig_num = line_elements[0][1:]
            contig_length = line_elements[1].split("=")[-1]
            contig_rel_depth = line_elements[2].split("=")[-1][:-1]
            contigs_dict[contig_num] = [contig_num, contig_length, contig_rel_depth, "deepskyblue", "blue"]

with open(args.input_contig_depths, "r") as infile2:
    for line in infile2:
        if line.startswith("/home"):
            contig_num = line.strip().split("_")[-1].split(".")[0]
            contig_depth = next(infile2).strip().split(" ")[-1]
            contigs_dict[contig_num].append(contig_depth)

for key in contigs_dict:
    if len(contigs_dict[key]) < 6:
        contigs_dict[key].append("0")
    contigs_dict[key] = "\t".join(contigs_dict[key])

output_handle = args.output_dir + args.isolate + "_parsed_output.txt"
to_write = "\n".join([contigs_dict[key] for key in contigs_dict])

with open(output_handle, "w") as outfile1:
    outfile1.write("Contig\tLength\tRel_Depth\tColor\tOutline\tDepth\n")
    outfile1.write(to_write)
