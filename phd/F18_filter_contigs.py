import argparse
from Bio import SeqIO

parser = argparse.ArgumentParser()

parser.add_argument("--mincontiglen", help="minimum contig length to keep")
parser.add_argument("--input_contigs", help="input contigs file to filter")
parser.add_argument("--select_blast_contigs", help="remove contigs that did not blast to target genus? 'yes' or 'no', default='no'", default="no")
parser.add_argument("--parsed_blast_output", help="parsed blast_output file")
parser.add_argument("--filtered_assembly", help="filtered assembly output file")
parser.add_argument("--assembler", help="unicycler or spades, used for parsing contig names, default=unicycler", default="unicycler")

args = parser.parse_args()

## Filter contigs based on length:
contigs = list(SeqIO.parse(args.input_contigs, "fasta"))
print("Original number of contigs: ", str(len(contigs)))

contigs_to_keep = []

# Identify contigs shorter than specified length and remove:
for contig in contigs:
    length = len(contig.seq)
    if length >= int(args.mincontiglen):
        contigs_to_keep.append(contig)

## Identify contigs that blast to target organism:
if args.select_blast_contigs == "yes":
    print("Parsing contigs to keep based on blast results...")

    if args.assembler == "unicycler": # if assembler is unicycler, parse filenames a certain way
        blast_contigs_to_keep = []
        with open(args.parsed_blast_output, 'r') as infile2:
            for line in infile2:
                line_elements = line.strip().split("\t")
                if line_elements[4] == "blue":
                    blast_contigs_to_keep.append(line_elements[0])
        blast_contigs_to_keep_set = set(blast_contigs_to_keep)

        final_contig_list = []
        for contig in contigs_to_keep:
            if contig.id in blast_contigs_to_keep_set:
                final_contig_list.append(contig)

    elif args.assembler == "spades":
        blast_contigs_to_keep = []
        with open(args.parsed_blast_output, 'r') as infile2:
            for line in infile2:
                line_elements = line.strip().split("\t")
                if line_elements[4] == "blue":
                    blast_contigs_to_keep.append(line_elements[0])
        blast_contigs_to_keep_set = set(blast_contigs_to_keep)

        final_contig_list = []
        for contig in contigs_to_keep:
            contig_num = contig.description.split("_")[1] # Only dif is that we need to parse contig number from contig description for SPAdes to match blast contig ids
            if contig_num in blast_contigs_to_keep_set:
                final_contig_list.append(contig)

elif args.select_blast_contigs == "no":
    print("Skipping parsing contigs based on blast results...")
    final_contig_list = contigs_to_keep

# Write output to file:
print("Final number of contigs: ", str(len(final_contig_list)))

SeqIO.write(final_contig_list, args.filtered_assembly, "fasta")


# Write good contigs to file:
# SeqIO.write(contigs_to_keep, args.out, "fasta")
