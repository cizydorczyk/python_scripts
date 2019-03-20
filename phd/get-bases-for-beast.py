import argparse
import subprocess
from Bio import SeqIO

parser = argparse.ArgumentParser()

parser.add_argument("--aln", help="recombination masked whole genome alignment from which to get # of bases for beast")
parser.add_argument("--ref", help="reference sequence")

args = parser.parse_args()

records = list(SeqIO.parse(args.aln, "fasta"))
reference = list(SeqIO.parse(args.ref, "fasta"))
ref_length = len(reference[0].seq)
print(ref_length)


# Check records contain only 'ATCGN-':
acceptable_bases = set(['A', 'T', 'C', 'G', 'N', '-'])
for record in records:
    bases_set = set(record.seq)
    if bases_set.issubset(acceptable_bases):
        continue
    else:
        print("Error: make sure sequences in input file only contain 'ATCGN-'.")
        print(bases_set)
        break

# Run snp-sites (must be in path!) to generate variant position VCF:

# Check if snp-sites is in path:
try:
    subprocess.run(['snp-sites', '-V'], check=True)
except subprocess.CalledProcessError as err:
    print('ERROR:', err)

print("Running snp-sites to generate variant position VCF...")
subprocess.run(['snp-sites', '-v', '-o', 'temp.vcf', args.aln])

variant_positions = []
with open("temp.vcf", "r") as infile:
    for line in infile:
        if not line.startswith("#"):
            variant_positions.append(int(line.strip().split("\t")[1]))

print("Number of variant positions: %s" % str(len(variant_positions)))

variant_positions_set = set(variant_positions)
bases_count_dict = {"A":0, "T":0, "C":0, "G":0, "N":0, "Other":0}

reference_seq = list(reference[0].seq)
excluded_pos_count = 0

print("Counting bases at invariant positions...")
for i in range(0, len(reference_seq)):
    if i+1 not in variant_positions_set:
        if reference_seq[i] in acceptable_bases:
            bases_count_dict[reference_seq[i]] += 1
        else:
            bases_count_dict["Other"] += 1
    else:
        excluded_pos_count += 1

print(bases_count_dict)
subprocess.run(["rm", "temp.vcf"])
print("Invariant bases sum: %s" % str(bases_count_dict["A"] + bases_count_dict["T"] + bases_count_dict["C"] + bases_count_dict["G"] + bases_count_dict["N"] + bases_count_dict["Other"]))
print("Done.")

# for index in range(0,seq_length):
#     # print(records[1].seq[index])
#     bases = []
#     for record in records:
#         # print(index, record.seq[index])
#         bases.append(record.seq[index])
#
#     bases_set = set(bases)
#     if len(bases_set) == 1:
#         print(index, bases_set)
#     elif len(bases_set) == 2 and 'N' in bases_set:
#         print(index, "N")
#     elif len(bases_set) == 2 and '-' in bases_set:
#         print(index, "-")
#     elif len(bases_set) > 2 and 'N' in bases_set and '-' in bases_set:
#         print(index, 'N,-')
#     else:
#         print(index, bases_set)
