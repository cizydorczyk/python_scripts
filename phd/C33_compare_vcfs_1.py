import argparse
import os.path
from Bio import SeqIO

parser = argparse.ArgumentParser()

parser.add_argument("--input_vcf1", help="input vcf file")
parser.add_argument("--input_vcf2", help="input vcf file 2")
parser.add_argument("--out", help="output file")
parser.add_argument("--isolate", help="isolate number")
parser.add_argument("--st", help="sequence type (just the number; 131, 73, or 1193)")

args = parser.parse_args()

# Create empty list to hold VCF variants:
snps1 = []
insertions1 = []
deletions1 = []

snps2 = []
insertions2 = []
deletions2 = []

# Define VcfVariants class;
class VcfVariant(object):
    def __init__(self, isolateid, chrom, pos, ref, alt, filename):
        self.isolateid = isolateid
        self.chrom = chrom
        self.pos = pos
        self.ref = ref
        self.alt = alt
        self.filename = filename

# Classify each variant position in VCF file by type of variant (SNP, insertion, or deletion):
with open(args.input_vcf1, 'r') as infile1:
    #isolate_id = os.path.basename(args.input_vcf)[0:-4]
    for line in infile1:
        if not line.startswith("#"):
            line_list = line.strip().split('\t')
            if len(line_list[3]) == 1 and len(line_list[4]) == 1:
                vcf_variant_obj = VcfVariant(args.isolate, line_list[0], line_list[1], line_list[3], line_list[4], args.input_vcf1)
                snps1.append(vcf_variant_obj)
            elif len(line_list[3]) > 1:
                vcf_variant_obj = VcfVariant(args.isolate, line_list[0], line_list[1], line_list[3], line_list[4], args.input_vcf1)
                deletions1.append(vcf_variant_obj)
            elif len(line_list[4]) > 1:
                vcf_variant_obj = VcfVariant(args.isolate, line_list[0], line_list[1], line_list[3], line_list[4], args.input_vcf1)
                insertions1.append(vcf_variant_obj)

with open(args.input_vcf2, 'r') as infile2:
    #isolate_id = os.path.basename(args.input_vcf)[0:-4]
    for line in infile2:
        if not line.startswith("#"):
            line_list = line.strip().split('\t')
            if len(line_list[3]) == 1 and len(line_list[4]) == 1:
                vcf_variant_obj = VcfVariant(args.isolate, line_list[0], line_list[1], line_list[3], line_list[4], args.input_vcf2)
                snps2.append(vcf_variant_obj)
            elif len(line_list[3]) > 1:
                vcf_variant_obj = VcfVariant(args.isolate, line_list[0], line_list[1], line_list[3], line_list[4], args.input_vcf2)
                deletions2.append(vcf_variant_obj)
            elif len(line_list[4]) > 1:
                vcf_variant_obj = VcfVariant(args.isolate, line_list[0], line_list[1], line_list[3], line_list[4], args.input_vcf2)
                insertions2.append(vcf_variant_obj)

snp_pos_1 = []
snp_pos_2 = []

for i in snps1:
    snp_pos_1.append(i.pos)
for j in snps2:
    snp_pos_2.append(j.pos)

snp_pos_1_set = set(snp_pos_1)
snp_pos_2_set = set(snp_pos_2)

print("Number of identical SNPs: ", len(snp_pos_1_set & snp_pos_2_set))
print("Number of different SNPs: ", len(snp_pos_1_set - snp_pos_2_set))
