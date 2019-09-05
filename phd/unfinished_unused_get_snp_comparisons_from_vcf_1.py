import argparse
import os.path
from Bio import SeqIO

parser = argparse.ArgumentParser()

parser.add_argument("--input_vcfs", nargs='*', help="list of vcf files")

args = parser.parse_args()

# Create empty list to hold VCF variants:
# snps1 = []
# insertions1 = []
# deletions1 = []
#
# snps2 = []
# insertions2 = []
# deletions2 = []

# Define VcfVariants class;
class VcfVariant(object):
    def __init__(self, isolateid, chrom, pos, ref, alt, filename):
        self.isolateid = isolateid
        self.chrom = chrom
        self.pos = pos
        self.ref = ref
        self.alt = alt
        self.filename = filename

vcfs_dict = {}
isolates = []
for i in args.input_vcfs:
    isolate = i.strip().split('/')[-2]
    isolates.append(isolate)

    # initialize dict entry for each isolate:
    vcfs_dict[isolate] = []

    # Initialize SNP, indel, and other lists for isolate:
    snps = {}
    insertions = {}
    deletions = {}
    other = {}

    # Get SNPs, indels, and other variants from isolate VCF file:
    with open(i, 'r') as infile1:
        for line in infile1:
            if not line.startswith("#"):
                line_list = line.strip().split('\t')
                if len(line_list[3]) == 1 and len(line_list[4]) == 1:
                    vcf_variant_obj = VcfVariant(isolate, line_list[0], line_list[1], line_list[3], line_list[4], i)
                    snps[int(line_list[1])] = vcf_variant_obj
                    # snps.append((line_list[1], vcf_variant_obj)) # format is a tuple: (position, VcfVariant object)
                elif len(line_list[3]) > 1 and len(line_list[3]) > len(line_list[4]):
                    vcf_variant_obj = VcfVariant(isolate, line_list[0], line_list[1], line_list[3], line_list[4], i)
                    deletions[int(line_list[1])] = vcf_variant_obj
                    # deletions.append((line_list[1], vcf_variant_obj))
                elif len(line_list[4]) > 1 and len(line_list[4]) > len(line_list[3]):
                    vcf_variant_obj = VcfVariant(isolate, line_list[0], line_list[1], line_list[3], line_list[4], i)
                    insertions[int(line_list[1])] = vcf_variant_obj
                    # insertions.append((line_list[1], vcf_variant_obj))
                else:
                    vcf_variant_obj = VcfVariant(isolate, line_list[0], line_list[1], line_list[3], line_list[4], i)
                    other[int(line_list[1])] = vcf_variant_obj
                    # other.append((line_list[1], vcf_variant_obj))

    # Append SNP, indel, and other lists to isolate dict entry:
    vcfs_dict[isolate].append(snps)
    vcfs_dict[isolate].append(insertions)
    vcfs_dict[isolate].append(deletions)
    vcfs_dict[isolate].append(other)

# Sanity check to make sure SNPs, indels, and other (MNPs) are being classified correctly
# Ref < Alt for insertion, Ref > Alt for deletion, Ref == Alt but both > 1 for Other (MNP)
# for i in vcfs_dict["A013-E93-23-04-2014"][2]:
#     print(vcfs_dict["A013-E93-23-04-2014"][2][i].ref, vcfs_dict["A013-E93-23-04-2014"][2][i].alt)

# Create list of all SNP positions in all files being analyzed:
all_snps = []
iso_snp_pres_abs_dict = {}

for iso in isolates:
    snp_positions = list(vcfs_dict[iso][0].keys())
    all_snps = all_snps + snp_positions
    iso_snp_pres_abs_dict[iso] = [[],[]]

all_snps_set = sorted(list(set(all_snps)))

print(iso_snp_pres_abs_dict)

for snp in all_snps_set:
    for iso in isolates:
        isolate_snp_pos = set(vcfs_dict[iso][0].keys())

        if snp in isolate_snp_pos:
            iso_snp_pres_abs_dict[iso][0].append(1) # 0 element is presence absence recorded as 1/0
            iso_snp_pres_abs_dict[iso][1].append(vcfs_dict[iso][0][snp].alt) # 1 element is alt base (if snp present) or '-' (if snp absent)'
        elif snp not in isolate_snp_pos:
            iso_snp_pres_abs_dict[iso][0].append(0)
            iso_snp_pres_abs_dict[iso][1].append("-")

print(iso_snp_pres_abs_dict)
