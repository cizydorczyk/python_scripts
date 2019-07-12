import argparse
import os.path
from Bio import SeqIO

parser = argparse.ArgumentParser()

parser.add_argument("--input_vcf", help="input vcf file")
parser.add_argument("--out", help="output file")
parser.add_argument("--isolate", help="isolate number")
parser.add_argument("--st", help="sequence type (just the number; 131, 73, or 1193)")

args = parser.parse_args()

# Create empty list to hold VCF variants:
variants = []

# Define VcfVariants class;
class VcfVariant(object):
    def __init__(self, hypermutationgene, isolateid, chrom, pos, ref, alt, info, filename):
        self.hypermutationgene = hypermutationgene
        self.isolateid = isolateid
        self.chrom = chrom
        self.pos = pos
        self.ref = ref
        self.alt = alt
        self.info = info
        self.filename = filename

# Classify each variant position in VCF file by type of variant (SNP, insertion, or deletion):
with open(args.input_vcf, 'r') as infile1:
    #isolate_id = os.path.basename(args.input_vcf)[0:-4]
    for line in infile1:
        if not line.startswith("#"):
            line_list = line.strip().split('\t')
            vcf_variant_obj = VcfVariant("NA", args.isolate, line_list[0], line_list[1], line_list[3], line_list[4], line_list[7], args.input_vcf)
            variants.append(vcf_variant_obj)

# Create dictionaries with gene corodinates for each ST:
st131_genes = {"dnaQ":[256011,256751], "mutS":[2995859,2998420], "mutL":[4812158,4814008], "mutH":[3153039,3153728], "uvrD":[4355994,4358210], "dam":[1031787,1032551], "mutT":[119861,120259], "mutM":[4099146,4099955], "mutY":[3277431,3278513], "ung":[2888234,2888923], "sodA":[4472754,4473380], "sodB":[1796883,1797464], "oxyR":[4532948,4533865], "polA":[4425861,4428647]}

st73_genes = {"dnaQ":[2918122,2918862], "mutS":[2042831,2045392], "mutL":[2395062,2396909], "mutH":[1904725,1905414], "uvrD":[731123,733285], "dam":[1242933,1243769], "mutT":[2787313,2787711], "mutM":[942748,943557], "mutY":[1783892,1784944], "ung":[520724,521413], "sodA":[613028,613648], "sodB":[3399369,3399950], "oxyR":[544933,545850], "polA":[658874,661660]}

st1193_genes = {"dnaQ":[3805643,3806374], "mutS":[1141749,1144310], "mutL":[4312697,4314544], "mutH":[992459,993148], "uvrD":[4822097,4824259], "dam":[384530,385366], "mutT":[3934381,3934779], "mutM":[75616,76425], "mutY":[871746,872798], "ung":[1269340,1270029], "sodA":[4697362,4697982], "sodB":[2313822,2314403], "oxyR":[4629447,4630364], "polA":[4749699,4752485]}

# Check if any SNPs/insertions/deletions occur in hypermutation genes:

# Create empty lists to hold results:
hypermutation_variants = []

if args.st == "131":
    for variant in variants:
        for gene in st131_genes:
            if int(variant.pos) >= st131_genes[gene][0] and int(variant.pos) <= st131_genes[gene][1]:
                variant.hypermutationgene = gene
                hypermutation_variants.append(variant)
            else:
                continue

elif args.st == "73":
    for variant in variants:
        for gene in st73_genes:
            if int(variant.pos) >= st73_genes[gene][0] and int(variant.pos) <= st73_genes[gene][1]:
                variant.hypermutationgene = gene
                hypermutation_variants.append(variant)
            else:
                continue

elif args.st == "1193":
    for variant in variants:
        for gene in st1193_genes:
            if int(variant.pos) >= st1193_genes[gene][0] and int(variant.pos) <= st1193_genes[gene][1]:
                variant.hypermutationgene = gene
                hypermutation_variants.append(variant)
            else:
                continue

# Write results to file:

header = "hypermutation_gene\tisolate\tchrom\tpos\tref\talt\tinfo"

if len(hypermutation_variants) == 0:
    to_write = 'NA' + '\t' + args.isolate + '\tNA\tNA\tNA\tNA\tNA'

    if os.path.isfile(args.out):
        with open(args.out, 'a') as outfile:
            outfile.write('\n' + to_write)

    else:
        with open(args.out, 'w') as outfile:
            outfile.write(header + '\n' + to_write)

elif len(hypermutation_variants) != 0:
    to_write = []
    for i in hypermutation_variants:
        write_str = i.hypermutationgene + '\t' + i.isolateid + '\t' + i.chrom + '\t' + i.pos + '\t' + i.ref + '\t' + i.alt + '\t' + i.info
        to_write.append(write_str)

    if os.path.isfile(args.out):
        with open(args.out, 'a') as outfile:
            outfile.write('\n' + '\n'.join(to_write))
    else:
        with open(args.out, 'w') as outfile:
            outfile.write(header + '\n' + '\n'.join(to_write))
