import argparse
import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import sys
import os.path

parser = argparse.ArgumentParser()

parser.add_argument("--input_vcf", help="isolate-specific vcf file")
parser.add_argument("--output_vcf", help="output vcf file with snps only")


args = parser.parse_args()

def parse_vcf(vcf_file):
    """
    Parse a VCF file & return a dictionary with VCF entry positions as keys
    and the entire VCF line (i.e. the VCF entry) as the value.

    The function parses the VCF, keeping only entries that are NOT identical
    to the reference and do not contain ambiguous sites.

    -------------------
    Parameters
    -------------------
        vcf_file = single isolate VCF file with genotype indicated by 0/1/2/etc.
        where 0 = reference allele, 1 = first alt allele, 2 = second alt allele,
        etc.
    """
    vcf_pos_dict = []
    with open(vcf_file, 'r') as infile1:
        for line in infile1:
            if not line.startswith("#"):

                line_elements = line.strip().split("\t") # get vcf line elements
                pos = line_elements[1]

                ref_base = [line_elements[3]] # Get reference base
                alt_bases = line_elements[4].strip().split(",") # get all alt bases
                possible_bases = ref_base + alt_bases # create list of ref + alt bases
                genotype = int(line_elements[9]) # get genotype
                alt_base = possible_bases[genotype] # assign alt base based on the index of all bases genotype corresponds to

                # print(line_elements[1], possible_bases, genotype, alt_base)

                if alt_base != "*":
                    if alt_base == ref_base[0]:
                        continue
                    elif alt_base != ref_base[0]:
                        vcf_pos_dict.append(line.strip())
                else:
                    # print(line_elements[1], possible_bases, genotype, alt_base)
                    continue
    return vcf_pos_dict

vcf_snp_list = parse_vcf(args.input_vcf)

with open(args.output_vcf, 'w') as outfile1:
    outfile1.write("\n".join(vcf_snp_list))
