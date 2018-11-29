import argparse
from Bio.SeqIO.FastaIO import SimpleFastaParser
import re

parser = argparse.ArgumentParser()

parser.add_argument("input_mlst", help="input file from mlst program")
parser.add_argument("output_fasta", help="output file for MSA")
# parser.add_argument("adk_fasta", help="fasta with alleles for gene 1...")
# parser.add_argument("fumc_fasta", help="fasta with alleles for gene 2...")
# parser.add_argument("gyrb_fasta", help="fasta with alleles for gene 3...")
# parser.add_argument("icd_fasta", help="fasta with alleles for gene 4...")
# parser.add_argument("mdh_fasta", help="fasta with alleles for gene 5...")
# parser.add_argument("pura_fasta", help="fasta with alleles for gene 6...")
# parser.add_argument("reca_fasta", help="fasta with alleles for gene 7...")

args = parser.parse_args()

# adk_dict = {}
# fumc_dict = {}
# gyrb_dict = {}
# icd_dict = {}
# mdh_dict = {}
# pura_dict = {}
# reca_dict = {}

# with open(args.adk_fasta, 'r') as infile1:
#     for title, seq in SimpleFastaParser(infile1):
#         adk_dict[title] = seq
#
# with open(args.fumc_fasta, 'r') as infile1:
#     for title, seq in SimpleFastaParser(infile1):
#         fumc_dict[title] = seq
#
# with open(args.gyrb_fasta, 'r') as infile1:
#     for title, seq in SimpleFastaParser(infile1):
#         gyrb_dict[title] = seq
#
# with open(args.icd_fasta, 'r') as infile1:
#     for title, seq in SimpleFastaParser(infile1):
#         icd_dict[title] = seq
#
# with open(args.mdh_fasta, 'r') as infile1:
#     for title, seq in SimpleFastaParser(infile1):
#         mdh_dict[title] = seq
#
# with open(args.pura_fasta, 'r') as infile1:
#     for title, seq in SimpleFastaParser(infile1):
#         pura_dict[title] = seq
#
# with open(args.reca_fasta, 'r') as infile1:
#     for title, seq in SimpleFastaParser(infile1):
#         reca_dict[title] = seq

with open(args.input_mlst, 'r') as infile1:
    raw_mlst_data = list(infile1)

for line in raw_mlst_data:
    if "~" in line or "?" in line or "(-)" in line or "(n)" in line:
        pass
    else:
        line = line.strip().split("\t")
        sample = line[0]
        adk = line[3].split("(")[0].upper() + str(re.search('\(([^)]+)', line[3]).group(1))
        fumc = line[4].split("(")[0].upper() + str(re.search('\(([^)]+)', line[4]).group(1))
        icd = line[6].split("(")[0].upper() + str(re.search('\(([^)]+)', line[6]).group(1))
        mdh = line[7].split("(")[0].upper() + str(re.search('\(([^)]+)', line[7]).group(1))
        pura = line[8].split("(")[0].upper() + str(re.search('\(([^)]+)', line[8]).group(1))
        reca = line[9].split("(")[0].upper() + str(re.search('\(([^)]+)', line[9]).group(1))