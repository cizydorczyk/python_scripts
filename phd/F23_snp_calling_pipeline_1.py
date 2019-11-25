import argparse
import os
import random
import subprocess

parser = argparse.ArgumentParser()

parser.add_argument("--fastq_dir", help="dir with fastq files on Synergy")
parser.add_argument("--project_dir", help="output dir")
parser.add_argument("--isolate_list", help="isolate list, one isolate per line")
parser.add_argument("--gz", help="gzipped reads? yes or no, default = 'yes'", default='yes')
parser.add_argument("--reference", help="reference genome; should be ST-specific; use mash distances/mashtree to find closest ref")
parser.add_argument("--st", help="identifier for the SNP calling run; could be ST or other identifier")
parser.add_argument("--fastq_endings", help="ending indicating fastq forward/reverse reads, e.g. _1.fastq/_2.fastq Default = '_1.fastq,_2.fastq'.", default='_1.fastq,_2.fastq')
parser.add_argument("--nt", help="number of threads to use")
parser.add_argument("--minfrac", help="snippy minfrac option")

# Control options:
parser.add_argument("--call_snps", action="store_true")
parser.add_argument("--run_snippy_core", action="store_true")
parser.add_argument("--run_iqtree", action="store_true")
parser.add_argument("--run_cfml", action="store_true")
parser.add_argument("--clean_snippy_core_full_aln", action="store_true")
parser.add_argument("--remove_ref", action="store_true")
parser.add_argument("--create_raw_matrices", action="store_true")

args = parser.parse_args()

## Set Synergy memory requirements:
mem = int(args.nt) * 4000
mem_max = mem + 2000

## Read isolate list:
isolate_list = []
with open(args.isolate_list, 'r') as infile1:
    for line in infile1:
        isolate_list.append(line.strip())

# Create indicator as to whether reference genome should be kept in SNP tables, etc. or not
# CFML requires at least 3 isolates, so if only 2 are present, keep ref as 3rd to remove recombination
if len(isolate_list) <= 3:
    keep_ref = True
elif len(isolate_list) > 3:
    keep_ref = False

## Create directory structure:

# Check project_dir exists:
if os.path.isdir(args.project_dir):
    print("Project dir exists: ", args.project_dir)
else:
    print("Making project dir: ", args.project_dir)
    os.mkdir(args.project_dir)

if not os.path.exists(os.path.join(args.project_dir, "raw_snippy_output")):
    os.mkdir(os.path.join(args.project_dir, "raw_snippy_output"))
raw_snippy_output_dir = os.path.join(args.project_dir, "raw_snippy_output")

if not os.path.exists(os.path.join(args.project_dir, "snippy_core_all_isolates")):
    os.mkdir(os.path.join(args.project_dir, "snippy_core_all_isolates"))
snippy_core_all_isolates_dir = os.path.join(args.project_dir, "snippy_core_all_isolates")

if not os.path.exists(os.path.join(args.project_dir, "iqtree_output")):
    os.mkdir(os.path.join(args.project_dir, "iqtree_output"))
iqtree_output_dir = os.path.join(args.project_dir, "iqtree_output")

if not os.path.exists(os.path.join(args.project_dir, "cfml_output")):
    os.mkdir(os.path.join(args.project_dir, "cfml_output"))
cfml_output_dir = os.path.join(args.project_dir, "cfml_output")

## Run snippy:

fastq_endings = args.fastq_endings.split(",")
isolate_snippy_dir_list = []

for isolate in isolate_list:

    if args.gz == 'yes':
        r1 = args.fastq_dir + isolate + fastq_endings[0] + ".gz"
        r2 = args.fastq_dir + isolate + fastq_endings[1] + ".gz"

    elif args.gz == 'no':
        r1 = args.fastq_dir + isolate + fastq_endings[0]
        r2 = args.fastq_dir + isolate + fastq_endings[1]

    isolate_snippy_dir = os.path.join(raw_snippy_output_dir, isolate)
    isolate_snippy_dir_list.append(isolate_snippy_dir)

    snippy_cmd = "snippy --cleanup --cpus " + args.nt + " --outdir " + isolate_snippy_dir + " --ref " + args.reference + " --R1 " + r1 + " --R2 " + r2 + " --minfrac " + str(args.minfrac)
    #print(snippy_cmd)

    if args.call_snps:
        subprocess.run(snippy_cmd, shell=True)

## Run snippy core:
snippycore_cmd = "cd " + snippy_core_all_isolates_dir + "; snippy-core --prefix " + args.st + " --ref " + args.reference + " " + " ".join(isolate_snippy_dir_list)
#print(snippycore_cmd)

if args.run_snippy_core:
    subprocess.run(snippycore_cmd, shell=True)

## Clean up snippy core & remove reference (if required):
snippy_clean_cmd = "snippy-clean_full_aln " + snippy_core_all_isolates_dir + "/" + args.st + ".full.aln > " + snippy_core_all_isolates_dir + "/clean." + args.st + ".full.aln"

#print(snippy_clean_cmd)

if args.clean_snippy_core_full_aln:
    subprocess.run(snippy_clean_cmd, shell=True)

import F23_remove_ref

clean_fasta = snippy_core_all_isolates_dir + "/clean." + args.st + ".full.aln"

# Set no ref/clean fasta file variable depending on whether removing ref or not.
# Doing so now simplifies future calls, always using "no_ref_clean_fasta" as
# fasta file to work with, regardless of whether ref was removed or not.
if keep_ref == False:
    no_ref_clean_fasta = snippy_core_all_isolates_dir + "/no_ref.clean." + args.st + ".full.aln"

elif keep_ref == True:
    no_ref_clean_fasta = snippy_core_all_isolates_dir + "/clean." + args.st + ".full.aln"

# Remove ref if >3 isolates:
if args.remove_ref:
    if keep_ref == False:
        F23_remove_ref.RemoveReference(clean_fasta, no_ref_clean_fasta)
    elif keep_ref == True:
        print("Keeping reference to ensure <=3 sequences present for CFML. Skipping removing reference...")

## Create SNP distances matrix for all & core SNPs before removing recombinant sites:

# Reference included for core SNPs (easier to run w/out removing...):
snpdists_core_cmd = "snp-dists " + snippy_core_all_isolates_dir + "/" + args.st + ".aln > " + snippy_core_all_isolates_dir + "/core_snps." + args.st + ".aln.matrix.txt"

snpdists_all_cmd = "snp-dists " + no_ref_clean_fasta + " > " + snippy_core_all_isolates_dir + "/all_snps." + args.st + ".aln.matrix.txt"

if args.create_raw_matrices:
    subprocess.run(snpdists_core_cmd, shell=True)
    subprocess.run(snpdists_all_cmd, shell=True)

## Run iqtree:
if args.run_iqtree:
    cp_cmd = "cp " + no_ref_clean_fasta + " " + iqtree_output_dir + "/copy." + args.st + ".full.aln"
    subprocess.run(cp_cmd, shell=True)

    if keep_ref == False:
        iqtree_no_ref_clean_fasta = iqtree_output_dir + "/copy." + args.st + ".full.aln"
        iqtree_cmd = "iqtree -pre " + iqtree_output_dir + "/" + args.st + " -m MFP -nt AUTO -s " + iqtree_no_ref_clean_fasta + " -bb 10000"
    elif keep_ref == True:

        if len(isolate_list) == 3: # bootstrap if by including ref we have 4 isolates total
            iqtree_no_ref_clean_fasta = iqtree_output_dir + "/copy." + args.st + ".full.aln"
            iqtree_cmd = "iqtree -pre " + iqtree_output_dir + "/" + args.st + " -m MFP -nt AUTO -s " + iqtree_no_ref_clean_fasta + " -bb 10000"
        elif len(isolate_list) < 3:
            iqtree_no_ref_clean_fasta = iqtree_output_dir + "/copy." + args.st + ".full.aln"
            iqtree_cmd = "iqtree -pre " + iqtree_output_dir + "/" + args.st + " -m MFP -nt AUTO -s " + iqtree_no_ref_clean_fasta # bootstrap makes no sense for <4 seq!

    subprocess.run(iqtree_cmd, shell=True)

## Run CFML:
if keep_ref == False:
    iqtree_tree = iqtree_output_dir + "/" + args.st + ".contree"
elif keep_ref == True:

    if len(isolate_list) == 3: # there will be a contree file if 3 isolate + ref = 4 (bootstrap runs)
        iqtree_tree = iqtree_output_dir + "/" + args.st + ".contree"
    elif len(isolate_list) < 3:
        iqtree_tree = iqtree_output_dir + "/" + args.st + ".treefile"

cfml_cmd = "cd " + cfml_output_dir + "; ClonalFrameML " + iqtree_tree + " " + no_ref_clean_fasta + " " + args.st + ".cfml" + " -em true -show_progress true"

if args.run_cfml:
    subprocess.run(cfml_cmd, shell=True)
