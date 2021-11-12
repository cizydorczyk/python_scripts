## SNP-Calling Pipeline
## Written by: Conrad Izydorczyk
## Last Updated: November 3, 2021

import argparse
import os
import random
import subprocess
import i18_remove_ref
import shlex # for converting string commands to lists for subprocess.run to avoid use of shell=True
import sys

parser = argparse.ArgumentParser()

# General script options:
parser.add_argument("--fastq_dir", help="dir with fastq files on local; all fastq must be in single directory, isolates & outgroups included")
parser.add_argument("--project_dir", help="output dir")
parser.add_argument("--isolate_list", help="isolate list, one isolate per line")
parser.add_argument("--outgroup_id", help=("outgroup id; "
    "name of outgroup used to name reads; comma separated list possible"), default="")
parser.add_argument("--gz", help=("gzipped reads? yes or no, default = 'yes'"),
    default='yes')
parser.add_argument("--reference", help=("reference genome; should be "
    "ST-specific; use mash distances/mashtree to find closest ref"))
parser.add_argument("--st", help=("identifier for the SNP calling run; "
    "could be ST or other identifier"))
parser.add_argument("--fastq_endings", help=("ending indicating fastq "
    "forward/reverse reads, e.g. _1.fastq/_2.fastq Default = "
    "'_1.fastq,_2.fastq'."), default='_1.fastq,_2.fastq')
parser.add_argument("--nt", help="number of threads to use")

# Snippy options:
parser.add_argument("--minfrac", help="snippy minfrac option", type=float, required=True)
parser.add_argument("--cleanup", help="include Snippy --cleanup flag? Include this flag to run (default not run).", action="store_true")
parser.add_argument("--unmapped", help="Snippy --unmapped flag. Include to run (default = no).", action="store_true")

# Control options:
parser.add_argument("--call_snps", action="store_true")
parser.add_argument("--snippy_core", action="store_true")
parser.add_argument("--remove_ref", action="store_true")
parser.add_argument("--remove_outgroup", action="store_true")
parser.add_argument("--create_raw_matrices", action="store_true")
parser.add_argument("--iqtree", action="store_true")
parser.add_argument("--cfml", action="store_true")

args = parser.parse_args()

## Can't set --cleanup and --unmapped together (it just deletes the bam/read files...):
if args.cleanup and args.unmapped:
    print("Can't set --cleanup and --unmapped together. Exiting.")
    sys.exit()

## Read isolate list:
isolate_list = []
with open(args.isolate_list, 'r') as infile1:
    for line in infile1:
        if line.strip() != "": # necessary to deal with blank lines at end of file..they break the entire script past SNP calling...
            isolate_list.append(line.strip())

## Add outgroup to isolate list; outgroup fastq must be in same dir & named same as for rest of isolates:
isolate_list_w_outgroup = isolate_list.copy()
#isolate_list_w_outgroup.append(args.outgroup_id)
if args.outgroup_id == "":
    pass
elif len(args.outgroup_id.split(",")) == 1:
    isolate_list_w_outgroup.append(args.outgroup_id)
elif len(args.outgroup_id.split(",")) > 1:
    isolate_list_w_outgroup = isolate_list_w_outgroup + args.outgroup_id.split(",")


## Create directory structure:
raw_snippy_output_dir = os.path.join(args.project_dir, "raw_snippy_output")
snippy_core_dir = os.path.join(args.project_dir, "snippy_core")
iqtree_dir = os.path.join(args.project_dir, "iqtree")
cfml_dir = os.path.join(args.project_dir, "cfml")

if not os.path.isdir(args.project_dir):
    os.mkdir(args.project_dir)
if not os.path.isdir(raw_snippy_output_dir):
    os.mkdir(raw_snippy_output_dir)
if not os.path.isdir(snippy_core_dir):
    os.mkdir(snippy_core_dir)
if not os.path.isdir(iqtree_dir):
    os.mkdir(iqtree_dir)
if not os.path.isdir(cfml_dir):
    os.mkdir(cfml_dir)

## Run snippy:

# Check appropriate Snippy control options:
if 0.0 < float(args.minfrac) <= 1.0:
    pass
else:
    print("--minfrac must be a decimal between 0 and 1. Exiting.")
    sys.exit()

if args.cleanup:
    cleanup_var = "--cleanup "
else:
    cleanup_var = ""

if args.unmapped:
    unmapped_var = "--unmapped "
else:
    unmapped_var = ""


fastq_endings = args.fastq_endings.split(",")
snippy_dir_list = []

for isolate in isolate_list_w_outgroup:
    if args.gz == 'yes':
        r1 = args.fastq_dir + isolate + fastq_endings[0] + ".gz"
        r2 = args.fastq_dir + isolate + fastq_endings[1] + ".gz"
    elif args.gz == 'no':
        r1 = args.fastq_dir + isolate + fastq_endings[0]
        r2 = args.fastq_dir + isolate + fastq_endings[1]

    isolate_snippy_dir = os.path.join(raw_snippy_output_dir, isolate)
    snippy_dir_list.append(isolate_snippy_dir) # append isolate raw snippy dir to facilitate running snippy-core later

    snippy_cmd = (
        f"snippy "
        f"{cleanup_var}"
        f"{unmapped_var}"
        f"--cpus {args.nt} "
        f"--outdir {isolate_snippy_dir} "
        f"--ref {args.reference} "
        f"--R1 {r1} "
        f"--R2 {r2} "
        f"--minfrac {args.minfrac}"
    )

    # run snippy cmd:
    print(f"Running Snippy command: {snippy_cmd}")
    if args.call_snps: # allows for control of pipeline steps
        subprocess.run(shlex.split(snippy_cmd), shell=False)

## Run snippy core:
core_aln = (f"{snippy_core_dir}/{args.st}.aln")
full_aln = (f"{snippy_core_dir}/{args.st}.full.aln")
clean_full_aln = (f"{snippy_core_dir}/clean.{args.st}.full.aln") # output for snippy-clean_full_aln cmd

if args.snippy_core:
    snippy_core_output_prefix = os.path.join(snippy_core_dir, args.st)
    snippy_core_cmd = (f"snippy-core "
        f"--prefix {snippy_core_output_prefix} "
        f"--ref {args.reference} {' '.join(snippy_dir_list)}"
    )
    print(snippy_core_cmd)
    subprocess.run(shlex.split(snippy_core_cmd), shell=False)

    snippy_clean_cmd = (f"snippy-clean_full_aln {full_aln}")
    with open(clean_full_aln, 'w') as outfile:
        subprocess.run(shlex.split(snippy_clean_cmd), shell=False, stdout=outfile)

## Create fasta files with/without ref, outgroup:
# NOTE: that clean.{args.st}.full.aln is the full aln WITH ref & outgroup
no_ref_fasta = (f"{snippy_core_dir}/no_ref.clean.{args.st}.full.aln")
no_outgroup_fasta = (f"{snippy_core_dir}/no_outgroup.clean.{args.st}.full.aln")
no_ref_no_outgroup_fasta = (f"{snippy_core_dir}/no_ref_no_outgroup.clean."
    f"{args.st}.full.aln")

if args.remove_ref:
    i18_remove_ref.RemoveReference("Reference", clean_full_aln, no_ref_fasta)
if args.remove_outgroup:
    i18_remove_ref.RemoveReference(args.outgroup_id, clean_full_aln, no_outgroup_fasta)
if args.remove_ref and args.remove_outgroup:
    tmp_no_ref_aln = (f"{snippy_core_dir}/tmp.no_ref.aln")
    i18_remove_ref.RemoveReference("Reference", clean_full_aln, tmp_no_ref_aln)
    i18_remove_ref.RemoveReference(args.outgroup_id, tmp_no_ref_aln, no_ref_no_outgroup_fasta)
    subprocess.run(["rm", f"{tmp_no_ref_aln}"], shell=False)

## Create raw distance matrices with ref, with outgroup, and without both
raw_all_dist_matrix = (f"{clean_full_aln}.raw.all_snps.matrix.txt")
raw_core_dist_matrix = (f"{core_aln}.raw.core_snps.matrix.txt")

if args.create_raw_matrices:
    create_raw_all_dist_matrix_cmd = (f"snp-dists {clean_full_aln}")
    with open(raw_all_dist_matrix, 'w') as outfile:
        subprocess.run(shlex.split(create_raw_all_dist_matrix_cmd), shell=False, stdout=outfile)

    create_raw_core_dist_matrix_cmd = (f"snp-dists {core_aln}")
    with open(raw_core_dist_matrix, 'w') as outfile:
        subprocess.run(shlex.split(create_raw_core_dist_matrix_cmd), shell=False, stdout=outfile)

## Run IQ-Tree & CFML & Create additional matrices for no_ref and no_ref_no_outgroup; options for others exist but not default

# IQ-Tree with/without outgroup
iqtree_no_ref_prefix = (f"{iqtree_dir}/with_outgroup.{args.st}")
if len(isolate_list_w_outgroup) <= 3:
    no_ref_tree = (f"{iqtree_no_ref_prefix}.treefile") # .treefile if no bootstrap run
    iqtree_cmd = (f"iqtree "
        f"-pre {iqtree_no_ref_prefix} "
        f"-m MFP "
        f"-s {no_ref_fasta} " # with outgroup
        f"-nt AUTO"
    )
elif len(isolate_list_w_outgroup) > 3:
    no_ref_tree = (f"{iqtree_no_ref_prefix}.contree") # .contree if bootstrap run
    iqtree_cmd = (f"iqtree "
        f"-pre {iqtree_no_ref_prefix} "
        f"-m MFP "
        f"-s {no_ref_fasta} " # with outgroup
        f"-nt AUTO "
        f"-bb 10000"
    )

iqtree_no_ref_no_outgroup_prefix = (f"{iqtree_dir}/without_outgroup.{args.st}")
if len(isolate_list) <= 3:
    no_ref_no_outgroup_tree = (f"{iqtree_no_ref_no_outgroup_prefix}.treefile") # .treefile if no bootstrap run
    iqtree_cmd2 = (f"iqtree "
        f"-pre {iqtree_no_ref_no_outgroup_prefix} "
        f"-m MFP "
        f"-s {no_ref_no_outgroup_fasta} " # without outgroup
        f"-nt AUTO"
    )
elif len(isolate_list) > 3:
    no_ref_no_outgroup_tree = (f"{iqtree_no_ref_no_outgroup_prefix}.contree") # .contree if bootstrap run
    iqtree_cmd2 = (f"iqtree "
        f"-pre {iqtree_no_ref_no_outgroup_prefix} "
        f"-m MFP "
        f"-s {no_ref_no_outgroup_fasta} " # without outgroup
        f"-nt AUTO "
        f"-bb 10000"
    )

if args.iqtree:
    subprocess.run(shlex.split(iqtree_cmd), shell=False)
    subprocess.run(shlex.split(iqtree_cmd2), shell=False)

## Run CFML with and without outgroup

# CFML with outgroup, no ref
cfml_with_outgroup_prefix = (f"{cfml_dir}/with_outgroup.{args.st}")
cfml_cmd = (f"ClonalFrameML "
    f"{no_ref_tree} "
    f"{no_ref_fasta} "
    f"{cfml_with_outgroup_prefix} "
    f"-em true -emsim 100 -show_progress true"
)

# CFML without outgroup, no ref
cfml_without_outgroup_prefix = (f"{cfml_dir}/without_outgroup.{args.st}")
cfml_cmd2 = (f"ClonalFrameML "
    f"{no_ref_no_outgroup_tree} "
    f"{no_ref_no_outgroup_fasta} "
    f"{cfml_without_outgroup_prefix} "
    f"-em true -emsim 100 -show_progress true"
)

if args.cfml:
    subprocess.run(shlex.split(cfml_cmd), shell=False)
    subprocess.run(shlex.split(cfml_cmd2), shell=False)

### NOTE: RC masking & final matrices must be created using separate script.
### maskrc-svg requires downgrading of packages I do not want to downgrade in
### case they mess up the rest of this pipleine (i.e. Snippy, etc.).
