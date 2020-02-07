import argparse
import os
import random

parser = argparse.ArgumentParser()

# General required arguments:
parser.add_argument("--fastq_dir", help="dir with fastq files on Synergy")
parser.add_argument("--project_dir", help="output dir")
parser.add_argument("--isolate_list", help="isolate list, one isolate per line")
parser.add_argument("--gz", help="gzipped reads? yes or no, default = 'yes'", default='yes')
parser.add_argument("--fastq_endings", help="ending indicating fastq forward/reverse reads, e.g. _1.fastq/_2.fastq Default = '_1.fastq,_2.fastq'.", default='_1.fastq,_2.fastq')
parser.add_argument("--assembler", help="assembler to use; currently only 'unicycler' or 'spades' (default) is supported.", default = "spades")

# Synergy options (required):
parser.add_argument("--nt", help="number of threads to use (<= 56)")
parser.add_argument("--walle", help="walltime estimate HH:MM")
parser.add_argument("--wallm", help="walltime max HH:MM")

# Unicycler-specific options:
parser.add_argument("--depth_filter", help="unicycler --depth_filter option, default = 0.25", default=0.25)

# BLAST-specific options:
parser.add_argument("--evalue", help="evalue to use for blasting contigs, default = 0.00001", default=0.00001)

# Prokka-specific options:
parser.add_argument("--proteins", help="gbk file from which to first annotate; prokka --proteins option. Must provide full path to file found on Synergy. Default=''; if not given, no file will be used.", default="")
parser.add_argument("--kingdom", help="prokka --kingdom flag, default='Bacteria'")
parser.add_argument("--genus", help="prokka --genus flag")
parser.add_argument("--species", help="prokka --species flag")
parser.add_argument("--use_genus", help="True or False, prokka --use_genus flag, default = 'False'; no longer recommended by Prokka", default="False")
parser.add_argument("--gcode", help="prokka --gcode flag, default=11 (for bacteria)", default=11)
parser.add_argument("--rnammer", help="prokka --rnammer flag, more accurate than using Barrnap (default in prokka), default=True", default="True")
parser.add_argument("--mincontiglen", help="prokka --mincontiglen, default = 200 (prokka default = 1)", default=200)

args = parser.parse_args()

# Set Synergy memory requirements:
mem = int(args.nt) * 4000
mem_max = mem + 2000

# Read isolate list:
isolate_list = []
with open(args.isolate_list, 'r') as infile1:
    for line in infile1:
        isolate_list.append(line.strip())

### Create directory structure:
# Check project_dir exists:
if os.path.isdir(args.project_dir):
    print("Project dir exists: ", args.project_dir)
else:
    print("Making project dir: ", args.project_dir)
    os.mkdir(args.project_dir)

## Create general subdirectories & set directory variables:

synergy_prefix = "/home/cizydorczyk"
synergy_project_dir = os.path.join(synergy_prefix, args.project_dir.split("/")[-2])

# Unicycler or SPAdes directories:
if args.assembler == "unicycler":
    if not os.path.exists(os.path.join(args.project_dir, "unicycler_assemblies")):
        os.mkdir(os.path.join(args.project_dir, "unicycler_assemblies"))
    s_unicycler_assemblies_dir = os.path.join(synergy_project_dir, "unicycler_assemblies")

elif args.assembler == "spades":
    if not os.path.exists(os.path.join(args.project_dir, "spades_assemblies")):
        os.mkdir(os.path.join(args.project_dir, "spades_assemblies"))
    s_spades_assemblies_dir = os.path.join(synergy_project_dir, "spades_assemblies")

# QUAST directories:
if not os.path.exists(os.path.join(args.project_dir, "quast")):
    os.mkdir(os.path.join(args.project_dir, "quast"))
s_quast_dir = os.path.join(synergy_project_dir, "quast")

if not os.path.exists(os.path.join(args.project_dir, "quast", "raw_assemblies")):
    os.mkdir(os.path.join(args.project_dir, "quast", "raw_assemblies"))
s_quast_raw_assemblies_dir = os.path.join(synergy_project_dir, "quast", "raw_assemblies")

if not os.path.exists(os.path.join(args.project_dir, "quast", "final_assemblies")):
    os.mkdir(os.path.join(args.project_dir, "quast", "final_assemblies"))
s_quast_final_assemblies_dir = os.path.join(synergy_project_dir, "quast", "final_assemblies")

# BLAST directories:
if not os.path.exists(os.path.join(args.project_dir, "blast")):
    os.mkdir(os.path.join(args.project_dir, "blast"))
s_blast_dir = os.path.join(synergy_project_dir, "blast")

if not os.path.exists(os.path.join(args.project_dir, "blast", "raw_blast_output")):
    os.mkdir(os.path.join(args.project_dir, "blast", "raw_blast_output"))
s_blast_raw_output_dir = os.path.join(synergy_project_dir, "blast", "raw_blast_output")

if not os.path.exists(os.path.join(args.project_dir, "blast", "parsed_blast_output")):
    os.mkdir(os.path.join(args.project_dir, "blast", "parsed_blast_output"))
s_blast_parsed_output_dir = os.path.join(synergy_project_dir, "blast", "parsed_blast_output")

# Filtered assembly directories:
if not os.path.exists(os.path.join(args.project_dir, "filtered_assemblies")):
    os.mkdir(os.path.join(args.project_dir, "filtered_assemblies"))
s_filtered_assemblies_dir = os.path.join(synergy_project_dir, "filtered_assemblies")

if not os.path.exists(os.path.join(args.project_dir, "filtered_assemblies", "contig_len_and_blast_filtered_assemblies")):
    os.mkdir(os.path.join(args.project_dir, "filtered_assemblies", "contig_len_and_blast_filtered_assemblies"))
s_contig_len_and_blast_filtered_assemblies_dir = os.path.join(synergy_project_dir, "filtered_assemblies", "contig_len_and_blast_filtered_assemblies")

if not os.path.exists(os.path.join(args.project_dir, "filtered_assemblies", "mauve_ordered_assemblies")):
    os.mkdir(os.path.join(args.project_dir, "filtered_assemblies", "mauve_ordered_assemblies"))
s_contig_len_and_blast_filtered_assemblies_dir = os.path.join(synergy_project_dir, "filtered_assemblies", "mauve_ordered_assemblies")

if not os.path.exists(os.path.join(args.project_dir, "filtered_assemblies", "final_assemblies")): # Final assemblies directory:
    os.mkdir(os.path.join(args.project_dir, "filtered_assemblies", "final_assemblies"))
s_final_assemblies_dir = os.path.join(synergy_project_dir, "filtered_assemblies", "final_assemblies")

# Prokka directories:
if not os.path.exists(os.path.join(args.project_dir, "prokka_annotations")):
    os.mkdir(os.path.join(args.project_dir, "prokka_annotations"))
s_prokka_dir = os.path.join(synergy_project_dir, "prokka_annotations")

# Contig coverage directories:
# if not os.path.exists(os.path.join(args.project_dir, "contig_coverage")):
#     os.mkdir(os.path.join(args.project_dir, "contig_coverage"))
# s_contig_coverage_dir = os.path.join(synergy_project_dir, "contig_coverage")

# Jobs directories:
if not os.path.exists(os.path.join(args.project_dir, "jobfiles_dir")):
    os.mkdir(os.path.join(args.project_dir, "jobfiles_dir"))
jobfiles_dir = os.path.join(args.project_dir, "jobfiles_dir")

if args.assembler == "unicycler":
    if not os.path.exists(os.path.join(jobfiles_dir, "unicycler_jobs")):
        os.mkdir(os.path.join(jobfiles_dir, "unicycler_jobs"))
    unicycler_jobs_dir = os.path.join(jobfiles_dir, "unicycler_jobs")

elif args.assembler == "spades":
    if not os.path.exists(os.path.join(jobfiles_dir, "spades_jobs")):
        os.mkdir(os.path.join(jobfiles_dir, "spades_jobs"))
    spades_jobs_dir = os.path.join(jobfiles_dir, "spades_jobs")

if not os.path.exists(os.path.join(jobfiles_dir, "raw_quast_job")):
    os.mkdir(os.path.join(jobfiles_dir, "raw_quast_job"))
quast_raw_jobs_dir = os.path.join(jobfiles_dir, "raw_quast_job")

if not os.path.exists(os.path.join(jobfiles_dir, "blast_jobs")):
    os.mkdir(os.path.join(jobfiles_dir, "blast_jobs"))
blast_jobs_dir = os.path.join(jobfiles_dir, "blast_jobs")

if not os.path.exists(os.path.join(jobfiles_dir, "prokka_jobs")):
    os.mkdir(os.path.join(jobfiles_dir, "prokka_jobs"))
prokka_jobs_dir = os.path.join(jobfiles_dir, "prokka_jobs")

# if not os.path.exists(os.path.join(jobfiles_dir, "contig_coverage_jobs")):
#     os.mkdir(os.path.join(jobfiles_dir, "contig_coverage_jobs"))
# contig_coverage_jobs_dir = os.path.join(jobfiles_dir, "contig_coverage_jobs")

### Create Unicycler jobs ###
if args.assembler == "unicycler":

    import F18_unicycler_jobs

    for isolate in isolate_list:
        F18_unicycler_jobs.CreateUnicyclerJobs(s_unicycler_assemblies_dir,
        args.fastq_dir, unicycler_jobs_dir, isolate, args.nt, args.walle,
        args.wallm, args.depth_filter, args.fastq_endings, args.gz, mem, mem_max)

elif args.assembler == "spades":

    import F18_spades_jobs

    for isolate in isolate_list:
        F18_spades_jobs.CreateSpadesJobs(s_spades_assemblies_dir, args.fastq_dir,
        spades_jobs_dir, isolate, args.nt, args.walle, args.wallm, args.fastq_endings,
        args.gz, mem, mem_max)

### Create QUAST job for raw assemblies ###
import F18_quast_jobs

if args.assembler == "unicycler":
    F18_quast_jobs.CreateQuastJobs(s_quast_raw_assemblies_dir,
    s_unicycler_assemblies_dir, quast_raw_jobs_dir, isolate_list, args.nt,
    args.walle, args.wallm, mem, mem_max, args.assembler)

elif args.assembler == "spades":
    F18_quast_jobs.CreateQuastJobs(s_quast_raw_assemblies_dir,
    s_spades_assemblies_dir, quast_raw_jobs_dir, isolate_list, args.nt,
    args.walle, args.wallm, mem, mem_max, args.assembler)

### Create contig coverage jobs ###
# import F18_contig_coverage_jobs
#
# for isolate in isolate_list:
#     F18_contig_coverage_jobs.CreateContigCoverageJobs(s_contig_coverage_dir, s_unicycler_assemblies_dir, args.fastq_dir,
#     contig_coverage_jobs_dir, isolate, args.nt, args.walle, args.wallm, mem, mem_max, args.gz, args.fastq_endings, args.assembler)

### Create BLAST jobs ###
import F18_blast_jobs

for isolate in isolate_list:
    if args.assembler == "unicycler":
        F18_blast_jobs.CreateBlastJob(s_blast_raw_output_dir, s_unicycler_assemblies_dir, blast_jobs_dir, isolate, args.nt, args.walle, args.wallm, mem, mem_max, args.evalue)

    elif args.assembler == "spades":
        F18_blast_jobs.CreateBlastJob(s_blast_raw_output_dir, s_spades_assemblies_dir, blast_jobs_dir, isolate, args.nt, args.walle, args.wallm, mem, mem_max, args.evalue)

### Create Prokka jobs ###
import F18_prokka_jobs

for isolate in isolate_list:
    F18_prokka_jobs.CreateProkkaJob(s_prokka_dir, s_final_assemblies_dir, prokka_jobs_dir, isolate, args.nt, args.walle, args.wallm, mem, mem_max, args.proteins, args.kingdom, args.genus, args.use_genus, args.species, args.gcode, args.rnammer, args.mincontiglen)
