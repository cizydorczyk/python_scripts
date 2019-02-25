import argparse
import subprocess
import os
from Bio import SeqIO

parser = argparse.ArgumentParser()

# Trimmomatic arguments:
# parser.add_argument('--trimmomatic_env', action='store', required=True, help='trimmomatic conda env name')
parser.add_argument("--trim_qual", action='store', type=int, required=True, help="quality threshold for trimmomatic SLIDINGWINDOW")
parser.add_argument("--read_length", action='store', type=int, required=True, help='read length. trimmomatic will crop to this length after removing adapters.')
parser.add_argument("--adapters", action='store', required=True, help='full path to fasta with adapters to trim')
parser.add_argument('--trim_min_len', action='store', required=True, help='minimum read length to keep after all trimming')

# FastQC arguments:
# parser.add_argument("--fastqc_env", action='store', required=True, help='fastqc conda environment name')

# Kraken arguments:
# parser.add_argument("--kraken_env", action='store', required=True, help='Kraken2 conda environment name')
parser.add_argument("--kraken_db", action='store', required=True, help='Kraken2 db to use (minikraken db)')

# SRST2 arguments:
# parser.add_argument("--srst2_env", action='store', required=True, help='SRST2 conda environment name')
# parser.add_argument("--srst2_mlst_db", action='store', required=True, help='MLST scheme/db to use')

# Snippy arguments:
# parser.add_argument("--snippy_env", action='store', required=True, help='Snippy conda env name')
parser.add_argument("--reference", action='store', required=True, help='Reference fasta for SNP calling')
# parser.add_argument("--cpus", action='store', type=int, required=True, help='number threads for Snippy')
parser.add_argument("--minfrac", action='store', type=float, required=True, help='minfrac Snippy value')

# IQ-Tree arguments:
# parser.add_argument("--iqtree_env", action='store', required=True, help='IQ-Tree conda env name')
# parser.add_argument("--iqtree_model", action='store', required=True, help='IQ-Tree model. Recomended is "GTR+I+R"')
parser.add_argument("--rapid_bs", action='store', type=int, required=True, help='# rapid bootstraps (rec. 10000)')

# Gubbins arguments:
# parser.add_argument("--gubbins_env", action='store', required=True, help='Gubbins conda env name')
parser.add_argument("--gubbins_iter", action='store', type=int, required=True, help='# Gubbins iterations')

# General arguments:
parser.add_argument("--project_dir", required=True, help='directory where all output will be generated. does not have to already exist')
parser.add_argument("--input_fastq_dir", required=True, help='directory with input (untrimmed) fastq files')
parser.add_argument("--isolate_list", required=True, help='isolate list')

args = parser.parse_args()

### Start of script

## Make directory structure:

# Check if project dir exists and create if not:
if os.path.isdir(args.project_dir):
    print("Project directory: %s" % args.project_dir)
elif not os.path.isdir(args.project_dir):
    os.mkdir(args.project_dir)
    print("Project directory created: %s" % args.project_dir)

# Trimmomatic directories:

if os.path.isdir(os.path.join(args.project_dir, 'fastq_files')):
    print("fastq directory: %s" % os.path.join(args.project_dir, 'fastq_files'))
else:
    os.mkdir(os.path.join(args.project_dir, 'fastq_files'))
    print("fastq directory: %s" % os.path.join(args.project_dir, 'fastq_files'))

if os.path.isdir(os.path.join(args.project_dir, 'fastq_files', 'trimmed_fastq')):
    print("fastq directory (trimmed): %s" % os.path.join(args.project_dir, 'fastq_files', 'trimmed_fastq'))
else:
    os.mkdir(os.path.join(args.project_dir, 'fastq_files', 'trimmed_fastq'))
    print("fastq directory (trimmed): %s" % os.path.join(args.project_dir, 'fastq_files', 'trimmed_fastq'))

if os.path.isdir(os.path.join(args.project_dir, 'fastq_files', 'trimmed_fastq', 'paired')):
    print("fastq directory (trimmed, paired) %s" % os.path.join(args.project_dir, 'fastq_files', 'trimmed_fastq', 'paired'))
else:
    os.mkdir(os.path.join(args.project_dir, 'fastq_files', 'trimmed_fastq', 'paired'))
    print("fastq directory (trimmed, paired): %s" % os.path.join(args.project_dir, 'fastq_files', 'trimmed_fastq', 'paired'))

if os.path.isdir(os.path.join(args.project_dir, 'fastq_files', 'trimmed_fastq', 'unpaired')):
    print("fastq directory (trimmed, unpaired): %s" % os.path.join(args.project_dir, 'fastq_files', 'trimmed_fastq', 'unpaired'))
else:
    os.mkdir(os.path.join(args.project_dir, 'fastq_files', 'trimmed_fastq', 'unpaired'))
    print("fastq directory (trimmed, unpaired): %s" % os.path.join(args.project_dir, 'fastq_files', 'trimmed_fastq', 'unpaired'))

if os.path.isdir(os.path.join(args.project_dir, 'fastq_files', 'trimmed_fastq', 'summary_files')):
    print("trimming summary files directory: %s" % os.path.join(args.project_dir, 'fastq_files', 'trimmed_fastq', 'summary_files'))
else:
    os.mkdir(os.path.join(args.project_dir, 'fastq_files', 'trimmed_fastq', 'summary_files'))
    print("trimming summary files directory: %s" % os.path.join(args.project_dir, 'fastq_files', 'trimmed_fastq', 'summary_files'))

# Fastqc directories:
if os.path.isdir(os.path.join(args.project_dir, "fastqc")):
    print("fastqc directory: %s" % os.path.join(args.project_dir, "fastqc"))
else:
    os.mkdir(os.path.join(args.project_dir, "fastqc"))
    print("fastqc directory: %s" % os.path.join(args.project_dir, "fastqc"))

if os.path.isdir(os.path.join(args.project_dir, "fastqc", "raw_fastq")):
    print("fastqc directory (raw fastq): %s" % os.path.join(args.project_dir, "fastqc", "raw_fastq"))
else:
    os.mkdir(os.path.join(args.project_dir, "fastqc", "raw_fastq"))
    print("fastqc directory (raw fastq): %s" % os.path.join(args.project_dir, "fastqc", "raw_fastq"))

if os.path.isdir(os.path.join(args.project_dir, "fastqc", "trimmed_fastq")):
    print("fastqc directory (trimmed fastq): %s" % os.path.join(args.project_dir, "fastqc", "trimmed_fastq"))
else:
    os.mkdir(os.path.join(args.project_dir, "fastqc", "trimmed_fastq"))
    print("fastqc directory (trimmed fastq): %s" % os.path.join(args.project_dir, "fastqc", "trimmed_fastq"))

# Kraken directories:
if os.path.isdir(os.path.join(args.project_dir, "kraken")):
    print("kraken directory: %s" % os.path.join(args.project_dir, "kraken"))
else:
    os.mkdir(os.path.join(args.project_dir, "kraken"))
    print("kraken directory: %s" % os.path.join(args.project_dir, "kraken"))

if os.path.isdir(os.path.join(args.project_dir, "kraken", "kraken_trimmed_fastq")):
    print("kraken directory (trimmed fastq): %s" % os.path.join(args.project_dir, "kraken", "kraken_trimmed_fastq"))
else:
    os.mkdir(os.path.join(args.project_dir, "kraken", "kraken_trimmed_fastq"))
    print("kraken directory (trimmed fastq): %s" % os.path.join(args.project_dir, "kraken", "kraken_trimmed_fastq"))

if os.path.isdir(os.path.join(args.project_dir, "kraken", "kraken_trimmed_fastq", "kraken_output")):
    print("kraken directory (kraken output): %s" % os.path.join(args.project_dir, "kraken", "kraken_trimmed_fastq", "kraken_output"))
else:
    os.mkdir(os.path.join(args.project_dir, "kraken", "kraken_trimmed_fastq", "kraken_output"))
    print("kraken directory (kraken output): %s" % os.path.join(args.project_dir, "kraken", "kraken_trimmed_fastq", "raw_kraken_output"))

if os.path.isdir(os.path.join(args.project_dir, "kraken", "kraken_trimmed_fastq", "kraken_reports")):
    print("kraken directory (reports): %s" % os.path.join(args.project_dir, "kraken", "kraken_trimmed_fastq", "kraken_reports"))
else:
    os.mkdir(os.path.join(args.project_dir, "kraken", "kraken_trimmed_fastq", "kraken_reports"))
    print("kraken directory (reports): %s" % os.path.join(args.project_dir, "kraken", "kraken_trimmed_fastq", "kraken_reports"))

# SRST2 directories:
# if os.path.isdir(os.path.join(args.project_dir, "mlst")):
#     print("mlst directory: %s" % os.path.join(args.project_dir, "mlst"))
# else:
#     os.mkdir(os.path.join(args.project_dir, "mlst"))
#     print("mlst directory: %s" % os.path.join(args.project_dir, "mlst"))
#
# if os.path.isdir(os.path.join(args.project_dir, "mlst", "srst2_mlst")):
#     print("mlst directory(srst2): %s" % os.path.join(args.project_dir, "mlst", "srst2_mlst"))
# else:
#     os.mkdir(os.path.join(args.project_dir, "mlst", "srst2_mlst"))
#     print("mlst directory(srst2): %s" % os.path.join(args.project_dir, "mlst", "srst2_mlst"))
#
# if os.path.isdir(os.path.join(args.project_dir, "mlst", "srst2_mlst", "output")):
#     print("mlst directory(srst2 output): %s" % os.path.join(args.project_dir, "mlst", "srst2_mlst", "output"))
# else:
#     os.mkdir(os.path.join(args.project_dir, "mlst", "srst2_mlst", "output"))
#     print("mlst directory(srst2 output): %s" % os.path.join(args.project_dir, "mlst", "srst2_mlst", "output"))

# Snippy directories:
if os.path.isdir(os.path.join(args.project_dir, "snp_calling")):
    print("SNP calling directory: %s" % os.path.join(args.project_dir, "snp_calling"))
else:
    os.mkdir(os.path.join(args.project_dir, "snp_calling"))
    print("SNP calling directory: %s" % os.path.join(args.project_dir, "snp_calling"))

if os.path.isdir(os.path.join(args.project_dir, "snp_calling", "all_isolates_snp_calling")):
    print("SNP calling directory (all isolates output): %s" % os.path.join(args.project_dir, "snp_calling", "all_isolates_snp_calling"))
else:
    os.mkdir(os.path.join(args.project_dir, "snp_calling", "all_isolates_snp_calling"))
    print("SNP calling directory (all isolates output): %s" % os.path.join(args.project_dir, "snp_calling", "all_isolates_snp_calling"))

# IQ-Tree directories:
if os.path.isdir(os.path.join(args.project_dir, "phylogenetic_analyses")):
    print("Phylo trees directory: %s" % os.path.join(args.project_dir, "phylogenetic_analyses"))
else:
    os.mkdir(os.path.join(args.project_dir, "phylogenetic_analyses"))
    print("Phylo trees directory: %s" % os.path.join(args.project_dir, "phylogenetic_analyses"))

if os.path.isdir(os.path.join(args.project_dir, "phylogenetic_analyses", "all_isolates_iqtree")):
    print("Phylo trees directory (all isolates): %s" % os.path.join(args.project_dir, "phylogenetic_analyses", "all_isolates_iqtree"))
else:
    os.mkdir(os.path.join(args.project_dir, "phylogenetic_analyses", "all_isolates_iqtree"))
    print("Phylo trees directory (all isolates): %s" % os.path.join(args.project_dir, "phylogenetic_analyses", "all_isolates_iqtree"))

# Gubbins directories:
if os.path.isdir(os.path.join(args.project_dir, "phylogenetic_analyses", "gubbins_all_isolates")):
    print("Phylo trees directory (gubbins): %s" % os.path.join(args.project_dir, "phylogenetic_analyses", "gubbins_all_isolates"))
else:
    os.mkdir(os.path.join(args.project_dir, "phylogenetic_analyses", "gubbins_all_isolates"))
    print("Phylo trees directory (gubbins): %s" % os.path.join(args.project_dir, "phylogenetic_analyses", "gubbins_all_isolates"))

## Analyses start:

# Get isolates list:
isolate_list = []
with open(args.isolate_list, "r") as isolatelist:
    for line in isolatelist:
        isolate_list.append(line.strip())

# Get list of fastq files:
input_fastq_list = sorted(os.listdir(args.input_fastq_dir))

prefix_list = []
for i in input_fastq_list:
    prefix = '_'.join(i.split('_')[:-1])
    if prefix not in prefix_list:
        prefix_list.append(prefix)

# Get fastq read 1/read 2 format:
suffix_ = input_fastq_list[0].split('_')[-1].split('.')[0]
if suffix_ == 'R1' or suffix_ == 'R2':
    suffix1 = 'R1'
    suffix2 = 'R2'
elif suffix_ == '1' or suffix_ == '2':
    suffix1 = '1'
    suffix2 = '2'

# Get file ending (.fastq or .gz):
file_type = input_fastq_list[0].split('_')[-1].split('.')[-1]
if file_type == 'gz':
    file_ending = '.fastq.gz'
elif file_type == 'fastq':
    file_ending = '.fastq'

## Trim reads:
print("Running Trimmomatic, FastQC, Kraken2, and Snippy...")
for prefix in prefix_list:
    isolate_num = prefix.split('_')[0]
    fastq1 = os.path.join(args.input_fastq_dir + prefix + '_' + suffix1 + file_ending)
    fastq2 = os.path.join(args.input_fastq_dir + prefix + '_' + suffix2 + file_ending)

    # Create variables for output files names for simplicity:
    summary_file = args.project_dir + "/fastq_files/trimmed_fastq/summary_files/" + isolate_num + "_summary_log.txt"
    paired_R1_file = args.project_dir + "/fastq_files/trimmed_fastq/paired/" + isolate_num + "_1.fastq"
    paired_R2_file = args.project_dir + "/fastq_files/trimmed_fastq/paired/" + isolate_num + "_2.fastq"
    unpaired_R1_file = args.project_dir + "/fastq_files/trimmed_fastq/unpaired/" + isolate_num + "_u_1.fastq"
    unpaired_R2_file = args.project_dir + "/fastq_files/trimmed_fastq/unpaired/" + isolate_num + "_u_2.fastq"

    # Run Trimmomatic (the print statement can be used to print the trimmomatic command that will be run):
    # print(' '.join(['trimmomatic', 'PE', '-threads', '8', '-summary', summary_file, fastq1, fastq2, paired_R1_file, unpaired_R1_file, paired_R2_file, unpaired_R2_file, 'ILLUMINACLIP:' + args.adapters + ':2:30:10:8:true', 'CROP:' + str(args.read_length), 'SLIDINGWINDOW:4:' + str(args.trim_qual), 'MINLEN:' + str(args.trim_min_len)]))
    subprocess.run(['trimmomatic', 'PE', '-threads', '8', '-summary', summary_file, fastq1, fastq2, paired_R1_file, unpaired_R1_file, paired_R2_file, unpaired_R2_file, 'ILLUMINACLIP:' + args.adapters + ':2:30:10:8:true', 'CROP:' + str(args.read_length), 'SLIDINGWINDOW:4:' + str(args.trim_qual), 'MINLEN:' + str(args.trim_min_len)])

    ## Run FastQC:

    raw_fastq_dir = os.path.join(args.project_dir, "fastqc", "raw_fastq")
    trimmed_fastq_dir = os.path.join(args.project_dir, "fastqc", "trimmed_fastq")

    subprocess.run(['fastqc', '-t', '8', '-o', raw_fastq_dir, fastq1])
    subprocess.run(['fastqc', '-t', '8', '-o', raw_fastq_dir, fastq2])

    subprocess.run(['fastqc', '-t', '8', '-o', trimmed_fastq_dir, paired_R1_file])
    subprocess.run(['fastqc', '-t', '8', '-o', trimmed_fastq_dir, paired_R2_file])

    ## Run Kraken2:

    kraken_output_file = os.path.join(args.project_dir, "kraken", "kraken_trimmed_fastq", "kraken_output", isolate_num + "_raw_output.txt")
    kraken_report_file = os.path.join(args.project_dir, "kraken", "kraken_trimmed_fastq", "kraken_reports", isolate_num + "_report.txt")

    subprocess.run(['kraken2', '--db', args.kraken_db, '--threads', '8', '--output', kraken_output_file, '--report', kraken_report_file, '--paired', '--use-names', paired_R1_file, paired_R2_file])

    # Run Snippy:

    snippy_output_dir = os.path.join(args.project_dir, "snp_calling", "all_isolates_snp_calling", isolate_num)
    subprocess.run(['snippy', '--outdir', snippy_output_dir, '--R1', paired_R1_file, '--R2', paired_R2_file, '--cpus', '8', '--reference', args.reference, '--minfrac', str(args.minfrac)])

## Run snippy-core/snippy-clean:
snp_calling_dir_list = os.listdir(os.path.join(args.project_dir, "snp_calling", "all_isolates_snp_calling"))
snp_directory_list = []
for i in snp_calling_dir_list:
    directory = os.path.join(args.project_dir, "snp_calling", "all_isolates_snp_calling", i)
    snp_directory_list.append(directory)
snp_calling_directory_list = ' '.join(snp_directory_list)

snippy_core_cmd = 'cd ' + os.path.join(args.project_dir, "snp_calling", "all_isolates_snp_calling") + ' && snippy-core --ref ' + args.reference + ' ' + snp_calling_directory_list
subprocess.run(snippy_core_cmd, shell=True)
snippy_clean_cmd = 'cd ' + os.path.join(args.project_dir, "snp_calling", "all_isolates_snp_calling") + ' && snippy-clean_full_aln core.full.aln > clean.core.full.aln'
subprocess.run(snippy_clean_cmd, shell=True)

# Remove reference sequence from clean full genome alignment:

records = list(SeqIO.parse(os.path.join(args.project_dir, "snp_calling", "all_isolates_snp_calling", "clean.core.full.aln"), "fasta"))
snp_records = []
for i in records:
    if i.id == 'Reference':
        continue
    else:
        snp_records.append(i)

with open(os.path.join(args.project_dir, "phylogenetic_analyses", "all_isolates_iqtree", "no_ref.clean.core.full.aln"), 'w') as outfile:
    SeqIO.write(snp_records, outfile, "fasta")

## Run IQ-Tree:
iqtree_cmd = 'iqtree -s ' + os.path.join(args.project_dir, "phylogenetic_analyses", "all_isolates_iqtree", "no_ref.clean.core.full.aln") + ' -bb ' + str(args.rapid_bs) + ' -nt AUTO -m GTR+I+R'
subprocess.run(iqtree_cmd, shell=True)

## Run Gubbins:
starting_tree = os.path.join(args.project_dir, "phylogenetic_analyses", "all_isolates_iqtree", "no_ref.clean.core.full.aln.contree")
gubbins_cmd = 'cd ' + os.path.join(args.project_dir, "phylogenetic_analyses", "gubbins_all_isolates") + ' && run_gubbins.py --starting_tree ' + starting_tree + ' --tree_builder raxml --raxml_model GTRGAMMA --prefix gubbins.all_isolates --verbose --iterations ' + str(args.gubbins_iter) + ' --threads 8 ' + os.path.join(args.project_dir, "phylogenetic_analyses", "all_isolates_iqtree", "no_ref.clean.core.full.aln")
subprocess.run(gubbins_cmd, shell=True)


