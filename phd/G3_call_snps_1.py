import argparse
import subprocess
import os

parser = argparse.ArgumentParser()

"""
Script to help with per-ST SNP calling. Note that all paths must be given in full.
"""

parser.add_argument("--project_dir", help="directory where all output files will be generated")
parser.add_argument("--output_prefix", help="prefix for all output files")
parser.add_argument("--isolate_list", help="list of isolate names, one per line")
parser.add_argument("--fastq_dir", help="directory with fastq files, one pair of files for each isolate in isolate list. Note that these must end in '.fastq' or '.fastq.gz'.")
parser.add_argument("--reference_fasta", help="reference sequence fasta file")
parser.add_argument("--threads", help="number of threads to use, default=8", default="8")
parser.add_argument("--snippy_minfrac", help="snippy --minfrac option, default=0.9", default="0.9")

args = parser.parse_args()

# Check to make sure all options are present (all are required):
if args.project_dir is None:
    parser.error("--project_dir is required")
if args.output_prefix is None:
    parser.error("--output_prefix is required")
if args.isolate_list is None:
    parser.error("--isolate_list is required")
if args.fastq_dir is None:
    parser.error("--fastq_dir is required")
if args.reference_fasta is None:
    parser.error("--reference_fasta is required")

# Check if required tools are in the path:
def is_tool(name):
    """Check whether `name` is on PATH and marked as executable."""

    from shutil import which

    return which(name) is not None

if is_tool("snippy") == False:
    raise Exception("Snippy not in path; make sure you are in an env with snippy in the path.")

if is_tool("snp-sites") == False:
    raise Exception("snp-sites not in path; make sure you are in an env with snp-sites in the path.")

if is_tool("snp-dists") == False:
    raise Exception("snp-dists not in path; make sure you are in an env with snp-dists in the path.")

if is_tool("run_gubbins.py") == False:
    raise Exception("run_gubbins.py not in path; make sure you are in an env with run_gubbins.py in the path.")

# Change to project directory:
os.chdir(args.project_dir)

# Read isolate list:
isolate_list = []
with open(args.isolate_list, 'r') as infile1:
    for line in infile1:
        isolate_list.append(line.strip())

# Set up and run Snippy:

#os.mkdir("raw_snippy_output")

for isolate in isolate_list:
    print("Running snippy with cmd: snippy --cpus " + args.threads + " --outdir " + os.path.join(args.project_dir, "raw_snippy_output", isolate) + " --ref " + args.reference_fasta + " --R1 " + args.fastq_dir + isolate + "_1.fastq.gz --R2 " + args.fastq_dir + isolate + "_2.fastq.gz --minfrac " + args.snippy_minfrac)

    subprocess.run()
