import argparse
import os

parser = argparse.ArgumentParser()

parser.add_argument("--nt", help="number of threads to use (<= 56)")
parser.add_argument("--walle", help="walltime estimate HH:MM")
parser.add_argument("--wallm", help="walltime max HH:MM")
parser.add_argument("--fastq_dir", help="Synergy dir with fastq files, must end with '/'")
parser.add_argument("--output_dir", help="Synergy output dir, must end with '/'")
parser.add_argument("--isolate", help="isolate name/number")
parser.add_argument("--job_file", help="output job file")
parser.add_argument("--crop", help="length to crop reads to (250, 300, etc.)", default=300)
parser.add_argument("--minlen", help="minimum length for reads after trimming", default=30)
parser.add_argument("--minqual", help="minimum quality for 4 bp sliding window", default=5)

args = parser.parse_args()

# Set memory requirements
mem = int(args.nt) * 4000
mem_max = mem + 2000
shovill_ram = str(int(((int(args.nt) * 4000) - 1000)/1000))

# Set directory names
project_dir = args.output_dir
assembly_dir = args.output_dir + args.isolate

# Header:
header1_l1 = "#! /usr/bin/env bash\n\n"
header1_l2 = "#BSUB -J " + args.isolate + "_trm_job" + '\n'
header1_l3 = "#BSUB -n " + str(args.nt) + '\n'
header1_l4 = '#BSUB -R "span[hosts=1]"' + '\n'
header1_l5 = '#BSUB -R "rusage[mem=' + str(mem) + ']"' + '\n'
# header1_l6 = '#BSUB -R "select[hname!=node013]"' + '\n'
header1_l7 = '#BSUB -M ' + str(mem_max) + '\n'
header1_l8 = '#BSUB -We ' + args.walle + '\n'
header1_l9 = '#BSUB -W ' + args.wallm + '\n'
header1_l10 = '#BSUB -o ' + project_dir + args.isolate + '.out' + '\n'
header1_l11 = '#BSUB -e ' + project_dir + args.isolate + '.err' + '\n'

header1 = header1_l1 + header1_l2 + header1_l3 + header1_l4 + header1_l5 + header1_l7 + header1_l8 + header1_l9 + header1_l10 + header1_l11 + '\n'

trim_cmd = "trimmomatic PE -threads " + args.nt + " " + args.fastq_dir + args.isolate + "_1.fq" + " " + args.fastq_dir + args.isolate + "_2.fq" + " " + args.output_dir + args.isolate + "_1.fastq.gz" + " " + args.output_dir + args.isolate + "_u_1.fastq.gz" + " " + args.output_dir + args.isolate + "_2.fastq.gz" + " " + args.output_dir + args.isolate + "_u_2.fastq.gz ILLUMINACLIP:/home/cizydorczyk/NexteraPE-PE.fa:2:30:10:8:true CROP:" + args.crop +" SLIDINGWINDOW:4:" + args.minqual + " MINLEN:" + args.minlen

with open(args.job_file, 'w') as outfile:
    outfile.write(header1 + trim_cmd)
