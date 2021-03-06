import argparse
import os.path
import random

parser = argparse.ArgumentParser()

parser.add_argument("--nt", help="number of threads to use (<= 56)")
parser.add_argument("--walle", help="walltime estimate HH:MM")
parser.add_argument("--wallm", help="walltime max HH:MM")
parser.add_argument("--fastq_dir", help="dir with fastq files")
parser.add_argument("--output_dir", help="output dir")
parser.add_argument("--isolate", help="isolate name/number")
parser.add_argument("--job_file", help="output job file")
parser.add_argument("--k", help="comma separated list of kmer sizes between 21 \
                    and 127")
parser.add_argument("--subsampled", help="subsampled reads (_1.sub.fastq.gz)? 'yes' or 'no',", default='no')
parser.add_argument("--covcutoff", help="coverage cutoff", default=10)

args = parser.parse_args()

# Set memory requirements
mem = int(args.nt) * 4000
mem_max = mem + 2000

# Set directory names
project_dir = args.output_dir
assembly_dir = args.output_dir + args.isolate

# Header:
header1_l1 = "#! /usr/bin/env bash\n\n"
header1_l2 = "#BSUB -J " + args.isolate + "_spds_job" + '\n'
header1_l3 = "#BSUB -n " + str(args.nt) + '\n'
header1_l4 = '#BSUB -R "span[hosts=1]"' + '\n'
header1_l5 = '#BSUB -R "rusage[mem=' + str(mem) + ']"' + '\n'
#header1_l6 = '#BSUB -R "select[hname!=node013]"' + '\n'
header1_l7 = '#BSUB -M ' + str(mem_max) + '\n'
header1_l8 = '#BSUB -We ' + args.walle + '\n'
header1_l9 = '#BSUB -W ' + args.wallm + '\n'
header1_l10 = '#BSUB -o ' + project_dir + args.isolate + '.out' + '\n'
header1_l11 = '#BSUB -e ' + project_dir + args.isolate + '.err' + '\n'

header1 = header1_l1 + header1_l2 + header1_l3 + header1_l4 + header1_l5 + \
          header1_l7 + header1_l8 + header1_l9 + header1_l10 +\
          header1_l11 + '\n'

if args.subsampled == 'no':
    spades_cmd = "spades.py -o " + assembly_dir + " -1 " + args.fastq_dir + \
                  args.isolate + "_1.fastq.gz -2 " + args.fastq_dir + \
                  args.isolate + '_2.fastq.gz --careful --threads ' + args.nt +\
                  ' --cov-cutoff ' + args.covcutoff + ' -k ' + args.k
elif args.subsampled == 'yes':
    spades_cmd = "spades.py -o " + assembly_dir + " -1 " + args.fastq_dir + \
                  args.isolate + "_1.sub.fastq.gz -2 " + args.fastq_dir + \
                  args.isolate + '_2.sub.fastq.gz --careful --threads ' + \
                  args.nt + ' --cov-cutoff ' + args.covcutoff + ' -k ' + args.k

with open(args.job_file, 'w') as outfile:
    outfile.write(header1 + spades_cmd)
