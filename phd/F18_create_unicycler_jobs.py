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
parser.add_argument("--depth_filter", help="unicycler --depth_filter option, default = 0.25", default=0.25)
parser.add_argument("--fastq_ending", help="ending indicating fastq forward/reverse reads, e.g. _1.fastq/_2.fastq Default = '_1,_2'.", default='_1.fastq,_2.fastq')
parser.add_argument("--gz", help="gzipped reads? yes or no, default = 'yes'", default='yes')

args = parser.parse_args()

# Set memory requirements
mem = int(args.nt) * 4000
mem_max = mem + 2000

# Output dir is the parent directory of all assemblies, which will be produced
# in their own subdirectories within this parent directory.
# Set individual assembly directory:
assembly_dir = args.output_dir + args.isolate

# Header:
header1_l1 = "#! /usr/bin/env bash\n\n"
header1_l2 = "#BSUB -J " + args.isolate + "_unic_job" + '\n'
header1_l3 = "#BSUB -n " + str(args.nt) + '\n'
header1_l4 = '#BSUB -R "span[hosts=1]"' + '\n'
header1_l5 = '#BSUB -R "rusage[mem=' + str(mem) + ']"' + '\n'
#header1_l6 = '#BSUB -R "select[hname!=node013]"' + '\n'
header1_l7 = '#BSUB -M ' + str(mem_max) + '\n'
header1_l8 = '#BSUB -We ' + args.walle + '\n'
header1_l9 = '#BSUB -W ' + args.wallm + '\n'
header1_l10 = '#BSUB -o ' + args.output_dir + args.isolate + '.out' + '\n'
header1_l11 = '#BSUB -e ' + args.output_dir + args.isolate + '.err' + '\n'

header1 = header1_l1 + header1_l2 + header1_l3 + header1_l4 + header1_l5 + \
          header1_l7 + header1_l8 + header1_l9 + header1_l10 +\
          header1_l11 + '\n'

# Fastq file names:
fastq_endings = args.fastq_ending.split(",")

if args.gz == 'yes':
    r1 = args.fastq_dir + args.isolate + fastq_endings[0] + ".gz"
    r2 = args.fastq_dir + args.isolate + fastq_endings[1] + ".gz"

elif args.gz == 'no':
    r1 = args.fastq_dir + args.isolate + fastq_endings[0]
    r2 = args.fastq_dir + args.isolate + fastq_endings[1]

unicycler_cmd = "unicycler -1 " + r1 + " -2 " + r2 + " -o " + assembly_dir + " -t " + str(args.nt)

with open(args.job_file, 'w') as outfile:
    outfile.write(header1 + unicycler_cmd)
