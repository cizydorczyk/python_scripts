import argparse
import os.path
import random

parser = argparse.ArgumentParser()

## End all directories with a "/"!

parser.add_argument("--nt", help="number of threads to use (<= 56)")
parser.add_argument("--walle", help="walltime estimate HH:MM")
parser.add_argument("--wallm", help="walltime max HH:MM")
parser.add_argument("--assembly", help="path to assembly on synergy")
parser.add_argument("--output_dir", help="output dir on synergy")
parser.add_argument("--isolate", help="isolate name/number")
parser.add_argument("--job_file", help="output job file")

args = parser.parse_args()

# Set memory requirements
mem = int(args.nt) * 4000
mem_max = mem + 2000

# Header:
header1_l1 = "#! /usr/bin/env bash\n\n"
header1_l2 = "#BSUB -J " + args.isolate + "_art_job" + '\n'
header1_l3 = "#BSUB -n " + str(args.nt) + '\n'
header1_l4 = '#BSUB -R "span[hosts=1]"' + '\n'
header1_l5 = '#BSUB -R "rusage[mem=' + str(mem) + ']"' + '\n'
# header1_l6 = '#BSUB -R "select[hname!=node013]"' + '\n'
header1_l7 = '#BSUB -M ' + str(mem_max) + '\n'
header1_l8 = '#BSUB -We ' + args.walle + '\n'
header1_l9 = '#BSUB -W ' + args.wallm + '\n'
header1_l10 = '#BSUB -o ' + args.output_dir + args.isolate + '.out' + '\n'
header1_l11 = '#BSUB -e ' + args.output_dir + args.isolate + '.err' + '\n'

header1 = header1_l1 + header1_l2 + header1_l3 + header1_l4 + header1_l5 + header1_l7 + header1_l8 + header1_l9 + header1_l10 + header1_l11 + '\n'

art_cmd = "/home/cizydorczyk/art_bin_MountRainier/art_illumina -ss MSv3 -i " + args.assembly + " -l 250 -p -m 400 -s 110 -f 100 -o " + args.output_dir + args.isolate + "_ -na -rs 42"

# shovill_cmd = "shovill --outdir " + assembly_dir + " --R1 " + args.fastq_dir + args.isolate + "_1.fastq.gz --R2 " + args.fastq_dir + args.isolate + '_2.fastq.gz --minlen 200 --mincov 10 --ram ' + shovill_ram + ' --cpus 0 --assembler spades --opts "--careful"'

with open(args.job_file, 'w') as outfile:
    outfile.write(header1 + art_cmd)
