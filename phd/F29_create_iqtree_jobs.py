import argparse
import os

parser = argparse.ArgumentParser()

parser.add_argument("--input_fa", help="path to fasta on Synergy for CFML")
parser.add_argument("--output_dir", help="path to output dir on synergy (per isolate)")
parser.add_argument("--isolate", help="sample/isolate #")
parser.add_argument("--bb", help="number rapid bootstraps")
parser.add_argument("--model", help="iqtree model; default = 'MFP'", default='MFP')
parser.add_argument("--nt", help="max # threads to request on Synergy")
parser.add_argument("--walle", help="wall time estimate")
parser.add_argument("--wallm", help="wall time max")
parser.add_argument("--job_file", help="jobfile on local")

args = parser.parse_args()

# Set memory requirements
mem = int(args.nt) * 4000
mem_max = mem + 1000

h1 = '#! /usr/bin/env bash'
h2 = '#BSUB -J iqtree' + "_" + args.isolate
h3 = '#BSUB -n ' + args.nt
h4 = '#BSUB -R "span[hosts=1]"'
h5 = '#BSUB -R "rusage[mem=' + str(mem) + ']"'
h6 = '#BSUB -M ' + str(mem_max)
h7 = '#BSUB -We ' + args.walle
h8 = '#BSUB -W ' + args.wallm
h9 = '#BSUB -o ' + args.output_dir + "/" + args.isolate + ".out"
h10 = '#BSUB -e ' + args.output_dir + "/" + args.isolate + ".err"

header = h1+'\n\n'+h2+'\n'+h3+'\n'+h4+'\n'+h5+'\n'+h6+'\n'+h7+'\n'+h8+'\n'+h9+'\n'+h10+'\n\n'

iqtree_cmd = "iqtree -nt AUTO -s " + args.input_fa + " -bb " + args.bb + " -m " + args.model + " -pre " + args.output_dir + "/" + args.isolate + ".iqtree -ntmax " + args.nt

with open(args.job_file, 'w') as outfile:
    outfile.write(header + iqtree_cmd)
