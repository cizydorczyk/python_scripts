import argparse

parser = argparse.ArgumentParser()

parser.add_argument("--input_fa", help="path to fasta on Synergy for CFML")
parser.add_argument("--output_dir", help="path to main project dir on synergy")
parser.add_argument("--sample", help="sample/isolate #")
parser.add_argument("--tree", help="full path to tree file on synergy")
parser.add_argument("--job_file", help="local job file")

args = parser.parse_args()

h1 = '#! /usr/bin/env bash'
h2 = '#BSUB -J cfml1'
h3 = '#BSUB -n 1'
h4 = '#BSUB -R "span[hosts=1]"'
h5 = '#BSUB -R "rusage[mem=4000]"'
h6 = '#BSUB -M 6000'
h7 = '#BSUB -We 96:00'
h8 = '#BSUB -W 144:00'
h9 = '#BSUB -o ' + args.output_dir + "/" + args.sample + ".cfml.out"
h10 = '#BSUB -e ' + args.output_dir + "/" + args.sample + ".cfml.err"

header = h1+'\n\n'+h2+'\n'+h3+'\n'+h4+'\n'+h5+'\n'+h6+'\n'+h7+'\n'+h8+'\n'+h9+'\n'+h10+'\n\n'

cfml_cmd = "ClonalFrameML " + args.tree + " " + args.input_fa + " " + args.output_dir + "/" + args.sample + ".cfml" + " -em true -show_progress true"

with open(args.job_file, 'w') as outfile:
    outfile.write(header + cfml_cmd)
