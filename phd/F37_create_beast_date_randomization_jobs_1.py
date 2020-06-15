import argparse
import os.path
import random
random.seed()

parser = argparse.ArgumentParser()

parser.add_argument("--input_xml", help="full path to xml file") # must be JUST the filename, so must run this script from the directory with xml files...
parser.add_argument("--nt", help="number of threads to use (<= 56)")
parser.add_argument("--project_dir", help="specify local project directory")
parser.add_argument("--job_name_prefix", help="prefix for job names")
parser.add_argument("--walle", help="walltime estimate HH:MM")
parser.add_argument("--wallm", help="walltime max HH:MM")
parser.add_argument("--synergy_project_dir", help="project directory in synergy")
parser.add_argument("--rep_num", help="rep number; integer")

args = parser.parse_args()

# Set seeds
seed1 = random.randint(1,1000000000)

# Create trial directories
model = args.input_xml.strip().split('.')[-2]

os.mkdir(args.project_dir + model)

os.mkdir(args.project_dir + model + '/trial_1/')

# Set synergy project directory:
synergy_project_dir = args.synergy_project_dir + model + '/'

# Create job names
jobname1 = "r" + args.rep_num

# Set memory requirements
mem = int(args.nt) * 4000
mem_max = mem + 1000

## Trial 1
# Header 1:
header1_l1 = "#! /usr/bin/env bash\n\n"
header1_l2 = "#BSUB -J " + jobname1 + '\n'
header1_l3 = "#BSUB -n " + str(args.nt) + '\n'
header1_l4 = '#BSUB -R "span[hosts=1]"' + '\n'
header1_l5 = '#BSUB -R "rusage[mem=' + str(mem) + ']"' + '\n'
header1_l6 = '#BSUB -R "select[hname!=node007]"' + '\n'
header1_l7 = '#BSUB -M ' + str(mem_max) + '\n'
header1_l8 = '#BSUB -We ' + args.walle + '\n'
header1_l9 = '#BSUB -W ' + args.wallm + '\n'
header1_l10 = '#BSUB -o ' + synergy_project_dir + 'trial_1/' + jobname1 + '.out' + '\n'
header1_l11 = '#BSUB -e ' + synergy_project_dir + 'trial_1/' + jobname1 + '.err' + '\n'

# header1 = header1_l1 + header1_l2 + header1_l3 + header1_l4 + header1_l5 + header1_l6 + header1_l7 + header1_l8 + header1_l9 + header1_l10 + header1_l11 + '\n'
header1 = header1_l1 + header1_l2 + header1_l3 + header1_l4 + header1_l5 + header1_l7 + header1_l8 + header1_l9 + header1_l10 + header1_l11 + '\n'


# BEAST command
beast_cmd1 = '/home/cizydorczyk/BEASTv1.10.4/bin/beast -seed ' + str(seed1) + ' -working ' + synergy_project_dir + 'trial_1/' + args.input_xml.split('/')[-1]

# Write job script:
outfile_1 = args.project_dir + model + '/' + 'trial_1/' + jobname1 + '.sh'
with open(outfile_1, 'w') as outfile1:
    to_write1 = header1 + beast_cmd1
    outfile1.write(to_write1)

# Copy XML to trial directories
import subprocess
subprocess.run(['cp', args.input_xml, args.project_dir + model + '/' + 'trial_1/'])
