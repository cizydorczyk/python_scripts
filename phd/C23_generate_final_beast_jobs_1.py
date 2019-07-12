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

args = parser.parse_args()

# Set seeds
seed1 = random.randint(1,1000000000)
seed2 = random.randint(1,1000000000)
seed3 = random.randint(1,1000000000)
seed4 = random.randint(1,1000000000)
seed5 = random.randint(1,1000000000)

# Create trial directories
model = args.input_xml.strip().split('.')[-2]

os.mkdir(args.project_dir + model)

os.mkdir(args.project_dir + model + '/trial_1/')
os.mkdir(args.project_dir + model + '/trial_2/')
os.mkdir(args.project_dir + model + '/trial_3/')
os.mkdir(args.project_dir + model + '/trial_4/')
os.mkdir(args.project_dir + model + '/trial_5/')

# Set synergy project directory:
synergy_project_dir = args.synergy_project_dir + model + '/'

# Create job names
jobname1 = args.job_name_prefix + "_trial_1"
jobname2 = args.job_name_prefix + "_trial_2"
jobname3 = args.job_name_prefix + "_trial_3"
jobname4 = args.job_name_prefix + "_trial_4"
jobname5 = args.job_name_prefix + "_trial_5"

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
#header1_l6 = '#BSUB -R "select[hname!=node013]"' + '\n'
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

## Trial 2
# Header 2:
header2_l1 = "#! /usr/bin/env bash\n\n"
header2_l2 = "#BSUB -J " + jobname2 + '\n'
header2_l3 = "#BSUB -n " + str(args.nt) + '\n'
header2_l4 = '#BSUB -R "span[hosts=1]"' + '\n'
header2_l5 = '#BSUB -R "rusage[mem=' + str(mem) + ']"' + '\n'
#header2_l6 = '#BSUB -R "select[hname!=node013]"' + '\n'
header2_l7 = '#BSUB -M ' + str(mem_max) + '\n'
header2_l8 = '#BSUB -We ' + args.walle + '\n'
header2_l9 = '#BSUB -W ' + args.wallm + '\n'
header2_l10 = '#BSUB -o ' + synergy_project_dir + 'trial_2/' + jobname2 + '.out' + '\n'
header2_l11 = '#BSUB -e ' + synergy_project_dir + 'trial_2/' + jobname2 + '.err' + '\n'

# header2 = header2_l1 + header2_l2 + header2_l3 + header2_l4 + header2_l5 + header2_l6 + header2_l7 + header2_l8 + header2_l9 + header2_l10 + header2_l11 + '\n'
header2 = header2_l1 + header2_l2 + header2_l3 + header2_l4 + header2_l5 + header2_l7 + header2_l8 + header2_l9 + header2_l10 + header2_l11 + '\n'


# BEAST command
beast_cmd2 = '/home/cizydorczyk/BEASTv1.10.4/bin/beast -seed ' + str(seed2) + ' -working ' + synergy_project_dir + 'trial_2/' + args.input_xml.split('/')[-1]

# Write job script:
outfile_2 = args.project_dir + model + '/' + 'trial_2/' + jobname2 + '.sh'
with open(outfile_2, 'w') as outfile2:
    to_write2 = header2 + beast_cmd2
    outfile2.write(to_write2)

## Trial 3
# Header 3:
header3_l1 = "#! /usr/bin/env bash\n\n"
header3_l2 = "#BSUB -J " + jobname3 + '\n'
header3_l3 = "#BSUB -n " + str(args.nt) + '\n'
header3_l4 = '#BSUB -R "span[hosts=1]"' + '\n'
header3_l5 = '#BSUB -R "rusage[mem=' + str(mem) + ']"' + '\n'
#header3_l6 = '#BSUB -R "select[hname!=node013]"' + '\n'
header3_l7 = '#BSUB -M ' + str(mem_max) + '\n'
header3_l8 = '#BSUB -We ' + args.walle + '\n'
header3_l9 = '#BSUB -W ' + args.wallm + '\n'
header3_l10 = '#BSUB -o ' + synergy_project_dir + 'trial_3/' + jobname3 + '.out' + '\n'
header3_l11 = '#BSUB -e ' + synergy_project_dir + 'trial_3/' + jobname3 + '.err' + '\n'

# header3 = header3_l1 + header3_l2 + header3_l3 + header3_l4 + header3_l5 + header3_l6 + header3_l7 + header3_l8 + header3_l9 + header3_l10 + header3_l11 + '\n'
header3 = header3_l1 + header3_l2 + header3_l3 + header3_l4 + header3_l5 + header3_l7 + header3_l8 + header3_l9 + header3_l10 + header3_l11 + '\n'


# BEAST command
beast_cmd3 = '/home/cizydorczyk/BEASTv1.10.4/bin/beast -seed ' + str(seed3) + ' -working ' + synergy_project_dir + 'trial_3/' + args.input_xml.split('/')[-1]

# Write job script:
outfile_3 = args.project_dir + model + '/trial_3/' + jobname3 + '.sh'
with open(outfile_3, 'w') as outfile3:
    to_write3 = header3 + beast_cmd3
    outfile3.write(to_write3)

## Trial 4
# Header 4:
header4_l1 = "#! /usr/bin/env bash\n\n"
header4_l2 = "#BSUB -J " + jobname4 + '\n'
header4_l3 = "#BSUB -n " + str(args.nt) + '\n'
header4_l4 = '#BSUB -R "span[hosts=1]"' + '\n'
header4_l5 = '#BSUB -R "rusage[mem=' + str(mem) + ']"' + '\n'
#header3_l6 = '#BSUB -R "select[hname!=node013]"' + '\n'
header4_l7 = '#BSUB -M ' + str(mem_max) + '\n'
header4_l8 = '#BSUB -We ' + args.walle + '\n'
header4_l9 = '#BSUB -W ' + args.wallm + '\n'
header4_l10 = '#BSUB -o ' + synergy_project_dir + 'trial_4/' + jobname4 + '.out' + '\n'
header4_l11 = '#BSUB -e ' + synergy_project_dir + 'trial_4/' + jobname4 + '.err' + '\n'

# header3 = header3_l1 + header3_l2 + header3_l3 + header3_l4 + header3_l5 + header3_l6 + header3_l7 + header3_l8 + header3_l9 + header3_l10 + header3_l11 + '\n'
header4 = header4_l1 + header4_l2 + header4_l3 + header4_l4 + header4_l5 + header4_l7 + header4_l8 + header4_l9 + header4_l10 + header4_l11 + '\n'


# BEAST command
beast_cmd4 = '/home/cizydorczyk/BEASTv1.10.4/bin/beast -seed ' + str(seed4) + ' -working ' + synergy_project_dir + 'trial_4/' + args.input_xml.split('/')[-1]

# Write job script:
outfile_4 = args.project_dir + model + '/trial_4/' + jobname4 + '.sh'
with open(outfile_4, 'w') as outfile4:
    to_write4 = header4 + beast_cmd4
    outfile4.write(to_write4)

## Trial 5
# Header 5:
header5_l1 = "#! /usr/bin/env bash\n\n"
header5_l2 = "#BSUB -J " + jobname5 + '\n'
header5_l3 = "#BSUB -n " + str(args.nt) + '\n'
header5_l4 = '#BSUB -R "span[hosts=1]"' + '\n'
header5_l5 = '#BSUB -R "rusage[mem=' + str(mem) + ']"' + '\n'
#header3_l6 = '#BSUB -R "select[hname!=node013]"' + '\n'
header5_l7 = '#BSUB -M ' + str(mem_max) + '\n'
header5_l8 = '#BSUB -We ' + args.walle + '\n'
header5_l9 = '#BSUB -W ' + args.wallm + '\n'
header5_l10 = '#BSUB -o ' + synergy_project_dir + 'trial_5/' + jobname5 + '.out' + '\n'
header5_l11 = '#BSUB -e ' + synergy_project_dir + 'trial_5/' + jobname5 + '.err' + '\n'

# header3 = header3_l1 + header3_l2 + header3_l3 + header3_l4 + header3_l5 + header3_l6 + header3_l7 + header3_l8 + header3_l9 + header3_l10 + header3_l11 + '\n'
header5 = header5_l1 + header5_l2 + header5_l3 + header5_l4 + header5_l5 + header5_l7 + header5_l8 + header5_l9 + header5_l10 + header5_l11 + '\n'


# BEAST command
beast_cmd5 = '/home/cizydorczyk/BEASTv1.10.4/bin/beast -seed ' + str(seed5) + ' -working ' + synergy_project_dir + 'trial_5/' + args.input_xml.split('/')[-1]

# Write job script:
outfile_5 = args.project_dir + model + '/trial_5/' + jobname5 + '.sh'
with open(outfile_5, 'w') as outfile5:
    to_write5 = header5 + beast_cmd5
    outfile5.write(to_write5)

# Copy XML to trial directories
import subprocess
subprocess.run(['cp', args.input_xml, args.project_dir + model + '/' + 'trial_1/'])
subprocess.run(['cp', args.input_xml, args.project_dir + model + '/' + 'trial_2/'])
subprocess.run(['cp', args.input_xml, args.project_dir + model + '/' + 'trial_3/'])
subprocess.run(['cp', args.input_xml, args.project_dir + model + '/' + 'trial_4/'])
subprocess.run(['cp', args.input_xml, args.project_dir + model + '/' + 'trial_5/'])
