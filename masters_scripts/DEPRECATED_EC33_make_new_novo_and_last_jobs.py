from sys import argv
import os.path

script, input_job, output_novoalign_jobs, output_samtools_fixmate_jobs = argv

lines = []

with open(input_job, 'r') as infile1:
    for line in infile1:
        lines.append(line.strip())

if os.path.isfile(output_novoalign_jobs):
    with open(output_novoalign_jobs, 'a+') as outfile1:
        outfile1.write('\n' + lines[1])
elif not os.path.isfile(output_novoalign_jobs):
    with open(output_novoalign_jobs, 'w') as outfile1:
        outfile1.write(lines[1])

if os.path.isfile(output_samtools_fixmate_jobs):
    with open(output_samtools_fixmate_jobs, 'a+') as outfile2:
        outfile2.write('\n' + lines[3])
elif not os.path.isfile(output_samtools_fixmate_jobs):
    with open(output_samtools_fixmate_jobs, 'w') as outfile2:
        outfile2.write(lines[3])