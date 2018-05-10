from sys import argv
import itertools

script, isolate_list, outputdir = argv

isolates = []
with open(isolate_list, 'r') as infile1:
	for line in infile1:
		isolates.append(line.strip())


for i in isolates:
    header = "#!/bin/bash\n#SBATCH --nodes=1\n#SBATCH --cpus-per-task=40\n#SBATCH --time=1:00:00\n#SBATCH --job-name " + i + ".samtools.sort.job\n"
    output_file = outputdir + i + ".new.pseudo.bwa.job"
    with open(output_file, 'w') as outfile1:
        outfile1.write(header + "\n/home/d/dguttman/cizydor/Software/samtools-1.8/samtools sort -@ 40 -o /scratch/d/dguttman/cizydor/J10_sorted_bam/" + i + "_sorted.bam -O bam /scratch/d/dguttman/cizydor/J10_bwa_output/" + i + "_aln.sam")