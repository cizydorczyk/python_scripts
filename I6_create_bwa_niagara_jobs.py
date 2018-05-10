from sys import argv
import itertools

script, isolate_list, outputdir = argv

isolates = []
with open(isolate_list, 'r') as infile1:
	for line in infile1:
		isolates.append(line.strip())


for i in isolates:
    header = "#!/bin/bash\n#SBATCH --nodes=1\n#SBATCH --cpus-per-task=40\n#SBATCH --time=1:00:00\n#SBATCH --job-name " + i + "_assembly_bwa_job\n"
    output_file = outputdir + i + ".orig.pseudo.bwa.job"
    with open(output_file, 'w') as outfile1:
        outfile1.write(header + "\n/home/d/dguttman/cizydor/Software/bwa-0.7.17/bwa index /scratch/d/dguttman/cizydor/J10_good_contigs/" + i + "_good_assembly.fa\n\n/home/d/dguttman/cizydor/Software/bwa-0.7.17/bwa mem -t 40 /scratch/d/dguttman/cizydor/J10_good_contigs/" + i + "_good_assembly.fa /scratch/d/dguttman/cizydor/J10_reads/" + i + "_R1.fq.gz /scratch/d/dguttman/cizydor/J10_reads/" + i + "_R2.fq.gz > /scratch/d/dguttman/cizydor/J10_bwa_output/" + i + "_aln.sam")