from sys import argv
import itertools

script, isolate_list, outputdir = argv

isolates = []
with open(isolate_list, 'r') as infile1:
	for line in infile1:
		isolates.append(line.strip())


for i in isolates:
	header = "#!/bin/bash\n#SBATCH --nodes=1\n#SBATCH --cpus-per-task=40\n#SBATCH --time=4:00:00\n#SBATCH --job-name " + i + "_new_pseudo_assembly_blast\n"
	output_file = outputdir + i + ".spades.pseudo.job"
	with open(output_file, 'w') as outfile1:
		outfile1.write(header + "\n/home/d/dguttman/cizydor/Software/SPAdes-3.11.1-Linux/bin/spades.py --careful -1 /scratch/d/dguttman/cizydor/I14_good_fastq/" + i + "_1.fq.gz -2 /scratch/d/dguttman/cizydor/I14_good_fastq/" + i + "_2.fq.gz -o /scratch/d/dguttman/cizydor/I14_spades_assemblies/" + i + " -t 40")