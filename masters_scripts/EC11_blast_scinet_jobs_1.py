from sys import argv
import itertools

script, isolate_list, outputdir = argv

isolates = []
with open(isolate_list, 'r') as infile1:
	for line in infile1:
		isolates.append(line.strip())


for i in isolates:
	header = "#!/bin/bash\n#PBS -l nodes=1:ppn=8\n#PBS -l walltime=34:00:00\n#PBS -N " + i + "_ecoli\n#PBS -A zke-503-ab"
	output_file = outputdir + i + ".ecoli.clc.blast.job"
	with open(output_file, 'w') as outfile1:
		outfile1.write(header)
		outfile1.write("\n\nexport BLASTDB=$SCRATCH/blastdb\n\n")
		outfile1.write("\n" + "\n" + "module load blast" + "\n")
		outfile1.write("\n" + "blastn -db /scratch/d/dguttman/cizydor/blastdb/nt -query /scratch/d/dguttman/cizydor/ecoli_project/EC14_clc_assemblies/" + i + "_assembly.fa -evalue 0.00001 -out /scratch/d/dguttman/cizydor/ecoli_project/EC15_blast_results/" + i + '_raw_blast_output.txt -outfmt "7 qseqid sseqid pident length evalue sscinames" -num_threads 8')



#!/bin/bash
#PBS -l nodes=1:ppn=8
#PBS -l walltime=48:00:00
#PBS -N 36_pae
#PBS -A zke-503-ab