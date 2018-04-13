from sys import argv
import itertools

script, isolate_list, outputdir = argv

isolates = []
with open(isolate_list, 'r') as infile1:
	for line in infile1:
		isolates.append(line.strip())


for i in isolates:
	header = "#!/bin/bash\n#SBATCH --nodes=1\n#SBATCH --cpus-per-task=40\n#SBATCH --time=12:00:00\n#SBATCH --job-name " + i + "_assembly_blast_job\n"
	output_file = outputdir + i + ".orig.pseudo.job"
	with open(output_file, 'w') as outfile1:
		outfile1.write(header + "\nexport BLASTDB=/home/d/dguttman/cizydor/blastdb\n" + "\n/home/d/dguttman/cizydor/Software/ncbi-blast-2.7.1+/bin/blastn -db /home/d/dguttman/cizydor/blastdb/nt -gilist /home/d/dguttman/cizydor/pae_sequence.gi -query /scratch/d/dguttman/cizydor/I4_assemblies/" + i + "_assembly.fa -evalue 0.00001 -out /scratch/d/dguttman/cizydor/I4_blast_results/" + i + '_raw_blast_output.txt -outfmt "7 qseqid sseqid pident length evalue" -num_threads 40')