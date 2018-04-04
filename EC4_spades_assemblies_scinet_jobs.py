from sys import argv
import itertools

script, isolate_list, outputdir = argv

isolates = []
with open(isolate_list, 'r') as infile1:
	for line in infile1:
		isolates.append(line.strip())


for i in isolates:
	header = "#!/bin/bash\n#PBS -l nodes=1:ppn=8\n#PBS -l walltime=00:30:00\n#PBS -N " + i + "_ecoli\n#PBS -A zke-503-ab"
	output_file = outputdir + i + ".ecoli.spades.job"
	with open(output_file, 'w') as outfile1:
		outfile1.write(header)
		# outfile1.write("\n" + "\n" + "cd /scratch/d/dguttman/cizydor/ecoli_assemblies/" + "\n")
		# outfile1.write("\n" + "mkdir " + i + "_assembly" + "\n" + "cd " + i + "_assembly" + "\n")
		outfile1.write("\n" + "\n" + "module load spades/spades3.11.1" + "\n")
		outfile1.write("\n" + "spades.py --careful -1 /scratch/d/dguttman/cizydor/ecoli_fastq/" + i + "_R1.fastq -2 /scratch/d/dguttman/cizydor/ecoli_fastq/" + i + "_R2.fastq -o /scratch/d/dguttman/cizydor/ecoli_assemblies/" + i + " -t 8" +  "\n")



#!/bin/bash
#PBS -l nodes=1:ppn=8
#PBS -l walltime=48:00:00
#PBS -N 36_pae
#PBS -A zke-503-ab

# do spades.py --careful -1 ~/Data/primary_project_3/fastq_files/G37_good_reads_rarefied/$i"_R1.fastq" -2 ~/Data/primary_project_3/fastq_files/G37_good_reads_rarefied/$i"_R2.fastq" -o ~/Data/primary_project_3/de_novo_assemblies/H51_spades_155_rarefied/$i/ -t 8 -m 50; done