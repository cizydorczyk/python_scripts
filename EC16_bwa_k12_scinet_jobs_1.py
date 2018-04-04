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
		outfile1.write("\n" + "\n" + "module load bwakit/0.7.13" + "\n")
		outfile1.write('module load samtools/1.3.1\n')
		outfile1.write("\n" + "bwa mem -t 8 /scratch/d/dguttman/cizydor/ecoli_project/EC16_k12_reference/k12mg1655.referece.fasta /scratch/d/dguttman/cizydor/ecoli_project/ecoli_trim5/paired/" + i + "_R1.fastq /scratch/d/dguttman/cizydor/ecoli_project/ecoli_trim5/paired/" + i + "_R2.fastq > /scratch/d/dguttman/cizydor/ecoli_project/EC16_k12_aln/" + i + "_bwa_aln.sam")
		outfile1.write('\n\n' + 'samtools sort -@ 8 -o /scratch/d/dguttman/cizydor/ecoli_project/EC16_k12_aln/' + i + '_sorted_aln.bam -O bam /scratch/d/dguttman/cizydor/ecoli_project/EC16_k12_aln/' + i + '_bwa_aln.sam')
		outfile1.write('\n\ncd /scratch/d/dguttman/cizydor/ecoli_project/EC16_k12_aln/\n\nsamtools index /scratch/d/dguttman/cizydor/ecoli_project/EC16_k12_aln/' + i + '_sorted_aln.bam')

