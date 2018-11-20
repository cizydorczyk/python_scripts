from sys import argv

script, isolate_list, outputdir = argv

isolates = []
with open(isolate_list, 'r') as infile1:
	for line in infile1:
		isolates.append(line.strip())

for i in isolates:
	header = "#!/bin/bash\n"#SBATCH --nodes=1\n#SBATCH --cpus-per-task=80\n#SBATCH --time=1:00:00\n#SBATCH --job-name " + i + "_last.novo.mpileup.call\n"
	output_file = outputdir + i + ".last_novo_mpileup_call.job.sh"
	with open(output_file, 'w') as outfile1:
		outfile1.write(header + "/home/d/dguttman/diazcaba/Software/bcftools-1.8/bcftools mpileup --threads 80 -Ob -o /scratch/d/dguttman/cizydor/K7_snp_calling/alignments/LAST-" + i + "/r.bcf -f /scratch/d/dguttman/cizydor/K7_snp_calling/references/pao1.fasta /scratch/d/dguttman/cizydor/K7_snp_calling/alignments/LAST-" + i + "/r.sorted.bam" + '\n' + "/home/d/dguttman/diazcaba/Software/bcftools-1.8/bcftools call --threads 80 --ploidy 1 -vmO v -o /scratch/d/dguttman/cizydor/K7_snp_calling/alignments/LAST-" + i + "/r.vcf /scratch/d/dguttman/cizydor/K7_snp_calling/alignments/LAST-" + i + "/r.bcf" + '\n' + "/home/d/dguttman/diazcaba/Software/bcftools-1.8/bcftools mpileup --threads 80 -Ob -o /scratch/d/dguttman/cizydor/K7_snp_calling/alignments/NOVOALIGN-" + i + "/r.bcf -f /scratch/d/dguttman/cizydor/K7_snp_calling/references/pao1.fasta /scratch/d/dguttman/cizydor/K7_snp_calling/alignments/NOVOALIGN-" + i + "/r.sorted.bam" + '\n' + "/home/d/dguttman/diazcaba/Software/bcftools-1.8/bcftools call --threads 80 --ploidy 1 -vmO v -o /scratch/d/dguttman/cizydor/K7_snp_calling/alignments/NOVOALIGN-" + i + "/r.vcf /scratch/d/dguttman/cizydor/K7_snp_calling/alignments/NOVOALIGN-" + i + "/r.bcf")