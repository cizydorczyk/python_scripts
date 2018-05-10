from sys import argv

script, isolate_list, outputdir = argv

isolates = []
with open(isolate_list, 'r') as infile1:
	for line in infile1:
		isolates.append(line.strip())

for i in isolates:
	header = "#!/bin/bash\n#SBATCH --nodes=1\n#SBATCH --cpus-per-task=40\n#SBATCH --time=4:00:00\n#SBATCH --job-name " + i + "_novoalign\n"
	output_file = outputdir + i + ".novoaln.job"
	with open(output_file, 'w') as outfile1:
		outfile1.write(header + "/home/d/dguttmancizydor/Software/novocraft/novoalign -o SAM -t 15,4 --hlimit 8 -p 5,20 --trim3hp -c 40 -f /scratch/d/dguttman/cizydor/ecoli/EC33_trimmed_fastq/" + i + "_1.fq /scratch/d/dguttman/cizydor/ecoli/EC33_trimmed_fastq/" + i + "_2.fq -d /scratch/d/dguttman/cizydor/ecoli/EC33_snp_calling/references/k12mg1655_NOVOALIGN > /scratch/d/dguttman/cizydor/ecoli/EC33_snp_calling/alignments/NOVOALIGN-101/r.sam")
