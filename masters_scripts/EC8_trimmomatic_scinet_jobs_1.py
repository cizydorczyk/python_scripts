from sys import argv

script, isolate_list, outputdir = argv

isolatelist = []
with open(isolate_list, 'r') as infile1:
	for line in infile1:
		isolatelist.append(line.strip())

for i in isolatelist:
	outhandle = outputdir + i + "_trimmomatic_scinet.job"
	with open(outhandle, 'w') as outfile1:
		header = "#!/bin/bash\n#PBS -l nodes=1:ppn=8\n#PBS -l walltime=00:30:00\n#PBS -N " + i + "_ecoli\n#PBS -A zke-503-ab"
		outfile1.write(header + '\n\nmodule load java\n\n')
		outfile1.write("java -jar /home/d/dguttman/cizydor/Software/Trimmomatic-0.36/trimmomatic-0.36.jar PE /scratch/d/dguttman/cizydor/ecoli_raw_fastq/" + i + "_R1.fastq.gz /scratch/d/dguttman/cizydor/ecoli_raw_fastq/" + i + "_R2.fastq.gz /scratch/d/dguttman/cizydor/ecoli_trim5/paired/" + i + "_R1.fastq /scratch/d/dguttman/cizydor/ecoli_trim5/unpaired/" + i + "_unpaired_R1.fastq /scratch/d/dguttman/cizydor/ecoli_trim5/paired/" + i + "_R2.fastq /scratch/d/dguttman/cizydor/ecoli_trim5/unpaired/" + i + "_unpaired_R2.fastq CROP:150 SLIDINGWINDOW:4:5 ILLUMINACLIP:/home/d/dguttman/cizydor/Software/Trimmomatic-0.36/adapters/NexteraPE-PE.fa:2:30:10 MINLEN:100")