from sys import argv

script, isolate_list, outputdir = argv

isolates = []
with open(isolate_list, 'r') as infile1:
	for line in infile1:
		isolates.append(line.strip())

for i in isolates:
	header = "#!/bin/bash\n"
	output_file = outputdir + i + ".novoaln.job"
	with open(output_file, 'w') as outfile1:
		outfile1.write(header + "/home/conrad/Software/novocraft/novoalign -c 8 -o SAM -f /home/conrad/Data/ecoli/fastq_files/EC8_trimmomatic_paired/" + i + "_1.fq /home/conrad/Data/ecoli/fastq_files/EC8_trimmomatic_paired/" + i + "_2.fq -d /home/conrad/Data/ecoli/EC33_snp_calling/references/k12mg1655_NOVOALIGN > /home/conrad/Data/ecoli/EC33_snp_calling/alignments/NOVOALIGN-" + i + "/r.sam")
