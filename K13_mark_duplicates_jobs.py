from sys import argv

script, input_isolate_list, output_scinet_script = argv

isolates = []
with open(input_isolate_list, 'r') as infile1:
    for line in infile1:
        isolates.append(line.strip())

to_write = []
for i in isolates:
    to_write.append("java -jar $EBROOTPICARD/picard.jar MarkDuplicates I=/scratch/d/dguttman/cizydor/K13_indel_calling/bam_files/" + i + "_sorted_rg.bam O=/scratch/d/dguttman/cizydor/K13_indel_calling/bam_files/" + i + "_sorted_rg_mrkdup.bam M=/scratch/d/dguttman/cizydor/K13_indel_calling/bam_files/" + i + "_mrkdup_metrics.txt &")

with open(output_scinet_script, 'w') as outfile1:
    to_write2 = '\n'.join(to_write)
    outfile1.write(to_write2)