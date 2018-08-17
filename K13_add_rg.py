from sys import argv

script, input_isolate_list, output_scinet_script = argv

isolates = []
with open(input_isolate_list, 'r') as infile1:
    for line in infile1:
        isolates.append(line.strip())

to_write = []
for i in isolates:
    to_write.append("java -jar $EBROOTPICARD/picard.jar AddOrReplaceReadGroups I=/scratch/d/dguttman/cizydor/K7_snp_calling/alignments/BWA-" + i + "/r.sorted.bam O=/scratch/d/dguttman/cizydor/K13_indel_calling/bam_files/" + i + "_sorted_rg.bam RGID=HHMW3BGXY.1 RGPU=HHMW3BGXY.1." + i + " RGSM=" + i + " RGPL=ILLUMINA RGLB=" + i + " CREATE_INDEX=True &")

with open(output_scinet_script, 'w') as outfile1:
    to_write2 = '\n'.join(to_write)
    outfile1.write(to_write2)