from sys import argv

script, input_isolate_list, output_dir = argv

isolates = []
with open(input_isolate_list, 'r') as infile1:
    for line in infile1:
        isolates.append(line.strip())

to_write = []
for i in isolates:
    header = "#!/bin/bash\n#SBATCH --nodes=1\n#SBATCH --cpus-per-task=80\n#SBATCH --time=4:00:00\n#SBATCH --job-name " + i + "_gatk\n\n"
    modules = "module load CCEnv\nmodule load nixpkgs/16.09\nmodule load java\nmodule load gatk/3.8\nmodule load gcc/5.4.0\nmodule load intel/2016.4\nmodule load bcftools/1.5\n\n"
    haplotype_command = "java -jar $EBROOTGATK/GenomeAnalysisTK.jar -T HaplotypeCaller -R /scratch/d/dguttman/cizydor/K13_indel_calling/reference/pao1.fasta -I /scratch/d/dguttman/cizydor/K13_indel_calling/bam_files/" + i + "_sorted_rg_mrkdup.bam -o /scratch/d/dguttman/cizydor/K13_indel_calling/raw_vcf_files/" + i + "_raw.vcf -ploidy 1 -nct 80 -stand_call_conf 30\n\n"
    selectvar_command = "java -jar $EBROOTGATK/GenomeAnalysisTK.jar -T SelectVariants -R /scratch/d/dguttman/cizydor/K13_indel_calling/reference/pao1.fasta -nt 80 -V /scratch/d/dguttman/cizydor/K13_indel_calling/raw_vcf_files/" + i + "_raw.vcf -selectType INDEL -o /scratch/d/dguttman/cizydor/K13_indel_calling/raw_vcf_files/" + i +"_raw_indels.vcf\n\n"
    filter_command = "bcftools filter -i 'QUAL>=30 & DP>=20' -o /scratch/d/dguttman/cizydor/K13_indel_calling/filtered_indels_vcf/" + i + "_filtered_indels.vcf --threads 80 /scratch/d/dguttman/cizydor/K13_indel_calling/raw_vcf_files/" + i + "_raw_indels.vcf"

    output_file = output_dir + i + "_gatk.job"
    with open(output_file, 'w') as outfile1:
        outfile1.write(header + modules + haplotype_command + selectvar_command + filter_command)