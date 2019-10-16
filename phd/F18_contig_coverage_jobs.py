import os

def CreateContigCoverageJobs(contig_dir, unicycler_dir, fastq_dir, jobs_dir, isolate, nt, walle, wallm, mem, mem_max, gz, fastq_endings, assembler):
    """Function to create a job to get average contig coverage for a de novo assembly."""

    # Header:
    header1_l1 = "#! /usr/bin/env bash\n\n"
    header1_l2 = "#BSUB -J " + "cc_job" + '\n'
    header1_l3 = "#BSUB -n " + str(nt) + '\n'
    header1_l4 = '#BSUB -R "span[hosts=1]"' + '\n'
    header1_l5 = '#BSUB -R "rusage[mem=' + str(mem) + ']"' + '\n'
    header1_l7 = '#BSUB -M ' + str(mem_max) + '\n'
    header1_l8 = '#BSUB -We ' + walle + '\n'
    header1_l9 = '#BSUB -W ' + wallm + '\n'
    header1_l10 = '#BSUB -o ' + contig_dir + "/" + isolate + '_contigcov.out' + '\n'
    header1_l11 = '#BSUB -e ' + contig_dir + "/" + isolate + '_contigcov.err' + '\n'

    header1 = header1_l1 + header1_l2 + header1_l3 + header1_l4 + header1_l5 + \
              header1_l7 + header1_l8 + header1_l9 + header1_l10 +\
              header1_l11 + '\n'

    # contig coverage commands:

    isolate_assembly_dir = os.path.join(unicycler_dir, isolate)
    bwa_unsorted_sam = contig_dir + "/" + isolate + "_unsorted.sam"
    bwa_unsorted_bam = contig_dir + "/" + isolate + "_unsorted.bam"
    bwa_sorted_bam = contig_dir + "/" + isolate + "_sorted.bam"
    samtools_depth_file = contig_dir + "/" + isolate + "_depths.txt"

    fq_endings = fastq_endings.split(",")

    if gz == 'yes':
        r1 = fastq_dir + isolate + fq_endings[0] + ".gz"
        r2 = fastq_dir + isolate + fq_endings[1] + ".gz"

    elif gz == 'no':
        r1 = fastq_dir + isolate + fq_endings[0]
        r2 = fastq_dir + isolate + fq_endings[1]

    cd_cmd = "cd " + isolate_assembly_dir

    if assembler == "unicycler":
        bwa_index_cmd = "bwa index " + isolate + "_assembly.fasta"
        bwa_mem_cmd = "bwa mem -t " + str(nt) + " " + isolate + "_assembly.fasta " + r1 + " " + r2 + " > " + bwa_unsorted_sam
    else:
        raise Exception("Assembler must be 'unicycler' at this time.")

    samtools_view_cmd = "samtools view -@ " + str(nt) + " -O BAM -o " + bwa_unsorted_bam + " " + bwa_unsorted_sam
    samtools_sort_cmd = "samtools sort -@ " + str(nt) + " " + bwa_unsorted_bam + " > " + bwa_sorted_bam
    samtools_depth_cmd = "samtools depth -a " + bwa_sorted_bam + " > " + samtools_depth_file

    cleanup_cmd = "rm " + bwa_unsorted_bam + " " + bwa_unsorted_sam + " " + bwa_sorted_bam

    to_write = header1 + "\n" + cd_cmd + "\n\n" + bwa_index_cmd + "\n\n" + bwa_mem_cmd + "\n\n" + samtools_view_cmd + "\n\n" + samtools_sort_cmd + "\n\n" + samtools_depth_cmd + "\n\n" + cleanup_cmd

    output_job_file = jobs_dir + "/" + isolate + "_contig_coverage_job.sh"
    with open(output_job_file, 'w') as outfile:
        outfile.write(to_write)

# Sample function call:
#CreateContigCoverageJobs("contig_dir", "unicycler_dir", "fastq_dir", "jobs_dir", "isolate", "nt", "walle", "wallm", "mem", "mem_max", "yes", "_1.fastq,_2.fastq", "unicycler")
