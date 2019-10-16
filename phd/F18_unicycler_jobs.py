import os

def CreateUnicyclerJobs(unicycler_dir, fastq_dir, jobs_dir, isolate, nt, walle, wallm, depth_filter, fastq_endings, gz, mem, mem_max):
    """Function to create Unicycler jobs for Synergy."""

    # Create directory for isolate output:
    assembly_dir = os.path.join(unicycler_dir, isolate)

    # Header:
    header1_l1 = "#! /usr/bin/env bash\n\n"
    header1_l2 = "#BSUB -J " + isolate + "_unic_job" + '\n'
    header1_l3 = "#BSUB -n " + str(nt) + '\n'
    header1_l4 = '#BSUB -R "span[hosts=1]"' + '\n'
    header1_l5 = '#BSUB -R "rusage[mem=' + str(mem) + ']"' + '\n'
    header1_l7 = '#BSUB -M ' + str(mem_max) + '\n'
    header1_l8 = '#BSUB -We ' + walle + '\n'
    header1_l9 = '#BSUB -W ' + wallm + '\n'
    header1_l10 = '#BSUB -o ' + assembly_dir + "/" + isolate + '.out' + '\n'
    header1_l11 = '#BSUB -e ' + assembly_dir + "/" + isolate + '.err' + '\n'

    header1 = header1_l1 + header1_l2 + header1_l3 + header1_l4 + header1_l5 + \
              header1_l7 + header1_l8 + header1_l9 + header1_l10 +\
              header1_l11 + '\n'

    # Fastq file names:
    fq_endings = fastq_endings.split(",")

    if gz == 'yes':
        r1 = fastq_dir + isolate + fq_endings[0] + ".gz"
        r2 = fastq_dir + isolate + fq_endings[1] + ".gz"

    elif gz == 'no':
        r1 = fastq_dir + isolate + fq_endings[0]
        r2 = fastq_dir + isolate + fq_endings[1]

    # Unicycler command:
    unicycler_cmd = "unicycler -1 " + r1 + " -2 " + r2 + " -o " + assembly_dir + " -t " + str(nt) + " --depth_filter " + str(depth_filter) + "\n\n"
    #print(header1 + unicycler_cmd)

    # Rename assembly to include isolate name cmd:
    rename_cmd = "mv " + assembly_dir + "/assembly.fasta " + assembly_dir + "/" + isolate + "_assembly.fasta"

    # Create job file:
    job_file = jobs_dir + "/" + isolate + ".sh"
    with open(job_file, 'w') as outfile:
        outfile.write(header1 + unicycler_cmd + rename_cmd)

# Sample function call:
#create_unicycler_jobs("unicycler_dir/", "fastq_dir/", "jobs_dir/", "isolate", "nt", "walle", "wallm", "depth_filter", "_1.fastq,_2.fastq", "yes", "mem", "mem_max")
