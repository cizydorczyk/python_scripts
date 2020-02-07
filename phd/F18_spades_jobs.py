import os

def CreateSpadesJobs(spades_dir, fastq_dir, jobs_dir, isolate, nt, walle, wallm, fastq_endings, gz, mem, mem_max):
    """Function to create SPAdes jobs for Synergy."""

    # Create directory for isolate output:
    assembly_dir = os.path.join(spades_dir, isolate)

    # Header:
    header1_l1 = "#! /usr/bin/env bash\n\n"
    header1_l2 = "#BSUB -J " + isolate + "_spd_job" + '\n'
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

    # SPAdes command:
    spades_cmd = "spades.py --careful --cov-cutoff 10 -k 21,33,55,77,99,127 -o " + assembly_dir + " -1 " + r1 + " -2 " + r2 + " --threads " + str(nt) + "\n\n"

    # Rename assembly to include isolate name cmd:
    rename_cmd = "mv " + assembly_dir + "/scaffolds.fasta " + assembly_dir + "/" + isolate + "_assembly.fasta"

    # Create job file:
    job_file = jobs_dir + "/" + isolate + ".sh"
    with open(job_file, 'w') as outfile:
        outfile.write(header1 + spades_cmd + rename_cmd)

# Sample function call:
#create_unicycler_jobs("unicycler_dir/", "fastq_dir/", "jobs_dir/", "isolate", "nt", "walle", "wallm", "depth_filter", "_1.fastq,_2.fastq", "yes", "mem", "mem_max")
