import os

def CreateBlastJob(s_blast_raw_output_dir, s_assemblies_dir, blast_jobs_dir, isolate, nt, walle, wallm, mem, mem_max, evalue):
    """Function to create a job to run BLAST on a set of de novo assemblies."""

    # Header:
    header1_l1 = "#! /usr/bin/env bash\n\n"
    header1_l2 = "#BSUB -J " + "blast_job" + '\n'
    header1_l3 = "#BSUB -n " + str(nt) + '\n'
    header1_l4 = '#BSUB -R "span[hosts=1]"' + '\n'
    header1_l5 = '#BSUB -R "rusage[mem=' + str(mem) + ']"' + '\n'
    header1_l7 = '#BSUB -M ' + str(mem_max) + '\n'
    header1_l8 = '#BSUB -We ' + walle + '\n'
    header1_l9 = '#BSUB -W ' + wallm + '\n'
    header1_l10 = '#BSUB -o ' + s_blast_raw_output_dir + "/" + isolate + '_blast.out' + '\n'
    header1_l11 = '#BSUB -e ' + s_blast_raw_output_dir + "/" + isolate + '_blast.err' + '\n'

    header1 = header1_l1 + header1_l2 + header1_l3 + header1_l4 + header1_l5 + \
              header1_l7 + header1_l8 + header1_l9 + header1_l10 +\
              header1_l11 + '\n'

    # Blast cmd:
    export_db_cmd = "export BLASTDB=/home/cizydorczyk/blastdb"

    isolate_assembly = s_assemblies_dir + "/" + isolate + "/" + isolate + "_assembly.fasta"
    isolate_raw_blast_output = s_blast_raw_output_dir + "/" + isolate + '_raw_blast_output.txt'

    blast_cmd = "blastn -db /home/cizydorczyk/blastdb/nt -query " + isolate_assembly + " -evalue " + str(evalue) + " -out " + isolate_raw_blast_output + ' -outfmt "7 qseqid sseqid pident length evalue sscinames" -num_threads ' + str(nt)

    # Write jobs file:
    to_write = header1 + "\n" + export_db_cmd + "\n\n" + blast_cmd
    blast_job_file = blast_jobs_dir + "/" + isolate + "_blast_job.sh"

    # Uncomment line below to run sample function call:
    #print(to_write)

    # Comment out 2 lines below to run sample function call:
    with open(blast_job_file, "w") as outfile:
        outfile.write(to_write)

# Sample function call:
# CreateBlastJob("s_blast_raw_output_dir", "s_unicycler_assemblies_dir", "blast_jobs_dir", "isolate", "nt", "walle", "wallm", "mem", "mem_max", "evalue")
