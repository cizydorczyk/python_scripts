import os

def CreateQuastJobs(quast_dir, unicycler_dir, jobs_dir, isolate_list, nt, walle, wallm, mem, mem_max, assembler):
    """Function to create a job to run QUAST on a set of de novo assemblies."""

    # Header:
    header1_l1 = "#! /usr/bin/env bash\n\n"
    header1_l2 = "#BSUB -J " + "raw_quast_job" + '\n'
    header1_l3 = "#BSUB -n " + str(nt) + '\n'
    header1_l4 = '#BSUB -R "span[hosts=1]"' + '\n'
    header1_l5 = '#BSUB -R "rusage[mem=' + str(mem) + ']"' + '\n'
    header1_l7 = '#BSUB -M ' + str(mem_max) + '\n'
    header1_l8 = '#BSUB -We ' + walle + '\n'
    header1_l9 = '#BSUB -W ' + wallm + '\n'
    header1_l10 = '#BSUB -o ' + quast_dir + "/" + "raw_quast" + '.out' + '\n'
    header1_l11 = '#BSUB -e ' + quast_dir + "/" + "raw_quast" + '.err' + '\n'

    header1 = header1_l1 + header1_l2 + header1_l3 + header1_l4 + header1_l5 + \
              header1_l7 + header1_l8 + header1_l9 + header1_l10 +\
              header1_l11 + '\n'

    isolate_assemblies = []
    for isolate in isolate_list:
        if assembler == "unicycler":
            isolate_assembly = unicycler_dir + "/" + isolate + "/" + isolate + "_assembly.fasta"
            isolate_assemblies.append(isolate_assembly)
        else:
            raise Exception("only unicycler is currently supported as an assembler.")

    # QUAST cmd:
    quast_cmd = "quast.py -o " + quast_dir + " -t " + str(nt) + " " + " ".join(isolate_assemblies)
    #print(quast_cmd)

    # Create QUAST job:
    raw_quast_job_file = jobs_dir + "/raw_quast_job.sh"
    with open(raw_quast_job_file, 'w') as outfile:
        outfile.write(header1 + quast_cmd)


# Sample function call:
#isolate_list = ["1", "2", "3"]
#CreateQuastJobs("quast_dir", "unicycler_dir", "jobs_dir", isolate_list, "nt", "walle", "wallm", "mem", "mem_max", "unicycler")
