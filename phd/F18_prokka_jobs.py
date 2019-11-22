import os

def CreateProkkaJob(s_prokka_dir, s_final_assemblies_dir, prokka_jobs_dir, isolate, nt, walle, wallm, mem, mem_max, proteins_file, kingdom, genus, use_genus, species, gcode, rnammer, mincontiglen):
    """Function to create a job to run BLAST on a set of de novo assemblies."""

    # Header:
    header1_l1 = "#! /usr/bin/env bash\n\n"
    header1_l2 = "#BSUB -J " + "prokka_job" + '\n'
    header1_l3 = "#BSUB -n " + str(nt) + '\n'
    header1_l4 = '#BSUB -R "span[hosts=1]"' + '\n'
    header1_l5 = '#BSUB -R "rusage[mem=' + str(mem) + ']"' + '\n'
    header1_l7 = '#BSUB -M ' + str(mem_max) + '\n'
    header1_l8 = '#BSUB -We ' + walle + '\n'
    header1_l9 = '#BSUB -W ' + wallm + '\n'
    header1_l10 = '#BSUB -o ' + s_prokka_dir + "/" + isolate + '_prokka.out' + '\n'
    header1_l11 = '#BSUB -e ' + s_prokka_dir + "/" + isolate + '_prokka.err' + '\n'

    header1 = header1_l1 + header1_l2 + header1_l3 + header1_l4 + header1_l5 + \
              header1_l7 + header1_l8 + header1_l9 + header1_l10 +\
              header1_l11 + '\n'

    # Prokka input/output files:
    isolate_assembly = s_final_assemblies_dir + "/" + isolate + "_filtered_ordered_assembly.fasta"
    output_dir = os.path.join(s_prokka_dir, isolate)

    # Set rnammer status:
    if rnammer == "True":
        rnammer_ = " --rnammer"
    elif rnammer == "False":
        rnammer_ = ""

    # Set use_genus status:
    if use_genus == "True":
        use_genus_ = " --use_genus"
    elif use_genus == "False":
        use_genus_ = ""

    # Set proteins file:
    if proteins_file != "":
        proteins_file_ = " --proteins " + proteins_file
    elif proteins_file == "":
        proteins_file_ = ""

    # Prokka cmd:
    prokka_cmd = "prokka --outdir " + output_dir + " --prefix " + isolate + " --kingdom " + kingdom + " --gcode " + str(gcode) + " --genus " + genus + " --species " + species + proteins_file_ + " --cpus " + str(nt) + rnammer_ + use_genus_ + " --rfam --addgenes --locustag " + isolate + " --mincontiglen " + str(mincontiglen) + " " + isolate_assembly


    # Write job file:
    output_file = prokka_jobs_dir + "/" + isolate + "_prokka_job.sh"

    to_write = header1 + prokka_cmd

    with open(output_file, 'w') as outfile:
        outfile.write(to_write)
