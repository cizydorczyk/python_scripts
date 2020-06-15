import os
import argparse

parser = argparse.ArgumentParser()

parser.add_argument("--nt", help="number of threads to use (<= 56)")
parser.add_argument("--walle", help="walltime estimate HH:MM")
parser.add_argument("--wallm", help="walltime max HH:MM")
parser.add_argument("--synergy_output_dir", help="synergy output directory; individual isolate directories will be "
                                                 "created in this directory")
parser.add_argument("--isolate", help="isolate name/number")
parser.add_argument("--input_assembly", help="full synergy path to input assembly")
parser.add_argument("--job_file", help="local job file; upload to synergy and used to run job")

parser.add_argument("--kingdom", help="prokka kingdom")
parser.add_argument("--genus", help="prokka genus")
parser.add_argument("--gcode", help="prokka gcode")
parser.add_argument("--species", help="prokka species")
parser.add_argument("--rfam", help="prokka rfam; include to use rfam", action="store_true")
parser.add_argument("--rnammer", help="prokka rnammer; include to use rnammer instead of barrnap", action="store_true")
parser.add_argument("--addgenes", help="prokka addgenes; include to use addgenes", action="store_true")
parser.add_argument("--mincontiglen", help="min contig len to annotate, default = 200", type=int, default=200)
parser.add_argument("--prodigaltf", help="prokka prodigaltf; provide file to use this argument")
parser.add_argument("--proteins", help="prokka proteins")

args = parser.parse_args()

# Set memory requirements
mem = int(args.nt) * 4000
mem_max = mem + 2000
snippy_ram = str(int(((int(args.nt) * 4000) - 1000)/1000))

# Header:
header1_l1 = "#! /usr/bin/env bash\n\n"
header1_l2 = "#BSUB -J " + "j" + args.isolate + '\n'
header1_l3 = "#BSUB -n " + str(args.nt) + '\n'
header1_l4 = '#BSUB -R "span[hosts=1]"' + '\n'
header1_l5 = '#BSUB -R "rusage[mem=' + str(mem) + ']"' + '\n'
header1_l7 = '#BSUB -M ' + str(mem_max) + '\n'
header1_l8 = '#BSUB -We ' + args.walle + '\n'
header1_l9 = '#BSUB -W ' + args.wallm + '\n'
header1_l10 = '#BSUB -o ' + args.synergy_output_dir + "/" + args.isolate + "/" + args.isolate + '_prokka.out' + '\n'
header1_l11 = '#BSUB -e ' + args.synergy_output_dir + "/" + args.isolate + "/" + args.isolate + '_prokka.err' + '\n'

header1 = header1_l1 + header1_l2 + header1_l3 + header1_l4 + header1_l5 + \
          header1_l7 + header1_l8 + header1_l9 + header1_l10 +\
          header1_l11 + '\n'

# Set optional prokka options:
if args.rfam:
    rfam = "--rfam"
else:
    rfam = ""

if args.rnammer:
    rnammer = "--rnammer"
else:
    rnammer = ""

if args.addgenes:
    addgenes = "--addgenes"
else:
    addgenes = ""

if args.proteins:
    proteins_arg = " --proteins " + args.proteins
else:
    proteins_arg = ""

if args.prodigaltf:
    prodigaltf_arg = " --prodigaltf " + args.prodigaltf
else:
    prodigaltf_arg = ""

# Create prokka cmd:
isolate_output_dir = os.path.join(args.synergy_output_dir, args.isolate)

prokka_cmd = "prokka --outdir " + isolate_output_dir + " --prefix " + args.isolate + " --kingdom " + args.kingdom + \
             " --gcode " + str(args.gcode) + " --genus " + args.genus + " --species " + args.species + " --cpus " + \
             str(args.nt) + " " + rnammer + " " + rfam + " " + addgenes + " --locustag " + args.isolate + " --mincontiglen "\
             + str(args.mincontiglen) + prodigaltf_arg + proteins_arg + " " + args.input_assembly

# print(prokka_cmd)

# Write output file:
to_write = header1 + prokka_cmd
with open(args.job_file, 'w') as outfile1:
    outfile1.write(to_write)

# def CreateProkkaJob(s_prokka_dir, s_final_assemblies_dir, prokka_jobs_dir, isolate, nt, walle, wallm, mem, mem_max, proteins_file, kingdom, genus, use_genus, species, gcode, rnammer, mincontiglen):
#     """Function to create a job to run BLAST on a set of de novo assemblies."""
#
#     # Header:
#     header1_l1 = "#! /usr/bin/env bash\n\n"
#     header1_l2 = "#BSUB -J " + "prokka_job" + '\n'
#     header1_l3 = "#BSUB -n " + str(nt) + '\n'
#     header1_l4 = '#BSUB -R "span[hosts=1]"' + '\n'
#     header1_l5 = '#BSUB -R "rusage[mem=' + str(mem) + ']"' + '\n'
#     header1_l7 = '#BSUB -M ' + str(mem_max) + '\n'
#     header1_l8 = '#BSUB -We ' + walle + '\n'
#     header1_l9 = '#BSUB -W ' + wallm + '\n'
#     header1_l10 = '#BSUB -o ' + s_prokka_dir + "/" + isolate + '_prokka.out' + '\n'
#     header1_l11 = '#BSUB -e ' + s_prokka_dir + "/" + isolate + '_prokka.err' + '\n'
#
#     header1 = header1_l1 + header1_l2 + header1_l3 + header1_l4 + header1_l5 + \
#               header1_l7 + header1_l8 + header1_l9 + header1_l10 +\
#               header1_l11 + '\n'
#
#     # Prokka input/output files:
#     isolate_assembly = s_final_assemblies_dir + "/" + isolate + "_filtered_ordered_assembly.fasta"
#     output_dir = os.path.join(s_prokka_dir, isolate)
#
#     # Set rnammer status:
#     if rnammer == "True":
#         rnammer_ = " --rnammer"
#     elif rnammer == "False":
#         rnammer_ = ""
#
#     # Set use_genus status:
#     if use_genus == "True":
#         use_genus_ = " --use_genus"
#     elif use_genus == "False":
#         use_genus_ = ""
#
#     # Set proteins file:
#     if proteins_file != "":
#         proteins_file_ = " --proteins " + proteins_file
#     elif proteins_file == "":
#         proteins_file_ = ""
#
#     # Prokka cmd:
#     prokka_cmd = "prokka --outdir " + output_dir + " --prefix " + isolate + " --kingdom " + kingdom + " --gcode " + \
#                  str(gcode) + " --genus " + genus + " --species " + species + proteins_file_ + " --cpus " + str(nt) + \
#                  rnammer_ + use_genus_ + " --rfam --addgenes --locustag " + isolate + " --mincontiglen " + \
#                  str(mincontiglen) + " " + isolate_assembly
#
#
#     # Write job file:
#     output_file = prokka_jobs_dir + "/" + isolate + "_prokka_job.sh"
#
#     to_write = header1 + prokka_cmd
#
#     with open(output_file, 'w') as outfile:
#         outfile.write(to_write)
