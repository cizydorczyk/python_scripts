import argparse

parser = argparse.ArgumentParser()

parser.add_argument("--nt", help="number of threads to use")
parser.add_argument("--maxtime", help="max time for computation")
parser.add_argument("--memmax", help="max memory to request")
parser.add_argument("--partitions", help="comma sep list of partitions")

parser.add_argument("--clean_mode", help="panaroo clean mode; sensitive, moderate, or strict")
parser.add_argument("--arc_output_dir", help="panaroo output dir on arc")
parser.add_argument("--arc_file_paths_file", help="file of file paths to gffs on arc")
parser.add_argument("--arc_alignment", default="", help="perform alignment? core or pan; if excluded, does not align")

parser.add_argument("--job_file", help="output job file")
parser.add_argument("--run", help="run name for output .out and .err files")

args = parser.parse_args()

header1 = "#!/bin/bash\n"
header7 = "#SBATCH --time=" + args.maxtime
header3 = "#SBATCH --nodes=1"
header2 = "#SBATCH --ntasks=1"
header4 = "#SBATCH --cpus-per-task=" + args.nt
header5 = "#SBATCH --partition=" + args.partitions
header6 = "#SBATCH --mem=" + args.memmax
header8 = "#SBATCH --output=" + args.run + ".out"
header9 = "#SBATCH --error=" + args.run + ".err"

header = "\n".join([header1, header7, header2, header3, header4, header5, header6, header8, header9])

# Set alignment parameters:
if args.arc_alignment == "":
    alignment = ""
elif args.arc_alignment == "core":
    alignment = " -a core --aligner mafft --core_threshold 0.95 "
elif args.arc_alignment == "pan":
    alignment = " -a pan --aligner mafft --core_threshold 0.95 "

# Create panaroo command
panaroo_cmd = "conda activate base\n\nconda activate panaroo-1.2.9\n\npanaroo -i " + args.arc_file_paths_file + " -o " + args.arc_output_dir + alignment + " -t " + args.nt + " --threshold 0.98 + --len_dif_percent 0.98 + --clean-mode " + args.clean_mode + " --remove-invalid-genes"

# Write job file:
to_write = header + "\n\n" + panaroo_cmd

print(to_write)

with open(args.job_file, "w") as outfile:
    outfile.write(to_write)