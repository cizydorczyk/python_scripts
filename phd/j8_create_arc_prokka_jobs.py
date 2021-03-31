import argparse
import os.path
import random

parser = argparse.ArgumentParser()

parser.add_argument("--nt", help="number of threads to use")
parser.add_argument("--maxtime", help="max time for computation")
parser.add_argument("--memmax", help="max memory to request")
parser.add_argument("--partitions", help="comma sep list of partitions")

parser.add_argument("--fasta", help="assembly fasta to annotate")
parser.add_argument("--output_dir", help="output dir")
parser.add_argument("--isolate", help="isolate name/number")
parser.add_argument("--job_file", help="output job file")
parser.add_argument("--proteins_file", help="prokka proteins file")
parser.add_argument("--genus")
parser.add_argument("--species")
parser.add_argument("--mincontiglen")

args = parser.parse_args()

header1 = "#!/bin/bash\n"
header7 = "#SBATCH --time=" + args.maxtime
header3 = "#SBATCH --nodes=1"
header2 = "#SBATCH --ntasks=1"
header4 = "#SBATCH --cpus-per-task=" + args.nt
header5 = "#SBATCH --partition=" + args.partitions
header6 = "#SBATCH --mem=" + args.memmax
header8 = "#SBATCH --output=" + args.isolate + ".out"
header9 = "#SBATCH --error=" + args.isolate + ".err"

header = "\n".join([header1, header7, header2, header3, header4, header5, header6, header8, header9])

unicycler_cmd = "module load biobuilds/conda\n\nsource activate prokka-env\n\nprokka --outdir " + args.output_dir + " --prefix " +\
    args.isolate + " --addgenes --locustag " + args.isolate + " --genus " + args.genus + " --species " + args.species + \
    " --kingdom Bacteria --gcode 11 --proteins " + args.proteins_file + " --cpus " + args.nt + " --mincontiglen " +\
    args.mincontiglen + " --rfam " + args.fasta

to_write = header + "\n\n" + unicycler_cmd

print(to_write)

with open(args.job_file, "w") as outfile:
    outfile.write(to_write)
