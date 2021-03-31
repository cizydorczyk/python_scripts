import argparse
import os.path
import random

parser = argparse.ArgumentParser()

parser.add_argument("--nt", help="number of threads to use")
parser.add_argument("--maxtime", help="max time for computation")
parser.add_argument("--memmax", help="max memory to request")
parser.add_argument("--partitions", help="comma sep list of partitions")

parser.add_argument("--f1", help="f1 fastq file")
parser.add_argument("--f2", help="f2 fastq file")
parser.add_argument("--output_dir", help="output dir")
parser.add_argument("--isolate", help="isolate name/number")
parser.add_argument("--job_file", help="output job file")

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

unicycler_cmd = "module load biobuilds/conda\n\nsource activate unicycler-0.4.8\n\nunicycler -1 " + args.f1 + " -2 " + args.f2 + " -o " + args.output_dir +\
    " -t " + args.nt + " --depth_filter 0.01"

to_write = header + "\n\n" + unicycler_cmd

with open(args.job_file, "w") as outfile:
    outfile.write(to_write)
