import argparse
import os.path
import random

parser = argparse.ArgumentParser()

parser.add_argument("--nt", help="number of threads to use")
parser.add_argument("--maxtime", help="max time for computation")
parser.add_argument("--memmax", help="max memory to request")
parser.add_argument("--partitions", help="comma sep list of partitions")

parser.add_argument("--r1")
parser.add_argument("--r2")

parser.add_argument("--ref", help="assembly fasta to annotate")
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

# snippy_cmd = "conda activate i18_snp_pipeline_env_testing\n\nsnippy --R1 " + args.r1 + " --R2 " + args.r2 + \
#     " --ref " + args.ref + " --outdir " + args.output_dir + "/" + args.isolate + " --cpus " + args.nt + " --minfrac 0.9 --cleanup"

snippy_cmd = f"source activate i18_snp_pipeline_env_testing\n\nsnippy --R1 {args.r1} --R2 {args.r2} --ref {args.ref} --outdir {args.output_dir}{args.isolate}/ --cpus {args.nt} --minfrac 0.9"

to_write = header + "\n\n" + snippy_cmd

with open(args.job_file, "w") as outfile:
    outfile.write(to_write)
