import argparse
import os.path
import random

parser = argparse.ArgumentParser()

parser.add_argument("--nt", help="number of threads to use")
parser.add_argument("--maxtime", help="max time for computation")
parser.add_argument("--memmax", help="max memory to request")
parser.add_argument("--partitions", help="comma sep list of partitions")

parser.add_argument("--isolate", help="isolate name/number")
parser.add_argument("--job_file", help="output job file")

parser.add_argument("--arc_project_dir", help="dir with xml file & where working directory will be on arc; output dir on arc")

parser.add_argument("--input_xml", help="full path to beauti xml on arc") 

args = parser.parse_args()

# Create header
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

# Set seeds
seed1 = random.randint(1,1000000000)

# Export beagle lib commands
beagle_cmd1 = "export LD_LIBRARY_PATH=$HOME/lib:$LD_LIBRARY_PATH"
beagle_cmd2 = "export PKG_CONFIG_PATH=$HOME/lib/pkgconfig:$PKG_CONFIG_PATH"

# BEAST command
beast_cmd1 = '/home/conrad.izydorczyk/software/BEASTv1.10.4/bin/beast -seed ' + str(seed1) + " " + args.input_xml

# Write job file
to_write = header + "\n\n" + beagle_cmd1 + "\n\n" + beagle_cmd2 + "\n\n" + beast_cmd1
with open(args.job_file, "w") as outfile1:
    outfile1.write(to_write)
