import argparse
import numpy

parser = argparse.ArgumentParser()

parser.add_argument("--nt", help="number of threads to use")
parser.add_argument("--maxtime", help="max time for computation")
parser.add_argument("--memmax", help="max memory to request")
parser.add_argument("--partitions", help="comma sep list of partitions")

parser.add_argument("--isolate_list")
parser.add_argument("--job_scripts_directory")
parser.add_argument("--chunk_size", help="number of individual trimmomatic commands to include in each job file. deafult=5", default=5, type=int)
parser.add_argument("--fastq_ending", help="fastq or fq or fastq.gz or fq.gz. default=fastq.gz", default="fastq.gz")

parser.add_argument("--arc_input_fastq_directory")
parser.add_argument("--arc_output_fastq_directory")

parser.add_argument("--crop_length")
parser.add_argument("--min_length")
parser.add_argument("--sliding_window_quality")
parser.add_argument("--arc_adapters_file")

args = parser.parse_args()

# Read isolate list:
with open(args.isolate_list, "r") as infile1:
    isolate_list = [line.strip() for line in infile1]

# Create trimmomatic commands:
trimmomatic_commands = {}
for isolate in isolate_list:
    input_r1 = args.arc_input_fastq_directory + isolate + "_1." + args.fastq_ending
    input_r2 = args.arc_input_fastq_directory + isolate + "_2." + args.fastq_ending
    output_r1 = args.arc_output_fastq_directory + isolate + "_1.fastq.gz"
    output_r2 = args.arc_output_fastq_directory + isolate + "_2.fastq.gz"
    output_r1_unpaired = args.arc_output_fastq_directory + isolate + "_u_1.fastq.gz"
    output_r2_unpaired = args.arc_output_fastq_directory + isolate + "_u_2.fastq.gz"

    trimmomatic_cmd = "trimmomatic PE -threads " + args.nt + " " + input_r1 + " " + input_r2 + " " + output_r1 + " " + output_r1_unpaired + " " + output_r2 + " " + output_r2_unpaired + " ILLUMINACLIP:" + args.arc_adapters_file + ":2:30:10:8:true CROP:" + args.crop_length + " SLIDINGWINDOW:4:" + args.sliding_window_quality + " MINLEN:" + args.min_length + "\n"

    trimmomatic_commands[isolate] = trimmomatic_cmd

# Create chunked job files:
chunks = numpy.array_split(isolate_list, args.chunk_size)

for chunk in chunks:
    job_script = args.job_scripts_directory + "_".join(chunk) + "_slurm.sh"

    # Create job header:
    header1 = "#!/bin/bash\n"
    header7 = "#SBATCH --time=" + args.maxtime
    header3 = "#SBATCH --nodes=1"
    header2 = "#SBATCH --ntasks=1"
    header4 = "#SBATCH --cpus-per-task=" + args.nt
    header5 = "#SBATCH --partition=" + args.partitions
    header6 = "#SBATCH --mem=" + args.memmax
    header8 = "#SBATCH --output=" + "_".join(chunk) + ".out"
    header9 = "#SBATCH --error=" + "_".join(chunk) + ".err"
    header10 = "\n\nsource activate trimmomatic-0.39\n\n"

    header = "\n".join([header1, header7, header2, header3, header4, header5, header6, header8, header9, header10])
    
    # Create commands:
    job_commands = []

    for isolate in chunk:
        job_cmd = trimmomatic_commands[isolate]
        job_commands.append(job_cmd)
    
    with open(job_script, "w") as outfile1:
        outfile1.write(header + "\n".join(job_commands))