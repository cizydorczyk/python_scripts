#!/usr/bin/env python3
"""
Author : Conrad Izydorczyk <conrad.izydorczyk@ucalgary.ca>
Date   : 2023-05-04
Purpose: create plasmidspades jobs for ARC
"""

import argparse
import os


# --------------------------------------------------
def get_args():
    """Get command-line arguments"""

    parser = argparse.ArgumentParser(
        description="Create plasmidspades jobs for ARC.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    # ----------------------------------------------
    # General options
    # ----------------------------------------------

    parser.add_argument(
        "-nt",
        "--num_threads",
        help="Number of threads to use for job. Set based on partitions desired. Default=8 appropriate for partitions single,lattice.",
        type=int,
        metavar="int",
    )

    parser.add_argument(
        "-mx",
        "--memmax",
        help="Max memory to use for job in MB; default=all available (0). Required.",
        type=int,
        metavar="int",
    )

    parser.add_argument(
        "-mu",
        "--mem_use",
        help="max memory to use for spades.py option -m (in GB); set based on ARC partition"
    )

    parser.add_argument(
        "-mt",
        "--maxtime",
        help="max time for job to run. Required.",
        type=str,
    )

    parser.add_argument(
        "-p",
        "--partitions",
        help='ARC partitions to use, comma-separated list; default="single,lattice". Required.',
        type=str,
        metavar="str",
    )

    parser.add_argument(
        "-i",
        "--isolate_list",
        help="path to isolate list on local computer; one isolate per line in list. Required.",
        metavar="file",
    )

    parser.add_argument(
        "-cs",
        "--chunk_size",
        help="number of isolates (chunk) to analyze in one job file. Required. Recommended size for Unicycler = 1.",
        type=int,
        metavar="int",
    )

    parser.add_argument(
        "-fe",
        "--fastq_ending",
        help="fastq file ending, one of: .fastq, .fq, .fastq.gz, or .fq.gz. Required.",
        type=str,
        metavar="str",
    )

    parser.add_argument(
        "-pe",
        "--pair_ending",
        help="how forward/reverse read pairs are specified. One of _1,_2 or _R1,_R2.",
        type=str,
        metavar="str",
    )

    parser.add_argument(
        "-d",
        "--input_dir",
        help="input directory on ARC containing fastq files to be trimmed. Default = current dir. Required.",
        type=str,
        metavar="path/to/dir",
    )

    parser.add_argument(
        "-o",
        "--output_dir",
        help="output directory on ARC where trimmed files will be placed. Default = current dir. Required.",
        type=str,
        metavar="path/to/dir",
    )

    parser.add_argument(
        "-j",
        "--jobs_dir",
        help="directory (local) to write job files to. default = current dir. Required.",
        type=str,
        metavar="path/to/dir",
    )

    parser.add_argument(
        "-pre",
        "--prefix",
        help="job name prefix for job scripts, .out, and .err files. Required.",
        type=str,
        metavar="str",
    )

    parser.add_argument(
        "-e",
        "--env",
        help="Conda environment on ARC that has desired analysis software installed. Required.",
        type=str,
        metavar="str",
    )

    return parser.parse_args()


# --------------------------------------------------
def main():

    args = get_args()
    
    # Read in isolate list:
    isolate_list = readIsolateList(args.isolate_list)

    # Create fastqc commands:
    pair_endings = args.pair_ending.strip().split(",")

    plasmidspades_commands = {}
    for isolate in isolate_list:
        input_r1 = os.path.join(
            args.input_dir, f"{isolate}{pair_endings[0]}{args.fastq_ending}"
        )
        input_r2 = os.path.join(
            args.input_dir, f"{isolate}{pair_endings[1]}{args.fastq_ending}"
        )
        output_dir = os.path.join(args.output_dir, isolate)

        plasmidspades_cmd = (
            f"spades.py --plasmid --careful -1 {input_r1} -2 {input_r2} -o {output_dir}/ "
            f"-t {args.num_threads} -m {args.mem_use}"
        )

        plasmidspades_commands[isolate] = plasmidspades_cmd

    # Create chunked job files:
    chunks = [
        isolate_list[i : i + args.chunk_size]
        for i in range(0, len(isolate_list), args.chunk_size)
    ]

    for idx, chunk in enumerate(chunks):
        job_script = os.path.join(args.jobs_dir, f"{idx}.{args.prefix}.slurm")

        # Create job header:
        header = (
            f"#!/bin/bash\n"
            f"#SBATCH --time={args.maxtime}\n"
            f"#SBATCH --nodes=1\n"
            f"#SBATCH --ntasks=1\n"
            f"#SBATCH --cpus-per-task={args.num_threads}\n"
            f"#SBATCH --partition={args.partitions}\n"
            f"#SBATCH --mem={args.memmax}\n"
            f"#SBATCH --output={idx}.{args.prefix}.out\n"
            f"#SBATCH --error={idx}.{args.prefix}.err\n"
            f"\nsource activate {args.env}\n\n"
        )

        # Create individual isolate commands:
        job_commands = []

        for isolate in chunk:
            job_cmd = plasmidspades_commands[isolate]
            job_commands.append(job_cmd)

        job_commands = "\n".join(job_commands)

        # Write job script:
        with open(job_script, "w") as job_file_1:
            job_file_1.write(f"{header}\n{job_commands}")

def readIsolateList(isolate_list_file):
    # Read in isolate list:
    with open(isolate_list_file, "r") as isolate_list_file:
        isolate_list = [isolate.strip() for isolate in isolate_list_file]

    return isolate_list

# --------------------------------------------------
if __name__ == '__main__':
    main()
