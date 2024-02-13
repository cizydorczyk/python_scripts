#!/usr/bin/env python3
"""
Author : Conrad Izydorczyk <conrad.izydorczyk@ucalgary.ca>
Date   : 2022-11-09
Purpose: Rock the Casbah
"""

import argparse
import os.path


# --------------------------------------------------
def get_args():
    """Get command-line arguments"""

    parser = argparse.ArgumentParser(
        description="Create FastQC, Trimmomatic, Unicycler, SKESA, or Bakta jobs for ARC.",
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
        "-d",
        "--input_dir",
        help="input directory on ARC containing fastq files to be trimmed. Default = current dir. Required.",
        type=str,
        metavar="path/to/dir",
    )

    parser.add_argument(
        "-o",
        "--output_file",
        help="output file on ARC where trimmed files will be placed. Default = current dir. Required.",
        type=str,
        metavar="path/to/file",
    )

    parser.add_argument(
        "-j",
        "--jobs_dir",
        help="directory (local) to write job files to. default = current dir. Required.",
        type=str,
        metavar="path/to/dir",
    )

    parser.add_argument(
        "-e",
        "--env",
        help="Conda environment on ARC that has desired analysis software installed. Required.",
        type=str,
        metavar="str",
    )

    parser.add_argument(
        "--stringmlst_db",
        help="directory on ARC with stringmlst database (must already be downloaded)",
        type=str,
        metavar="directory"
    )

    parser.add_argument(
        "--isolate_list",
        help="isolate list, one isolate per line",
        required=True,
    )

    parser.add_argument(
        "--chunk_size",
        help="number of isolates to process per job",
        required=True,
        type=int,
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
        "-pre",
        "--prefix",
        help="job name prefix for job scripts, .out, and .err files. Required.",
        type=str,
        metavar="str",
    )

    return parser.parse_args()


# --------------------------------------------------
def main():
    """Make a jazz noise here"""

    args = get_args()

    isolate_list = read_isolate_list(args.isolate_list)
    pair_ending = args.pair_ending.split(",")

    # Create trimmomatic commands:
    stringmlst_commands = {}
    for isolate in isolate_list:
        input_r1 = os.path.join(
            args.input_dir, f"{isolate}{pair_ending[0]}{args.fastq_ending}"
        )
        input_r2 = os.path.join(
            args.input_dir, f"{isolate}{pair_ending[1]}{args.fastq_ending}"
        )
        
        stringmlst_cmd = (f"stringMLST.py --predict --prefix {args.stringmlst_db} -o {args.output_file} "
                          f"-1 {input_r1} -2 {input_r2}")

        stringmlst_commands[isolate] = stringmlst_cmd

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
            job_cmd = stringmlst_commands[isolate]
            job_commands.append(job_cmd)

        job_commands = "\n".join(job_commands)

        # Write job script:
        with open(job_script, "w") as job_file_1:
            job_file_1.write(f"{header}\n{job_commands}")

# --------------------------------------------------
def read_isolate_list(isolate_list_handle):
    isolate_list = []
    with open(isolate_list_handle, "r") as infile1:
        for line in infile1:
            isolate_list.append(line.strip())
    return(isolate_list)

# def create_stringmlst_jobs(
#     num_threads_arg,
#     memmax_arg,
#     maxtime_arg,
#     partitions_arg,
#     input_dir_arg,
#     output_file_arg,
#     jobs_dir_arg,
#     env_arg,
#     stringmlst_db_arg,
#     ):

#     """Function to create stringMLST jobs for ARC."""
  
#     stringmlst_cmd = (
#         f"stringMLST.py --predict -d {input_dir_arg} --prefix {stringmlst_db_arg} -o {output_file_arg}"
#     )

#     # Create job header:
#     header = (
#         f"#!/bin/bash\n"
#         f"#SBATCH --time={maxtime_arg}\n"
#         f"#SBATCH --nodes=1\n"
#         f"#SBATCH --ntasks=1\n"
#         f"#SBATCH --cpus-per-task={num_threads_arg}\n"
#         f"#SBATCH --partition={partitions_arg}\n"
#         f"#SBATCH --mem={memmax_arg}\n"
#         f"#SBATCH --output=stringmlst.out\n"
#         f"#SBATCH --error=stringmlst.err\n"
#         f"\nsource activate {env_arg}\n\n"
#     )

#     # Write job script:        
#     with open(jobs_dir_arg, "w") as job_file_1:
#         job_file_1.write(f"{header}\n{stringmlst_cmd}")

# --------------------------------------------------
if __name__ == '__main__':
    main()
