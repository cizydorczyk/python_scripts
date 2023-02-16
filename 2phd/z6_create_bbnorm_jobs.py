#!/usr/bin/env python3
"""
Author : Conrad Izydorczyk <conrad.izydorczyk@ucalgary.ca>
Date   : 2023-02-02
Purpose: Create bbnorm.sh jobs for arc.
"""

import argparse
import os.path

# --------------------------------------------------
def get_args():
    """Get command-line arguments"""

    parser = argparse.ArgumentParser(
        description='Create bbnorm.sh jobs for arc.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

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

    # BBNorm-specific options
    parser.add_argument(
        '-t',
        '--target_kmer_depth',
        help='target kmer depth; ReadDepth = KmerDepth*(ReadLength/(ReadLength-KmerLength+1))',
        metavar='str',
        type=str,
        required=True
    )

    return parser.parse_args()


# --------------------------------------------------
def main():
    """Function to create bbnorm jobs on ARC."""

    args = get_args()
    num_threads_arg=args.num_threads
    memmax_arg=args.memmax
    maxtime_arg=args.maxtime
    partitions_arg=args.partitions
    isolate_list_arg=args.isolate_list
    chunk_size_arg=args.chunk_size
    fastq_ending_arg=args.fastq_ending
    input_dir_arg=args.input_dir
    output_dir_arg=args.output_dir
    job_prefix_arg=args.prefix
    jobs_dir_arg=args.jobs_dir
    pair_ending_arg=args.pair_ending.split(",")
    kmer_depth_arg = args.target_kmer_depth


    # Read in isolate list:
    isolate_list = []
    with open(isolate_list_arg, "r") as infile1:
        for line in infile1:
            isolate_list.append(line.strip())
    
    # Create trimmomatic commands:
    bbnorm_commands = {}
    for isolate in isolate_list:
        input_r1 = os.path.join(
            input_dir_arg, f"{isolate}{pair_ending_arg[0]}{fastq_ending_arg}"
        )
        input_r2 = os.path.join(
            input_dir_arg, f"{isolate}{pair_ending_arg[1]}{fastq_ending_arg}"
        )
        output_r1 = os.path.join(
            output_dir_arg, f"{isolate}{pair_ending_arg[0]}{fastq_ending_arg}"
        )
        output_r2 = os.path.join(
            output_dir_arg, f"{isolate}{pair_ending_arg[1]}{fastq_ending_arg}"
        )

        bbnorm_cmd = (
            f"/home/conrad.izydorczyk/software/bbmap/bbnorm.sh "
            f"in={input_r1} in2={input_r2} out1={output_r1} out2={output_r2} target={kmer_depth_arg}\n"
        )

        bbnorm_commands[isolate] = bbnorm_cmd

    # Create chunked job files:
    chunks = [
        isolate_list[i : i + chunk_size_arg]
        for i in range(0, len(isolate_list), chunk_size_arg)
    ]

    for idx, chunk in enumerate(chunks):
        job_script = os.path.join(jobs_dir_arg, f"{idx}.{job_prefix_arg}.slurm")

        # Create job header:
        header = (
            f"#!/bin/bash\n"
            f"#SBATCH --time={maxtime_arg}\n"
            f"#SBATCH --nodes=1\n"
            f"#SBATCH --ntasks=1\n"
            f"#SBATCH --cpus-per-task={num_threads_arg}\n"
            f"#SBATCH --partition={partitions_arg}\n"
            f"#SBATCH --mem={memmax_arg}\n"
            f"#SBATCH --output={idx}.{job_prefix_arg}.out\n"
            f"#SBATCH --error={idx}.{job_prefix_arg}.err\n"
        )

        # Create individual isolate commands:
        job_commands = []

        for isolate in chunk:
            job_cmd = bbnorm_commands[isolate]
            job_commands.append(job_cmd)

        job_commands = "\n".join(job_commands)

        # Write job script:
        with open(job_script, "w") as job_file_1:
            job_file_1.write(f"{header}\n{job_commands}")


# --------------------------------------------------
if __name__ == '__main__':
    main()
