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
        "--jobs_script",
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

    return parser.parse_args()


# --------------------------------------------------
def main():
    """Make a jazz noise here"""

    args = get_args()
    
    create_stringmlst_jobs(
        num_threads_arg=args.num_threads,
        memmax_arg=args.memmax,
        maxtime_arg=args.maxtime,
        partitions_arg=args.partitions,
        input_dir_arg=args.input_dir,
        output_file_arg=args.output_file,
        jobs_dir_arg=args.jobs_script,
        env_arg=args.env,
        stringmlst_db_arg=args.stringmlst_db,
    )

# --------------------------------------------------
def create_stringmlst_jobs(
    num_threads_arg,
    memmax_arg,
    maxtime_arg,
    partitions_arg,
    input_dir_arg,
    output_file_arg,
    jobs_dir_arg,
    env_arg,
    stringmlst_db_arg,
    ):

    """Function to create stringMLST jobs for ARC."""
  
    stringmlst_cmd = (
        f"stringMLST.py --predict -d {input_dir_arg} --prefix {stringmlst_db_arg} -o {output_file_arg}"
    )

    # Create job header:
    header = (
        f"#!/bin/bash\n"
        f"#SBATCH --time={maxtime_arg}\n"
        f"#SBATCH --nodes=1\n"
        f"#SBATCH --ntasks=1\n"
        f"#SBATCH --cpus-per-task={num_threads_arg}\n"
        f"#SBATCH --partition={partitions_arg}\n"
        f"#SBATCH --mem={memmax_arg}\n"
        f"#SBATCH --output=stringmlst.out\n"
        f"#SBATCH --error=stringmlst.err\n"
        f"\nsource activate {env_arg}\n\n"
    )

    # Write job script:        
    with open(jobs_dir_arg, "w") as job_file_1:
        job_file_1.write(f"{header}\n{stringmlst_cmd}")

# --------------------------------------------------
if __name__ == '__main__':
    main()
