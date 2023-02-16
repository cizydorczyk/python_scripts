#!/usr/bin/env python3
"""
Author : Conrad Izydorczyk <conrad.izydorczyk@ucalgary.ca>
Date   : 2022-12-02
Purpose: Create cfml jobs arc jobs.
"""

import argparse


# --------------------------------------------------
def get_args():
    """Get command-line arguments"""

    parser = argparse.ArgumentParser(
        description="Create CFML jobs for ARC.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    # ----------------------------------------------
    # General options
    # ----------------------------------------------
    parser.add_argument(
        "-nt",
        "--num_threads",
        help="Max # of threads IQ-Tree can use; set based on partitions selected",
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
        "--input_alignment",
        help="input alignment file on ARC",
        type=str,
        metavar="path/to/dir",
    )

    parser.add_argument(
        "--input_tree",
        help="input tree for CFML",
        type=str,
        metavar="str",
    )

    parser.add_argument(
        "-o",
        "--output_dir",
        help="output dir on ARC",
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

    return parser.parse_args()


# --------------------------------------------------
def main():
    """Create CFML jobs."""

    args = get_args()
    
    create_cfml_jobs(
        num_threads_arg = args.num_threads,

        memmax_arg=args.memmax,
        maxtime_arg=args.maxtime,
        partitions_arg=args.partitions,
        input_alignment_arg=args.input_alignment,
        input_tree_arg=args.input_tree,
        output_dir_arg=args.output_dir,
        jobs_dir_arg=args.jobs_script,
        env_arg=args.env,
    )

def create_cfml_jobs(
    num_threads_arg,

    memmax_arg,
    maxtime_arg,
    partitions_arg,
    input_alignment_arg,
    input_tree_arg,
    output_dir_arg,
    jobs_dir_arg,
    env_arg,
    ):

    """Function to create CFML jobs for ARC."""
  
    cfml_cmd = (
        f"ClonalFrameML {input_tree_arg} {input_alignment_arg} {output_dir_arg} -em true"
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
        f"#SBATCH --output=cfml.out\n"
        f"#SBATCH --error=cfml.err\n"
        f"\nsource activate {env_arg}\n\n"
    )

    # Write job script:        
    with open(jobs_dir_arg, "w") as job_file_1:
        job_file_1.write(f"{header}\n{cfml_cmd}")

# --------------------------------------------------
if __name__ == '__main__':
    main()