#!/usr/bin/env python3
"""
Author : Conrad Izydorczyk <conrad.izydorczyk@ucalgary.ca>
Date   : 2022-12-02
Purpose: Create arc iqtree jobs.
"""

import argparse
import os.path


# --------------------------------------------------
def get_args():
    """Get command-line arguments"""

    parser = argparse.ArgumentParser(
        description="Create IQ-Tree jobs for ARC.",
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

    parser.add_argument(
        "--bb",
        help="iqtree -bb; number of rapid bootstraps",
        type=int,
        metavar="int",
    )

    parser.add_argument(
        "--model",
        help="iqtree -model; model to use",
        metavar="str",
        type=str,
    )

    parser.add_argument(
        "--mem",
        help="memory for iqtree to use; should not exceed max available on partition in gb, ex. 64G",
        type=str,
        metavar="str",
    )

    return parser.parse_args()


# --------------------------------------------------
def main():
    """Create IQ-Tree jobs."""

    args = get_args()
    
    create_iqtree_jobs(
        num_threads_arg = args.num_threads,
        mem_arg = args.mem,
        memmax_arg=args.memmax,
        maxtime_arg=args.maxtime,
        partitions_arg=args.partitions,
        input_alignment_arg=args.input_alignment,
        output_dir_arg=args.output_dir,
        jobs_dir_arg=args.jobs_script,
        env_arg=args.env,
        bb_arg=args.bb,
        model_arg=args.model,
    )

def create_iqtree_jobs(
    num_threads_arg,
    mem_arg,
    memmax_arg,
    maxtime_arg,
    partitions_arg,
    input_alignment_arg,
    output_dir_arg,
    jobs_dir_arg,
    env_arg,
    bb_arg,
    model_arg,
    ):

    """Function to create IQ-Tree jobs for ARC."""
  
    iqtree_cmd = (
        f"iqtree -s {input_alignment_arg} -bb {bb_arg} -m {model_arg} -pre {output_dir_arg} -nt auto -mem {mem_arg} -ntmax {num_threads_arg}"
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
        f"#SBATCH --output=iqtree.out\n"
        f"#SBATCH --error=iqtree.err\n"
        f"\nsource activate {env_arg}\n\n"
    )

    # Write job script:        
    with open(jobs_dir_arg, "w") as job_file_1:
        job_file_1.write(f"{header}\n{iqtree_cmd}")

# --------------------------------------------------
if __name__ == '__main__':
    main()

