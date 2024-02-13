#!/usr/bin/env python3
"""
Author : Conrad Izydorczyk <conrad.izydorczyk@ucalgary.ca>
Date   : 2023-06-20
Purpose: Create snakemake snp calling pipeline config files.
"""

import argparse


# --------------------------------------------------
def get_args():
    """Get command-line arguments"""

    parser = argparse.ArgumentParser(
        description='Create snakemake snp calling pipeline config files.',
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
        help="path to isolate list on ARC; one isolate per line in list. Required.",
        metavar="file",
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
        help="input directory on ARC containing fastq files. Required.",
        type=str,
        metavar="path/to/dir",
    )

    parser.add_argument(
        "-o",
        "--output_dir",
        help="Project directory on ARC. Required.",
        type=str,
        metavar="path/to/dir",
    )

    parser.add_argument(
        "--st",
        help="mlst sequence type",
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

    parser.add_argument(
        "-e",
        "--env",
        help="Conda environment on ARC that has desired analysis software installed. Required.",
        type=str,
        metavar="str",
    )

    parser.add_argument(
        "--reference",
        help="reference gbk on ARC",
        type=str,
        metavar="str",
    )

    parser.add_argument(
        "--snippy_minfrac",
        help="snippy minfrac option",
        type=str,
        metavar="str",
    )

    parser.add_argument(
        "--snippy_cleanup",
        help="Snippy cleanup option. Required.",
        type=str,
        metavar="str",
        default="yes"
    )

    parser.add_argument(
        "--snippy_unmapped",
        help="Snippy unmapped option. Required.",
        type=str,
        metavar="str",
        default="no"
    )

    parser.add_argument(
        "--snippy_reference",
        help="Should be 'Reference' or empty string '' if want to keep reference. Required.",
        type=str,
        metavar="str",
        default="Reference"
    )

    parser.add_argument(
        "--outgroups",
        help="Comma separated list of outgroup isolates. Required.",
        type=str,
        metavar="str",
        default=""
    )

    parser.add_argument(
        "--bb",
        help="IQ-Tree -bb option. Required.",
        type=str,
        metavar="str"
    )

    parser.add_argument(
        "--iqtree_model",
        help="Phylogenetic model to use in iqtree. Required.",
        type=str,
        metavar="str"
    )

    parser.add_argument(
        "--snakefile",
        help="Snakefile to use on ARC. Required.",
        type=str,
        metavar="str"
    )

    parser.add_argument(
        "--arc_configfile",
        help="Config file to use on ARC for snakemake; made by this script, but put path to where it will be uploaded on ARC. Required.",
        type=str,
        metavar="str"
    )

    parser.add_argument(
        "-j",
        "--job_script",
        help="Local job script to produce. Required.",
        type=str,
        metavar="path/to/dir",
    )

    parser.add_argument(
        "--output_config",
        help="Local config file to produce. Required.",
        type=str,
        metavar="path/to/dir",
    )

    return parser.parse_args()


# --------------------------------------------------
def main():
    """Make a jazz noise here"""

    args = get_args()

    pair_endings = args.pair_ending.split(",")
    
    fastq_endings = [f"{pair_endings[0]}{args.fastq_ending}", f"{pair_endings[1]}{args.fastq_ending}"]

    
    config_text = f'''#snippy options
fastq_dir: "{args.input_dir}"
project_dir: "{args.output_dir}"
ref: "{args.reference}"
st: "{args.st}"
isolate_list: "{args.isolate_list}"
fastq_endings: "{fastq_endings[0]},{fastq_endings[1]}"
    
#snippy options
snippy_minfrac: "{args.snippy_minfrac}"
snippy_cleanup: "{args.snippy_cleanup}" # makes output dirs smaller; yes or no
snippy_unmapped: "{args.snippy_unmapped}" # for getting unmapped reads; do not set to yes yet    
    
# remove ref and outgroups options:
remove_ref: "{args.snippy_reference}" # snippy always names the reference sequence 'Reference'; can be an empty string
outgroups: "{args.outgroups}" # as many or as few outgroups as needed; comma separated list for >1; can be an empty string

# iqtree options:
bb: "{args.bb}" # change to set different number of rapid bootstraps
model: "{args.iqtree_model}" # change to specify a specific model
'''

    header = (
            f"#!/bin/bash\n"
            f"#SBATCH --time={args.maxtime}\n"
            f"#SBATCH --nodes=1\n"
            f"#SBATCH --ntasks=1\n"
            f"#SBATCH --cpus-per-task={args.num_threads}\n"
            f"#SBATCH --partition={args.partitions}\n"
            f"#SBATCH --mem={args.memmax}\n"
            f"#SBATCH --output={args.st}.{args.prefix}.out\n"
            f"#SBATCH --error={args.st}.{args.prefix}.err\n"
            f"\nsource activate {args.env}\n\n"
        )
    
    snakemake_cmd = f"snakemake -s {args.snakefile} --configfile {args.arc_configfile} --cores 1 -p"

    # Write config file (local):
    with open(args.output_config, "w") as outfile1:
        outfile1.write(config_text)
    
    # Write associated job script:
    with open(args.job_script, "w") as outfile2:
        outfile2.write(header + "\n\n" + snakemake_cmd)



# --------------------------------------------------
if __name__ == '__main__':
    main()
