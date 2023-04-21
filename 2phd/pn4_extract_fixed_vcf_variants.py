#!/usr/bin/env python3
"""
Author : Conrad Izydorczyk <conrad.izydorczyk@ucalgary.ca>
Date   : 2023-04-07
Purpose: Filter bcftools combined vcf file to get fixed variants only
"""

import argparse


# --------------------------------------------------
def get_args():
    """Get command-line arguments"""

    parser = argparse.ArgumentParser(
        description='Filter bcftools combined vcf file to get fixed variants only',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-i',
                        '--input_vcf',
                        help='A named string argument',
                        metavar='str',
                        type=str,
                        required=True)
    
    parser.add_argument('-r',
                        '--reference_genotype',
                        help='reference genotype (how the reference genotype is recorded in the vcf file, e.g. "./.")',
                        type=str,
                        metavar='str',
                        required=True)
    
    parser.add_argument('-of',
                        '--output_fixed',
                        help='output vcf with fixed variants only',
                        metavar='str',
                        type=str,
                        required=True)
    
    parser.add_argument('-os',
                        '--output_segregating',
                        help='output vcf with segregating variants only',
                        metavar='str',
                        type=str,
                        required=True)    


    return parser.parse_args()


# --------------------------------------------------
def main():
    """Filter bcftools combined vcf file to get fixed variants only"""

    args = get_args()
    
    header = []
    fixed_positions = []
    segregating_positions = []
    with open(args.input_vcf, "r") as infile1:
        for line in infile1:
            if not line.startswith("#"):
                genotypes = [i.strip().split(":")[0] for i in line.strip().split("\t")[9:]]
                if args.reference_genotype in genotypes:
                    segregating_positions.append(line.strip())
                else:
                    if len(set(genotypes)) == 1:
                        fixed_positions.append(line.strip())
            else:
                header.append(line.strip())
    
    to_write_fixed = "\n".join(header) + "\n" + "\n".join(fixed_positions)
    to_write_segregating = "\n".join(header) + "\n" + "\n".join(segregating_positions)

    with open(args.output_fixed, "w") as outfile1:
        outfile1.write(to_write_fixed)
    
    with open(args.output_segregating, "w") as outfile2:
        outfile2.write(to_write_segregating)


# --------------------------------------------------
if __name__ == '__main__':
    main()
