#!/usr/bin/env python3
"""
Author : Conrad Izydorczyk <conrad.izydorczyk@ucalgary.ca>
Date   : 2023-06-02
Purpose: filter VCFs to produce recombination masked versions
"""

import argparse
import shutil
import sys

# --------------------------------------------------
def get_args():
    """Get command-line arguments"""

    parser = argparse.ArgumentParser(
        description='filter VCFs to produce recombination masked versions',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-i',
                        '--input_vcf',
                        help='Unfiltered VCF',
                        metavar='str',
                        type=str,
                        required=True)
    
    parser.add_argument('-o',
                        '--output_vcf',
                        help='output filtered vcf',
                        type=str,
                        metavar='str',
                        required=True)
    
    parser.add_argument('-r',
                        '--rc_regions',
                        help='rc regions file from maskrc-svg',
                        type=str,
                        metavar='str',
                        required=True)
    
    parser.add_argument('-s',
                        '--isolate',
                        help='sample name. must be in regions file',
                        type=str,
                        metavar='str',
                        required=True)

    parser.add_argument('--summary',
                        help='summary of # of snps removed for each isolate',
                        type=str,
                        metavar='str',
                        required=True)

    return parser.parse_args()


# --------------------------------------------------
def main():
    """Make a jazz noise here"""

    args = get_args()

    print(f"\nProcessing isolate {args.isolate} to remove recombinant regions as identified in the rc regions file: {args.rc_regions}")
    
    # Create dict from regions file
    regions_dict = {}
    with open(args.rc_regions, "r") as infile1:
        for line in infile1:
            line_elements = line.strip().split("\t")
            line_sample = line_elements[0]
            line_start = int(line_elements[1])
            line_end = int(line_elements[2])

            if line_sample not in regions_dict:
                regions_dict[line_sample] = [(line_start, line_end)]
            elif line_sample in regions_dict:
                regions_dict[line_sample].append((line_start, line_end))

    # Remove VCF entries if position falls in any corresponding RC region for that sample
    output_vcf_entries = []
    orig_vcf_len = 0

    with open(args.input_vcf, "r") as infile2:
        for line in infile2:
            orig_vcf_len += 1
            if line.startswith("#"):
                output_vcf_entries.append(line)
            elif not line.startswith("#"):
                line_elements = line.strip().split("\t")
                try:
                    if any(int(line_elements[1]) in range(interval[0], interval[1]+1) for interval in regions_dict[args.isolate]):
                        continue
                    else:
                        output_vcf_entries.append(line)
                except KeyError:
                    print(f"\nIsolate {args.isolate} has no recombinant regions. Creating copy of input file as 'filtered output'.")
                    shutil.copyfile(args.input_vcf, args.output_vcf)
                    sys.exit()
    
    # Write filtered VCF:
    with open(args.output_vcf, "w") as outfile1:
        outfile1.write("".join(output_vcf_entries))

    # Write summary file:
    with open(args.summary, "a") as outfile2:
        outfile2.write(f"{args.isolate}\t{orig_vcf_len}\t{len(output_vcf_entries)}\n")
                        





# --------------------------------------------------
if __name__ == '__main__':
    main()
