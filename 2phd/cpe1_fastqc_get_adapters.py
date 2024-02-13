#!/usr/bin/env python3
"""
Author : Conrad Izydorczyk <conradizydorczk@outlook.com>
Date   : 2024-01-30
Purpose: get fastq adapter contamination
"""

import argparse


# --------------------------------------------------
def get_args():
    """Get command-line arguments"""

    parser = argparse.ArgumentParser(
        description='get fastq adapter contamination',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-i',
                        '--multiqc_adapters_file',
                        help='multiqc adapter contamination file',
                        metavar='str',
                        type=str,
                        required=True)
    
    parser.add_argument('-o',
                        '--output_file',
                        help='output tsv with read lengths',
                        metavar='str',
                        type=str,
                        required=True)

    return parser.parse_args()


# --------------------------------------------------
def main():
    """Make a jazz noise here"""

    args = get_args()
    
    output_dict = {}
    with open(args.multiqc_adapters_file, "r") as infile1:
        for line_number, line in enumerate(infile1):
            line_elements = line.strip().split("\t")
            
            if line_number % 2 == 0:
                continue
            else:
                isolate = line_elements[0].split(" - ")[0].strip()
                adapter = line_elements[0].split(" - ")[1].strip()
                
                isolate_ = isolate.split("_")[0]
                
                if isolate_ not in output_dict:
                    output_dict[isolate_] = adapter

    to_write = []
    for key, value in output_dict.items():
        to_write.append(f"{key}\t{value}")
    
    with open(args.output_file, "w") as outfile1:
        outfile1.write("\n".join(to_write))



# --------------------------------------------------
if __name__ == '__main__':
    main()
