#!/usr/bin/env python3
"""
Author : Conrad Izydorczyk <conrad.izydorczyk@ucalgary.ca>
Date   : 2022-10-12
Purpose: Get blaZ types from nuc sequences.
"""

import argparse
from Bio import SeqIO
import os.path


# --------------------------------------------------
def get_args():
    """Get command-line arguments"""

    parser = argparse.ArgumentParser(
        description='Get blaZ types from nuc sequences.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)


    parser.add_argument('-i',
                        '--input',
                        help='Input kraken report',
                        metavar='str',
                        type=str,
                        required=True)
    
    parser.add_argument('-o',
                        '--output',
                        help='Output file',
                        metavar='str',
                        type=str,
                        required=True)
    
    parser.add_argument('--isolate',
                        help='isolate',
                        required=True,
                        type=str,
                        metavar="str")


    return parser.parse_args()


# --------------------------------------------------
def main():
    """Make a jazz noise here"""

    args = get_args()

    report_lines = []
    with open(args.input, "r") as infile1:
        for line in infile1:
            line_elements = line.strip().split("\t")
            if line_elements[3] == "G":
                # print(line_elements[0], line_elements[-1])
                gid = line_elements[-1].strip()
                gpercent = line_elements[0].strip()
                break
        for line in infile1:
            line_elements = line.strip().split("\t")
            if line_elements[3] == "S":
                # print(line_elements[0], line_elements[-1])
                sid = line_elements[-1].strip()
                spercent = line_elements[0].strip()
                break
            
    line_to_output = f'{args.isolate}\t{gpercent}\t{gid}\t{spercent}\t{sid}\n'

    if os.path.isfile(args.output):
        with open(args.output, "a") as outfile1:
            outfile1.write(line_to_output)
    else:
        header = f'Isolate\tGenus_Percent\tGenus_ID\tSpecies_Percent\tSpecies_ID\n'
        with open(args.output, "w") as outfile1:
            outfile1.write(header + line_to_output)
    
    

# --------------------------------------------------
if __name__ == '__main__':
    main()