#!/usr/bin/env python3
"""
Author : conrad <conrad@localhost>
Date   : 2023-06-22
Purpose: Fix headers to all have two ## symbols (bakta adds one in some places)
"""

import argparse


# --------------------------------------------------
def get_args():
    """Get command-line arguments"""

    parser = argparse.ArgumentParser(
        description='Rock the Casbah',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-ig',
                        '--input_gff',
                        help='input gff to fix',
                        metavar='FILE',
                        default=None)

    parser.add_argument('-og',
                        '--output_gff',
                        help='output gff to write',
                        metavar='FILE',
                        default=None)

    return parser.parse_args()


# --------------------------------------------------
def main():
    """Make a jazz noise here"""

    args = get_args()
    
    output_lines = []

    with open(args.input_gff, "r") as infile1:
        for line in infile1:
            if line.startswith("# "):
                newline = f"#{line}"
                output_lines.append(newline)
            else:
                output_lines.append(line)
    
    with open(args.output_gff, "w") as outfile1:
        outfile1.write("".join(output_lines))


# --------------------------------------------------
if __name__ == '__main__':
    main()
