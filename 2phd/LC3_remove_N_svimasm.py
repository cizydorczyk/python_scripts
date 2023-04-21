#!/usr/bin/env python3
"""
Author : Conrad Izydorczyk <conrad.izydorczyk@ucalgary.ca>
Date   : 2023-04-11
Purpose: remove insertions with Ns from ssvim-asm output
"""

import argparse


# --------------------------------------------------
def get_args():
    """Get command-line arguments"""

    parser = argparse.ArgumentParser(
        description='remove insertions with Ns from ssvim-asm output',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    parser.add_argument('-i',
                        '--input_vcf',
                        help='svim-asm vcf file',
                        metavar='str',
                        type=str,
                        required=True)
    
    parser.add_argument('-o',
                        '--output_vcf',
                        help='output vcf with insertions with Ns removed',
                        metavar='str',
                        type=str,
                        required=True)
    
    parser.add_argument('--isolate',
                        help='isolate/genome name/number',
                        metavar='str',
                        type=str,
                        required=True)

    return parser.parse_args()


# --------------------------------------------------
def main():
    """Remove insertion variants with Ns in them (even just a single N)"""

    args = get_args()
    
    header = []
    pos_to_keep = []
    with open(args.input_vcf, "r") as infile1:
        for line in infile1:
            if not line.startswith("#"):
                if line.strip().split("\t")[2].split(".")[1] == "INS":
                    if 'N' in line.strip().split("\t")[4]:
                        continue
                    else:
                        pos_to_keep.append(line.strip())
                else:
                    pos_to_keep.append(line.strip())
            else:
                header.append(line.strip())

    pos_to_keep_edited = []
    for i in pos_to_keep:
        var_elements = i.strip().split("\t")
        new_var_id = args.isolate + "." + var_elements[2]
        var_elements[2] = new_var_id
        pos_to_keep_edited.append("\t".join(var_elements))
        
    to_write = "\n".join(header) + "\n" + "\n".join(pos_to_keep_edited)

    with open(args.output_vcf, "w") as outfile1:
        outfile1.write(to_write)


# --------------------------------------------------
if __name__ == '__main__':
    main()
