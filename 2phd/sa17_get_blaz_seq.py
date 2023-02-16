#!/usr/bin/env python3
"""
Author : Conrad Izydorczyk <conrad.izydorczyk@ucalgary.ca>
Date   : 2022-10-29
Purpose: Get blaZ sequences from annotated bakta .ffn or .faa files.
"""

import argparse
import os.path
from Bio import SeqIO


# --------------------------------------------------
def get_args():
    """Get command-line arguments"""

    parser = argparse.ArgumentParser(
        description='Rock the Casbah',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-a',
                        '--amrfinderplus',
                        help='file with amrfinderplus output',
                        metavar='FILE',)

    parser.add_argument('-f',
                        '--fasta',
                        help='nucl (bakta .ffn) or prot (bakta .faa) fasta with isolate gene sequences',
                        metavar='int')

    parser.add_argument('-g',
                        '--gene',
                        help='amrfinderplus gene symbol(s) for target genes; must match exactly to what amrfinderplus uses!',
                        metavar='FILE',
                        nargs = "+")

    parser.add_argument('-o',
                        '--output',
                        help='output fasta',
                        metavar='FILE')

    parser.add_argument("-i",
                        "--isolate",
                        help="isolate name")

    return parser.parse_args()


# --------------------------------------------------
def main():
    """Get target gene sequence"""

    args = get_args()

    # Read in fasta seqs:
    fasta_seqs = list(SeqIO.parse(args.fasta, "fasta"))

    # Get locus tag for target seq from amrfinderplus file:
    target_loci = []
    with open(args.amrfinderplus, "r") as infile1:
        for line in infile1:
            line_elements = line.strip().split("\t")
            if line_elements[6] in args.gene:
                target_locus = line_elements[1]
                target_loci.append(target_locus)
    
    target_seqs = []
    for locus in target_loci:
        for seq in fasta_seqs:
            if locus in seq.id:
                seq.id = args.isolate + "_" + seq.id
                target_seqs.append(seq)
    
    SeqIO.write(target_seqs, args.output, "fasta")
                
            

    



def check_isfile(handle):
    if os.path.isfile(handle):
        return(True)
    else:
        return(False)


# --------------------------------------------------
if __name__ == '__main__':
    main()
